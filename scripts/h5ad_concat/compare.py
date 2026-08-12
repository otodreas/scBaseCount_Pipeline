import hashlib
import json
from collections.abc import Iterator
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from pydantic import BaseModel, Field
from shared.repo import REPO_ROOT, rel_to_repo

from h5ad_concat.models import H5adConcatResult


class AtlasCompareConfig(BaseModel):
    baselinePath: Path
    candidatePath: Path
    baselineManifestPath: Path | None = None
    candidateManifestPath: Path | None = None
    baselineStatusCsvPath: Path | None = None
    reportPath: Path
    chunkSize: int = Field(default=2048, ge=1)


class AtlasCompareReport(BaseModel):
    baselinePath: str
    candidatePath: str
    byteIdentical: bool
    fullLogicalIdentical: bool
    nBaselineCells: int
    nCandidateCells: int
    nCommonCells: int
    nCandidateOnlyCells: int
    nBaselineOnlyCells: int
    commonMatrixEqual: bool
    commonObsEqual: bool
    commonVarEqual: bool
    baselineOnlyDirectRiboRemovals: int = 0
    baselineOnlyWholeAccessionRemovals: int = 0
    baselineOnlyOther: int = 0
    fileAdmissionDiffs: dict[str, list[str]] = Field(default_factory=dict)
    notes: list[str] = Field(default_factory=list)


def _chunk_ranges(n_rows: int, chunk_size: int) -> Iterator[tuple[int, int]]:
    for start in range(0, n_rows, chunk_size):
        yield start, min(start + chunk_size, n_rows)


def _byte_identical(left: Path, right: Path, chunk_size: int = 1024 * 1024) -> bool:
    if left.stat().st_size != right.stat().st_size:
        return False
    with left.open("rb") as left_handle, right.open("rb") as right_handle:
        while True:
            left_chunk = left_handle.read(chunk_size)
            right_chunk = right_handle.read(chunk_size)
            if left_chunk != right_chunk:
                return False
            if not left_chunk:
                return True


def _to_dense(block: object) -> np.ndarray:
    toarray = getattr(block, "toarray", None)
    if callable(toarray):
        return np.asarray(toarray())
    return np.asarray(block)


def _matrix_rows_equal(left: ad.AnnData, right: ad.AnnData, obs_names: list[str], chunk_size: int) -> bool:
    if not obs_names:
        return True
    if left.X is None or right.X is None:
        return False
    left_idx = left.obs_names.get_indexer(obs_names)
    right_idx = right.obs_names.get_indexer(obs_names)
    if (left_idx < 0).any() or (right_idx < 0).any():
        return False
    for start, stop in _chunk_ranges(len(obs_names), chunk_size):
        left_block = _to_dense(left.X[left_idx[start:stop]])
        right_block = _to_dense(right.X[right_idx[start:stop]])
        if not np.array_equal(left_block, right_block):
            return False
    return True


def _obs_subset_equal(left: pd.DataFrame, right: pd.DataFrame, obs_names: list[str]) -> bool:
    if not obs_names:
        return True
    left_obs = left.loc[obs_names]
    right_obs = right.loc[obs_names]
    if list(left_obs.columns) != list(right_obs.columns):
        return False
    return left_obs.equals(right_obs)


def _accession_from_obs_name(obs_name: str) -> str:
    if "_" not in obs_name:
        return obs_name
    return obs_name.rsplit("_", 1)[-1]


def _load_status_map_from_csv(path: Path) -> dict[str, str]:
    frame = pd.read_csv(path, dtype="string")
    out: dict[str, str] = {}
    for accession, status, reason in frame[["accession", "status", "reason"]].itertuples(index=False, name=None):
        status_text = str(status)
        reason_text = str(reason or "")
        out[str(accession)] = status_text if status_text == "success" else f"skip:{reason_text}"
    return out


def _load_status_map_from_manifest(path: Path) -> dict[str, str]:
    result = H5adConcatResult.model_validate_json(path.read_text())
    out: dict[str, str] = {}
    for record in result.files:
        if record.status == "success":
            out[record.accession] = "success"
        else:
            reason = record.skipReason.value if record.skipReason is not None else ""
            out[record.accession] = f"skip:{reason}"
    return out


def _ribo_removed_cells(baseline: ad.AnnData, baseline_only: list[str], max_pct_ribo: float) -> set[str]:
    obs = pd.DataFrame(baseline.obs)
    if not baseline_only or "pct_counts_ribo" not in obs.columns:
        return set()
    threshold = max_pct_ribo * 100
    ribo = obs.loc[baseline_only, "pct_counts_ribo"]
    return set(ribo.index[ribo >= threshold].astype(str))


def compare_atlases(cfg: AtlasCompareConfig) -> AtlasCompareReport:
    """Compare baseline and candidate atlases with bounded-memory matrix checks."""
    notes: list[str] = []
    byte_identical = _byte_identical(cfg.baselinePath, cfg.candidatePath)

    baseline = ad.read_h5ad(cfg.baselinePath, backed="r")
    candidate = ad.read_h5ad(cfg.candidatePath, backed="r")
    try:
        baseline_names = set(map(str, baseline.obs_names))
        candidate_names = set(map(str, candidate.obs_names))
        common = sorted(baseline_names & candidate_names)
        candidate_only = sorted(candidate_names - baseline_names)
        baseline_only = sorted(baseline_names - candidate_names)

        common_var_equal = list(baseline.var_names) == list(candidate.var_names)
        common_matrix_equal = common_var_equal and _matrix_rows_equal(baseline, candidate, common, cfg.chunkSize)
        common_obs_equal = _obs_subset_equal(pd.DataFrame(baseline.obs), pd.DataFrame(candidate.obs), common)

        full_logical_identical = (
            not candidate_only and not baseline_only and common_matrix_equal and common_obs_equal and common_var_equal
        )

        file_admission_diffs: dict[str, list[str]] = {
            "onlyInBaseline": [],
            "onlyInCandidate": [],
            "statusChanged": [],
        }

        baseline_status: dict[str, str] = {}
        candidate_status: dict[str, str] = {}
        if cfg.baselineStatusCsvPath is not None and cfg.baselineStatusCsvPath.exists():
            baseline_status = _load_status_map_from_csv(cfg.baselineStatusCsvPath)
        elif cfg.baselineManifestPath is not None and cfg.baselineManifestPath.exists():
            try:
                baseline_status = _load_status_map_from_manifest(cfg.baselineManifestPath)
            except Exception as exc:
                notes.append(f"Could not parse baseline manifest as post-cutover result: {exc}")
                sibling_csv = cfg.baselinePath.with_suffix(".csv")
                if sibling_csv.exists():
                    baseline_status = _load_status_map_from_csv(sibling_csv)

        if cfg.candidateManifestPath is not None and cfg.candidateManifestPath.exists():
            candidate_status = _load_status_map_from_manifest(cfg.candidateManifestPath)

        if baseline_status or candidate_status:
            all_accessions = sorted(set(baseline_status) | set(candidate_status))
            for accession in all_accessions:
                left = baseline_status.get(accession)
                right = candidate_status.get(accession)
                if left is None:
                    file_admission_diffs["onlyInCandidate"].append(accession)
                elif right is None:
                    file_admission_diffs["onlyInBaseline"].append(accession)
                elif left != right:
                    file_admission_diffs["statusChanged"].append(f"{accession}:{left}->{right}")

        max_pct_ribo = 0.5
        if cfg.candidateManifestPath is not None and cfg.candidateManifestPath.exists():
            candidate_payload = json.loads(cfg.candidateManifestPath.read_text())
            config_path = candidate_payload.get("configPath")
            if config_path:
                resolved = Path(config_path)
                if not resolved.is_absolute():
                    resolved = REPO_ROOT / resolved
                if resolved.exists():
                    max_pct_ribo = float(json.loads(resolved.read_text()).get("maxPctRibo", 0.5))

        ribo_removed = _ribo_removed_cells(baseline, baseline_only, max_pct_ribo)
        baseline_only_by_accession: dict[str, list[str]] = {}
        for obs_name in baseline_only:
            baseline_only_by_accession.setdefault(_accession_from_obs_name(obs_name), []).append(obs_name)

        candidate_accessions = {_accession_from_obs_name(name) for name in map(str, candidate.obs_names)}
        direct_ribo = 0
        whole_accession = 0
        other = 0
        for accession, cells in baseline_only_by_accession.items():
            if accession not in candidate_accessions:
                whole_accession += len(cells)
            else:
                for cell in cells:
                    if cell in ribo_removed:
                        direct_ribo += 1
                    else:
                        other += 1

        if candidate_only:
            notes.append("Candidate-only cells are unexpected for a stricter ribosomal filter")

        report = AtlasCompareReport(
            baselinePath=rel_to_repo(cfg.baselinePath),
            candidatePath=rel_to_repo(cfg.candidatePath),
            byteIdentical=byte_identical,
            fullLogicalIdentical=full_logical_identical,
            nBaselineCells=len(baseline_names),
            nCandidateCells=len(candidate_names),
            nCommonCells=len(common),
            nCandidateOnlyCells=len(candidate_only),
            nBaselineOnlyCells=len(baseline_only),
            commonMatrixEqual=common_matrix_equal,
            commonObsEqual=common_obs_equal,
            commonVarEqual=common_var_equal,
            baselineOnlyDirectRiboRemovals=direct_ribo,
            baselineOnlyWholeAccessionRemovals=whole_accession,
            baselineOnlyOther=other,
            fileAdmissionDiffs=file_admission_diffs,
            notes=notes,
        )
    finally:
        baseline.file.close()
        candidate.file.close()

    cfg.reportPath.parent.mkdir(parents=True, exist_ok=True)
    cfg.reportPath.write_text(json.dumps(report.model_dump(mode="json"), indent=2) + "\n")
    return report


def file_digest(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()
