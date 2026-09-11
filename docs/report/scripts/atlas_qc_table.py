import argparse
import json
from collections import Counter
from pathlib import Path
from typing import Any

_FILE_LEVEL_REASONS = (
    "md5_mismatch",
    "cell_type_all_missing",
    "too_few_cells",
    "excessive_cell_dropout",
)
_PER_CELL_FILTERS = ("minGenesPerCell", "maxPctMito", "maxPctRibo")


def main() -> None:
    parser = argparse.ArgumentParser(description="Write the atlas QC table LaTeX fragment.")
    parser.add_argument("-i", "--input", type=Path, required=True, help="atlas_result.json from concatenation")
    parser.add_argument("-o", "--output", type=Path, required=True, help="output LaTeX path")
    args = parser.parse_args()
    table = build_table(load_result(args.input))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(table)
    print(f"Wrote {args.output}")


def load_result(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(f"atlas result not found: {path}")
    payload = json.loads(path.read_text())
    required = ("nObs", "nFilesConcatenated", "nFilesSkipped", "studiesSeen", "qcSummary", "skipped")
    missing = [key for key in required if key not in payload]
    if missing:
        raise ValueError(f"{path} is missing required keys: {', '.join(missing)}")
    return payload


def build_table(result: dict[str, Any]) -> str:
    n_concat = int(result["nFilesConcatenated"])
    n_skipped = int(result["nFilesSkipped"])
    n_obs = int(result["nObs"])
    n_studies = len(result["studiesSeen"])
    qc_summary = result["qcSummary"]
    all_qc = _cohort(qc_summary, "allQcProcessedFiles")
    concat_qc = _cohort(qc_summary, "concatenatedFiles")
    n_eligible = int(all_qc["nFiles"])
    n_barcodes_eligible = int(all_qc["nCellsBefore"])
    if n_eligible != n_concat + n_skipped:
        raise ValueError(f"file counts disagree: eligible={n_eligible} concatenated={n_concat} skipped={n_skipped}")
    if int(concat_qc["nFiles"]) != n_concat:
        raise ValueError("concatenatedFiles.nFiles does not match nFilesConcatenated")
    if int(concat_qc["nCellsAfter"]) != n_obs:
        raise ValueError("concatenatedFiles.nCellsAfter does not match nObs")

    skip_counts, skip_barcodes = _skip_totals(result["skipped"])
    if sum(skip_counts.values()) != n_skipped:
        raise ValueError("sum of skip-reason counts does not match nFilesSkipped")

    dropped_datasets = n_eligible - n_concat
    dropped_barcodes = n_barcodes_eligible - n_obs
    per_cell = all_qc["nCellsDroppedByFilter"]

    rows = [
        _row(
            "MD5 hash (mismatch)",
            skip_counts["md5_mismatch"],
            skip_barcodes["md5_mismatch"],
            n_eligible,
            n_barcodes_eligible,
        ),
        _row(
            r"CZ CELLxGENE \texttt{cell_type} (labels absent)",
            skip_counts["cell_type_all_missing"],
            skip_barcodes["cell_type_all_missing"],
            n_eligible,
            n_barcodes_eligible,
        ),
        _row(
            r"Min.\ genes per cell ($<$ 200)",
            None,
            int(per_cell["minGenesPerCell"]),
            n_eligible,
            n_barcodes_eligible,
        ),
        _row(
            r"Max.\ mitochondrial fraction ($\geq$ 20\%)",
            None,
            int(per_cell["maxPctMito"]),
            n_eligible,
            n_barcodes_eligible,
        ),
        _row(
            r"Max.\ ribosomal fraction ($\geq$ 50\%)",
            None,
            int(per_cell["maxPctRibo"]),
            n_eligible,
            n_barcodes_eligible,
        ),
        _row(
            r"Min.\ cells remaining ($<$ 100)",
            skip_counts["too_few_cells"],
            skip_barcodes["too_few_cells"],
            n_eligible,
            n_barcodes_eligible,
        ),
        _row(
            r"Min.\ fraction of cells remaining ($<$ 50\%)",
            skip_counts["excessive_cell_dropout"],
            skip_barcodes["excessive_cell_dropout"],
            n_eligible,
            n_barcodes_eligible,
        ),
    ]

    caption = (
        rf"\textbf{{Atlas contained {_fmt_int(n_concat)} datasets ({_fmt_int(n_obs)} barcodes) "
        rf"across {_fmt_int(n_studies)} BioProject accessions, of {_fmt_int(n_eligible)} eligible "
        rf"datasets ({_fmt_int(n_barcodes_eligible)} barcodes) following filtering.}}"
    )
    footer = (
        f"Dataset percentages use the {_fmt_int(n_eligible)} eligible datasets as the denominator. "
        f"Barcode percentages use the {_fmt_int(n_barcodes_eligible)} barcodes present in those datasets before QC. "
        "Per-cell filter counts include every file that reached QC, including files later rejected by file-level gates; "
        f"they are not restricted to concatenated files. A total of {_fmt_int(dropped_datasets)} datasets were dropped. "
        f"Of those, {_fmt_int(skip_counts['md5_mismatch'])} were dropped due to MD5 hash mismatch, "
        f"{_fmt_int(skip_counts['cell_type_all_missing'])} "
        f"({_pct(skip_counts['cell_type_all_missing'], n_eligible)}\\% total datasets) were dropped due to fully absent "
        r"CZ CELLxGENE \texttt{cell_type} labels. "
        f"Of the {_fmt_int(dropped_barcodes)} barcodes that were dropped from the eligible set, "
        f"{_fmt_int(int(per_cell['minGenesPerCell']))} "
        f"({_pct(int(per_cell['minGenesPerCell']), n_barcodes_eligible)}\\% total barcodes) were dropped for gene counts "
        f"of less than 200, {_fmt_int(int(per_cell['maxPctMito']))} "
        f"({_pct(int(per_cell['maxPctMito']), n_barcodes_eligible)}\\% total barcodes) were dropped for having 20\\% or more "
        f"of total gene count being mitochondrial, and {_fmt_int(int(per_cell['maxPctRibo']))} "
        f"({_pct(int(per_cell['maxPctRibo']), n_barcodes_eligible)}\\% total barcodes) for having 50\\% or more of total "
        f"gene counts being ribosomal. Following cell filtering, {_fmt_int(skip_counts['too_few_cells'])} "
        f"({_pct(skip_counts['too_few_cells'], n_eligible)}\\% total datasets) were dropped for retaining fewer than 100 cells, "
        f"and {_fmt_int(skip_counts['excessive_cell_dropout'])} "
        f"({_pct(skip_counts['excessive_cell_dropout'], n_eligible)}\\% total datasets) were dropped for retaining less than "
        f"50\\% of cells. Among concatenated files only, cell filters dropped {_fmt_int(int(concat_qc['nCellsDropped']))} barcodes."
    )

    body = "\n".join(rows)
    return (
        "\\begin{table}[!ht]\n"
        "\\centering\n"
        f"\\caption{{{caption}}}\n"
        "\\begin{tabular}{|l|l|l|}\n"
        "\\hline\n"
        "Filter & Datasets dropped, $n$ (\\%) & Barcodes dropped, $n$ (\\%)\\\\ \\thickhline\n"
        f"{body}\n"
        "\\end{tabular}\n"
        f"\\begin{{flushleft}} {footer}\n"
        "\\end{flushleft}\n"
        "\\label{table1}\n"
        "\\end{table}\n"
    )


def _cohort(qc_summary: dict[str, Any], key: str) -> dict[str, Any]:
    if key not in qc_summary:
        raise ValueError(f"qcSummary is missing {key}")
    cohort = qc_summary[key]
    required = ("nFiles", "nCellsBefore", "nCellsAfter", "nCellsDropped", "nCellsDroppedByFilter")
    missing = [name for name in required if name not in cohort]
    if missing:
        raise ValueError(f"qcSummary.{key} is missing: {', '.join(missing)}")
    filters = cohort["nCellsDroppedByFilter"]
    missing_filters = [name for name in _PER_CELL_FILTERS if name not in filters]
    if missing_filters:
        raise ValueError(f"qcSummary.{key}.nCellsDroppedByFilter is missing: {', '.join(missing_filters)}")
    return cohort


def _skip_totals(skipped: list[dict[str, Any]]) -> tuple[Counter[str], dict[str, int]]:
    counts: Counter[str] = Counter({reason: 0 for reason in _FILE_LEVEL_REASONS})
    barcodes = {reason: 0 for reason in _FILE_LEVEL_REASONS}
    unknown: set[str] = set()
    for record in skipped:
        reason = str(record["reason"])
        if reason not in _FILE_LEVEL_REASONS:
            unknown.add(reason)
            continue
        counts[reason] += 1
        qc = record.get("qc")
        if reason == "md5_mismatch":
            if qc is not None:
                raise ValueError("md5_mismatch skip unexpectedly includes QC stats")
            continue
        if qc is None:
            raise ValueError(f"{reason} skip is missing QC stats")
        barcodes[reason] += int(qc["nCellsAfter"])
    if unknown:
        raise ValueError(f"unhandled skip reasons: {', '.join(sorted(unknown))}")
    return counts, barcodes


def _row(label: str, n_datasets: int | None, n_barcodes: int, n_eligible: int, n_barcodes_eligible: int) -> str:
    if n_datasets is None:
        datasets = "---"
    else:
        datasets = f"{_fmt_int(n_datasets)} ({_pct(n_datasets, n_eligible)})"
    barcodes = f"{_fmt_int(n_barcodes)} ({_pct(n_barcodes, n_barcodes_eligible)})"
    return f"{label} & {datasets} & {barcodes}\\\\ \\hline"


def _fmt_int(n: int) -> str:
    return f"{n:,}"


def _pct(n: int, denom: int) -> str:
    if n == 0:
        return "0"
    value = 100.0 * n / denom
    if value >= 1:
        return f"{value:.1f}"
    return f"{value:.2f}"


if __name__ == "__main__":
    main()
