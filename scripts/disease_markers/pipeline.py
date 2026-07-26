import datetime
import time
from collections.abc import Callable

from shared.logger import add_stdout_handler, configure_file_logger, log_run_separator
from shared.repo import rel_to_repo

from disease_markers.annotate import annotate_atlas
from disease_markers.config import DiseaseMarkersConfig
from disease_markers.consistency import flag_study_consistency
from disease_markers.de import de_areas_in_cluster
from disease_markers.labels import build_sample_label_table
from disease_markers.outputs import (
    write_area_cluster_counts,
    write_eligibility_labels,
    write_marker_tables,
    write_summary_json,
)
from disease_markers.pseudobulk import iter_cluster_pseudobulks
from disease_markers.transfer import load_full_atlas_transfer_clusters

_LOG_FILENAME = "disease_markers.log"
log = configure_file_logger(_LOG_FILENAME, __name__)
add_stdout_handler()


def _timed[T](name: str, action: Callable[[], T]) -> T:
    """Run action, logging START/DONE with elapsed seconds, and return its result."""
    started = time.perf_counter()
    log.info("START %s", name)
    result = action()
    log.info("DONE %s in %.1fs", name, time.perf_counter() - started)
    return result


def run_disease_markers(cfg: DiseaseMarkersConfig, *, runDe: bool = True, labelsOnly: bool = False) -> None:
    """Run cluster transfer, annotation, optional pseudobulk DE, and write outputs."""
    log_run_separator(log)
    log.info("disease markers run started")
    log.info("config: %s", cfg.model_dump_json())

    label_table = _timed(
        "sample labels",
        lambda: build_sample_label_table(cfg.contextsPath, cfg.atlasCsvPath),
    )
    elig_path = write_eligibility_labels(label_table, cfg.outputDir)
    log.info("Wrote %s", rel_to_repo(elig_path))

    if labelsOnly:
        log.info("labels-only run complete")
        return

    adata = _timed(
        "cluster transfer",
        lambda: load_full_atlas_transfer_clusters(
            cfg.inputAtlasH5ad,
            cfg.harmonyAtlasH5ad,
            clusterKey=cfg.clusterKey,
        ),
    )

    if cfg.writeTransferredAtlas and cfg.transferredAtlasH5ad is not None:
        cfg.transferredAtlasH5ad.parent.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(cfg.transferredAtlasH5ad, compression=cfg.compression)
        log.info("Wrote %s", rel_to_repo(cfg.transferredAtlasH5ad))

    adata = _timed("annotate + eligibility subset", lambda: annotate_atlas(adata, label_table, cfg))
    counts_path = write_area_cluster_counts(adata, cfg, cfg.outputDir)
    log.info("Wrote %s", rel_to_repo(counts_path))

    marker_results: dict[str, dict[str, object]] = {}
    marker_paths: list = []

    if runDe:

        def _de_pass() -> dict[str, dict[str, object]]:
            results: dict[str, dict[str, object]] = {}
            for cluster, pdata in iter_cluster_pseudobulks(adata, cfg):
                area_tables = de_areas_in_cluster(pdata, cfg)
                for area, hits in list(area_tables.items()):
                    area_tables[area] = flag_study_consistency(pdata, area, hits, cfg)
                if area_tables:
                    results[cluster] = area_tables
            return results

        marker_results = _timed("pseudobulk + DE", _de_pass)
        marker_paths = write_marker_tables(marker_results, cfg.outputDir)
        log.info("Wrote %d marker tables under %s", len(marker_paths), rel_to_repo(cfg.outputDir / "markers"))

    summary_path = write_summary_json(cfg, label_table, adata, marker_paths, cfg.outputDir)
    log.info("Wrote %s", rel_to_repo(summary_path))
    log.info("disease markers run complete in %s", datetime.datetime.now())
