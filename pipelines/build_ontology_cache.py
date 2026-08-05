"""Build release-pinned MONDO/UBERON term caches from OLS for atlas disease labeling."""

import argparse
from pathlib import Path

import pandas as pd
from ontology_lookup import OntologyCache, OntologyLookupConfig, ontology_tokens
from shared.repo import REPO_ROOT


def _collect_ids(path: Path, column: str) -> list[str]:
    frame = pd.read_csv(path, dtype="string", usecols=[column])
    seen: set[str] = set()
    out: list[str] = []
    for value in frame[column].tolist():
        for token in ontology_tokens(None if value is None else str(value)):
            if token not in seen:
                seen.add(token)
                out.append(token)
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--datasets-csv",
        type=Path,
        default=REPO_ROOT / "output" / "metadata" / "datasets_v2.csv",
    )
    args = parser.parse_args()
    if not args.datasets_csv.exists():
        raise FileNotFoundError(f"datasets CSV not found: {args.datasets_csv}")

    cfg = OntologyLookupConfig()
    mondo_ids = [c for c in _collect_ids(args.datasets_csv, "disease_ontology_term_id") if c.startswith("MONDO:")]
    uberon_ids = [c for c in _collect_ids(args.datasets_csv, "tissue_ontology_term_id") if c.startswith("UBERON:")]
    # Area and respiratory roots used by classifiers, even if absent from the CSV.
    mondo_ids = list(
        dict.fromkeys(
            [
                *mondo_ids,
                "MONDO:0800504",
                "MONDO:0002771",
                "MONDO:0002429",
                "MONDO:0100096",
                "MONDO:0100320",
                "MONDO:0008903",
                "MONDO:0005002",
                "MONDO:0009061",
                "MONDO:0005149",
                "MONDO:0015925",
            ]
        )
    )
    uberon_ids = list(dict.fromkeys([*uberon_ids, "UBERON:0001004", "UBERON:0002048"]))

    mondo = OntologyCache.for_mondo(cfg)
    uberon = OntologyCache.for_uberon(cfg)
    print(f"Fetching {len(mondo_ids)} MONDO terms for release {cfg.mondoRelease}", flush=True)
    mondo.ensure(mondo_ids, allowNetwork=True)
    print(f"Fetching {len(uberon_ids)} UBERON terms for release {cfg.uberonRelease}", flush=True)
    uberon.ensure(uberon_ids, allowNetwork=True)
    print(f"Wrote {mondo.termsPath}", flush=True)
    print(f"Wrote {uberon.termsPath}", flush=True)


if __name__ == "__main__":
    main()
