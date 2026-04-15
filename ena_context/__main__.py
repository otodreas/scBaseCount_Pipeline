from __future__ import annotations

import argparse
import json
import sys

from ena_context.fetch import fetch_experiment_context


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Fetch experiment context from ENA for one or more SRX/ERX accessions."
    )
    parser.add_argument(
        "accessions",
        nargs="+",
        help="ENA experiment accessions (e.g. SRX12345678 or ERX12345678)",
    )
    args = parser.parse_args()

    results = [fetch_experiment_context(acc) for acc in args.accessions]
    payload = [r.model_dump(mode="json") for r in results]
    json.dump(payload, sys.stdout, indent=2)
    sys.stdout.write("\n")


if __name__ == "__main__":
    main()
