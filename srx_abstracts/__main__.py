from __future__ import annotations

import argparse
import json
import sys

from srx_abstracts.entrez_client import EntrezClient
from srx_abstracts.pipeline import pipeline_for_srx_list


def main() -> None:
    parser = argparse.ArgumentParser(description="Resolve SRX accessions to PubMed abstracts via Entrez.")
    parser.add_argument(
        "srx",
        nargs="+",
        help="SRA experiment accessions (e.g. SRX12345678)",
    )
    args = parser.parse_args()

    with EntrezClient() as client:
        bundles = pipeline_for_srx_list(client, args.srx)

    payload = [b.model_dump(mode="json") for b in bundles]
    json.dump(payload, sys.stdout, indent=2)
    sys.stdout.write("\n")


if __name__ == "__main__":
    main()
