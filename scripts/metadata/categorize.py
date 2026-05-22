from __future__ import annotations

import json
from pathlib import Path
from typing import TypedDict

import pandas as pd

from metadata.config import MetadataConfig
from metadata.regexes import DISEASE_MAP


class AccessionDiseaseCategory(TypedDict):
    disease: str
    categories: list[str]


def disease_categories_for(diseaseValue: str) -> list[str]:
    """Return every DISEASE_MAP label whose pattern matches diseaseValue.

    DISEASE_MAP is a hierarchy where parents and children can both match (for example
    a LUAD string matches Lung Cancer, NSCLC, and LUAD). Labels are returned in
    DISEASE_MAP order (parent before child). Empty list when no label matches.
    """
    text = str(diseaseValue)
    return [label for label, pat in DISEASE_MAP if pat.search(text)]


def build_accession_disease_categories(samplesDf: pd.DataFrame) -> dict[str, AccessionDiseaseCategory]:
    """Build a {srx_accession: {disease, categories}} mapping from a samples frame.

    samplesDf must contain srx_accession and disease columns. The categories list
    contains every matching DISEASE_MAP label (empty list when nothing matches).
    """
    required = {"srx_accession", "disease"}
    missing = required - set(samplesDf.columns)
    if missing:
        raise KeyError(f"samplesDf missing required columns: {sorted(missing)}")

    out: dict[str, AccessionDiseaseCategory] = {}
    for _, row in samplesDf[["srx_accession", "disease"]].iterrows():
        srx = str(row["srx_accession"])
        disease = "" if pd.isna(row["disease"]) else str(row["disease"])
        out[srx] = {"disease": disease, "categories": disease_categories_for(disease)}
    return out


def export_accession_disease_categories(
    samplesDf: pd.DataFrame,
    cfg: MetadataConfig,
    filename: str = "accession_disease_categories.json",
    outputPath: Path | None = None,
) -> Path:
    """Write a JSON file mapping each accession to its disease string and matching DISEASE_MAP labels.

    Writes to cfg.outputDir / filename unless outputPath is provided. Accessions are
    sorted alphabetically in the output so diffs stay stable. Returns the written path.
    """
    mapping = build_accession_disease_categories(samplesDf)
    sortedMapping = {srx: mapping[srx] for srx in sorted(mapping)}

    if outputPath is None:
        cfg.outputDir.mkdir(parents=True, exist_ok=True)
        target = cfg.outputDir / filename
    else:
        outputPath.parent.mkdir(parents=True, exist_ok=True)
        target = outputPath

    with target.open("w") as f:
        json.dump(sortedMapping, f, indent=2)
    return target
