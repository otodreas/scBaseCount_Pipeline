from __future__ import annotations

from pydantic import BaseModel


class CellTypistRunnerConfig(BaseModel):
    modelName: str = "Nuclei_Lung_Airway.pkl"
    # CellTypist expects normalization to 10,000 total counts per cell.
    targetSum: int = 10_000
    geneSymbolCol: str = "gene_symbols"
    predictedLabelKey: str = "predicted_labels"
    downloadIfMissing: bool = True
