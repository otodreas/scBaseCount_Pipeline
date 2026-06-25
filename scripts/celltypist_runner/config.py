from __future__ import annotations

from pydantic import BaseModel


class CellTypistRunnerConfig(BaseModel):
    modelName: str = "Nuclei_Lung_Airway.pkl"
    targetSum: int = 10_000
    majorityVoting: bool = False
    geneSymbolCol: str = "gene_symbols"
    predictedLabelKey: str = "predicted_labels"
    downloadIfMissing: bool = True
