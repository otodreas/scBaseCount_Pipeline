from __future__ import annotations

import re

NORMAL_HEALTHY_RE = re.compile(
    (
        r"\b(?:"
        r"normal|healthy|control|none|unknown|unstimulated|naive|"
        r"uninvolved|unaffected|unexposed|vehicle|wild[-_\s]?type|"
        r"wt|no treatment|baseline|"
        r"unsure|not specified|not stated|not reported|not available"
        r")\b"
    ),
    re.IGNORECASE,
)

LUNG_DISEASE_RE = re.compile(
    (
        r"\b(?:"
        r"lung|pulmonary|NSCLC|SCLC|adenocarcinoma of the lung|mesothelioma|pleura|"
        r"squamous cell carcinoma of the lung|large[- ]cell|bronch|EGFR|(?<!\w)ALK(?!\w)|ROS1|"
        r"thymoma|thymic|"
        r"cancer|carcinoma|adenocarcinoma|sarcoma|tumou?r|malignancy|"
        r"malignant|benign|neoplasm|mass|adenoma|adenomatous|lesion|"
        r"disease|fibrosis|emphysema|asthma|COPD|bronchitis|pneumonia|"
        r"infection|metastasis|injury"
        r")\b"
        r"|"
        r"\b(?:"
        r"emphysema|asthma|COPD|bronchitis|pneumonia|pulmonary fibrosis|"
        r"interstitial lung disease|"
        r"COVID-19|COVID19|COVID 19|sars-cov-2"
        r")\b"
    ),
    re.IGNORECASE,
)

LUNG_TISSUE_RE = re.compile(
    (
        r"\b(?:"
        r"lung|pulmonary|alveolus|alveoli|bronchus|bronchi|bronchiole|"
        r"airway|respiratory tract|bronchiolar|thoracic|trachea|pleura|"
        r"interstitium|parenchyma"
        r")\b"
    ),
    re.IGNORECASE,
)

CANCER_RE = re.compile(
    (
        r"\b(?:"
        r"cancer|carcinoma|adenocarcinoma|sarcoma|tumou?r|malignancy|malignant|"
        r"neoplasm|neoplastic|leukemia|leukaemia|lymphoma|myeloma|melanoma|"
        r"NSCLC|SCLC|mesothelioma|glioblastoma|glioma|blastoma|"
        r"metastasis|metastatic|oncolog"
        r")\b"
    ),
    re.IGNORECASE,
)

DISEASE_MAP: list[tuple[str, re.Pattern[str]]] = [
    ("IPF / Pulmonary Fibrosis",    re.compile(r"pulmonary fibrosis|IPF|idiopathic pulmonary fibrosis", re.IGNORECASE)),
    ("COVID-19 / SARS-CoV-2",       re.compile(r"COVID|SARS.CoV", re.IGNORECASE)),
    ("Lung Adenocarcinoma (LUAD)",   re.compile(r"lung adenocarcinoma|LUAD", re.IGNORECASE)),
    ("NSCLC",                        re.compile(r"NSCLC|non.small.cell lung|non.small.cell carcinoma|carcinoma non.small", re.IGNORECASE)),
    ("Lung Squamous Cell Carcinoma", re.compile(r"squamous cell carcinoma of the lung|lung squamous", re.IGNORECASE)),
    ("COPD",                         re.compile(r"\bCOPD\b", re.IGNORECASE)),
    ("Lung Cancer (general)",        re.compile(r"lung cancer|lung carcinoma|MPLC|KRAS.mutant lung|SCLC", re.IGNORECASE)),
    ("Cystic Fibrosis",              re.compile(r"cystic fibrosis", re.IGNORECASE)),
    ("Interstitial Lung Disease",    re.compile(r"interstitial lung|ILD|SSc", re.IGNORECASE)),
]
