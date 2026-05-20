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

# Lung-cancer regex components, composed below into nested labels.
# Hierarchy: Lung Cancer > {SCLC, NSCLC > {LUAD, LUSC, LCC}}.
_LUAD_RE = r"\bLUAD\b|lung adenocarcinoma|adenocarcinoma of the lung"
_LUSC_RE = r"\bLUSC\b|squamous cell carcinoma of the lung|lung squamous|squamous cell lung"
_LCC_RE = r"\bLCC\b|large[\s-]cell carcinoma|large[\s-]cell lung|lung large[\s-]cell"
_NSCLC_GENERIC_RE = r"\bNSCLC\b|non[\s-]small[\s-]cell lung|non[\s-]small[\s-]cell carcinoma|carcinoma non[\s-]small"
_NSCLC_RE = "|".join([_NSCLC_GENERIC_RE, _LUAD_RE, _LUSC_RE, _LCC_RE])
_SCLC_RE = r"\bSCLC\b|small[\s-]cell lung"
_LUNG_CANCER_GENERIC_RE = r"lung cancer|lung carcinoma|\bMPLC\b|KRAS.mutant lung"
_LUNG_CANCER_RE = "|".join([_LUNG_CANCER_GENERIC_RE, _SCLC_RE, _NSCLC_RE])

DISEASE_MAP: list[tuple[str, re.Pattern[str]]] = [
    ("IPF / Pulmonary Fibrosis", re.compile(r"pulmonary fibrosis|\bIPF\b|idiopathic pulmonary fibrosis", re.IGNORECASE)),
    ("COVID-19 / SARS-CoV-2", re.compile(r"\bCOVID\b|SARS.CoV", re.IGNORECASE)),
    ("Lung Cancer", re.compile(_LUNG_CANCER_RE, re.IGNORECASE)),
    ("Small Cell Lung Cancer (SCLC)", re.compile(_SCLC_RE, re.IGNORECASE)),
    ("Non-small Cell Lung Cancer (NSCLC)", re.compile(_NSCLC_RE, re.IGNORECASE)),
    ("Lung Adenocarcinoma (LUAD)", re.compile(_LUAD_RE, re.IGNORECASE)),
    ("Lung Squamous Cell Carcinoma (LUSC)", re.compile(_LUSC_RE, re.IGNORECASE)),
    ("Lung Large Cell Carcinoma (LCC)", re.compile(_LCC_RE, re.IGNORECASE)),
    ("COPD", re.compile(r"\bCOPD\b", re.IGNORECASE)),
    ("Cystic Fibrosis", re.compile(r"cystic fibrosis", re.IGNORECASE)),
    ("Interstitial Lung Disease", re.compile(r"interstitial lung|\bILD\b|\bSSc\b", re.IGNORECASE)),
]
