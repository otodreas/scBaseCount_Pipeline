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
        r"|"
        r"\bno\s+(?:"
        r"disease|COPD|"
        r"diagnosed\s+disease|specific\s+disease|overt\s+disease|donor\s+disease|"
        r"disease\s+diagnosis|record\s+of\s+lung\s+disease"
        r")\b"
        r"|"
        r"\bnon[-\s]?(?:disease|COPD)\b"
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

# Lung-cancer regex components, composed below into nested labels.
# Hierarchy: Lung Cancer > {SCLC, NSCLC > {LUAD, LUSC, LCC}}.
_LUAD_RE = r"\bLUAD\b|lung adenocarcinoma|adenocarcinoma of the lung"
_LUSC_RE = (
    r"\bLUSC\b|squamous cell carcinoma of the lung|lung squamous|squamous cell lung"
    r"|non[\s-]small[\s-]cell\s+squamous"
)
_LCC_RE = r"\bLCC\b|large[\s-]cell carcinoma|large[\s-]cell lung|lung large[\s-]cell"
_NSCLC_GENERIC_RE = r"\bNSCLC\b|non[\s-]small[\s-]cell lung|non[\s-]small[\s-]cell carcinoma|carcinoma non[\s-]small"
_NSCLC_RE = "|".join([_NSCLC_GENERIC_RE, _LUAD_RE, _LUSC_RE, _LCC_RE])
_SCLC_RE = r"\bSCLC\b|(?<!non[\s-])(?<!non)small[\s-]cell lung"
_CANCER_TERMS = r"(?:cancer|carcinoma|tumou?r|malignancy|neoplasm)"
_LUNG_CANCER_GENERIC_RE = (
    rf"lung {_CANCER_TERMS}"
    rf"|{_CANCER_TERMS} of the lung"
    r"|\bMPLC\b|KRAS.mutant lung"
)
_LUNG_CANCER_RE = "|".join([_LUNG_CANCER_GENERIC_RE, _SCLC_RE, _NSCLC_RE])
LUNG_CANCER_RE = re.compile(_LUNG_CANCER_RE, re.IGNORECASE)

NON_CF_RE = re.compile(r"\bnon[-\s]?CF\b|non[-\s]?cystic\s+fibrosis", re.IGNORECASE)

DISEASE_MAP: list[tuple[str, re.Pattern[str]]] = [
    ("IPF / Pulmonary Fibrosis", re.compile(r"pulmonary fibrosis|\bIPF\b|idiopathic pulmonary fibrosis", re.IGNORECASE)),
    ("COVID-19 / SARS-CoV-2", re.compile(r"\bCOVID\b|SARS.CoV", re.IGNORECASE)),
    ("Lung Cancer", LUNG_CANCER_RE),
    ("Small Cell Lung Cancer (SCLC)", re.compile(_SCLC_RE, re.IGNORECASE)),
    ("Non-small Cell Lung Cancer (NSCLC)", re.compile(_NSCLC_RE, re.IGNORECASE)),
    ("Lung Adenocarcinoma (LUAD)", re.compile(_LUAD_RE, re.IGNORECASE)),
    ("Lung Squamous Cell Carcinoma (LUSC)", re.compile(_LUSC_RE, re.IGNORECASE)),
    ("Lung Large Cell Carcinoma (LCC)", re.compile(_LCC_RE, re.IGNORECASE)),
    ("COPD", re.compile(r"\bCOPD\b", re.IGNORECASE)),
    ("Cystic Fibrosis", re.compile(r"cystic fibrosis", re.IGNORECASE)),
    ("Interstitial Lung Disease", re.compile(r"interstitial lung|\bILD\b|\bSSc\b", re.IGNORECASE)),
    (
        "Pulmonary Hypertension",
        re.compile(
            r"\bPAH\b|\bIPAH\b|\bCTEPH\b|"
            r"pulmonary arterial hypertension|"
            r"pulmonary hypertension|"
            r"thromboembolic pulmonary",
            re.IGNORECASE,
        ),
    ),
]
