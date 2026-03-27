# Lung Metadata Analysis Report

**Dataset:** scBaseCount 2026-01-12 · Homo sapiens · `GeneFull`

---

## 1. Sample overview

Samples are filtered to those with **known, non-healthy labels** in both the `disease` and `tissue` fields (i.e. entries labelled normal, healthy, control, unknown, etc. are discarded). The remaining samples are split based on whether both the disease *and* tissue label are lung-related (**lung intersection**, 858 SRX) or either one is (**lung union**, 7,986 SRX). The intersection is the primary analysis set.

Of the 35,266 total unique SRX accessions, 16,716 (47.4%) are discarded as healthy/unknown. Of the remaining known samples, 858 carry a lung-specific disease *and* tissue label — 361 of which (42.1%) involve a cancer diagnosis.

### 1.1 Regex at a glance

Four regular expressions drive the filtering pipeline (all case-insensitive):

| Name | Purpose | Key terms matched |
|---|---|---|
| `NORMAL_HEALTHY_REGEX` | Discard healthy / control samples | `normal`, `healthy`, `control`, `unknown`, `wild-type`, `baseline`, … |
| `LUNG_DISEASE_RE` | Keep samples with a lung-related disease label | `lung`, `pulmonary`, `NSCLC`, `SCLC`, `fibrosis`, `COPD`, `COVID-19`, `cancer`, `carcinoma`, … |
| `LUNG_TISSUE_RE` | Keep samples with a lung-related tissue label | `lung`, `pulmonary`, `alveol*`, `bronch*`, `pleura`, `trachea`, `parenchyma`, … |
| `CANCER_RE` | Flag samples as cancer vs. non-cancer | `cancer`, `carcinoma`, `tumour`, `malignant`, `NSCLC`, `SCLC`, `mesothelioma`, … |



---

## 2. Lung disease breakdown

### 2.1 Raw label breakdown

The tables below show the top 10 free-text labels exactly as they appear in the dataset, before any normalisation is applied.

**Top 10 tissue labels**


| Tissue label                 | % of samples |
| ---------------------------- | ------------ |
| lung                         | 50.70        |
| lung tumor                   | 3.26         |
| lung tumour                  | 1.52         |
| lung adenocarcinoma          | 1.28         |
| lung tissue                  | 1.28         |
| lung tumor central margin    | 1.17         |
| lung (primary basal cells)   | 1.05         |
| lung biopsy                  | 1.05         |
| lung tumor subpleural margin | 0.93         |
| parietal pleura              | 0.70         |


**Top 10 disease labels**


| Disease label                       | % of samples |
| ----------------------------------- | ------------ |
| lung adenocarcinoma                 | 12.82        |
| idiopathic pulmonary fibrosis (IPF) | 8.51         |
| Idiopathic Pulmonary Fibrosis (IPF) | 7.46         |
| SARS-CoV-2 infection                | 4.66         |
| idiopathic pulmonary fibrosis       | 3.96         |
| pulmonary fibrosis                  | 2.91         |
| COPD                                | 2.80         |
| non-small cell lung cancer (NSCLC)  | 2.56         |
| lung adenocarcinoma (LUAD)          | 2.56         |
| carcinoma non-small cell            | 1.86         |


### 2.2 Normalised disease breakdown

The table below shows the regex rules used to map free-text disease labels into broad categories. Rules are applied in order — the first match wins. Labels that match none of the rules are placed in **Other**.

| Category | Pattern |
|---|---|
| IPF / Pulmonary Fibrosis | `pulmonary fibrosis \| IPF \| idiopathic pulmonary fibrosis` |
| COVID-19 / SARS-CoV-2 | `COVID \| SARS.CoV` |
| Lung Adenocarcinoma (LUAD) | `lung adenocarcinoma \| LUAD` |
| NSCLC | `NSCLC \| non.small.cell lung \| non.small.cell carcinoma \| carcinoma non.small` |
| Lung Squamous Cell Carcinoma | `squamous cell carcinoma of the lung \| lung squamous` |
| COPD | `\bCOPD\b` |
| Lung Cancer (general) | `lung cancer \| lung carcinoma \| MPLC \| KRAS.mutant lung \| SCLC` |
| Cystic Fibrosis | `cystic fibrosis` |
| Interstitial Lung Disease | `interstitial lung \| ILD \| SSc` |

Within the 858 lung-intersection samples, free-text disease labels are normalised into broad categories via regex. Categories representing < 2% are collapsed into **Other**. IPF / Pulmonary Fibrosis and Lung Adenocarcinoma together account for roughly 40% of the set.



---

## 3. Cell count distribution

The 858 lung-intersection SRX samples contain ~6.1M cells in total (mean 7,126 / median 5,125 per SRX). The distribution is right-skewed, with most samples in the 3,000–9,000 cell range and a small number of outliers above 40,000.

