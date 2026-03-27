# 1. Sample overview

## 1.1 Samples by disease

Samples are filtered to those with **known, non-healthy labels** in both the `disease` and `tissue` fields (i.e. entries labelled normal, healthy, control, unknown, etc. are discarded). The remaining samples are split based on whether both the disease *and* tissue label are lung-related (**lung intersection**, 858 SRX) or either one is (**lung union**, 7,986 SRX). The intersection is the primary analysis set.


| Stage                                                           | SRX     | Notes                                      |
| --------------------------------------------------------------- | ------- | ------------------------------------------ |
| Total unique SRX                                                | 35,266  | Raw dataset                                |
| Healthy / unknown (discarded)                                   | 16,716  | Matched by `NORMAL_HEALTHY_REGEX` (47.4 %) |
| Known, non-healthy (`sample_known`)                             | 18,550  | Remainder after discarding                 |
| Lung union (disease **or** tissue is lung-related)              | 7,986   | At least one lung-related label            |
| **Lung intersection** (disease **and** tissue are lung-related) | **858** | Primary analysis set                       |
| └─ Cancer subset                                                | 361     | 42.1 % of intersection                     |
| Not in lung intersection                                        | 17,692  | Non-lung known samples (see §2.3)          |


Of the 35,266 total unique SRX accessions, 16,716 (47.4%) are discarded as healthy/unknown. Of the remaining known samples, 858 carry a lung-specific disease *and* tissue label, 361 of which (42.1%) involve a cancer diagnosis. The broader lung union (either label is lung-related) contains 7,986 samples, but the stricter intersection is used as the primary analysis set.

## 1.2 Regex at a glance

Four regular expressions drive the filtering pipeline (all case-insensitive):


| Name                   | Purpose                                        | Key terms matched                                                                              |
| ---------------------- | ---------------------------------------------- | ---------------------------------------------------------------------------------------------- |
| `NORMAL_HEALTHY_REGEX` | Discard healthy / control samples              | `normal`, `healthy`, `control`, `unknown`, `wild-type`, `baseline`, …                          |
| `LUNG_DISEASE_RE`      | Keep samples with a lung-related disease label | `lung`, `pulmonary`, `NSCLC`, `SCLC`, `fibrosis`, `COPD`, `COVID-19`, `cancer`, `carcinoma`, … |
| `LUNG_TISSUE_RE`       | Keep samples with a lung-related tissue label  | `lung`, `pulmonary`, `alveol`*, `bronch`*, `pleura`, `trachea`, `parenchyma`, …                |
| `CANCER_RE`            | Flag samples as cancer vs. non-cancer          | `cancer`, `carcinoma`, `tumour`, `malignant`, `NSCLC`, `SCLC`, `mesothelioma`, …               |


---

# 2. Lung disease breakdown

## 2.1 Raw label breakdown

### 2.1.1 Top excluded labels

The 17,692 samples in `sample_known` that did **not** pass the lung-intersection filter (i.e. they lacked a lung-specific disease label, a lung-specific tissue label, or both). The tables below show the most common free-text labels in that excluded set, confirming the filter is working as expected.

Notice that there are a good number of COVID-19 disease labeled data that get discarded due to tissue filtering. These are likely not labeled `lung` or similar in the tissue field. 

TODO: consider adding rule to let through specific rows where labels only match for disease

**Top 10 disease labels (excluded)**


| Disease label                | % of excluded samples |
| ---------------------------- | --------------------- |
| chronic myelogenous leukemia | 4.47                  |
| COVID-19                     | 2.49                  |
| acute T cell leukemia        | 2.33                  |
| Chronic Myelogenous Leukemia | 1.55                  |
| other                        | 1.54                  |
| *(blank)*                    | 1.45                  |
| glioblastoma                 | 1.03                  |
| multiple myeloma             | 1.02                  |
| breast cancer                | 0.94                  |
| melanoma                     | 0.86                  |


**Top 10 tissue labels (excluded)**


| Tissue label                              | % of excluded samples |
| ----------------------------------------- | --------------------- |
| CML                                       | 4.30                  |
| blood                                     | 4.25                  |
| K562 cell line                            | 3.28                  |
| bone marrow                               | 3.13                  |
| acute T cell leukemia                     | 2.08                  |
| Peripheral Blood Mononuclear Cells (PBMC) | 1.97                  |
| peripheral blood                          | 1.74                  |
| liver                                     | 1.46                  |
| skin                                      | 1.46                  |
| breast                                    | 1.36                  |


### 2.1.2 Top lung labels

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


## 2.2 Normalised disease breakdown

The table below shows the regex rules used to map free-text disease labels into broad categories. Rules are applied in order — the first match wins. Labels that match none of the rules are placed in **Other**.


| Category                     | Pattern                              |
| ---------------------------- | ------------------------------------ |
| IPF / Pulmonary Fibrosis     | `pulmonary fibrosis`                  |
| COVID-19 / SARS-CoV-2        | `COVID`                               |
| Lung Adenocarcinoma (LUAD)   | `lung adenocarcinoma`                 |
| NSCLC                        | `NSCLC`                               |
| Lung Squamous Cell Carcinoma | `squamous cell carcinoma of the lung` |
| COPD                         | `\bCOPD\b`                           |
| Lung Cancer (general)        | `lung cancer`                         |
| Cystic Fibrosis              | `cystic fibrosis`                    |
| Interstitial Lung Disease    | `interstitial lung`                   |


Within the 858 lung-intersection samples, free-text disease labels are normalised into broad categories via regex. Categories representing < 2% are collapsed into **Other**. IPF / Pulmonary Fibrosis and Lung Adenocarcinoma together account for roughly 40% of the set.

---

# 3. Cell count distribution

The 858 lung-intersection SRX samples contain ~6.1M cells in total (mean 7,126 / median 5,125 per SRX). The distribution is right-skewed, with most samples in the 3,000–9,000 cell range and a small number of outliers above 40,000.


# Appendix

## Top lung union labels

The 7,986 samples in the **lung union** (disease *or* tissue is lung-related). Because the union casts a wider net, many non-lung diseases appear, the disease column reflects whatever was annotated regardless of whether the tissue label is lung-specific.

**Top 10 tissue labels (union)**


| Tissue label                               | % of union samples |
| ------------------------------------------ | ------------------ |
| lung                                       | 6.30               |
| blood                                      | 3.69               |
| breast                                     | 2.55               |
| Peripheral Blood Mononuclear Cells (PBMC)  | 2.23               |
| liver                                      | 1.89               |
| Peripheral blood mononuclear cells (PBMCs) | 1.78               |
| peripheral blood mononuclear cells         | 1.68               |
| blood monocyte                             | 1.38               |
| esophageal carcinoma                       | 1.36               |
| brain                                      | 1.35               |


**Top 10 disease labels (union)**


| Disease label                       | % of union samples |
| ----------------------------------- | ------------------ |
| COVID-19                            | 5.71               |
| breast cancer                       | 2.09               |
| Crohn's disease                     | 1.78               |
| colorectal cancer                   | 1.48               |
| lung adenocarcinoma                 | 1.48               |
| high-grade serous ovarian cancer    | 1.35               |
| cardiovascular disease risk factors | 1.15               |
| hepatocellular carcinoma            | 1.13               |
| Hepatocellular carcinoma            | 1.05               |
| Alzheimer's Disease                 | 1.01               |