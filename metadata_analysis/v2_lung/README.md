# V2 Metadata Analysis

## Rationale for starting a new metadata analysis

On 2026-03-26, Parashar and Oliver decided that the previous approach to data collection was unfesable. We had overestimated the value of scBaseCount. We determined it lacks rich enough metadata for our previous goals.

## Steps forward

We will focus on lung data. We suspect it is a manageably small yet diverse dataset. Informed by the previous round of metadata analysis, we suspect there to be < 1,000,000 cells with a lung-like tissue label in scBaseCount.

## In this directory

```
metadata_analysis/v2_lung
├── images
├── README.md
├── analysis.ipynb
├── report.md
└── scbasecount_metadata -> ../../data/scbasecount/2026-01-12/metadata/GeneFull/Homo_sapiens
```