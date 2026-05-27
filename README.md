# Lactate, Logistic Regression, and ROC Analysis

This repository contains reproducible R code to analyze admission lactate (T1/DIA1), 24-hour lactate (T2/DIA2), in-hospital mortality, logistic regression models, and ROC curves.

## Important Note About the Data

The original dataset is not included in this repository because it contains sensitive or non-public information. To run the analysis locally, place the private dataset at:

```text
data/SPSSlactato.csv
```

Use `data/SPSSlactato_template.csv` as a reference for the expected data structure.

## Repository Structure

```text
.
├── R/
│   └── analysis_lactate.R
├── data/
│   └── SPSSlactato_template.csv
├── outputs/
│   └── .gitkeep
├── README.md
└── .gitignore
```

## Required Variables

The script expects the following columns:

| Variable | Description | Example |
|---|---|---|
| `Edad` | Age in years | `65` |
| `SEXO` | Biological sex as recorded in the dataset | `MASCULINO` or `FEMENINO` |
| `MORTALIDAD` | In-hospital mortality | `SI` or `NO` |
| `DIA1` | Admission lactate / T1, mmol/L | `2,4` or `2.4` |
| `DIA2` | 24-hour lactate / T2, mmol/L | `3,2` or `3.2` |

Additional columns may be present in the CSV file; they will be ignored by the analysis script.

## How to Run

From RStudio or from a terminal in the repository root:

```r
source("R/analysis_lactate.R")
```

You can also run:

```bash
Rscript R/analysis_lactate.R
```

To use a custom input or output path:

```bash
Rscript R/analysis_lactate.R path/to/data.csv path/to/outputs
```

## R Packages

The script automatically installs missing packages:

- `dplyr`
- `broom`
- `pROC`
- `ggplot2`

## Generated Outputs

The following files are saved in `outputs/`:

- `Logistic_regression_results.csv`
- `Logistic_regression_table_latex.txt`
- `ROC_summary.csv`
- `DIA2_cutoff_3_confusion_matrix.csv`
- `ROC_lactate_T1_T2_publication.pdf`
- `ROC_lactate_T1_T2_publication.tiff`
- `ROC_lactate_T1_T2_publication.png`

## Statistical Models

**Model 1**

```text
Mortality ~ DIA1 + DIA2
```

**Model 2**

```text
Mortality ~ Age + male sex + DIA2 > 3 mmol/L
```

The outcome is coded as:

```text
1 = in-hospital death
0 = survival
```

## Reproducibility Note

Before pushing changes to GitHub, verify that `data/SPSSlactato.csv` is not included in the commit. The `.gitignore` file is configured to exclude real datasets inside `data/` while keeping the template file.

## Data Availability

The original data are not publicly available due to CEISH confidentiality restrictions. Data may be requested from the corresponding authors, subject to applicable ethical and institutional approvals.
