# lifelines-personality-stress-correlations

Analysis scripts and reproducibility materials for the EJP manuscript *Familial Correlation between Personality and Stress Exposure: Insights from the Multigenerational Lifelines Cohort Study* (EJP-26-4766).

## Contents

- [Overview](#overview)
- [Repository tree](#repository-tree)
- [How to run (script by script)](#how-to-run-script-by-script)
- [Upload to GitHub](#upload-to-github)
- [ACE models with ASReml-R](#ace-models-with-asreml-r)
- [Software](#software)
- [Data availability](#data-availability)
- [Open Science](#open-science)
- [Citation](#citation)

## Overview

**Background:** Personality traits and stress exposures co-occur in the population, but the relative roles of shared genetic and environmental factors in their overlap are not fully understood in large, multigenerational samples. We analysed eight phenotypes (five stress-related measures and three personality domains) in the Lifelines cohort using pedigree-based ACE models in ASReml-R.

**Results:** Univariate ACE models quantified the relative contributions of additive genetic (A), shared household environment (C), and unique environment (E) variance for each phenotype. Bivariate ACE models estimated genetic (r_g), shared-environmental (r_c), unique-environmental (r_e), and phenotypic (r_P) correlations between all stress–stress and personality–stress pairings, with per-component 1-df likelihood-ratio tests and Bonferroni correction (8 × 2 = 16 univariate tests; 25 × 3 = 75 bivariate tests; Ani et al., 2026 protocol).

**Conclusions:** Individual-level Lifelines data cannot be shared publicly. This repository provides the analysis code for approved researchers working in the Lifelines secure HPC environment.

## Repository tree

Folders are numbered so they appear in analysis order on GitHub.

```
lifelines-personality-stress-correlations/
│
├── README.md
├── LICENSE
├── CITATION.cff
├── OPEN_SCIENCE_STATEMENT.md
├── .gitignore
├── lifelines-personality-stress-correlations.Rproj
│
├── 01_Config/
│   ├── 01_setup.R
│   └── 01_helpers.R
│
├── 02_Data_Preparation/
│   ├── 02a_score_construction.R
│   └── 02_data_preparation.R
│
├── 03_Descriptives/
│   ├── 03_descriptives.R
│   ├── 03b_anova_sex_diff.R
│   ├── 03c_age_brackets.R
│   └── 03d_cross_wave_correlations.R
│
├── 04_Main_Analysis/
│   ├── 04_univariate_heritability.R
│   ├── 04b_bivariate_stress_stress.R
│   ├── 04c_bivariate_personality_stress.R
│   └── 04d_significance_audit.R
│
├── 05_Sensitivity/
│   ├── 05_sensitivity_int.R
│   └── 05b_perwave_univariate.R
│
├── 06_CT_Liability/
│   ├── 06_rrr_liability_heritability.R
│   └── 06b_ct_subscale_heritability.R
│
├── 07_Figures/
│   ├── 07_figures_main.R
│   └── 07b_figures_supplementary.R
│
├── docs/
│   └── Lifelines_Logbook.Rmd
│
├── data/
│   └── README.md          ← expected Lifelines files (not in repo)
│
└── output/                ← tables + figures (gitignored, created when you run scripts)
    ├── tables/
    └── figures/
```

Each numbered folder contains a short `README.md` listing its scripts.

## How to run (script by script)

**Requirements:** R ≥ 4.3; [ASReml-R 4.2](https://vsni.co.uk/software/asreml-r/) for scripts in folders 04–05; CRAN packages listed in `01_Config/01_setup.R`.

Open the repository root in RStudio (`lifelines-personality-stress-correlations.Rproj`) or run `setwd()` to the folder that contains this `README.md`.

### Step 0 — setup (always first)

```r
setwd("/path/to/lifelines-personality-stress-correlations")
source("01_Config/01_setup.R")
```

### Step 1 — data preparation

```r
source("02_Data_Preparation/02a_score_construction.R")
source("02_Data_Preparation/02_data_preparation.R")
```

### Step 2 — descriptives

```r
source("03_Descriptives/03_descriptives.R")
source("03_Descriptives/03b_anova_sex_diff.R")
source("03_Descriptives/03c_age_brackets.R")
source("03_Descriptives/03d_cross_wave_correlations.R")
```

### Step 3 — main analysis (ASReml-R)

```r
source("04_Main_Analysis/04_univariate_heritability.R")
source("04_Main_Analysis/04b_bivariate_stress_stress.R")
source("04_Main_Analysis/04c_bivariate_personality_stress.R")
source("04_Main_Analysis/04d_significance_audit.R")
```

### Step 4 — sensitivity

```r
source("04_Main_Analysis/04_univariate_heritability.R")   # defines fit_ace()
source("04_Main_Analysis/04b_bivariate_stress_stress.R")  # defines fit_biv()
source("05_Sensitivity/05_sensitivity_int.R")
source("05_Sensitivity/05b_perwave_univariate.R")
```

### Step 5 — CT liability (supplementary)

```r
source("06_CT_Liability/06_rrr_liability_heritability.R")
source("06_CT_Liability/06b_ct_subscale_heritability.R")
```

### Step 6 — figures

```r
source("07_Figures/07_figures_main.R")
source("07_Figures/07b_figures_supplementary.R")
```

**HPC note:** one ASReml licence seat — do not run parallel ASReml jobs.

## Upload to GitHub

Create a **new private** repository named **`lifelines-personality-stress-correlations`**, then upload from:

`C:\Users\ReicheltF\OneDrive - UMCG\1_lifelines_new\scripts\github_short_2026-07-15`

### Option A — upload via browser (simplest)

1. On GitHub: **New repository** → name `lifelines-personality-stress-correlations` → Private → **Create repository** (no README — you have one locally).
2. Click **Add file** → **Upload files**.
3. Drag the **entire contents** of `github_short_2026-07-15`:
   - Root: `README.md`, `LICENSE`, `CITATION.cff`, `.gitignore`, `lifelines-personality-stress-correlations.Rproj`, …
   - Folders: `01_Config` … `07_Figures`, `docs`, `data`
4. Commit message: `Initial upload — numbered analysis folders`
5. Refresh — you should see folders `01_Config` … `07_Figures`, not loose `.R` files at the root.

**Do not upload:** `.Rdata`, `.csv`, or `output/`; only `data/README.md` in `data/`.

### Option B — git from your PC

```powershell
cd "C:\Users\ReicheltF\OneDrive - UMCG\1_lifelines_new\scripts\github_short_2026-07-15"
git init
git add .
git commit -m "Reorganise into numbered analysis folders"
git branch -M main
git remote add origin https://github.com/FelixReichelt13/lifelines-personality-stress-correlations.git
git push -u origin main --force
```

Use `--force` only if you are replacing the first upload entirely and no one else has cloned the repo.

### After upload — quick check

- [ ] Root shows `README.md` + folders `01_Config` … `07_Figures`, not 20 loose `.R` files
- [ ] `data/` contains only `README.md`
- [ ] No `.Rdata` or `.csv` files visible in the repo

## ACE models with ASReml-R

### Univariate (script 04)

```
y = mu + b1*sex + b2*age + b3*age^2 + a + c + e
```

with `a ~ N(0, A * sigma_a^2)` (additive genetic), `c ~ N(0, I * sigma_c^2)` (shared household) and `e ~ N(0, I * sigma_e^2)`. Heritability, shared environment and unique environment:

```r
heritability <- vpredict(model, h2 ~ V2 / (V1 + V2 + V3))
shared_env   <- vpredict(model, c2 ~ V1 / (V1 + V2 + V3))
unique_env   <- vpredict(model, e2 ~ V3 / (V1 + V2 + V3))
```

### Model comparison (scripts 04 and 05)

| Model | Random terms | Tests via LRT |
|-------|--------------|---------------|
| ACE | `vm(ID, Ainv) + idv(HouseID1)` | reference |
| AE | `vm(ID, Ainv)` | ACE vs AE → c² = 0 |
| CE | `idv(HouseID1)` | ACE vs CE → h² = 0 |

Bonferroni: 8 phenotypes × 2 tests = **16** (α = .05/16).

### Bivariate (scripts 04b, 04c, 05)

Each pair: `rp` (estimate + SE + CI, no LRT); `rg`, `rc`, `re` each with a separate 1-df LRT. Bonferroni: 25 pairings × 3 = **75** tests (α = .05/75).

## Software

| Purpose | Package |
|---------|---------|
| Variance components (ACE) | **`asreml-R`** *(commercial)* |
| Latent factors for LTE / LDI | `lavaan` |
| RRR / liability *h²* | `marginaleffects`, `sandwich` |
| Data handling | `data.table`, `dplyr`, `tidyr` |
| Figures | `ggplot2`, `patchwork`, `scales` |

## Data availability

Individual-level Lifelines data cannot be shared because of participant consent. Bona-fide researchers can apply via the [Lifelines Research Office](https://www.lifelines.nl/researcher). See [`data/README.md`](data/README.md) for the file list expected on the HPC.

## Open Science

See [OPEN_SCIENCE_STATEMENT.md](OPEN_SCIENCE_STATEMENT.md). Zenodo DOI will be added in `CITATION.cff` at acceptance.

## Citation

> Reichelt, F., [co-authors]. *Familial Correlation between Personality and Stress Exposure: Insights from the Multigenerational Lifelines Cohort Study.* European Journal of Personality (in revision).

> Reichelt, F. (2026). *lifelines-personality-stress-correlations*. GitHub: https://github.com/FelixReichelt13/lifelines-personality-stress-correlations

**Contact:** Felix Reichelt — f.reichelt@umcg.nl
