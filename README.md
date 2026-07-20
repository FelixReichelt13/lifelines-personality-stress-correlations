# lifelines-personality-stress-correlations

Analysis code for EJP-26-4766 — *Familial Correlation between Personality and Stress Exposure: Insights from the Multigenerational Lifelines Cohort Study* (Lifelines cohort; ASReml-R ACE models).

Individual-level data are **not** included. See [Data availability](#data-availability) and [`data/README.md`](data/README.md).

## Contents

- [Reviewer quick path](#reviewer-quick-path)
- [Primary vs exploratory analyses](#primary-vs-exploratory-analyses)
- [Phenotypes](#phenotypes)
- [CT coding](#ct-coding)
- [ACE models with ASReml-R](#ace-models-with-asreml-r)
- [Bivariate pairings](#bivariate-pairings)
- [Transformations and scoring](#transformations-and-scoring)
- [Software](#software)
- [Manuscript output map](#manuscript-output-map)
- [How to run](#how-to-run)
- [Data availability](#data-availability)
- [Open Science](#open-science)
- [Affiliation](#affiliation)
- [Citation](#citation)

## Reviewer quick path

You do **not** need an ASReml licence to audit most of the pipeline; several steps need no Lifelines data at all:

1. **[`docs/Lifelines_Logbook.Rmd`](docs/Lifelines_Logbook.Rmd)** — item-level scoring, exclusions, and analysis narrative (knit to HTML locally, or read the `.Rmd` in GitHub).
2. **Scripts 4, 14–16** — sex-difference ANOVA from published means, CT subscale table, and all figures; numbers are embedded or assembled without individual-level data.
3. **Scripts 3, 5–6** — descriptives and cross-wave diagnostics; CRAN only (require prepared Lifelines extracts, no ASReml).
4. **Script 10** (`04d_significance_audit.R`) — Bonferroni / FDR audit from `*_INTERNAL.tsv` files produced by scripts 7–9; inspectable once those tables exist from a prior run.

Full reproduction requires Lifelines HPC access, the files listed in [`data/README.md`](data/README.md), and ASReml-R for scripts 7–12.

## Primary vs exploratory analyses

**Primary (confirmatory):** univariate and bivariate ACE on the **raw** scale for all eight phenotypes; 25 bivariate pairings with Bonferroni correction (/16 univariate, /75 bivariate). CT enters as continuous `Childhoodtrauma`.

**Exploratory** (reported separately; not part of the primary multiplicity correction):

- Inverse-normal sensitivity — script 11; Suppl. Tables S10–S10e
- Per-wave LTE/LDI univariate ACE — script 12; Suppl. Table S9
- Age-bracket descriptives — script 5; Suppl. Table S7d; Suppl. Figs S3–S16
- CTQ-SF subscale heritability — script 14; Suppl. Table S12b
- Liability-scale CT (*h²* from RRR) — folder `06_CT_Liability`; Suppl. Tables S11–S12

Details: [OPEN_SCIENCE_STATEMENT.md](OPEN_SCIENCE_STATEMENT.md) §1 (Preregistration).

## Phenotypes

Eight phenotypes enter the primary ACE models (scripts 7–9). Analytic sample: adults ≥ 18 years (N ≈ 141,751) from the Lifelines extraction with three-generation pedigree data. Main exclusions: participants < 18 years; pedigree members with unknown or inconsistent parental IDs; LTE/LDI sum scores with any missing item within a wave (complete-case scoring). Further sampling details are in [OPEN_SCIENCE_STATEMENT.md](OPEN_SCIENCE_STATEMENT.md).

| Phenotype | Instrument | Scale in primary ACE |
|-----------|------------|----------------------|
| `Childhoodtrauma` | CTQ-SF total (25 clinical items) | Continuous (raw) |
| `Loneliness` | De Jong Gierveld Loneliness Scale (6 items) | Continuous (raw) |
| `Lack_of_social_support` | SPF-IL social well-being (9 items, reversed) | Continuous (raw) |
| `Stressful_life_events` | LTE latent factor across waves (`lavaan`, FIML) | Continuous (SD = 1 after prep) |
| `Life_difficulties` | LDI latent factor across waves (`lavaan`, FIML) | Continuous (SD = 1 after prep) |
| `Neuroticism` | NEO-PI-R six-facet domain sum | z-scored (Lifelines extract) |
| `Extraversion` | NEO-FFI domain SUM3 | z-scored (Lifelines extract) |
| `Conscientiousness` | NEO-FFI domain SUM3 | z-scored (Lifelines extract) |

Inverse-normal-transformed stress phenotypes are used in the sensitivity analysis (script 11) only.

## CT coding

Two childhood-trauma definitions are used — do not conflate them:

- **Primary ACE models (scripts 7–9):** continuous `Childhoodtrauma` (CTQ-SF total score).
- **Liability analyses (folder `06_CT_Liability`):** dichotomous CT when **any** CTQ-SF subscale meets its moderate/severe threshold (Bernstein et al., 2003): emotional abuse ≥ 12, physical abuse ≥ 9, sexual abuse ≥ 7, emotional neglect ≥ 14, physical neglect ≥ 9. Coded 0 only when all five subscales are observed and none meet threshold (see `02_data_preparation.R`).

## ACE models with ASReml-R

### Univariate (script 07)

```
y = mu + b1*sex + b2*age + b3*age^2 + a + c + e
```

with `a ~ N(0, A * sigma_a^2)` (additive genetic), `c ~ N(0, I * sigma_c^2)` (shared household, one level per residential HouseID) and `e ~ N(0, I * sigma_e^2)`. *A* is the additive relationship matrix from the three-generation Lifelines pedigree via `asreml::ainverse()`. Heritability, shared environment and unique environment are extracted with `vpredict()` (ASReml-R Reference Manual 4.2, §4.6):

```r
heritability <- vpredict(model, h2 ~ V2 / (V1 + V2 + V3))
shared_env   <- vpredict(model, c2 ~ V1 / (V1 + V2 + V3))
unique_env   <- vpredict(model, e2 ~ V3 / (V1 + V2 + V3))
```

### Model comparison (scripts 07 and 11)

Scripts 07 and 11 fit **three nested submodels** per phenotype (ACE, AE, CE) for a full decomposition table with boundary-corrected LRTs and information criteria:

| Model | Random terms | Tests via LRT | df (boundary) |
|-------|--------------|---------------|---------------|
| ACE | `vm(ID, Ainv) + idv(HouseID1)` | (reference model) | – |
| AE | `vm(ID, Ainv)` | ACE vs AE → c² = 0 | 1 (0.5·χ²₀ + 0.5·χ²₁) |
| CE | `idv(HouseID1)` | ACE vs CE → h² = 0 | 1 (0.5·χ²₀ + 0.5·χ²₁) |

e² is not tested univariately (always present / on the boundary). Bonferroni correction spans 8 phenotypes × 2 tests = **16** (α = .05/16 = .003125).

Each fit contributes to Suppl. Tables **S8b** (LRTs) and **S8c** (logLik / AIC / BIC). Raw ASReml `varcomp` and fixed effects go to **S8d** / **S8e**.

### Bivariate (scripts 08, 09, 11)

Two-trait REML uses wide `cbind(P1, P2)` with `random = ~ us(trait):vm(ID, Ainv) + us(trait):HouseID1` and `rcov = ~ us(trait):units`. ASReml orders variance components alphabetically by random-term name — always confirm V-indexing with:

```r
cbind(summary(model)$varcomp, n = seq_along(model$vparameters))
```

Each pair yields `rp` (estimate + SE + 95 % CI; **no LRT**), plus `rg`, `rc`, and `re` each tested with a **separate 1-df LRT** (full model vs sub-model with that cross-trait covariance fixed to zero). Bonferroni correction: 25 pairings × 3 correlations = **75** tests (α = .05/75 ≈ .000667).

Tables **S13/S14** hold estimates and LRT *p* values; **S13c/S14c** hold model comparisons; **S15/S16** hold A/C/E decomposition. `*_INTERNAL.tsv` files (Wald *p*, Bonferroni flags) are for audit only. Script **10** (`04d_significance_audit.R`) summarises those INTERNAL tables independently (Bonferroni / FDR counts and decomposition warnings).
## Bivariate pairings

All pairings are fit on the raw scale (script **8** or **9**). Each pairing contributes three Bonferroni-tested correlations (*r*<sub>g</sub>, *r*<sub>c</sub>, *r*<sub>e</sub>); phenotypic *r*<sub>P</sub> is estimate + 95 % CI only.

**Stress × stress (10 pairs, Table S13 / Figure 2):**

| | |
|---|---|
| `Childhoodtrauma` × `Loneliness` | `Childhoodtrauma` × `Lack_of_social_support` |
| `Childhoodtrauma` × `Stressful_life_events` | `Childhoodtrauma` × `Life_difficulties` |
| `Loneliness` × `Lack_of_social_support` | `Loneliness` × `Stressful_life_events` |
| `Loneliness` × `Life_difficulties` | `Lack_of_social_support` × `Stressful_life_events` |
| `Lack_of_social_support` × `Life_difficulties` | `Stressful_life_events` × `Life_difficulties` |

**Personality × stress (15 pairs, Table S14 / Figure 3):** every combination of `Neuroticism`, `Extraversion`, or `Conscientiousness` with each of the five stress phenotypes above (3 × 5).

## Transformations and scoring

| Transformation | Script | Notes |
|----------------|--------|-------|
| Rank-based inverse normal (stress phenotypes) | `05_sensitivity_int.R` | Sensitivity only; primary ACE uses raw scale |
| SPF-IL sign flip → `Lack_of_social_support` | `02a_score_construction.R` | Higher = less support |
| LTE / LDI latent factors (FIML) | `02_data_preparation.R` | Cross-wave one-factor SEM; factor scores scaled to SD = 1 |
| CT dichotomisation (liability) | `02_data_preparation.R` | Subscale Bernstein thresholds; continuous total used in ACE |
| Personality z-scores | Lifelines NEO release | Consumed as-is in scripts 7–9 |

Non-standard scoring conventions (CTQ mean imputation when ≤ 2 items missing per subscale; LTE/LDI complete-case per wave; NEO `_SUM3` rescaling) are documented in [`docs/Lifelines_Logbook.Rmd`](docs/Lifelines_Logbook.Rmd) and [OPEN_SCIENCE_STATEMENT.md](OPEN_SCIENCE_STATEMENT.md).

## Software

| Purpose | Package |
|---------|---------|
| Variance components (ACE) | **`asreml-R`** *(commercial)* |
| Latent factors for LTE / LDI | `lavaan` |
| RRR / liability *h²* | `marginaleffects`, `sandwich` |
| Data handling | `data.table`, `dplyr`, `tidyr` |
| Inverse normal transformation | base R (`qnorm`, `rank`) |
| Figures | `ggplot2`, `patchwork`, `scales` |

`asreml-R` is required for scripts 7–12; obtain an academic licence from [vsni.co.uk/software/asreml-r](https://vsni.co.uk/software/asreml-r/).

## Manuscript output map

| Manuscript output | Script(s) | Key repo outputs |
|-------------------|-----------|------------------|
| **Table 1** | 3 | `output/tables/table_1_descriptives_overall.tsv` |
| **Figure 1** | 7, 15 | `table_S8_*`; `output/figures/figure_1_univariate_ace.*` |
| **Figure 2** | 8, 15 | `table_S13_*`, `table_S15_*`; `figure_2_stress_stress.*` |
| **Figure 3** | 9, 15 | `table_S14_*`, `table_S16_*`; `figure_3_personality_stress.*` |
| **S5–S6** | 6 | `table_S5_lte_ldi_cross_wave_correlations.tsv`, `table_S6_lte_ldi_factor_loadings.tsv` |
| **S7a–S7b** | 3 | `table_S7a_descriptives_by_sex.tsv`, `table_S7b_sex_differences_age_adjusted.tsv` |
| **S7c** | 4 | `table_S7c_sex_differences_anova.tsv` |
| **S7d** | 5 | `table_S7d_age_bracket_descriptives.tsv` |
| **S7e** | 1 | `table_S7e_internal_consistency.tsv` |
| **S8–S8e** | 7 | `table_S8_univariate_ace_full.tsv`, `table_S8b`–`S8e` |
| **S9** | 12 | `table_S9_perwave_lte_ldi_variance.tsv` |
| **S10–S10e** | 11 | `table_S10_int_*` (inverse-normal sensitivity) |
| **S11–S12** | 13 | `table_S11_ct_dichotomization_prevalence.tsv`, `table_S12_ct_rrr_liability_h2.tsv` |
| **S12b** | 14 | `table_S12b_ct_subscale_heritability.tsv` |
| **S13–S15** | 8 | `table_S13_bivariate_stress_stress.tsv`, `S13c`, `S15` |
| **S14–S16** | 9 | `table_S14_bivariate_personality_stress.tsv`, `S14c`, `S16` |
| **Figs S3–S21** | 16 | `output/figures/` (see header of `07b_figures_supplementary.R`) |
| **Audit** | 10 | `table_audit_*.tsv`, `audit_summary.txt` (independent post-hoc summary) |
| **Logbook** | — | [`docs/Lifelines_Logbook.Rmd`](docs/Lifelines_Logbook.Rmd) |

## How to run

Run from the repository root. Always start with `01_Config/01_setup.R`.

```r
setwd("/path/to/lifelines-personality-stress-correlations")
source("01_Config/01_setup.R")

source("02_Data_Preparation/02a_score_construction.R")
source("02_Data_Preparation/02_data_preparation.R")
# … continue in table order below
```

Script **11** (`05_sensitivity_int.R`) needs `fit_ace()` and `fit_biv()` from scripts **7** and **8** in the same R session. Set the working directory to the folder that contains this `README.md` before sourcing.

Generated tables and figures go to `output/tables/` and `output/figures/` (gitignored).

**Requirements:** R ≥ 4.3; ASReml-R 4.2 for scripts 7–12; CRAN packages in `01_Config/01_setup.R`. One ASReml licence seat on HPC — no parallel ASReml jobs.

### Compute expectations

| Stage | Scripts | ASReml? | Typical wall time (serial, Lifelines HPC)* |
|-------|---------|---------|---------------------------------------------|
| Item scoring + merge | 1–2 | No | Minutes–1 h (depends on raw CSVs) |
| Descriptives + SEM diagnostics | 3–6 | No | Minutes |
| Univariate ACE (8 phenotypes × ACE/AE/CE) | 7 | Yes | Hours |
| Bivariate ACE (25 pairs × nested LRTs) | 8–9 | Yes | Dominant cost; hours–days |
| Significance audit | 10 | No | Seconds |
| INT sensitivity | 11 | Yes | Hours (subset of 7–8) |
| Per-wave LTE/LDI | 12 | Yes | ~1 h |
| CT liability RRR | 13 | No | Minutes |
| Figures | 15–16 | No | Minutes |

\*Order-of-magnitude from a full production run (~141k adults, three-generation pedigree). Do not run multiple ASReml jobs in parallel on a single licence.

### Script index

| # | Script | What it does | Key output |
|---|--------|--------------|-------------------|
| — | [`01_Config/01_setup.R`](01_Config/01_setup.R) | Paths, packages, phenotypes, Bonferroni thresholds (/16 uni, /75 bivar), plotting palette | — |
| 1 | [`01_Config/01_helpers.R`](01_Config/01_helpers.R) | Shared utilities (`inverse_normal`, `save_table`, `save_figure`, `ct_dichotomous_moderate`) | — |
| 2 | [`02_Data_Preparation/02_data_preparation.R`](02_Data_Preparation/02_data_preparation.R) | Merge pedigree, household, phenotypes, personality; adult sample; LTE/LDI latent factors (`lavaan`); dichotomous CT for liability (subscale Bernstein thresholds) → `data/data_prepared.Rdata` | — |
| 3 | [`03_Descriptives/03_descriptives.R`](03_Descriptives/03_descriptives.R) | Table 1; sex-stratified descriptives | Table 1; S7a–S7b |
| 4 | [`03_Descriptives/03b_anova_sex_diff.R`](03_Descriptives/03b_anova_sex_diff.R) | Sex-difference ANOVA from published means (no Lifelines data) | S7c |
| 5 | [`03_Descriptives/03c_age_brackets.R`](03_Descriptives/03c_age_brackets.R) | Phenotype descriptives by age bracket | S7d |
| 6 | [`03_Descriptives/03d_cross_wave_correlations.R`](03_Descriptives/03d_cross_wave_correlations.R) | LTE/LDI cross-wave correlations and lavaan loadings | S5–S6 |
| 7 | [`04_Main_Analysis/04_univariate_heritability.R`](04_Main_Analysis/04_univariate_heritability.R) | Univariate ACE + AE + CE; boundary 1-df LRTs for h² and c²; Bonferroni /16; 95 % CIs | Figure 1; S8, S8b–S8e |
| 8 | [`04_Main_Analysis/04b_bivariate_stress_stress.R`](04_Main_Analysis/04b_bivariate_stress_stress.R) | Bivariate ACE, 10 stress × stress pairs; per-component 1-df LRTs (r_g, r_c, r_e); Bonferroni /75 | Figure 2; S13, S13c, S15 |
| 9 | [`04_Main_Analysis/04c_bivariate_personality_stress.R`](04_Main_Analysis/04c_bivariate_personality_stress.R) | Bivariate ACE, 15 personality × stress pairs; same LRT structure as 04b | Figure 3; S14, S14c, S16 |
| 10 | [`04_Main_Analysis/04d_significance_audit.R`](04_Main_Analysis/04d_significance_audit.R) | Independent Bonferroni / FDR audit and `audit_summary.txt` from INTERNAL tables | `table_audit_*`, `audit_summary.txt` |
| 11 | [`05_Sensitivity/05_sensitivity_int.R`](05_Sensitivity/05_sensitivity_int.R) | Inverse-normal sensitivity (univariate + stress × stress bivariate); sources `fit_ace` / `fit_biv` from 04 + 04b | S10–S10e |
| 12 | [`05_Sensitivity/05b_perwave_univariate.R`](05_Sensitivity/05b_perwave_univariate.R) | Univariate ACE per LTE/LDI wave | S9 (per-wave) |
| 13 | [`06_CT_Liability/06_rrr_liability_heritability.R`](06_CT_Liability/06_rrr_liability_heritability.R) | First-degree-relative RRR → Wray–Gottesman liability h² for dichotomous CT | S11–S12 |
| 14 | [`06_CT_Liability/06b_ct_subscale_heritability.R`](06_CT_Liability/06b_ct_subscale_heritability.R) | Assemble CT subscale heritability table (no Lifelines data) | S12b |
| 15 | [`07_Figures/07_figures_main.R`](07_Figures/07_figures_main.R) | Main Figures 1–3 (embedded summary numbers) | Figures 1–3 |
| 16 | [`07_Figures/07b_figures_supplementary.R`](07_Figures/07b_figures_supplementary.R) | Supplementary Figures S3–S21 | S3–S21 |
| — | [`docs/Lifelines_Logbook.Rmd`](docs/Lifelines_Logbook.Rmd) | Narrative walk-through for reviewers | HTML logbook |

## Data availability

Individual-level Lifelines data cannot be shared because of the participant consent. Bona-fide researchers can apply for access to the same extraction via the Lifelines Research Office (<https://www.lifelines.nl/researcher>); the analytic dataset can then be opened inside the Lifelines secure HPC environment using the scripts in this repository. See [`data/README.md`](data/README.md) for the list of files the scripts expect.

## Open Science

EJP Open Materials compliance: [OPEN_SCIENCE_STATEMENT.md](OPEN_SCIENCE_STATEMENT.md). Repository archived on Zenodo (DOI in `CITATION.cff`) at acceptance.

## Affiliation

This study uses data from the [Lifelines Cohort Study](https://www.lifelines.nl), a multidisciplinary research project on health and health-related behaviours in a large representative sample of the population of the north of the Netherlands.

This work was carried out in the context of the [Stress in Action](https://stress-in-action.nl) consortium (NWO Gravitation programme), which studies daily-life stress and its impact on mental and cardiometabolic health.

## Citation

> Reichelt, F., [co-authors]. *Familial Correlation between Personality and Stress Exposure: Insights from the Multigenerational Lifelines Cohort Study.* European Journal of Personality (in revision).

> Reichelt, F. (2026). *lifelines-personality-stress-correlations*. GitHub: https://github.com/FelixReichelt13/lifelines-personality-stress-correlations

See [CITATION.cff](CITATION.cff) for machine-readable metadata.

**Contact:** Felix Reichelt — University Medical Center Groningen — f.reichelt@umcg.nl — ORCID: *(to be added)*
