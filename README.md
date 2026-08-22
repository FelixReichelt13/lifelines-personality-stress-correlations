# Familial Correlation between Personality and Stress Exposure
### Insights from the Multigenerational Lifelines Cohort Study

**Authors:** Reichelt, F., Wang, R., Triatin, R. D., Bos, M., Kavousi, M., de Geus, E. J. C., Snieder, H.  
**Journal:** *European Journal of Personality* (manuscript EJP-26-4766, in revision)  
**Institution:** University Medical Center Groningen · [Stress in Action](https://stress-in-action.nl) · [Lifelines](https://www.lifelines.nl)

### About this study

People differ widely in how frequently and powerfully they report stressful circumstances. In this study we tested to what extent these differences are influenced by two factors that family members partly share: the family environment and genetic resemblance. We also tested how far personality mediates the link between family background and self-reported stress.

We used data from more than 140,000 adults in the Lifelines Cohort Study in the Netherlands. The stress measures were childhood trauma, stressful life events, ongoing life difficulties, loneliness, and lack of social support. The personality traits were neuroticism, extraversion, and conscientiousness.

Individuals with higher extraversion and conscientiousness typically reported less stress, especially for loneliness and lack of social support. Individuals with higher neuroticism reported more stress across all domains. With pedigree-based models we estimated genetic, shared household, and person-specific environmental contributions to stress exposure, personality, and their association. Person-specific environmental influences and the shared household environment explained a large part of individual differences in self-reported stress. Self-reported stress exposures were also partly heritable, with the largest genetic influence on retrospectively reported childhood trauma. There was clear overlap in the genetic factors influencing personality and stress exposures. This fits with personality shaping the environments people evoke or select, and with personality influencing appraisal, memory, and disclosure in self-report.

Self-reported stress exposure does not reflect external events alone. Environments, perception, and recall are also linked to personality through overlapping genetic and shared environmental factors.

*Keywords: stress exposure, personality traits, genetic correlation, heritability, stressful life events*

This repository holds the R code for the analyses, tables, and figures in the manuscript.

**Start here:** [`docs/Lifelines_Logbook.Rmd`](docs/Lifelines_Logbook.Rmd)  
**License:** MIT. **Data:** no individual-level Lifelines data included ([details](#data-availability))

Analysis scripts are in numbered folders (`02_Data_Preparation` through `07_Figures`). Flat copies for browsing are in [`00_All_Scripts/`](00_All_Scripts/). Open-science documentation: [`OPEN_SCIENCE_STATEMENT.md`](OPEN_SCIENCE_STATEMENT.md).

## Contents

- [About this study](#about-this-study)
- [Phenotypes](#phenotypes)
- [ACE models with ASReml-R](#ace-models-with-asreml-r)
- [Bivariate pairings](#bivariate-pairings)
- [Transformations and scoring](#transformations-and-scoring)
- [Software](#software)
- [Manuscript output map](#manuscript-output-map)
- [Data availability](#data-availability)
- [Open Science](#open-science)
- [Affiliation](#affiliation)
- [Citation](#citation)

## Phenotypes

Eight phenotypes were analysed in adults aged ≥ 18 years (N ≈ 141,751) from the Lifelines extraction with three-generation pedigree data. Further sampling details are in [OPEN_SCIENCE_STATEMENT.md](OPEN_SCIENCE_STATEMENT.md).

| Phenotype | Instrument | Scale in primary analysis |
|-----------|------------|---------------------------|
| `Childhoodtrauma` | CTQ-SF total (25 clinical items) | Continuous (raw) |
| `Loneliness` | De Jong Gierveld Loneliness Scale (6 items) | Continuous (raw) |
| `Lack_of_social_support` | SPF-IL social well-being (9 items, reversed) | Continuous (raw) |
| `Stressful_life_events` | LTE latent factor across waves | Continuous (SD = 1 after prep) |
| `Life_difficulties` | LDI latent factor across waves | Continuous (SD = 1 after prep) |
| `Neuroticism` | NEO-PI-R six-facet domain sum | z-scored (Lifelines extract) |
| `Extraversion` | NEO-FFI domain SUM3 | z-scored (Lifelines extract) |
| `Conscientiousness` | NEO-FFI domain SUM3 | z-scored (Lifelines extract) |

Item-level scoring, exclusions, and the cross-wave latent-factor step for LTE and LDI are documented in [`docs/Lifelines_Logbook.Rmd`](docs/Lifelines_Logbook.Rmd).

## ACE models with ASReml-R

Univariate and bivariate ACE models were fit with ASReml-R (scripts 7–9; inverse-normal sensitivity in script 11). All models adjust for sex and a quadratic age trend. Additive genetic effects use the three-generation pedigree; shared environment is modelled at the residential household level. Full model specification and variance-component definitions are in the script headers and the logbook.

### Univariate (script 7)

For each phenotype, scripts 7 and 11 fit **three nested submodels** (ACE, AE, CE) with boundary-corrected likelihood-ratio tests and information criteria:

| Model | Random terms | Tests via LRT | df (boundary) |
|-------|--------------|---------------|---------------|
| ACE | `vm(ID, Ainv) + idv(HouseID1)` | (reference model) | – |
| AE | `vm(ID, Ainv)` | ACE vs AE → c² = 0 | 1 (0.5·χ²₀ + 0.5·χ²₁) |
| CE | `idv(HouseID1)` | ACE vs CE → h² = 0 | 1 (0.5·χ²₀ + 0.5·χ²₁) |

*e²* is not tested univariately (always present / on the boundary). Bonferroni correction spans 8 phenotypes × 2 tests = **16** (α = .05/16 = .003125).

Wide outputs **S8b** / **S8c** (from script **04**) feed script **04a**, which writes manuscript **Table S9** (long ACE model comparison, one row per nested model). After the inverse-normal sensitivity run (**05**), **04a** reads **S11b** / **S11c** and writes **Table S12** (INT long model comparison). Raw ASReml `varcomp` and fixed effects go to **S8d** / **S8e**.

### Bivariate (scripts 8, 9, 11)

Two-trait REML is used for all 25 pairings. Each pair yields phenotypic correlation (*r*<sub>P</sub>: estimate + SE + 95 % CI; **no LRT**), plus *r*<sub>g</sub>, *r*<sub>c</sub>, and *r*<sub>e</sub> each tested with a **separate 1-df LRT** (full model vs sub-model with that cross-trait covariance fixed to zero). Bonferroni correction: 25 pairings × 3 correlations = **75** tests (α = .05/75 ≈ .000667).

Manuscript tables **S15/S18** report estimates and LRT *p* values; **S16/S19** show nested-model comparisons; **S17/S20** show A/C/E decomposition. `*_INTERNAL.tsv` files (Wald *p*, Bonferroni flags) stay in the repository for audit only. See script **04d** (`04d_significance_audit.R`).

## Bivariate pairings

All pairings are fit on the raw scale (script **8** or **9**). Each pairing contributes three Bonferroni-tested correlations (*r*<sub>g</sub>, *r*<sub>c</sub>, *r*<sub>e</sub>); phenotypic *r*<sub>P</sub> is estimate + 95 % CI only.

**Stress × stress (10 pairs, Table S15 / Figure 2):**

| | |
|---|---|
| `Childhoodtrauma` × `Loneliness` | `Childhoodtrauma` × `Lack_of_social_support` |
| `Childhoodtrauma` × `Stressful_life_events` | `Childhoodtrauma` × `Life_difficulties` |
| `Loneliness` × `Lack_of_social_support` | `Loneliness` × `Stressful_life_events` |
| `Loneliness` × `Life_difficulties` | `Lack_of_social_support` × `Stressful_life_events` |
| `Lack_of_social_support` × `Life_difficulties` | `Stressful_life_events` × `Life_difficulties` |

**Personality × stress (15 pairs, Table S18 / Figure 3):** every combination of `Neuroticism`, `Extraversion`, or `Conscientiousness` with each of the five stress phenotypes above (3 × 5).

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

`asreml-R` is required for scripts 7–11; obtain an academic licence from [vsni.co.uk/software/asreml-r](https://vsni.co.uk/software/asreml-r/).

## Manuscript output map

| Manuscript output | Script(s) | Key repo outputs |
|-------------------|-----------|------------------|
| **Table 1** | 3 | `output/tables/table_1_descriptives_overall.tsv` |
| **Figure 1** | 04, 07 | `table_S8_*`; `output/figures/figure_1_univariate_ace.*` |
| **Figure 2** | 04b, 07 | `table_S15_*`, `table_S17_*`; `figure_2_stress_stress.*` |
| **Figure 3** | 04c, 07 | `table_S18_*`, `table_S20_*`; `figure_3_personality_stress.*` |
| **Suppl. S5–S6** | 03d | `table_S5_lte_ldi_cross_wave_correlations.tsv`, `table_S6_lte_ldi_factor_loadings.tsv` |
| **Suppl. S7a–S7b** | 03 | `table_S7a_descriptives_by_sex.tsv`, `table_S7b_sex_differences_age_adjusted.tsv` |
| **Suppl. S7c** | 03b | `table_S7c_sex_differences_anova.tsv` |
| **Suppl. S7d** | 03c | `table_S7d_age_bracket_descriptives.tsv` |
| **Suppl. S8–S8e** | 04 | `table_S8_univariate_ace_full.tsv`, `table_S8b`–`S8e` |
| **Suppl. S9** | 04a | `table_S9_univariate_model_comparison_long.tsv` |
| **Suppl. S11–S11e, S12** | 05, 04a | `table_S11_int_*`; `table_S12_int_model_comparison_long.tsv` |
| **Suppl. S13–S14** | 06 | `table_S13_ct_dichotomization_prevalence.tsv`, `table_S14_ct_rrr_liability_h2.tsv` |
| **Suppl. S12b** | 06b | `table_S12b_ct_subscale_heritability.tsv` (exploratory) |
| **Suppl. S15–S17** | 04b | `table_S15_bivariate_stress_stress.tsv`, `S16`, `S17` |
| **Suppl. S18–S20** | 04c | `table_S18_bivariate_personality_stress.tsv`, `S19`, `S20` |
| **Suppl. Figs S3–S21** | 15 | `output/figures/` (see header of `07b_figures_supplementary.R`) |
| **Audit tables** | 10 | `table_audit_*.tsv`, `audit_summary.txt` (repo only) |
| **Logbook** | n/a | [`docs/Lifelines_Logbook.Rmd`](docs/Lifelines_Logbook.Rmd) |

## Data availability

Individual-level Lifelines data cannot be shared because of the participant consent. Bona-fide researchers can apply for access to the same extraction via the Lifelines Research Office (<https://www.lifelines.nl/researcher>); the analytic dataset can then be opened inside the Lifelines secure HPC environment using the scripts in this repository. See [`data/README.md`](data/README.md) for the list of files the scripts expect.

## Open Science

EJP Open Materials compliance: [OPEN_SCIENCE_STATEMENT.md](OPEN_SCIENCE_STATEMENT.md). Repository archived on Zenodo (DOI in `CITATION.cff`) at acceptance.

## Affiliation

This study uses data from the [Lifelines Cohort Study](https://www.lifelines.nl), a multidisciplinary research project on health and health-related behaviours in a large representative sample of the population of the north of the Netherlands.

This work was carried out in the context of the [Stress in Action](https://stress-in-action.nl) consortium (NWO Gravitation programme), which studies daily-life stress and its impact on mental and cardiometabolic health.

## Citation

> Reichelt, F., Wang, R., Triatin, R. D., Bos, M., Kavousi, M., de Geus, E. J. C., Snieder, H. *Familial Correlation between Personality and Stress Exposure: Insights from the Multigenerational Lifelines Cohort Study.* European Journal of Personality (in revision).

> Reichelt, F. (2026). *lifelines-personality-stress-correlations*. GitHub: https://github.com/FelixReichelt13/lifelines-personality-stress-correlations

See [CITATION.cff](CITATION.cff) for machine-readable metadata.

**Contact:** Felix Reichelt, University Medical Center Groningen, f.reichelt@umcg.nl
