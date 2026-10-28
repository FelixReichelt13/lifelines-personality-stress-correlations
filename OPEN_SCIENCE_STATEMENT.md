# Open Science Statement

This document describes the availability of the study materials, data, and analysis scripts, 
together with the reporting status of the eight European Journal of Personality Open Science standards. 
Where a requirement cannot be fully met, the relevant limitations are explained.

## Permanent links

| Item | URL |
|---|---|
| Public repository and documentation: | <https://github.com/FelixReichelt13/lifelines-personality-stress-correlations> |
| Archived release | Not yet available; DOI to be added after archiving. |
| Lifelines cohort information and data access | <https://www.lifelines.nl/researcher> |
| Lifelines cohort overview paper (Sijtsma et al., 2022) | <https://doi.org/10.1093/ije/dyab257> |
| ASReml-R 4.2 (commercial software for ACE models) | <https://vsni.co.uk/software/asreml-r> |

## 1. Preregistration

The hypotheses and analysis plan were not preregistered. 
This study involved secondary analysis of existing data from the Lifelines Cohort Study. 
The primary hypotheses and sensitivity analyses are described in the manuscript and supplementary information.

The following analyses are explicitly *exploratory*:

- CTQ-SF subscale heritability (Suppl. Table S12b, script `07b_ct_subscale_heritability.R`)
- Inverse-normal-transformation sensitivity (Suppl. Table S10, script `06_sensitivity_int.R`)
- Age-bracket descriptives (Suppl. Table S7d / Suppl. Figs S3 – S16, script `09_age_brackets.R`)
- Liability-scale heritability at the low Bernstein cut-off (sensitivity row in Suppl. Table S12b)

## 2. Sampling

Sample size was *not* set by an a-priori power calculation; we use the full
available Lifelines extraction.

| Quantity | Value |
|---|---|
| Total Lifelines cohort | 152,722 |
| Adults aged ≥ 18 yr (analytic sample) |Adult eligibility was defined as age ≥ 18 years. The number of participants included in each analysis varied according to phenotype availability and model requirements. Phenotype-specific sample sizes and bivariate sample sizes are reported in the manuscript and supplementary information.|
| Multigenerational families | 30,914 |
| Singletons | 40,496 |
| Demographic composition | Northern Netherlands, ~58 % female, mean age ~46 yr (range 18 – 93) |

Exclusions:

| Exclusion | Rationale | Script |
|---|---|---|
| Participants < 18 yr at the relevant wave | Removes parent-proxy responses on SPF-IL Social Support that bias estimates | `01_data_preparation.R` |
| Pedigree members with unknown / inconsistent parental IDs | `asreml::ainverse()` requires consistent factor IDs with explicit `NA` for unknown parents | `01_data_preparation.R` |
| Wave-level rows with any missing LTE / LDI item | Discrete events are not theoretically suitable for mean-imputation | `01a_score_construction.R` |

## 3. Open material

Every questionnaire item, recoding rule, imputation rule and scale-construction
step is documented in `docs/Lifelines_Logbook.Rmd` §4 ("Item-level scoring") and in
the executable code in `01a_score_construction.R`.

| Measure | Items used | Primary citation | Public source (Lifelines Wiki) |
|---|---|---|---|
| CTQ-SF (Childhood Trauma Questionnaire — Short Form) | 25 clinical + 3 minimisation/denial (wave 2B only) | Bernstein et al., 1998; Thombs et al., 2009 | <https://wiki.lifelines.nl/doku.php?id=childhood_trauma_ctq> |
| De Jong Gierveld Loneliness Scale (6-item short form) | 6 (wave 2B) | de Jong Gierveld & van Tilburg, 2006 | <https://wiki.lifelines.nl/doku.php?id=loneliness_ls> |
| SPF-IL Social Production Function Instrument | 9 social-well-being items (of 15 total) | Nieboer et al., 2005; Elsman et al., 2022 | <https://wiki.lifelines.nl/doku.php?id=social_functioning_spf-il> |
| List of Threatening Experiences (LTE) | 12 binary events × 4 waves (1A, 1B, 1C, 2A) | Brugha & Cragg, 1990; Rosmalen et al., 2012 | <https://wiki.lifelines.nl/doku.php?id=stress_lte> |
| Long-term Difficulties Inventory (LDI) | The full instrument comprises 12 items in waves (1A, 1B, 1C, 2A) | Rosmalen et al., 2012 |
| NEO-PI-R Neuroticism + NEO-FFI Extraversion / Conscientiousness | Neuroticism: 48-item NEO-PI-R scale covering six facets (range 48–240); Extraversion and Conscientiousness: 12-item NEO-FFI scales (range 12–60 each). All three personality scores were standardised before analysis. | Costa & McCrae, 1992; Hoekstra et al., 2007 | <https://wiki.lifelines.nl/doku.php?id=personality_neo>, <https://wiki.lifelines.nl/doku.php?id=personality_sum_scores_neo> |

## 4. Open data

Individual-level Lifelines data **cannot** be made openly accessible. Lifelines
is a closed-access research infrastructure governed by Dutch data-protection
regulation and a participant-consent framework that does not permit
redistribution of individual-level data outside the Lifelines secure HPC
environment.

Bona-fide researchers can apply for access at
<https://www.lifelines.nl/researcher>. The exact phenotype / pedigree
extraction used here is documented in [`data/README.md`](data/README.md) so
that an approved researcher can recreate it.

A synthetic dataset is not provided. Reproducing the primary variance decomposition analyses 
requires access to the individual-level Lifelines data, including pedigree and household information. 
The publicly available scripts and documentation support reproduction by researchers 
who obtain authorised access to the relevant data and software.

## 5. Reproducible script

All scripts are openly accessible in the public repository above. At
acceptance the commit producing the print version of the manuscript will be
tagged `v1.0` and archived on Zenodo for a permanent DOI.

Reproducibility tiers:

| Tier | Scripts | What you need |
|---|---|---|
| Fully standalone (no data, no ASReml) | `02b_anova_sex_diff.R`, `07b_ct_subscale_heritability.R`, `08_figures.R`, `10_figures_supplementary.R` | CRAN R only |
| Needs Lifelines data, no ASReml | `01a_score_construction.R`, `01_data_preparation.R`, `02_descriptives.R`, `09_age_brackets.R`, `11_cross_wave_correlations.R` | Lifelines secure environment |
| Needs ASReml-R 4.2 and Lifelines data | `03_univariate_heritability.R`, `04_bivariate_stress_stress.R`, `05_bivariate_personality_stress.R`, `06_sensitivity_int.R`, `07_rrr_heritability_childhoodtrauma.R` | Lifelines secure environment + ASReml-R licence |

The manuscript reports analyses conducted using R version 4.2 and ASReml-R. 
The exact software environment and package versions should be taken from the final analysis session information. 
ASReml-R requires a licence.

## 6. Transformations and non-standard scoring

| Transformation | Where | Standard procedure also reported |
|---|---|---|
| Rank-based inverse normal transformation on stress phenotypes | `06_sensitivity_int.R` (sensitivity) | Raw-scale ACE in `03_univariate_heritability.R` (primary) |
| Sign flip for SPF-IL (`Lack_of_social_support = max(ss_total) − ss_total`) | `01a_score_construction.R` | Raw SPF-IL retained for verification |
| Dichotomisation of CT for liability h² (any subscale >= moderate threshold: ea ≥ 12, pa ≥ 9, sa ≥ 7, en ≥ 14, pn ≥ 9; see `01_data_preparation.R`) | `07_rrr_heritability_childhoodtrauma.R` | Continuous CTQ-SF total in `03_univariate_heritability.R` |
| Standardisation of LTE / LDI latent factor scores to SD = 1 | `01_data_preparation.R` | Unstandardised lavaan output in Suppl. Table S6 |

Non-standard scoring:

| Convention | Justification |
|---|---|
| NEO `_SUM3` (mean-imputed sum, rescaled by `(SUM2/NVALID) × k`) | The only variant available for all three traits (Lifelines does not distribute `SUM1_EXTRA`). Correlates ~0.93 with `_SUM1` where the latter exists. |
| CTQ-SF mean imputation when ≤ 2 of 5 subscale items missing | Bernstein et al. (1998) convention. |
| LTE / LDI complete-case (no imputation) | Discrete-event imputation is not theoretically defensible. |
| Bernstein moderate rule (any CTQ-SF subscale ≥ subscale-specific cut-off; `01_data_preparation.R`) for liability h² | Clinical convention (Bernstein et al., 2003). Subscale liability in Suppl. Table S12b. |

## 7. Covariates

All variance-decomposition models adjust for **sex** and a **quadratic age
trend**:

```r
fixed = phenotype ~ sex + age + I(age^2)                                # univariate
fixed = cbind(y1, y2) ~ trait + trait:sex + trait:age + trait:I(age^2)  # bivariate
```

Sex, age, and age squared were included to account for demographic differences in trait means before 
estimating residual genetic and environmental variance components. The quadratic age term allows for 
non-linear age associations. The primary analyses used these adjusted models. 
Parallel unadjusted variance-decomposition results are not presented.

## 8. Effects

Variance proportions and correlations are reported with corresponding standard errors, 
and confidence intervals are provided where available. Statistical significance of variance components 
and cross-trait covariances was evaluated using the likelihood-ratio tests described in the manuscript. 
The supplementary tables provide the corresponding test statistics and p-values, 
with very small p-values sometimes reported as p < .001.

Two separate Bonferroni families are applied for the primary analyses;
sensitivity and exploratory analyses are not corrected. The design follows
the per-component likelihood-ratio protocol used in our group's prior
Lifelines ACE work (Ani et al., 2026, *Diabetes Obes Metab*).

| Family | Tests | Threshold |
|---|---|---|
| Univariate ACE (1-df LRT for h² = 0 and c² = 0) | 8 phenotypes × 2 variance components = 16 | α = .05 / 16 = .003125 |
| Bivariate ACE (1-df LRT per cross-trait correlation) | 25 pairings × 3 correlations (rg, rc, re) = 75 | α = .05 / 75 ≈ .000667 |

**Univariate.** For each phenotype, h² is tested by ACE vs CE and c² by
ACE vs AE. Unique-environment variance (e²) is not tested (bounded,
always present). Each LRT is 1 df.

**Bivariate.** For each of the 25 pairings (10 stress × stress + 15
personality × stress), rg, rc and re are tested separately: the full
bivariate ACE model is compared with a sub-model in which only that one
cross-trait covariance is fixed to zero (1-df LRT each). All 75 p values
are Bonferroni-corrected with the same threshold α = .05 / 75. Phenotypic
correlation (rP) is reported with SE and 95 % CI but not LRT-tested.

**Wald caveat.** The formal bivariate tests are the per-component LRT
p values (`p_rg`, `p_rc`, `p_re`) in the model-comparison supplement
tables. Normal-approximation Wald p-values (`2 × pnorm(−|est / se|)`)
are retained in the internal audit file (`12_significance_audit.R`) only;
they can mislead near boundaries (rg ≈ ±1, c² ≈ 0).

Internal consistency (Cronbach's α) is computed in `01a_score_construction.R`
Step 6 for the multi-item stress scales. Published reliability estimates for 
the Dutch NEO instruments are distinguished from internal-consistency estimates 
calculated in the present study. 

## Compliance summary

| EJP standard | Compliant? | Where |
|---|---|---|
| 1. Preregistration | Not preregistered, disclosed | §1 |
| 2. Sampling | Yes — Sampling information reported | §2 |
| 3. Open material | Public documentation available | §3 |
| 4. Open data | Closed access | §4 |
| 5. Reproducible script |Publicly available | §5 |
| 6. Transformations | Yes — raw + INT sensitivity; standard imputation conventions | §6 |
| 7. Covariates | Yes — sex, age, age² with documented rationale | §7 |
| 8. Effects | Yes — 	Estimates and tests reported | §8 |

The authors request consideration for the Open Materials badge on the basis of the publicly 
accessible study documentation and analysis materials. Individual-level Lifelines data 
are not openly available, and the study was not preregistered. 
Badge eligibility is determined by the journal.
