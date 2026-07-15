#### =========================================================================
#### 06_rrr_liability_heritability.R
#### -------------------------------------------------------------------------
#### Paper     : Reichelt F. et al. -- Familial Correlation between Personality
####             and Stress Exposure: Insights from the Multigenerational
####             Lifelines Cohort Study. European Journal of Personality
####             (manuscript EJP-26-4766, in revision).
#### Author    : Felix Reichelt  <f.reichelt@umcg.nl>
####
#### What this script does
#### ---------------------
#### Estimates the liability-scale heritability of *dichotomous* childhood
#### trauma from CTQ-SF subscale moderate/severe thresholds (Bernstein
#### et al., 2003; see `02_data_preparation.R` and `01_helpers.R`) from the
#### first-degree-relative recurrence risk ratio (RRR) and the population
#### prevalence, following Wray & Gottesman (2012).
####
#### Two-step procedure
#### ------------------
#### 1. RRR is estimated from a logistic regression of CT (1/0) on a
####    "has a first-degree relative with CT" indicator, age, age^2 and
####    sex. We use `marginaleffects::avg_predictions()` to obtain the
####    average predicted probability in the general population and in
####    relatives of cases, and form the ratio
####
####        RRR = P(CT = 1 | FDR with CT) / P(CT = 1 | population).
####
####    SEs are robust to clustering on the family ID (`sandwich::vcovCL`).
####
#### 2. The Wray & Gottesman (2012) liability-threshold model converts
####    K (population prevalence) and the prevalence in relatives
####    (K_r = RRR * K) into a heritability on the liability scale:
####
####        t   = qnorm(1 - K)
####        t_r = qnorm(1 - K_r)
####        h^2_liab = (t - t_r) / (a * t),  with a = 0.5 for FDR.
####
#### Inputs    : data/data_prepared.Rdata (from script 02); the dataset
####             must contain `CT_dichotomous` (1/0 from subscale thresholds
####             in script 02) and a column flagging participants whose
####             first-degree relative has CT (the script searches for
####             `fdr_ct`, `FDR_CT` or `fdr_childhoodtrauma`).
#### Outputs   : output/tables/table_S11_ct_dichotomization_prevalence.tsv
####             output/tables/table_S12_ct_rrr_liability_h2.tsv
####
#### Manuscript: Methods 2.4 ("Sensitivity analyses"); Suppl. Table S11
####             (prevalence by sex), Suppl. Table S12 (RRR + liability
####             h^2). The reported manuscript values are
####             K = 14.2 %, K_male = 11.1 %, K_female = 16.3 %,
####             RRR = 1.84 (95 % CI 1.68-2.00), h^2_liab = 57.2 %
####             (95 % CI 48.7-65.8 %).
####
#### Method citation
#### ---------------
#### The recurrence-risk-ratio step follows Bos et al. (2025,
#### Psychological Medicine, doi:10.1017/S003329172500100X), who used
#### the same Lifelines pedigree to estimate familial coaggregation
#### of internalising disorders. The liability-threshold conversion
#### follows Wray & Gottesman (2012, Genome Research).
####
#### Implementation note
#### -------------------
#### The Wray-Gottesman liability-threshold conversion is implemented
#### here using `marginaleffects::avg_predictions()` for the RRR step
#### with cluster-robust standard errors (clustered at the family
#### level) via `sandwich::vcovCL()`. Subscale-level liability h² at
#### subscale-specific Bernstein thresholds is in `06b_ct_subscale_heritability.R`.
####
#### CT coding  : Same rule as script 02 — at least one subscale >= its
####              moderate threshold (ea >= 12, pa >= 9, sa >= 7, en >= 14,
####              pn >= 9); coded 0 only when all five observed and none met.
####
#### Pitfalls
#### --------
#### * The Wray method assumes that the disorder follows a normally
####   distributed liability and that relatives' liabilities correlate
####   exactly with the additive genetic coefficient `a` (0.5 for FDRs).
####   When dominance is large or shared environment is non-negligible
####   this can bias h^2 upward.
#### * `K * RRR` must be in (0, 1). The script clamps to (1e-4, 1 - 1e-4)
####   so very small samples don't blow up at the tails.
####
#### License   : This script does NOT require ASReml-R. It uses CRAN
####             packages only, so reviewers without a license can run
####             and inspect it.
#### Version history:
#### 2026-07-15  v1.2  Align with subscale-specific CT dichotomisation (01)
#### 2026-07-15  v1.1  Ani-aligned multiple testing (/16 uni, /75 bivar);
####                   per-component 1-df bivariate LRTs; INTERNAL tables
#### 2026-07-15  v1.0  Initial release
#### =========================================================================


## Required libraries -----------------------------------------------------
##
## Loaded centrally by `01_setup.R`; restated here for self-documentation.
## No ASReml dependency: the heritability is derived from the
## recurrence-risk ratio of a logistic regression, not from REML.

library(data.table)        # `baseline` is a `data.table`; needed for `[i, j]` indexing
library(marginaleffects)   # `avg_predictions()` for marginal probabilities P(CT=1 | ...)
library(sandwich)          # `vcovCL()` for cluster-robust standard errors (clustered on FamID)


load(file.path(DATA_DIR, "data_prepared.Rdata"))   # provides `baseline` with the `CT_dichotomous` column added by script 02


## Step 1 — Find the FDR indicator column ------------------------------

fdr_col <- intersect(c("fdr_ct", "FDR_CT", "fdr_childhoodtrauma"),
                     names(baseline))[1]                          # auto-detect column name across data versions
stopifnot(!is.na(fdr_col))                                        # fail fast if no FDR flag column is found

d <- baseline[!is.na(CT_dichotomous) & !is.na(get(fdr_col))]      # analytic sample: both CT status and FDR flag observed


## Step 2 — Prevalences in the population and in FDRs of cases --------
##
## We also compute K separately for males and females (Suppl. Table
## S11). The manuscript reports K_male = 11.1 % and K_female = 16.3 %.

K       <- mean(d$CT_dichotomous)
N       <- nrow(d)
K_fdr   <- mean(d$CT_dichotomous[d[[fdr_col]] == 1])
K_male  <- mean(d$CT_dichotomous[d$GENDER == 1L])
K_female<- mean(d$CT_dichotomous[d$GENDER == 2L])
N_male  <- sum(d$GENDER == 1L)
N_female<- sum(d$GENDER == 2L)

message(sprintf("Population prevalence K = %.4f  (N = %d)", K, N))
message(sprintf("Prevalence in FDR cases  = %.4f", K_fdr))
message(sprintf("Male prevalence   = %.4f  (N = %d)", K_male,   N_male))
message(sprintf("Female prevalence = %.4f  (N = %d)", K_female, N_female))


## Step 2b -- Write Suppl. Table S11 (prevalence by sex) -------------

prev_table <- data.frame(
  Group       = c("Overall", "Male", "Female"),
  N           = c(N, N_male, N_female),
  N_cases     = c(sum(d$CT_dichotomous),
                  sum(d$CT_dichotomous[d$GENDER == 1L]),
                  sum(d$CT_dichotomous[d$GENDER == 2L])),
  Prevalence  = c(K, K_male, K_female)
)
save_table(prev_table, "table_S11_ct_dichotomization_prevalence")


## Step 3 — RRR from logistic regression with cluster-robust SEs ------

fml   <- as.formula(paste0("CT_dichotomous ~ ", fdr_col,
                           " + AGE_c + AGE_c2 + GENDER"))         # CT (0/1) ~ FDR (0/1) + age + age^2 + sex
model <- glm(fml, data = d, family = binomial())                  # logistic regression on the dichotomous outcome

denom <- avg_predictions(
  model,
  type = "response",                                              # response-scale prediction = P(CT = 1 | covariates)
  vcov = vcovCL(model, cluster = d$FamID)                         # cluster-robust covariance: relatives share families
)                                                                 # `denom`: average P(CT = 1) in the full sample
numer <- avg_predictions(
  model,
  variables = fdr_col,                                            # condition on the FDR indicator (computes one row per level)
  type      = "response",
  vcov      = vcovCL(model, cluster = d$FamID)
)
numer <- numer[numer[[fdr_col]] == 1, ]                           # keep only the row for FDR == 1 (relatives of cases)

RRR    <- numer$estimate / denom$estimate                         # recurrence risk ratio = P(CT|FDR) / P(CT|population)

## Delta-method SE for the ratio of two predicted probabilities.
RRR_se <- sqrt(
  (numer$std.error / denom$estimate)^2 +
  (numer$estimate * denom$std.error / denom$estimate^2)^2
)
RRR_lo <- RRR - 1.96 * RRR_se                                     # symmetric Wald 95 % CI on the RRR scale
RRR_hi <- RRR + 1.96 * RRR_se

message(sprintf("RRR = %.3f  (95%% CI %.3f - %.3f)", RRR, RRR_lo, RRR_hi))


## Step 4 — Liability heritability (Wray & Gottesman, 2012) -----------

liability_h2 <- function(K, RRR, a = 0.5) {
  t   <- qnorm(1 - K)                                # liability threshold in the population (upper tail)
  Kr  <- pmin(pmax(K * RRR, 1e-4), 1 - 1e-4)         # prevalence in FDRs, clamped to (0, 1) to avoid Inf
  t_r <- qnorm(1 - Kr)                               # liability threshold in relatives of cases
  (t - t_r) / (a * t)                                # Wray & Gottesman (2012) eq. 1; a = 0.5 for first-degree relatives
}

h2    <- liability_h2(K, RRR)
h2_lo <- liability_h2(K, RRR_lo)
h2_hi <- liability_h2(K, RRR_hi)


## Step 5 — Save ------------------------------------------------------

out <- data.frame(
  Phenotype        = "Childhood_Trauma",
  N                = N,
  Prevalence_K     = K,
  Prevalence_K_FDR = K_fdr,
  RRR              = RRR,
  RRR_CI_lower    = RRR_lo,
  RRR_CI_upper    = RRR_hi,
  h2_liability    = h2,
  h2_CI_lower     = h2_lo,
  h2_CI_upper     = h2_hi
)
save_table(out, "table_S12_ct_rrr_liability_h2")

message("Step 06 complete. Tables S11 and S12 written.")
