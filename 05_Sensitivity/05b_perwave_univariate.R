#### =========================================================================
#### 05b_perwave_univariate.R
#### -------------------------------------------------------------------------
#### Paper     : Reichelt F. et al. -- Familial Correlation between Personality
####             and Stress Exposure: Insights from the Multigenerational
####             Lifelines Cohort Study. European Journal of Personality
####             (manuscript EJP-26-4766, in revision).
#### Author    : Felix Reichelt  <f.reichelt@umcg.nl>
####
#### What this script does
#### ---------------------
#### Fits the same univariate pedigree-based REML linear mixed model as
#### script 04 separately to *each individual wave* of LTE and LDI,
#### rather than to the latent factor that combines them. The output is
#### Suppl. Table S9 (per-wave univariate variance components), which
#### shows that the heritability of the wave-level sum scores is
#### substantially lower than the heritability of the *latent* trait
#### produced by the lavaan one-factor model (Suppl. Fig. S1 + S2,
#### Suppl. Table S6). The discrepancy is the empirical justification
#### for the latent-factor step used in the manuscript Methods.
####
#### Model (identical to script 04; see header there for the full
#### derivation):
####
####     y_iw = mu + b1*sex_i + b2*age_iw + b3*age_iw^2 + a_i + c_i + e_iw,
####
####     a ~ N(0, A * sigma_a^2)   additive genetic, A = numerator
####                                relationship matrix (`ainverse`).
####     c ~ N(0, I * sigma_c^2)   shared household environment.
####     e ~ N(0, I * sigma_e^2)   residual / unique environment.
####
#### One model is fitted per (instrument, wave) cell, eight cells total
#### (LTE x {1A, 1B, 1C, 2A} + LDI x {1A, 1B, 1C, 2A}).
####
#### Inputs    : data/data_prepared.Rdata (from script 02). The script
####             expects the per-wave sum-score columns
####               LTE_Baseline_1A, LTE_Followup_1B, LTE_Followup_1C,
####               LTE_Second_Assess_2A,
####               LDI_Baseline_1A, LDI_Followup_1B, LDI_Followup_1C,
####               LDI_Second_Assess_2A
####             produced by `02a_score_construction.R` (LDI uses items
####             1-9 only -- items 10-12 were added in later waves).
####
#### Outputs   : output/tables/table_S9_perwave_lte_ldi_variance.tsv
####
#### Manuscript: Methods 2.2 (latent factor justification);
####             Suppl. Table S9 (per-wave ACE);
####             Suppl. Fig. S1 + S2 (latent factor schematic).
####
#### Pitfalls
#### --------
#### * The per-wave sums have a lot more residual variance than the
####   latent factor (single-wave noise + acute events that don't
####   recur). Expect h^2 ~ 5-15 % for the LTE waves and h^2 ~ 15-25 %
####   for the LDI waves, well below the latent-factor h^2 reported in
####   Table 2.
#### * AGE_c and AGE_c2 in this script are *wave-specific* (i.e. the
####   participant's age at that wave), not the latent-factor age used
####   in script 04. `02_data_preparation.R` writes both: AGE_c at the
####   latent-factor anchor wave (Baseline 1A) is the one used in
####   scripts 03 - 07; per-wave ages live in AGE_c_<wave>.
#### * Wave 2A typically has the smallest analytic N. If the model fails
####   to converge for that wave, switch to AE (drop `idv(HouseID1)`) --
####   the household-level variance is usually tiny at 2A.
####
#### License   : ASReml-R is a commercial product of VSN International.
####             https://vsni.co.uk/software/asreml-r/
#### Version history:
#### 2026-07-15  v1.1  Ani-aligned multiple testing (/16 uni, /75 bivar);
####                   per-component 1-df bivariate LRTs; INTERNAL tables
#### 2026-07-15  v1.0  Initial release
#### =========================================================================


## Required libraries -----------------------------------------------------
##
## Loaded centrally by `01_setup.R`; restated here so the script reads
## as self-contained.

library(asreml)        # commercial; `asreml()`, `ainverse()`, `vpredict()`
library(data.table)    # `copy()`, in-place `:=`, `[, .(col1, col2)]` syntax


load(file.path(DATA_DIR, "data_prepared.Rdata"))   # provides `baseline` with the per-wave sum-score columns


## Step 1 -- Define the eight (instrument, wave) cells ---------------------

perwave_phenotypes <- data.frame(
  Inventory = c(rep("LTE", 4), rep("LDI", 4)),
  Wave      = rep(c("Baseline_1A", "Followup_1B",
                    "Followup_1C", "Second_Assess_2A"), 2),
  Column    = c("LTE_Baseline_1A", "LTE_Followup_1B",
                "LTE_Followup_1C", "LTE_Second_Assess_2A",
                "LDI_Baseline_1A", "LDI_Followup_1B",
                "LDI_Followup_1C", "LDI_Second_Assess_2A"),
  Age_col   = c("AGE_c_1A", "AGE_c_1B", "AGE_c_1C", "AGE_c_2A",
                "AGE_c_1A", "AGE_c_1B", "AGE_c_1C", "AGE_c_2A"),
  stringsAsFactors = FALSE
)


## Step 2 -- Fit ACE per (instrument, wave) cell ---------------------------

fit_ace_wave <- function(dt, column, age_col) {

  d <- copy(dt)                                           # data.table copy so we don't mutate `baseline`
  d[, PHENO := get(column)]                               # response column for this wave
  d[, AGE_c  := get(age_col)]                             # wave-specific centred age (AGE_c_1A, AGE_c_1B, ...)
  d[, AGE_c2 := AGE_c^2]                                  # quadratic age term, recomputed per wave
  d <- d[!is.na(PHENO) & !is.na(AGE_c)]                   # drop rows where the wave wasn't observed or age is missing

  ped  <- unique(d[, .(PSEUDOIDEXT, FatherID, MotherID)]) # pedigree built on the analytic frame (one row per individual)
  Ainv <- ainverse(ped)                                   # sparse inverse A matrix

  full <- asreml(
    fixed     = PHENO ~ GENDER + AGE_c + AGE_c2,          # same fixed effects as script 04
    random    = ~ vm(PSEUDOIDEXT, Ainv) + idv(HouseID1),  # additive genetic + shared household
    data      = d,
    na.action = na.method(x = "omit"),                    # complete-case fit (already filtered above)
    workspace = "2gb",                                    # raise if asreml asks for more memory
    ai.sing   = TRUE,                                     # tolerate AI singularities on huge pedigrees
    maxiter   = 50
  )

  h2 <- vp_row(vpredict(full, h2 ~ V2 / (V1 + V2 + V3)))  # V1 = c2, V2 = h2, V3 = e2  (see script 04 header)
  c2 <- vp_row(vpredict(full, c2 ~ V1 / (V1 + V2 + V3)))
  e2 <- vp_row(vpredict(full, e2 ~ V3 / (V1 + V2 + V3)))

  data.frame(
    N     = nrow(d),
    h2    = h2$estimate, h2_SE = h2$se,
    c2    = c2$estimate, c2_SE = c2$se,
    e2    = e2$estimate, e2_SE = e2$se
  )
}


## Step 3 -- Loop over the eight cells -------------------------------------

results <- vector("list", nrow(perwave_phenotypes))
for (i in seq_len(nrow(perwave_phenotypes))) {

  row     <- perwave_phenotypes[i, ]
  message("Per-wave ACE: ", row$Inventory, " - ", row$Wave)

  fit             <- fit_ace_wave(baseline, row$Column, row$Age_col)
  results[[i]]    <- cbind(row, fit)
}
results <- do.call(rbind, results)


## Step 4 -- Write Suppl. Table S9 -----------------------------------------

save_table(results, "table_S9_perwave_lte_ldi_variance")

message("Step 06b complete. Per-wave ACE table written to ",
        file.path(OUT_TABLES, "table_S9_perwave_lte_ldi_variance.tsv"))
