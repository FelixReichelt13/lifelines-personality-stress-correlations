#### =========================================================================
#### 04c_bivariate_personality_stress.R
#### -------------------------------------------------------------------------
#### Paper     : Reichelt F. et al. -- Familial Correlation between Personality
####             and Stress Exposure: Insights from the Multigenerational
####             Lifelines Cohort Study. European Journal of Personality
####             (manuscript EJP-26-4766, in revision).
#### Author    : Felix Reichelt  <f.reichelt@umcg.nl>
#### Cleaning  : item-level questionnaires were originally cleaned by
####             Dr. Rujia Wang (CTQ, Loneliness, Social Support) and
####             Lizanne JS Schweren (LTE, LDI) -- see 02a_score_construction.R.
####
#### What this script does
#### ---------------------
#### Fits a bivariate ACE linear mixed model for every combination of one
#### personality trait (Neuroticism, Extraversion, Conscientiousness) and
#### one stress phenotype (Childhood Trauma, Loneliness, Lack of Social
#### Support, Stressful Life Events, Life Difficulties). 3 x 5 = 15
#### models in total. The genetic (rg), shared-environmental (rc) and
#### unique-environmental (re) correlations feed **Figure 3** (left
#### panels) and **Suppl. Table S14**. The A/C/E decomposition of r_P
#### feeds the right panels of Figure 3 and **Suppl. Table S16**.
####
#### How ASReml-R is used here
#### -------------------------
#### Same bivariate setup as in script 04b -- see the header of
#### `04b_bivariate_stress_stress.R` for the full V1..V10 variance-component
#### mapping (shared env in V1-V3, genetic in V4-V6, residual scale at V7,
#### residual covariances in V8-V10) and the rg / rc / re / rp formulas
#### (ASReml-R Reference Manual 4.2, Section 4.6, "Estimates from
#### functions of variance components"). The only thing that changes
#### between scripts 04 and 05 is which two columns of the data are
#### plugged in as `P1` and `P2`.
####
#### Multiple testing -- per-component 1-df LRTs (Ani et al. design)
#### ----------------------------------------------------------------
#### Together with `04b_bivariate_stress_stress.R` (10 stress x stress
#### pairs), this script contributes 15 of the **25 bivariate pairings**
#### reported in the manuscript. Following Ani, Nolte, Vaez & Snieder
#### (2026, Diabetes Obes Metab 28:3344), each of the three cross-trait
#### correlations (rg, rc, re) is tested SEPARATELY with a 1-df likelihood-
#### ratio test comparing the full unstructured model against a sub-model
#### with only that one covariance fixed to zero (see the header of
#### `04b_bivariate_stress_stress.R` for the sub-model specifications and the
#### chi^2(1) / no-boundary-correction rationale). The Bonferroni-corrected
#### family-wise threshold is alpha = .05 / 75 = 6.67e-4 (25 pairings x 3
#### correlations = 75 tests), matching Ani et al.
####
#### With N > 140,000, statistical significance is expected even for
#### very small effects; the substantive interpretation in the
#### manuscript rests on effect sizes, confidence intervals and the
#### A/C/E decomposition of r_P rather than on p-values alone.
####
#### Wald p-value caveat
#### -------------------
#### Same Wald asymptotic caveat as in script 04b -- the r*_wald_p columns
#### (INTERNAL only) are normal-approximation two-sided p-values, which can
#### be misleading near correlation-space boundaries; the 1-df LRT is the
#### primary test. A footnote to S13 + S14 records this explicitly.
####
#### Trait-level expectations
#### ------------------------
#### Pre-specified expectations (manuscript Introduction):
####   * Neuroticism      : strongest positive overlap with all stress
####                        phenotypes (negative affect, stress
####                        reactivity, maladaptive coping).
####   * Conscientiousness: protective negative associations (self-
####                        regulation, goal-directed behaviour).
####   * Extraversion     : protective negative associations for the
####                        social-connectedness phenotypes, but
####                        potentially weaker or direction-dependent
####                        for stressful life events (stimulation
####                        seeking may also increase exposure).
####
#### Inputs    : data/data_prepared.Rdata (from script 02).
#### Outputs   : output/tables/table_S14_bivariate_personality_stress.tsv
####             output/tables/table_S16_decomp_personality_stress.tsv
####
#### Manuscript: Methods 2.3, Results 3.4, Figure 3,
####             Suppl. Tables S14 (correlations) + S16 (decomposition).
####
#### Pitfalls
#### --------
#### * The personality file is a smaller release than the stress file, so
####   the effective N for each pair is roughly N(personality) -- not the
####   full N from Table 1. Check `nrow(d)` after the NA filter.
#### * For trait pairs with very different scales (e.g. NEO sum scores
####   versus a 0/5 life-event count), starting values can matter. If a
####   model fails to converge, refit it with start values copied from a
####   univariate model.
####
#### License   : ASReml-R is commercial. https://vsni.co.uk/software/asreml-r/
#### Version history:
#### 2026-07-15  v1.1  Ani-aligned multiple testing (/16 uni, /75 bivar);
####                   per-component 1-df bivariate LRTs; INTERNAL tables
#### 2026-07-15  v1.0  Initial release
#### =========================================================================


## Required libraries -----------------------------------------------------
##
## Loaded centrally by `01_setup.R`. The actual fitting helper
## (`fit_biv()`) is defined in `04b_bivariate_stress_stress.R`, which
## must be sourced first in the same R session.

library(asreml)        # commercial; called inside `fit_biv()`
library(data.table)    # `fit_biv()` uses `data.table` internally


load(file.path(DATA_DIR, "data_prepared.Rdata"))   # provides `baseline`, `baseline_int`, `family_data`

traits <- intersect(PERSONALITY_PHENOTYPES, names(baseline))
stress <- intersect(STRESS_PHENOTYPES,       names(baseline))

pairs  <- expand.grid(Personality = traits, Stress = stress,
                      stringsAsFactors = FALSE)   # all 3 x 5 = 15 personality x stress combinations


## Step 1 — Reuse the bivariate function from script 04b ----------------
##
## Run `04b_bivariate_stress_stress.R` first in the same R session so
## that the function `fit_biv()` is in your global environment.


## Step 2 — Loop over the 15 personality x stress combinations --------

results <- lapply(seq_len(nrow(pairs)), function(i) {
  message("Bivariate (personality x stress): ",
          pairs$Personality[i], " x ", pairs$Stress[i])
  out <- fit_biv(baseline, pairs$Personality[i], pairs$Stress[i])
  names(out)[1:2] <- c("Personality", "Stress")
  out
})
results <- do.call(rbind, results)

## Per-component significance flags (Bonferroni over 75 tests; primary).
results$rg_sig_bonf <- results$p_rg < BONFERRONI_ALPHA
results$rc_sig_bonf <- results$p_rc < BONFERRONI_ALPHA
results$re_sig_bonf <- results$p_re < BONFERRONI_ALPHA

## Suppl. Table S14 (MANUSCRIPT): correlations, 95 % CIs, and per-component
## 1-df LRTs. Same table note as S13 (see 04b_bivariate_stress_stress.R).
save_table(results[, c("Personality", "Stress", "N",
                       "rp", "rp_SE", "rp_CI_lo", "rp_CI_hi",
                       "rg", "rg_SE", "rg_CI_lo", "rg_CI_hi", "LR_rg", "df_rg", "p_rg",
                       "rc", "rc_SE", "rc_CI_lo", "rc_CI_hi", "LR_rc", "df_rc", "p_rc",
                       "re", "re_SE", "re_CI_lo", "re_CI_hi", "LR_re", "df_re", "p_re")],
           "table_S14_bivariate_personality_stress")

## Suppl. Table S14c (MANUSCRIPT): per-component model comparison (full fit
## vs each single-covariance-dropped sub-model).
save_table(results[, c("Personality", "Stress", "N",
                       "logLik_full", "AIC_full", "BIC_full",
                       "logLik_noG", "AIC_noG", "BIC_noG", "LR_rg", "p_rg",
                       "logLik_noC", "AIC_noC", "BIC_noC", "LR_rc", "p_rc",
                       "logLik_noE", "AIC_noE", "BIC_noE", "LR_re", "p_re")],
           "table_S14c_personality_stress_model_comparison")

## Suppl. Table S16 (MANUSCRIPT): decomposition + unstable flag.
save_table(results[, c("Personality", "Stress", "N", "rp",
                       "A_pcor", "A_pcor_SE",
                       "C_pcor", "C_pcor_SE",
                       "E_pcor", "E_pcor_SE",
                       "prop_A", "prop_C", "prop_E",
                       "decomp_unstable")],
           "table_S16_decomp_personality_stress")

## INTERNAL audit file: every column including sig flags; NOT for upload.
save_table(results, "table_S14S16_bivariate_personality_stress_full_INTERNAL")


## Step 3 -- sessionInfo ------------------------------------------------

writeLines(capture.output(sessionInfo()),
           file.path(OUT_TABLES, "sessionInfo_bivariate_personality_stress.txt"))

message("Step 05 complete. ", nrow(results),
        " personality x stress bivariate models fitted; Tables S14, S14c, S16 (+ INTERNAL) written.")
