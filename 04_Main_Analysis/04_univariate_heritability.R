#### =========================================================================
#### 04_univariate_heritability.R
#### -------------------------------------------------------------------------
#### Paper     : Reichelt F. et al. -- Familial Correlation between Personality
####             and Stress Exposure: Insights from the Multigenerational
####             Lifelines Cohort Study. European Journal of Personality
####             (manuscript EJP-26-4766, in revision).
#### Author    : Felix Reichelt  <f.reichelt@umcg.nl>
####
#### What this script does
#### ---------------------
#### Fits the univariate pedigree-based REML linear mixed model for every
#### stress and personality phenotype and, for each phenotype, refits the
#### two nested submodels needed for the ACE model-comparison table:
####
####     ACE  =  vm(PSEUDOIDEXT, Ainv) + idv(HouseID1) + residual
####     AE   =  vm(PSEUDOIDEXT, Ainv)                 + residual
####     CE   =                          idv(HouseID1) + residual
####
#### From these three fits it extracts:
####
####   * h^2, c^2, e^2 with delta-method SEs from `vpredict()`;
####   * boundary-corrected LRT p-values for
####       - ACE vs CE  (tests A ; 1df mixture 0.5*chi^2_0 + 0.5*chi^2_1)
####       - ACE vs AE  (tests C ; same 1df mixture);
####   * logLik / AIC / BIC per model (ACE / AE / CE; Suppl. Table S8c);
####
#### NOTE (2026 revision): the joint ACE-vs-E test and the E-only model
#### have been removed. Under REML the E-only marginal likelihood is not
#### comparable to models with a different random-effects structure, and
#### for Stressful Life Events (a low-variance LTE factor score) the
#### E-only fit returned an invalid log-likelihood that exceeded its own
#### nested AE/CE parents. The two single-parameter tests (A via ACE-vs-CE
#### and C via ACE-vs-AE) each converge cleanly from the full ACE fit and
#### are the sole basis for component inference. See Methods and the
#### Table S8b footnote.
####   * the fixed-effect coefficients (age, age^2, sex; Suppl. Table S8d);
####   * the raw variance components (Suppl. Table S8e), with the
####     V-index annotated so the reader can cross-check the vpredict
####     formulas.
####
#### The linear mixed model is (ACE form)
####
####     y_i = mu + b1*sex_i + b2*age_i + b3*age_i^2 + a_i + c_i + e_i,
####
####     a ~ N(0, A * sigma_a^2)   additive genetic,
####                                A = numerator relationship matrix
####                                from the pedigree (via `ainverse`).
####     c ~ N(0, I * sigma_c^2)   shared household environment,
####                                modelled by `idv(HouseID1)`.
####     e ~ N(0, I * sigma_e^2)   residual / unique environment.
####
####     h^2 = sigma_a^2 / (sigma_a^2 + sigma_c^2 + sigma_e^2).
####
#### How ASReml-R is used here
#### -------------------------
#### * Variance components are extracted with `vpredict()`, which refers
####   to variance parameters by their position in the `vparameters`
####   slot of the asreml object using the strings "V1", "V2", ...
####   (ASReml-R Reference Manual 4.2, Section 4.6 -- "Estimates from
####    functions of variance components").
####
#### * For the univariate ACE model the order is:
####
####     V1 = HouseID1!HouseID1.var        (shared environment, c^2)
####     V2 = vm(PSEUDOIDEXT, Ainv)!ped    (additive genetic, h^2)
####     V3 = units!R                      (residual, e^2)
####
####   So:    h2 ~ V2 / (V1 + V2 + V3)
####          c2 ~ V1 / (V1 + V2 + V3)
####          e2 ~ V3 / (V1 + V2 + V3)
####
#### * The likelihood-ratio tests use boundary-corrected reference
####   distributions because sigma_a^2 = 0 and sigma_c^2 = 0 lie on the
####   boundary of the parameter space; see Self & Liang (1987) and
####   Visscher (2006) for the classical mixture-chi^2 result.
####
#### Inputs    : data/data_prepared.Rdata (from script 02).
#### Outputs   : output/tables/table_S8_univariate_ace_full.tsv        (headline; h^2 +/- SE, c^2 +/- SE, e^2 +/- SE, LRT p)
####             output/tables/table_S8b_univariate_LRT.tsv             (two boundary-corrected single-parameter LRTs per phenotype)
####             output/tables/table_S8c_univariate_model_comparison.tsv (logLik / AIC / BIC for ACE / AE / CE)
####             output/tables/table_S8d_univariate_fixed_effects.tsv    (fixed-effect coefficients per phenotype)
####             output/tables/table_S8e_univariate_variance_components_raw.tsv (raw ASReml varcomp with V-index)
####             output/tables/sessionInfo_univariate.txt                (R/ASReml session info)
####             output/figures/residuals_univariate_ACE.pdf             (residual diagnostic plots, one page per phenotype)
####
#### Related outputs: Figure 1; Table S8 (headline) with supporting
####             Tables S8b - S8e.
####
#### Pitfalls
#### --------
#### * Always confirm the V-indexing for *your* dataset by running
####   `cbind(summary(model)$varcomp, n = seq_along(model$vparameters))`.
####   Adding or removing a random term shifts every subsequent V index.
#### * Models that don't converge usually mean a variance component on
####   the boundary (often c^2 ~= 0). The nested AE, CE and E fits
####   already remove one or both of those terms, so the ACE model can
####   still be summarised even when its sub-parameter is at zero.
#### * `na.action = na.method(x = "omit")` drops rows with any missing
####   value used by the model -- check `nrow(d)` after the call.
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
## Loaded centrally by `01_setup.R`; restated here so the script is
## self-documenting. ASReml-R must be installed and licensed before
## sourcing -- the licence check happens inside `library(asreml)`.

library(asreml)        # commercial; provides `asreml()`, `ainverse()`,
                       # `vpredict()`, `lrt.asreml()`
library(data.table)    # used for the `copy()` + `:=` pattern inside `fit_ace()`


load(file.path(DATA_DIR, "data_prepared.Rdata"))   # provides `baseline`, `baseline_int`, `family_data`

## Primary univariate ACE on the raw (untransformed) scale — Figure 1 and
## Suppl. Table S8. Childhood trauma enters as the continuous CTQ-SF total
## score (same as other stress phenotypes). INT-scale sensitivity: script 04bc.
phenotypes <- intersect(c(STRESS_PHENOTYPES, PERSONALITY_PHENOTYPES),
                        names(baseline))


## Step 0 -- Boundary-corrected LRT p-values ----------------------------
##
## Under the null "sigma^2 = 0" the standard chi^2 reference distribution
## is wrong because zero is on the boundary of the parameter space. The
## corrected reference distribution is a mixture of chi^2 with different
## degrees of freedom. Only single-parameter boundary tests are used here
## (A alone via ACE-vs-CE, C alone via ACE-vs-AE):
##
##     LRT ~ 0.5 * chi^2_0 + 0.5 * chi^2_1  (Self & Liang, 1987).
##
## Because chi^2_0 places all mass at zero (LRT = 0 has p = 1 by
## construction), only the non-degenerate tail matters for the p-value.
## The joint two-parameter (A + C) test against the E-only model has been
## removed (see header): it relies on an E-only REML fit that is not
## comparable across random-effects structures and misbehaves for
## low-variance phenotypes such as Stressful Life Events.

p_bd_1df <- function(LR) 0.5 * (1 - pchisq(LR, df = 1))                                        # single boundary parameter


## Step 1 -- Univariate ACE (+ AE / CE / E) for one phenotype ----------
##
## `fit_ace()` fits four nested models -- the full ACE and the three
## submodels needed for the LRTs and AIC/BIC comparison -- and returns
## a *list* of tidy data.frames with everything downstream tables need.
##
## Returned structure
##   $main : one-row data.frame with h^2, c^2, e^2 and their SEs plus
##           the three LRT p-values (Suppl. Table S8).
##   $lrt  : one-row data.frame with the LRT statistic and boundary-
##           corrected p-value for each of the three tests (Suppl.
##           Table S8b).
##   $ic   : one-row data.frame with logLik / AIC / BIC and N per model
##           (Suppl. Table S8c).
##   $fx   : one-row-per-fixed-effect data.frame with the ACE model's
##           `summary(...)$coef.fixed` (Suppl. Table S8d).
##   $vc   : one-row-per-variance-component data.frame with the raw
##           ASReml varcomp table (Suppl. Table S8e) with the V-index
##           attached so the reader can match the vpredict formulas.

fit_ace <- function(dt, phenotype) {

  d <- copy(dt)                            # `copy()` because `data.table` modifies in place; we don't want to touch the caller's `baseline`
  d[, PHENO := get(phenotype)]             # copy the requested column into a generic name -> the asreml formula below works for any phenotype
  d <- d[!is.na(PHENO)]                    # complete-case for the response; fixed-effect NAs are dropped by `na.action` below
  N_used <- nrow(d)                        # analytic sample size for this phenotype

  ## Pedigree -> inverse numerator relationship matrix.
  ## This is the "A^{-1}" that ASReml needs for the `vm()` term.
  ## Use full pedigree (family_data) so ainverse() has all ancestors present —
  ## not just those with complete phenotype data. Subsetting to `d` would
  ## drop parents who lack phenotype measurements, causing chkPed to fail.
  ped  <- unique(family_data[, .(PSEUDOIDEXT, FatherID, MotherID)])
  ## ainverse() requires all three columns to be the same type (character).
  ## Converting from factor avoids "programming error in chkPed".
  ped[, PSEUDOIDEXT := as.character(PSEUDOIDEXT)]
  ped[, FatherID    := as.character(FatherID)]
  ped[, MotherID    := as.character(MotherID)]
  ## Ensure "0" or "" are treated as NA (unknown parents).
  ped[FatherID %in% c("0", "NA", ""), FatherID := NA_character_]
  ped[MotherID %in% c("0", "NA", ""), MotherID := NA_character_]
  Ainv <- ainverse(ped)                                    # sparse inverse of the additive relationship matrix

  ## Expose HouseID1 as a local object so ASReml's internal boundary-restart
  ## re-evaluation (eval in asremlEnv) can find it. Convert to data.frame to
  ## avoid data.table-specific eval() scoping issues.
  HouseID1 <- d$HouseID1
  .d <- as.data.frame(d)

  ## Full ACE model.
  mACE <- asreml(
    fixed     = PHENO ~ GENDER + AGE_c + AGE_c2,             # sex + linear/quadratic age; rationale in OPEN_SCIENCE_STATEMENT.md §7
    random    = ~ vm(PSEUDOIDEXT, Ainv) + idv(HouseID1),    # `vm(..., Ainv)` = additive genetic; `idv(HouseID1)` = shared household
    data      = .d,
    na.action = na.method(x = "omit"),                       # rows with any NA in the model columns are dropped here too
    workspace = "2gb",                                       # raise if asreml errors with "insufficient workspace"
    ai.sing   = TRUE,                                        # tolerate singularities in the AI matrix (common on huge pedigrees)
    maxiter   = 50                                           # REML iterations; the model usually converges in <10
  )

  ## AE  -- drops the shared-environment term (test C).
  mAE <- asreml(
    fixed     = PHENO ~ GENDER + AGE_c + AGE_c2,
    random    = ~ vm(PSEUDOIDEXT, Ainv),
    data      = .d,
    na.action = na.method(x = "omit"),
    workspace = "2gb",
    ai.sing   = TRUE,
    maxiter   = 50)

  ## CE  -- drops the genetic term (test A).
  mCE <- asreml(
    fixed     = PHENO ~ GENDER + AGE_c + AGE_c2,
    random    = ~ idv(HouseID1),
    data      = .d,
    na.action = na.method(x = "omit"),
    workspace = "2gb",
    ai.sing   = TRUE,
    maxiter   = 50)

  ## Residual diagnostic plot (one page per phenotype in the shared PDF).
  try(plot(mACE), silent = TRUE)

  ## Boundary-corrected single-parameter LRTs.
  ##   LR_A  = 2*(LL_ACE - LL_CE)  ~ test H0: sigma_a^2 = 0     (1df boundary)
  ##   LR_C  = 2*(LL_ACE - LL_AE)  ~ test H0: sigma_c^2 = 0     (1df boundary)
  ## The joint A+C test against an E-only model is intentionally not
  ## computed (see header): the E-only REML fit is not comparable across
  ## random-effects structures and returns an invalid log-likelihood for
  ## low-variance phenotypes such as Stressful Life Events.
  ll <- function(model) as.numeric(summary(model)$loglik)
  LR_A  <- 2 * (ll(mACE) - ll(mCE))
  LR_C  <- 2 * (ll(mACE) - ll(mAE))

  ## Variance components (see header for the V1/V2/V3 mapping).
  ##   V1 = sigma_c^2 (HouseID1)
  ##   V2 = sigma_a^2 (vm(PSEUDOIDEXT, Ainv))
  ##   V3 = sigma_e^2 (units!R, residual)
  h2 <- vp_row(vpredict(mACE, h2 ~ V2 / (V1 + V2 + V3)))     # delta-method SE of h^2 via `vpredict`
  c2 <- vp_row(vpredict(mACE, c2 ~ V1 / (V1 + V2 + V3)))
  e2 <- vp_row(vpredict(mACE, e2 ~ V3 / (V1 + V2 + V3)))

  ## Information-criterion table across all four models.
  gIC <- function(s, k) {
    ll_val <- as.numeric(s$loglik)
    aic    <- if ("aic" %in% names(s)) as.numeric(s$aic) else -2 * ll_val + 2 * k
    bic    <- if ("bic" %in% names(s)) as.numeric(s$bic) else -2 * ll_val + log(N_used) * k
    c(logLik = ll_val, AIC = aic, BIC = bic)
  }
  ic_ACE <- gIC(summary(mACE), k = 3)   # A, C, E variance components
  ic_AE  <- gIC(summary(mAE),  k = 2)   # A, E
  ic_CE  <- gIC(summary(mCE),  k = 2)   # C, E

  ## Fixed-effect coefficients from the ACE fit.
  fx <- NULL
  cf <- try(summary(mACE)$coef.fixed, silent = TRUE)
  if (!inherits(cf, "try-error") && !is.null(cf) && length(cf)) {
    fx <- data.frame(phenotype = phenotype,
                     term      = rownames(cf),
                     as.data.frame(cf),
                     row.names = NULL,
                     stringsAsFactors = FALSE)
  }

  ## Raw ASReml varcomp table (kept verbatim + phenotype + V-index).
  vc_df <- as.data.frame(summary(mACE)$varcomp)
  vc <- data.frame(phenotype      = phenotype,
                   component_name = rownames(vc_df),
                   V_index        = paste0("V", seq_len(nrow(vc_df))),
                   vc_df,
                   row.names = NULL,
                   stringsAsFactors = FALSE)

  ## 95 % CIs via normal approximation on the delta-method SE.
  ## These are approximate (the true posterior is not symmetric for h^2
  ## near 0 or 1), but they are the standard report for REML-estimated
  ## variance components in pedigree analyses.
  ci95 <- function(est, se) list(lo = est - 1.96 * se, hi = est + 1.96 * se)
  h2_ci <- ci95(h2$estimate, h2$se)
  c2_ci <- ci95(c2$estimate, c2$se)
  e2_ci <- ci95(e2$estimate, e2$se)

  ## Assemble the five tidy tables.
  list(
    main = data.frame(Phenotype = phenotype, N = N_used,
      h2          = h2$estimate, h2_SE = h2$se,
      h2_CI_lo    = h2_ci$lo,   h2_CI_hi = h2_ci$hi,
      c2          = c2$estimate, c2_SE = c2$se,
      c2_CI_lo    = c2_ci$lo,   c2_CI_hi = c2_ci$hi,
      e2          = e2$estimate, e2_SE = e2$se,
      e2_CI_lo    = e2_ci$lo,   e2_CI_hi = e2_ci$hi,
      p_A         = p_bd_1df(LR_A),
      p_C         = p_bd_1df(LR_C),
      stringsAsFactors = FALSE),
    lrt = data.frame(Phenotype = phenotype,
      LR_ACE_vs_CE = LR_A,  p_ACE_vs_CE = p_bd_1df(LR_A),
      LR_ACE_vs_AE = LR_C,  p_ACE_vs_AE = p_bd_1df(LR_C),
      stringsAsFactors = FALSE),
    ic  = data.frame(Phenotype  = phenotype, N = N_used,
      logLik_ACE = ic_ACE["logLik"], AIC_ACE = ic_ACE["AIC"], BIC_ACE = ic_ACE["BIC"],
      logLik_AE  = ic_AE["logLik"],  AIC_AE  = ic_AE["AIC"],  BIC_AE  = ic_AE["BIC"],
      logLik_CE  = ic_CE["logLik"],  AIC_CE  = ic_CE["AIC"],  BIC_CE  = ic_CE["BIC"],
      row.names = NULL, stringsAsFactors = FALSE),
    fx = fx, vc = vc)
}


## Step 2 -- Loop over phenotypes ---------------------------------------
##
## One ACE + AE + CE + E fit per phenotype. Residual diagnostic plots
## are stacked in one PDF so any single-page anomaly is easy to eyeball
## before trusting the corresponding h^2.

pdf(file.path(OUT_FIGURES, "residuals_univariate_ACE.pdf"), width = 8, height = 6)

results <- lapply(phenotypes, function(p) {
  message("Univariate ACE: ", p)
  fit_ace(baseline, p)
})

dev.off()

## Bind the five parallel tables together.
main_results <- do.call(rbind, lapply(results, `[[`, "main"))
lrt_results  <- do.call(rbind, lapply(results, `[[`, "lrt"))
ic_results   <- do.call(rbind, lapply(results, `[[`, "ic"))
fx_results   <- do.call(rbind, Filter(Negate(is.null), lapply(results, `[[`, "fx")))
vc_results   <- do.call(rbind, lapply(results, `[[`, "vc"))


## Step 3 -- Multiple-testing flag --------------------------------------
##
## Correction is parallel to Ani et al. (2026): the univariate family is
## the phenotypes x 2 tested variance components (h^2 via ACE-vs-CE, `p_A`;
## c^2 via ACE-vs-AE, `p_C`), i.e. 8 x 2 = 16 tests, so the Bonferroni
## threshold is alpha = .05 / 16 = .003125. (Ani et al. used the same
## phenotypes-x-2-components rule, there 16 x 2 = 32.) With N > 140,000
## every LRT is highly significant regardless of the exact threshold. The
## Bonferroni flags (h2_sig_bonf / c2_sig_bonf) are written only to the
## INTERNAL audit file; the public Table S8 keeps estimates, SEs, CIs and
## LRT p-values without star symbols. Threshold from 01_setup.R: α = .05 / 16.

main_results$h2_sig_bonf    <- main_results$p_A < UNIVARIATE_BONFERRONI_ALPHA
main_results$c2_sig_bonf    <- main_results$p_C < UNIVARIATE_BONFERRONI_ALPHA


## Step 4 -- Save the five tables + sessionInfo -------------------------
##
## Two versions of Table S8:
##   table_S8_univariate_ace_full.tsv           -- public table
##     (estimates, SEs, 95 % CIs, LRT p-values; NO sig flags)
##   table_S8_univariate_ace_full_INTERNAL.tsv  -- internal audit only
##     (same + h2_sig_bonf / c2_sig_bonf)

mscp_cols <- c("Phenotype", "N",
               "h2", "h2_SE", "h2_CI_lo", "h2_CI_hi",
               "c2", "c2_SE", "c2_CI_lo", "c2_CI_hi",
               "e2", "e2_SE", "e2_CI_lo", "e2_CI_hi",
               "p_A", "p_C")
save_table(main_results[, mscp_cols], "table_S8_univariate_ace_full")
save_table(main_results,              "table_S8_univariate_ace_full_INTERNAL")
save_table(lrt_results,  "table_S8b_univariate_LRT")
save_table(ic_results,   "table_S8c_univariate_model_comparison")
if (!is.null(fx_results) && nrow(fx_results) > 0) {
  save_table(fx_results, "table_S8d_univariate_fixed_effects")
}
save_table(vc_results,   "table_S8e_univariate_variance_components_raw")

writeLines(capture.output(sessionInfo()),
           file.path(OUT_TABLES, "sessionInfo_univariate.txt"))

message("Step 03 complete. Univariate ACE + AE + CE fitted for ",
        nrow(main_results), " phenotypes. Tables S8 and S8b - S8e written.")
