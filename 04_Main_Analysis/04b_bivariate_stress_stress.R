#### =========================================================================
#### 04b_bivariate_stress_stress.R
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
#### Fits a bivariate ACE linear mixed model for every pair of stress phenotypes
#### (10 pairs in total) and extracts the genetic, shared-environmental,
#### unique-environmental and phenotypic correlations between the two
#### traits. These correlations feed **Figure 2** (lower panels) of the
#### manuscript and **Suppl. Table S13**.
####
#### The script also computes the A/C/E **decomposition** of the
#### phenotypic correlation (upper panels of Figure 2; **Suppl. Table
#### S15**):
####
####     A_pcor = sigma_a1a2 / sqrt(Var_P1 * Var_P2)
####     C_pcor = sigma_c1c2 / sqrt(Var_P1 * Var_P2)
####     E_pcor = sigma_e1e2 / sqrt(Var_P1 * Var_P2)
####     prop_A = A_pcor / rp;   prop_C = C_pcor / rp;   prop_E = E_pcor / rp
####
#### How ASReml-R is used here (Section 4.6 of the Reference Manual 4.2)
#### -------------------------------------------------------------------
#### The two-trait REML linear mixed model stacks both phenotypes in
#### long format. The random part has an unstructured 2x2 covariance
#### for the additive genetic and for the shared-environmental factor
#### (see Section 4.6 of the ASReml-R Reference Manual 4.2):
####
####     random = ~ us(trait):vm(PSEUDOIDEXT, Ainv) +     # genetic
####                 us(trait):HouseID1                    # shared env
####     residual is the ASReml default us(trait) on units
####
#### KEY POINT -- the order in which the components appear in the
#### `vparameters` slot (and therefore the V1, V2, ... names used by
#### `vpredict()`) is **not the order of the random formula**. ASReml
#### reorders the random terms alphabetically: in our model "HouseID1"
#### precedes "vm(PSEUDOIDEXT, ...)" because H < v. The actual layout
#### that produced the manuscript estimates is:
####
####     V1   trait:HouseID1!t1:t1            sigma_c1^2     (shared env, trait 1)
####     V2   trait:HouseID1!t2:t1            sigma_c1c2     (shared env, covariance)
####     V3   trait:HouseID1!t2:t2            sigma_c2^2     (shared env, trait 2)
####     V4   trait:vm(...)!t1:t1             sigma_a1^2     (genetic, trait 1)
####     V5   trait:vm(...)!t2:t1             sigma_a1a2     (genetic covariance)
####     V6   trait:vm(...)!t2:t2             sigma_a2^2     (genetic, trait 2)
####     V7   units:trait!R                   1   (FIXED)    (residual scale; placeholder)
####     V8   units:trait!t1:t1               sigma_e1^2     (residual, trait 1)
####     V9   units:trait!t2:t1               sigma_e1e2     (residual covariance)
####     V10  units:trait!t2:t2               sigma_e2^2     (residual, trait 2)
####
#### Note that V7 is the residual scale parameter (`!R`) which ASReml
#### fixes at 1 (Constraint = F). It still occupies a slot in
#### `vparameters`, so the residual covariances live at V8/V9/V10.
####
#### The four correlations of interest (manual 4.6):
####
####     rg = V5 / sqrt(V4 * V6)                                # genetic
####     rc = V2 / sqrt(V1 * V3)                                # shared env
####     re = V9 / sqrt(V8 * V10)                               # unique env
####     rp = (V2 + V5 + V9) /
####          sqrt( (V1 + V4 + V8) * (V3 + V6 + V10) )          # phenotypic
####
#### Significance testing -- per-component 1-df LRTs (Ani et al. design)
#### --------------------------------------------------------------------
#### Together with `04c_bivariate_personality_stress.R` (15 personality
#### x stress pairs), this script contributes 10 of the **25 bivariate
#### pairings** reported in the manuscript.
####
#### Each of the three cross-trait correlations is tested SEPARATELY,
#### exactly as in Ani, Nolte, Vaez & Snieder (2026, Diabetes Obes Metab
#### 28:3344; the group's standard ASReml bivariate protocol). For a given
#### correlation the full unstructured model is compared with a sub-model
#### in which ONLY that one cross-trait covariance is fixed to zero (its
#### random term switched from `us(trait)` to `diag(trait)`), leaving the
#### other two covariances free:
####
####   r_g test : diag(trait):vm(...)  + us(trait):HouseID1  + us(trait):units
####   r_c test : us(trait):vm(...)    + diag(trait):HouseID1 + us(trait):units
####   r_e test : us(trait):vm(...)    + us(trait):HouseID1  + us(trait):units
####              with the residual cross-trait covariance fixed at 0 via R.param
####              (diag(trait):units does NOT constrain E cov in wide cbind models)
####
#### Each LR = 2*(logLik_full - logLik_sub) is referred to a chi^2 with
#### **1 df** (columns LR_rg/p_rg, LR_rc/p_rc, LR_re/p_re). A covariance is
#### a sign-unrestricted parameter, so zero is an INTERIOR point of the
#### parameter space and NO boundary (Self & Liang) correction is required
#### -- an ordinary chi^2(1) applies. (This is the opposite of the
#### univariate variance-component tests h2 = 0 and c2 = 0 in script 04,
#### where the parameter is on the boundary and the 1/2 chi^2_0 + 1/2 chi^2_1
#### mixture is used.)
####
#### Multiple testing: Bonferroni across all component tests, i.e.
#### 25 pairings x 3 correlations = 75 tests, alpha = .05 / 75 = 6.67e-4
#### (flags rg_sig_bonf, rc_sig_bonf, re_sig_bonf). This matches Ani et al.,
#### who corrected over (pairs x 3 correlations).
####
#### Because N > 140,000, statistical significance is expected even for
#### tiny effects, so the substantive interpretation rests on effect sizes,
#### confidence intervals, and the A/C/E decomposition of r_P rather than
#### on p-values alone. A Benjamini-Hochberg FDR sensitivity analysis is
#### still produced in the INTERNAL audit file (script 04bd).
####
#### Wald p-value caveat
#### -------------------
#### The formal test for each correlation is the 1-df LRT above. The
#### Wald-style normal-approximation p-values (2 * pnorm(-|est / se|)) are
#### retained in the INTERNAL audit file only; they can be misleading near
#### the boundary of the correlation space (r_g approx +-1) or when an
#### underlying variance component is near zero (c^2 ~ 0 gives a degenerate
#### r_c), whereas the LRT does not depend on that approximation. r_P has no
#### sub-model (it is a deterministic function of all six variances and
#### three covariances) and is reported descriptively with its 95 % CI.
####
#### Outputs   : output/tables/table_S13_bivariate_stress_stress.tsv
####             output/tables/table_S13c_stress_stress_model_comparison.tsv
####             output/tables/table_S15_decomp_stress_stress.tsv
####             output/tables/table_S13S15_bivariate_stress_stress_full_INTERNAL.tsv
####
#### Manuscript: Methods 2.3, Results 3.3, Figure 2,
####             Suppl. Tables S13 (estimates + LRTs), S13c (model comparison),
####             S15 (decomposition). LRT significance counts: script 04bd.
####
#### Pitfalls
#### --------
#### * ALWAYS re-derive the V-mapping for *your* model. The most common
####   bug in bivariate ASReml work is trusting the random formula order
####   instead of the actual `vparameters` order. Run
####       cbind(summary(model)$varcomp, n = seq_along(model$vparameters))
####   once before you trust the formulas above. If your component
####   names use a different alphabetical order (e.g. you renamed
####   HouseID to FamilyID -- F < v, but also F < H, so the order
####   shifts), update the formulas accordingly.
####
#### * `na.action = na.method(x = "include", y = "include")` keeps
####   participants with one trait observed and the other missing, which
####   is essential for an efficient bivariate analysis. Switching to
####   `na.omit` here would throw away thousands of relatives.
####
#### License   : ASReml-R is commercial. https://vsni.co.uk/software/asreml-r/
#### Version history:
#### 2026-07-15  v1.1  Ani-aligned multiple testing (/16 uni, /75 bivar);
####                   per-component 1-df bivariate LRTs; INTERNAL tables
#### 2026-07-15  v1.0  Initial release
#### =========================================================================


## Required libraries -----------------------------------------------------
##
## Loaded centrally by `01_setup.R`; restated here for self-documentation.

library(asreml)        # commercial; provides `asreml()`, `ainverse()`, `vpredict()`
library(data.table)    # `data.table()` constructor, `rbindlist()`, in-place `:=`


load(file.path(DATA_DIR, "data_prepared.Rdata"))   # provides `baseline`, `baseline_int`, `family_data`

## Primary bivariate ACE on the raw scale (continuous CTQ-SF total for CT).
stress <- intersect(STRESS_PHENOTYPES, names(baseline))
pairs  <- as.data.frame(t(combn(stress, 2)), stringsAsFactors = FALSE)
names(pairs) <- c("Var1", "Var2")


## Step 1 -- Helper that fits the bivariate ACE for one pair ------------

fit_biv <- function(dt, p1, p2) {

  d <- copy(dt)                                              # data.table is mutable: copy before in-place edits
  d[, P1 := get(p1)]                                         # generic alias for trait 1
  d[, P2 := get(p2)]                                         # generic alias for trait 2
  ## Keep individuals with at least one of the two traits (see header).
  d <- d[!(is.na(P1) & is.na(P2))]                           # drop rows missing BOTH; keep rows missing only one (bivariate FIML works)

  ## Wide-format bivariate model (one row per person, two response columns).
  ## ASReml's cbind() notation generates the 'trait' factor internally and
  ## handles individuals missing one trait via na.action = include.
  ## This avoids the long-format residual specification problem: with long
  ## format, ~us(trait):units requires units = individuals (N) but gets rows
  ## (2N), causing "residual model implies 254992" errors.
  ##
  ## With wide format, units = individuals = N rows, so ~us(trait):units
  ## correctly models the unstructured residual covariance between the two
  ## traits (V7 = scale fixed at 1; V8 = e1 var; V9 = e1e2 cov; V10 = e2 var).
  ## This matches the V1..V10 mapping documented in the script header.
  ##
  ## Reference: 2025 HPC analysis (bivar_loneliness_childhoodtrauma.R) which
  ## used this same cbind() + rcov = ~us(trait):units specification.

  ## Use full pedigree (family_data) so ainverse() has all ancestors present.
  ped <- unique(family_data[, .(PSEUDOIDEXT, FatherID, MotherID)])
  ped[, PSEUDOIDEXT := as.character(PSEUDOIDEXT)]
  ped[, FatherID    := as.character(FatherID)]
  ped[, MotherID    := as.character(MotherID)]
  ped[FatherID %in% c("0", "NA", ""), FatherID := NA_character_]
  ped[MotherID %in% c("0", "NA", ""), MotherID := NA_character_]
  Ainv <- ainverse(ped)

  ## Expose variables locally for ASReml's boundary-restart re-evaluation;
  ## convert to data.frame to avoid data.table eval() scoping issues.
  HouseID1    <- d$HouseID1
  PSEUDOIDEXT <- d$PSEUDOIDEXT
  .d <- as.data.frame(d)

  m <- asreml(
    fixed     = cbind(P1, P2) ~ trait + trait:GENDER + trait:AGE_c + trait:AGE_c2,
    random    = ~ us(trait):vm(PSEUDOIDEXT, Ainv) +                          # unstructured 2x2 additive-genetic covariance
                   us(trait):HouseID1,                                       # unstructured 2x2 shared-household covariance
    rcov      = ~ us(trait):units,                                           # unstructured residual covariance (E); wide format: units = individuals
    data      = .d,
    na.action = na.method(x = "include", y = "include"),
    workspace = "3gb",
    ai.sing   = TRUE,
    maxiter   = 50)

  ## Per-component 1-df LRT sub-models ----------------------------------
  ##
  ## Following Ani et al. (2026), each cross-trait correlation is tested by
  ## comparing the full model against a sub-model in which ONLY that one
  ## covariance is fixed to zero (its term switched from us() to diag()),
  ## leaving the other two covariances free. `fit_sub()` builds the three
  ## sub-models from a compact (genetic, household, residual) structure
  ## specification, where each entry is "us" (free covariance) or "diag"
  ## (covariance fixed to zero).
  fit_sub <- function(gen, hh, res) {
    try(asreml(
      fixed     = cbind(P1, P2) ~ trait + trait:GENDER + trait:AGE_c + trait:AGE_c2,
      random    = as.formula(sprintf("~ %s(trait):vm(PSEUDOIDEXT, Ainv) + %s(trait):HouseID1",
                                     gen, hh)),
      rcov      = as.formula(sprintf("~ %s(trait):units", res)),
      data      = .d,
      na.action = na.method(x = "include", y = "include"),
      workspace = "3gb",
      ai.sing   = TRUE,
      maxiter   = 50
    ), silent = TRUE)
  }
  m_noG <- fit_sub("diag", "us",   "us")     # genetic covariance = 0    -> tests r_g
  m_noC <- fit_sub("us",   "diag", "us")     # shared-household cov = 0  -> tests r_c

  ## Residual (E) covariance sub-model. diag(trait):units leaves logLik
  ## identical to the full model in wide cbind() format (verified Jul 2026:
  ## LR_re = 0, p_re = 1 for all pairs). Fix units:trait!t2:t1 at 0 via
  ## R.param on the fitted full model; fall back to dsum(~ units | trait).
  fit_residual_zero_cov <- function(m_full) {
    vc <- summary(m_full)$varcomp
    rn <- rownames(vc)
    cov_label <- grep("units:trait.*!t2:t1", rn, value = TRUE)
    if (!length(cov_label))
      cov_label <- grep("!t2:t1", rn, value = TRUE)
    if (!length(cov_label))
      stop("residual cross-trait covariance not found in varcomp")
    cov_label <- cov_label[1]

    patch_R_param <- function(Rdf) {
      if (!is.data.frame(Rdf) || !"Component" %in% names(Rdf))
        return(NULL)
      hit <- which(Rdf$Component == cov_label)
      if (length(hit) != 1)
        hit <- grep("!t2:t1", Rdf$Component)
      if (length(hit) != 1)
        return(NULL)
      Rdf <- Rdf
      Rdf$Value[hit] <- 0
      Rdf$Constraint[hit] <- "F"
      Rdf
    }

    ll_full <- as.numeric(summary(m_full)$loglik)
    R2 <- patch_R_param(m_full$R.param)
    if (!is.null(R2)) {
      m_fix <- update(m_full, R.param = R2, maxiter = 50, workspace = "3gb",
                      ai.sing = TRUE)
      ll_fix <- as.numeric(summary(m_fix)$loglik)
      if (abs(ll_full - ll_fix) > 1e-6)
        return(m_fix)
      message("fit_biv: R.param E-cov fix left logLik unchanged; trying dsum residual.")
    }

    update(m_full, residual. = ~ dsum(~ units | trait), maxiter = 50,
           workspace = "3gb", ai.sing = TRUE)
  }
  m_noE <- try(fit_residual_zero_cov(m), silent = TRUE)   # residual cov = 0 -> tests r_e

  ## Free (co)variance parameters for AIC/BIC. Each sub-model has exactly
  ## one fewer free covariance than the full model by construction, so
  ## npar_sub = npar_full - 1 (df of each 1-df LRT), independent of how
  ## ASReml reparameterises the residual scale in the diagonal fit.
  count_free_vc <- function(mod) {
    vc <- summary(mod)$varcomp
    if (!is.null(vc$bound)) sum(vc$bound != "F") else nrow(vc)
  }
  N_biv     <- nrow(d)
  ll_full   <- as.numeric(summary(m)$loglik)
  npar_full <- count_free_vc(m)
  AIC_full  <- -2 * ll_full + 2 * npar_full
  BIC_full  <- -2 * ll_full + log(N_biv) * npar_full

  ## One-parameter LRT + information criteria for a single sub-model.
  ## The tested covariance is sign-unrestricted (interior point), so the
  ## reference distribution is an ordinary chi^2(1); no boundary mixture.
  lrt_component <- function(msub) {
    if (inherits(msub, "try-error"))
      return(list(ll = NA_real_, LR = NA_real_, p = NA_real_,
                  AIC = NA_real_, BIC = NA_real_))
    ll_s   <- as.numeric(summary(msub)$loglik)
    LR     <- 2 * (ll_full - ll_s)
    if (!is.na(LR) && LR < 0) LR <- 0          # guard against convergence noise
    npar_s <- npar_full - 1L
    list(ll  = ll_s, LR = LR,
         p   = pchisq(LR, df = 1, lower.tail = FALSE),
         AIC = -2 * ll_s + 2 * npar_s,
         BIC = -2 * ll_s + log(N_biv) * npar_s)
  }
  lr_g <- lrt_component(m_noG)
  lr_c <- lrt_component(m_noC)
  lr_e <- lrt_component(m_noE)

  ## V-indices -- see header for the full mapping.
  ##   shared env  : V1, V2, V3       (sigma_c1^2, sigma_c1c2, sigma_c2^2)
  ##   genetic     : V4, V5, V6       (sigma_a1^2, sigma_a1a2, sigma_a2^2)
  ##   residual    : V8, V9, V10      (sigma_e1^2, sigma_e1e2, sigma_e2^2; V7 is the !R scale, fixed at 1)
  rc <- vp_row(vpredict(m, rc ~ V2 / sqrt(V1 * V3)))                            # shared-environment correlation
  rg <- vp_row(vpredict(m, rg ~ V5 / sqrt(V4 * V6)))                            # additive-genetic correlation
  re <- vp_row(vpredict(m, re ~ V9 / sqrt(V8 * V10)))                           # unique-environment correlation
  rp <- vp_row(vpredict(m, rp ~ (V2 + V5 + V9) /
                                 sqrt((V1 + V4 + V8) * (V3 + V6 + V10))))       # phenotypic correlation = sum of A, C, E covariances / sqrt(total variances)

  ## A / C / E decomposition of r_P (Suppl. Table S15).
  A_p <- vp_row(vpredict(m, A_p ~ V5 / sqrt((V1 + V4 + V8) *
                                             (V3 + V6 + V10))))
  C_p <- vp_row(vpredict(m, C_p ~ V2 / sqrt((V1 + V4 + V8) *
                                             (V3 + V6 + V10))))
  E_p <- vp_row(vpredict(m, E_p ~ V9 / sqrt((V1 + V4 + V8) *
                                             (V3 + V6 + V10))))

  ## Wald-style two-sided p-values for each correlation. The formal test is
  ## now the per-component 1-df LRT above (p_rg/p_rc/p_re); these Wald values
  ## are retained in the INTERNAL audit file only, as a cross-check.
  wald_p <- function(est, se) 2 * pnorm(-abs(est / se))

  ## 95 % CIs on all four correlations (normal approximation on delta-method SE).
  ci95 <- function(est, se) c(lo = est - 1.96 * se, hi = est + 1.96 * se)
  rp_ci <- ci95(rp$estimate, rp$se)
  rg_ci <- ci95(rg$estimate, rg$se)
  rc_ci <- ci95(rc$estimate, rc$se)
  re_ci <- ci95(re$estimate, re$se)

  ## Flag unstable decompositions: small rp (denominator near 0) or any
  ## proportion outside [−0.20, 1.20] indicate that the decomposition
  ## should not be emphasised in the Results text.
  prop_A_val <- A_p$estimate / rp$estimate
  prop_C_val <- C_p$estimate / rp$estimate
  prop_E_val <- E_p$estimate / rp$estimate
  decomp_unstable <- abs(rp$estimate) < 0.05 ||
                     any(c(prop_A_val, prop_C_val, prop_E_val) < -0.20) ||
                     any(c(prop_A_val, prop_C_val, prop_E_val) >  1.20)

  data.frame(Phenotype1 = p1, Phenotype2 = p2, N = nrow(d),
    logLik_full = ll_full, npar_full = npar_full, AIC_full = AIC_full, BIC_full = BIC_full,
    rp = rp$estimate, rp_SE = rp$se, rp_CI_lo = rp_ci["lo"], rp_CI_hi = rp_ci["hi"],
    rp_p = wald_p(rp$estimate, rp$se),
    ## Additive-genetic correlation + its 1-df LRT (sub-model: sigma_a1a2 = 0)
    rg = rg$estimate, rg_SE = rg$se, rg_CI_lo = rg_ci["lo"], rg_CI_hi = rg_ci["hi"],
    LR_rg = lr_g$LR, df_rg = 1L, p_rg = lr_g$p,
    logLik_noG = lr_g$ll, AIC_noG = lr_g$AIC, BIC_noG = lr_g$BIC,
    ## Shared-household correlation + its 1-df LRT (sub-model: sigma_c1c2 = 0)
    rc = rc$estimate, rc_SE = rc$se, rc_CI_lo = rc_ci["lo"], rc_CI_hi = rc_ci["hi"],
    LR_rc = lr_c$LR, df_rc = 1L, p_rc = lr_c$p,
    logLik_noC = lr_c$ll, AIC_noC = lr_c$AIC, BIC_noC = lr_c$BIC,
    ## Unique-environment correlation + its 1-df LRT (sub-model: sigma_e1e2 = 0)
    re = re$estimate, re_SE = re$se, re_CI_lo = re_ci["lo"], re_CI_hi = re_ci["hi"],
    LR_re = lr_e$LR, df_re = 1L, p_re = lr_e$p,
    logLik_noE = lr_e$ll, AIC_noE = lr_e$AIC, BIC_noE = lr_e$BIC,
    ## Wald normal-approximation p-values (INTERNAL audit only; LRT is primary)
    rg_wald_p = wald_p(rg$estimate, rg$se), rc_wald_p = wald_p(rc$estimate, rc$se), re_wald_p = wald_p(re$estimate, re$se),
    A_pcor = A_p$estimate, A_pcor_SE = A_p$se, C_pcor = C_p$estimate, C_pcor_SE = C_p$se, E_pcor = E_p$estimate, E_pcor_SE = E_p$se,
    prop_A = prop_A_val, prop_C = prop_C_val, prop_E = prop_E_val,
    decomp_unstable = decomp_unstable)
}

## Bonferroni-corrected threshold for the per-component bivariate tests
## (Ani et al. design): 25 pairings x 3 correlations (rg, rc, re) = 75
## tests -> alpha = .05 / 75 = 6.67e-4. The 25 pairings are the 10
## stress-stress pairs in this script + 15 personality-stress pairs in
## script 04bc. Used as the family-wise significance flag in S13 and S14.
BONFERRONI_ALPHA <- BIVARIATE_BONFERRONI_ALPHA


## Step 2 -- Loop over the 10 stress-stress pairs ----------------------

results <- lapply(seq_len(nrow(pairs)), function(i) {
  message("Bivariate (stress x stress): ",
          pairs$Var1[i], " x ", pairs$Var2[i])
  fit_biv(baseline, pairs$Var1[i], pairs$Var2[i])
})
results <- do.call(rbind, results)

## Per-component significance flags (Bonferroni over 75 tests). These are
## the primary family-wise flags for the genetic, shared-household and
## unique-environmental correlations, one 1-df LRT each.
results$rg_sig_bonf <- results$p_rg < BONFERRONI_ALPHA
results$rc_sig_bonf <- results$p_rc < BONFERRONI_ALPHA
results$re_sig_bonf <- results$p_re < BONFERRONI_ALPHA

## Suppl. Table S13 (MANUSCRIPT): correlations, 95 % CIs, and per-component
## 1-df LRTs. Table note (to be added in Word/LaTeX):
##   "The additive-genetic (rg), shared-household (rc) and unique-
##    environmental (re) correlations were each tested with a 1-df
##    likelihood-ratio test comparing the full bivariate ACE model against
##    a sub-model with that single cross-trait covariance fixed to zero
##    (chi^2(1); a covariance is sign-unrestricted, so no boundary
##    correction is applied). To account for multiple testing, p values
##    were Bonferroni-corrected across 25 pairings x 3 correlations = 75
##    tests (alpha = .05/75 = 6.67e-4). rp is the phenotypic correlation,
##    reported descriptively with its 95 % CI. Estimates are bivariate REML
##    with delta-method SEs and normal-approximation 95 % CIs."
save_table(results[, c("Phenotype1", "Phenotype2", "N",
                       "rp", "rp_SE", "rp_CI_lo", "rp_CI_hi",
                       "rg", "rg_SE", "rg_CI_lo", "rg_CI_hi", "LR_rg", "df_rg", "p_rg",
                       "rc", "rc_SE", "rc_CI_lo", "rc_CI_hi", "LR_rc", "df_rc", "p_rc",
                       "re", "re_SE", "re_CI_lo", "re_CI_hi", "LR_re", "df_re", "p_re")],
           "table_S13_bivariate_stress_stress")

## Suppl. Table S13c (MANUSCRIPT): per-component model comparison. Shows the
## full (unstructured ACE) fit and each single-covariance-dropped sub-model
## side by side so every 1-df LR test is fully reproducible.
save_table(results[, c("Phenotype1", "Phenotype2", "N",
                       "logLik_full", "AIC_full", "BIC_full",
                       "logLik_noG", "AIC_noG", "BIC_noG", "LR_rg", "p_rg",
                       "logLik_noC", "AIC_noC", "BIC_noC", "LR_rc", "p_rc",
                       "logLik_noE", "AIC_noE", "BIC_noE", "LR_re", "p_re")],
           "table_S13c_stress_stress_model_comparison")

## Suppl. Table S15 (MANUSCRIPT): decomposition; flag unstable rows.
## Rows with decomp_unstable = TRUE (very small rp or out-of-range
## proportions) should be footnoted rather than described in the text.
save_table(results[, c("Phenotype1", "Phenotype2", "N", "rp",
                       "A_pcor", "A_pcor_SE",
                       "C_pcor", "C_pcor_SE",
                       "E_pcor", "E_pcor_SE",
                       "prop_A", "prop_C", "prop_E",
                       "decomp_unstable")],
           "table_S15_decomp_stress_stress")

## INTERNAL audit file: every column including sig flags; NOT for upload.
save_table(results, "table_S13S15_bivariate_stress_stress_full_INTERNAL")


## Step 3 -- sessionInfo ------------------------------------------------

writeLines(capture.output(sessionInfo()),
           file.path(OUT_TABLES, "sessionInfo_bivariate_stress_stress.txt"))

message("Step 04 complete. ", nrow(results),
        " bivariate ACE models fitted; Tables S13, S13c, S15 (+ INTERNAL) written.")
