#### =========================================================================
#### 05_sensitivity_int.R
#### -------------------------------------------------------------------------
#### Paper     : Reichelt F. et al. -- Familial Correlation between Personality
####             and Stress Exposure: Insights from the Multigenerational
####             Lifelines Cohort Study. European Journal of Personality
####             (manuscript EJP-26-4766, in revision).
#### Author    : Felix Reichelt  <f.reichelt@umcg.nl>
####
#### What this script does
#### ---------------------
#### Repeats the univariate and bivariate ACE analyses (scripts 03 and
#### 04) after rank-based inverse normal transformation of all eight
#### phenotypes (five stress + three personality traits). Reported in
#### Supplementary Section S3 (Tables S10–S10c) as a sensitivity analysis.
####
#### Why bother
#### ----------
#### Several stress phenotypes are right-skewed (e.g. CTQ-SF total
#### counts of severe events; Suppl. Figure S2). Personality traits are
#### approximately normal (Suppl. Fig. S17) but are included here for
#### completeness so Supplementary Table S11 reports INT model
#### comparison for all phenotypes on the same footing.
####
#### How ASReml-R is used here
#### -------------------------
#### Exactly the same V-indexing as in scripts 04 and 04b -- the
#### transformation only changes the values of the phenotype, not the
#### structure of the model. See the headers of those two scripts for
#### the explicit `vpredict` formulas (ASReml-R Reference Manual 4.2,
#### Section 4.6). In the bivariate the order is shared env (V1-V3),
#### genetic (V4-V6), residual scale placeholder (V7), residual
#### covariances (V8-V10).
####
#### Inputs    : data/data_prepared.Rdata (from script 02).
#### Outputs   : output/tables/table_S10_int_univariate_ace.tsv         (INT: h^2 +/- SE for all 8 phenotypes)
####             output/tables/table_S10b_int_univariate_LRT.tsv         (boundary-corrected LRTs, all 8)
####             output/tables/table_S10c_int_univariate_model_comparison.tsv (logLik / AIC / BIC, all 8)
####             output/tables/table_S10d_int_bivariate_stress_stress.tsv (bivariate INT ACE, stress pairs only)
####             output/tables/table_S10e_int_univariate_variance_components_raw.tsv
####             output/tables/sessionInfo_int.txt
####
#### Manuscript: Suppl. Section S3, Table S10 (INT univariate) +
####             Table S10b/c/d (LRT, model comparison, INT bivariate).
####
#### Pitfalls
#### --------
#### * INT eliminates outliers by construction. If a "real" outlier was
####   driving an effect on the raw scale, INT will hide it. Always
####   report both scales side by side.
#### * Random tie-breaking is controlled by `set.seed(12345)` in
####   01_setup.R; do not skip the setup file or results will differ
####   slightly run-to-run.
####
#### License   : ASReml-R is commercial. https://vsni.co.uk/software/asreml-r/
#### Version history:
#### 2026-07-15  v1.1  Ani-aligned multiple testing (/16 uni, /75 bivar);
####                   per-component 1-df bivariate LRTs; INTERNAL tables
#### 2026-07-15  v1.0  Initial release
#### =========================================================================


## Required libraries -----------------------------------------------------
##
## Loaded centrally by `01_setup.R`. Source `04_univariate_heritability.R`
## and `04b_bivariate_stress_stress.R` first in the same R session to
## bring `fit_ace()` and `fit_biv()` into the global environment --
## this script only changes which columns are passed to those helpers.

library(asreml)        # commercial; used by `fit_ace()` / `fit_biv()`
library(data.table)    # `data.table` syntax used by the helpers


load(file.path(DATA_DIR, "data_prepared.Rdata"))   # provides `baseline_int` (already pre-INT-transformed by Lifelines)


## Step 1 — Apply inverse normal transformation (stress + personality)
##
## `inverse_normal()` is defined in `01_helpers.R`. Apply INT to RAW scores
## from `baseline` (not from `baseline_int`, which would double-INT stress).

for (p in c(STRESS_PHENOTYPES, PERSONALITY_PHENOTYPES)) {
  if (p %in% names(baseline)) {
    baseline_int[[paste0(p, "_INT")]] <- inverse_normal(baseline[[p]])
  }
}
int_phenotypes <- c(paste0(STRESS_PHENOTYPES, "_INT"),
                    paste0(PERSONALITY_PHENOTYPES, "_INT"))


## Step 2 — Reuse the fitting functions from scripts 04 and 04b --------
##
## Run `04_univariate_heritability.R` and `04b_bivariate_stress_stress.R`
## first in the same R session so `fit_ace()` and `fit_biv()` are in
## your global environment.


## Step 3 -- Univariate ACE (+ AE / CE) on INT phenotypes ------------
##
## The enriched `fit_ace()` from script 04 returns a list of five tidy
## tables per phenotype (main, lrt, ic, fx, vc). We collect them
## separately so the INT sensitivity produces the same suite of
## supplementary tables as the baseline analysis in script 04.

pdf(file.path(OUT_FIGURES, "residuals_univariate_ACE_INT.pdf"),
    width = 8, height = 6)

uni_int_results <- lapply(int_phenotypes, function(p) {
  message("Univariate ACE (INT): ", p)
  fit_ace(baseline_int, p)
})

dev.off()

uni_int_main <- do.call(rbind, lapply(uni_int_results, `[[`, "main"))
uni_int_lrt  <- do.call(rbind, lapply(uni_int_results, `[[`, "lrt"))
uni_int_ic   <- do.call(rbind, lapply(uni_int_results, `[[`, "ic"))
uni_int_vc   <- do.call(rbind, lapply(uni_int_results, `[[`, "vc"))

save_table(uni_int_main, "table_S10_int_univariate_ace")
save_table(uni_int_lrt,  "table_S10b_int_univariate_LRT")
save_table(uni_int_ic,   "table_S10c_int_univariate_model_comparison")
save_table(uni_int_vc,   "table_S10e_int_univariate_variance_components_raw")


## Step 4 -- Bivariate ACE on INT phenotypes --------------------------
##
## `fit_biv()` still returns one tidy row per pair (rp, rg, rc, re with
## SE + A_pcor / C_pcor / E_pcor + prop_A / C / E). See the header of
## `04b_bivariate_stress_stress.R` for the full V-index mapping.

## Step 4 -- Bivariate ACE on INT stress phenotypes only ----------------
##
## Bivariate INT restricted to stress x stress (10 pairs); personality
## bivariate INT is not part of the manuscript sensitivity section.

int_stress_phenotypes <- paste0(STRESS_PHENOTYPES, "_INT")
pairs <- as.data.frame(t(combn(int_stress_phenotypes, 2)), stringsAsFactors = FALSE)
names(pairs) <- c("Var1", "Var2")

biv_int <- lapply(seq_len(nrow(pairs)), function(i) {
  message("Bivariate ACE (INT): ", pairs$Var1[i], " x ", pairs$Var2[i])
  fit_biv(baseline_int, pairs$Var1[i], pairs$Var2[i])
})
biv_int <- do.call(rbind, biv_int)
save_table(biv_int, "table_S10d_int_bivariate_stress_stress")


## Step 5 -- sessionInfo ---------------------------------------------

writeLines(capture.output(sessionInfo()),
           file.path(OUT_TABLES, "sessionInfo_int.txt"))

message("Step 05 complete. INT tables (S10/S10b/c: ",
        length(int_phenotypes), " phenotypes; S10d: stress bivariate) written.")
