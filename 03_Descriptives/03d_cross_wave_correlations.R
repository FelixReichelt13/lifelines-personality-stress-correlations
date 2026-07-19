#### =========================================================================
#### 03d_cross_wave_correlations.R
#### -------------------------------------------------------------------------
#### Paper     : Reichelt F. et al. -- Familial Correlation between Personality
####             and Stress Exposure: Insights from the Multigenerational
####             Lifelines Cohort Study. European Journal of Personality
####             (manuscript EJP-26-4766, in revision).
#### Author    : Felix Reichelt  <f.reichelt@umcg.nl>
####
#### What this script does
#### ---------------------
#### Computes the cross-wave **Pearson** and **Spearman** correlations
#### for LTE (BRUGHA1-12 sum, 12 items per wave) and LDI (GLLM1-9
#### sum, 9 items per wave -- items 10-12 are not in 1A/1B/1C/2A)
#### across the four
#### Lifelines waves (1A baseline, 1B, 1C, 2A second assessment). The
#### numbers underpin the methodological argument for the latent-factor
#### approach (Methods §2.2, paragraph "We collapsed the four waves
#### into a single latent factor because the sum scores were strongly
#### correlated across waves [LDI Spearman 0.54-0.62] yet ...").
####
#### The script also writes the *factor loadings* of the LTE and LDI
#### one-factor `lavaan` models that produced the latent scores used
#### in scripts 04-05, so the supplement can quote them directly.
####
#### Inputs    : data/data_prepared.Rdata (from script 02) -- needed
####             for the cross-wave correlation. The lavaan loadings
####             at the bottom of the script are EMBEDDED (taken from
####             the latest lavaan run) so the supplementary table renders
####             without re-running the model.
####
#### Outputs   : output/tables/table_S5_lte_ldi_cross_wave_correlations.tsv
####             output/tables/table_S6_lte_ldi_factor_loadings.tsv
####
#### Related outputs: Tables S5–S6 (cross-wave correlations / loadings).
#### Version history:
#### 2026-07-15  v1.1  Ani-aligned multiple testing (/16 uni, /75 bivar);
####                   per-component 1-df bivariate LRTs; INTERNAL tables
#### 2026-07-15  v1.0  Initial release
#### =========================================================================


## Required libraries -----------------------------------------------------
##
## Loaded centrally by `01_setup.R`. This script does not call `lavaan`
## directly: the factor loadings in Step 2 are EMBEDDED from a prior
## `lavaan::sem()` run in `02_data_preparation.R` so the supplementary
## table renders without re-fitting the model.

library(data.table)   # for the `[, ..cols]` indexing pattern used in `wave_cor()`


## Step 1 -- Cross-wave correlations (Pearson + Spearman) -------------------
##
## The sum scores per wave (computed in `02a_score_construction.R`)
## are: LTE_Baseline, LTE_Followup_1B, LTE_Followup_1C,
## LTE_Second_Assess, LDI_Baseline, LDI_Followup_1B, LDI_Followup_1C,
## LDI_Second_Assess.

load(file.path(DATA_DIR, "data_prepared.Rdata"))

lte_waves <- c("LTE_Baseline", "LTE_Followup_1B",
                "LTE_Followup_1C", "LTE_Second_Assess")
ldi_waves <- c("LDI_Baseline", "LDI_Followup_1B",
                "LDI_Followup_1C", "LDI_Second_Assess")

wave_cor <- function(dt, cols, method) {
  m <- as.matrix(dt[, ..cols])
  cor(m, use = "pairwise.complete.obs", method = method)
}

cross_wave <- rbind(
  data.frame(inventory = "LTE", method = "Pearson",
              wave_a = combn(lte_waves, 2)[1, ],
              wave_b = combn(lte_waves, 2)[2, ],
              r = wave_cor(baseline, lte_waves, "pearson")[
                lower.tri(wave_cor(baseline, lte_waves, "pearson"))]),
  data.frame(inventory = "LTE", method = "Spearman",
              wave_a = combn(lte_waves, 2)[1, ],
              wave_b = combn(lte_waves, 2)[2, ],
              r = wave_cor(baseline, lte_waves, "spearman")[
                lower.tri(wave_cor(baseline, lte_waves, "spearman"))]),
  data.frame(inventory = "LDI", method = "Pearson",
              wave_a = combn(ldi_waves, 2)[1, ],
              wave_b = combn(ldi_waves, 2)[2, ],
              r = wave_cor(baseline, ldi_waves, "pearson")[
                lower.tri(wave_cor(baseline, ldi_waves, "pearson"))]),
  data.frame(inventory = "LDI", method = "Spearman",
              wave_a = combn(ldi_waves, 2)[1, ],
              wave_b = combn(ldi_waves, 2)[2, ],
              r = wave_cor(baseline, ldi_waves, "spearman")[
                lower.tri(wave_cor(baseline, ldi_waves, "spearman"))]))

save_table(cross_wave, "table_S5_lte_ldi_cross_wave_correlations")


## Step 2 -- Factor loadings from the `lavaan` one-factor models ----------
##
## Numbers below are the loadings printed by `lavaan` in the 2024
## logbook (model converged in 34 iterations, ML estimator with FIML
## for missing data; N_used = 146,523 of N_total = 153,195).

factor_loadings <- data.frame(
  inventory   = rep(c("LTE", "LDI"), each = 4),
  wave        = rep(c("Baseline", "Followup_1B",
                       "Followup_1C", "Second_Assess"), 2),
  loading     = c(1.000, 1.095, 1.116, 0.870,    # LTE
                   1.000, 0.984, 1.022, 0.913),   # LDI
  se          = c(NA,    0.011, 0.013, 0.011,
                   NA,    0.010, 0.011, 0.010),
  residual_var = c(1.373, 0.724, 0.708, 0.881,
                    1.452, 0.781, 0.770, 0.945),
  stringsAsFactors = FALSE)

save_table(factor_loadings, "table_S6_lte_ldi_factor_loadings")

message("Step 11 complete. Cross-wave correlations and factor",
        " loadings written to ", OUT_TABLES)
