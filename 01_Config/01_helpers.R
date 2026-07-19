#### =========================================================================
#### 01_helpers.R
#### -------------------------------------------------------------------------
#### Paper     : Reichelt F. et al. -- Familial Correlation between Personality
####             and Stress Exposure: Insights from the Multigenerational
####             Lifelines Cohort Study. European Journal of Personality
####             (manuscript EJP-26-4766, in revision).
#### Author    : Felix Reichelt  <f.reichelt@umcg.nl>
####
#### Purpose   : Tiny utility functions that are reused by more than one
####             analysis script. Nothing here is specific to a single step
####             of the pipeline.
####
#### Loaded automatically by 01_setup.R.
#### Version history:
#### 2026-07-15  v1.2  CT dichotomous flag from subscale-specific thresholds
#### 2026-07-15  v1.1  Ani-aligned multiple testing (/16 uni, /75 bivar);
####                   per-component 1-df bivariate LRTs; INTERNAL tables
#### 2026-07-15  v1.0  Initial release
#### =========================================================================


## inverse_normal() -------------------------------------------------------
##
## Rank-based inverse normal transformation. Used in script 04bc
## as a sensitivity analysis: ASReml ACE models assume residual normality,
## and several of the stress phenotypes are noticeably right-skewed (see
## Suppl. Fig. S18). The INT mapping used in this pipeline is
##
##     y_i  ->  qnorm( (rank(y_i) - 0.5) / N )
##
## where N is the number of non-missing observations. This is the
## (rank - 0.5) / N variant (sometimes called the "Bliss" or "Tukey"
## continuity-corrected INT). It brings every phenotype to an
## approximately standard normal scale without changing the ordering
## of observations.
##
## Pitfall: with ties, `rank(..., ties.method = "average")` would map
## several observations to the same quantile. We break ties at random
## (`ties.method = "random"`); reproducibility is guaranteed because
## 01_setup.R calls `set.seed(12345)` once before any data is loaded.

inverse_normal <- function(x) {
  n <- sum(!is.na(x))
  if (n == 0) return(rep(NA_real_, length(x)))
  ranks <- rank(x, na.last = "keep", ties.method = "random")
  qnorm((ranks - 0.5) / n)
}


## to_factor() ------------------------------------------------------------
##
## ASReml requires the pedigree IDs and the household / sex columns to be
## factors -- numeric IDs are silently treated as covariates. This wrapper
## converts every column in `cols` that exists in the data.table `dt` to
## a factor, in place, without printing warnings for missing ones.

to_factor <- function(dt, cols) {
  for (col in cols) {
    if (col %in% names(dt)) dt[, (col) := factor(get(col))]
  }
  dt
}


## vp_row() ---------------------------------------------------------------
##
## `vpredict()` returns a one-row data.frame with the columns "Estimate"
## and "SE". The two helpers below standardise the output across scripts
## so the result tables all use the same column names.

vp_row <- function(vp) {
  data.frame(estimate = vp[1, "Estimate"], se = vp[1, "SE"])
}


## save_table() / save_figure() -------------------------------------------
##
## Convenience writers so that every script saves to the same place with
## the same conventions (tab-separated tables, PNG + PDF figures at 300
## dpi).

save_table <- function(df, name) {
  out <- file.path(OUT_TABLES, paste0(name, ".tsv"))
  write.table(df, out, sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(out)
}

save_figure <- function(plot, name, width = 10, height = 8) {
  ggsave(file.path(OUT_FIGURES, paste0(name, ".png")),
         plot, width = width, height = height, dpi = 300)
  ggsave(file.path(OUT_FIGURES, paste0(name, ".pdf")),
         plot, width = width, height = height)
  invisible(NULL)
}


## ct_dichotomous_moderate() ----------------------------------------------
##
## Dichotomous childhood trauma for liability / RRR analyses (script 04bc).
## CT coding: present (1) when at least one CTQ-SF subscale meets the
## literature-based moderate/severe threshold (Bernstein et al., 2003):
##   emotional abuse >= 12, physical abuse >= 9, sexual abuse >= 7,
##   emotional neglect >= 14, physical neglect >= 9.
## Coded 0 only when all five subscales are observed and none meets its
## threshold; otherwise NA. Continuous `Childhoodtrauma` total is unchanged
## and is used in primary ACE models (03–05).

ct_dichotomous_moderate <- function(dt,
                                    subcols = CTQ_SUBSCALE_COLS,
                                    thresholds = CTQ_MODERATE_THRESHOLDS) {
  out <- rep(NA_integer_, nrow(dt))
  if (!all(subcols %in% names(dt))) return(out)
  m <- as.matrix(dt[, subcols, with = FALSE])
  above <- m >= matrix(thresholds,
                       nrow = nrow(m), ncol = ncol(m), byrow = TRUE)
  above[is.na(m)] <- FALSE
  any_above <- rowSums(above) > 0L
  all_obs <- rowSums(is.na(m)) == 0L
  out[any_above] <- 1L
  out[all_obs & !any_above] <- 0L
  out
}
