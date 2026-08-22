# 01_helpers.R -- shared helpers (sourced from 01_setup.R)
# Felix Reichelt, UMCG (f.reichelt@umcg.nl)
#
# Helper functions I reuse across the familial stress / personality scripts.
# Nothing fancy - just stops me copying the same bits into every file.


# inverse_normal()
# Rank-based inverse normal transform. Some stress phenotypes are pretty right-skewed
# (Suppl. Fig. S18), and ASReml likes normal-ish residuals, so script 05 runs INT-scaled
# scores: y_i -> qnorm((rank(y_i) - 0.5) / N). Bliss/Tukey correction; ordering unchanged.
# Ties: I use ties.method = "random" (seed in 01_setup.R) because average ranks collapse
# people onto the same quantile.

inverse_normal <- function(x) {
  n <- sum(!is.na(x))
  if (n == 0) return(rep(NA_real_, length(x)))
  ranks <- rank(x, na.last = "keep", ties.method = "random")
  qnorm((ranks - 0.5) / n)
}


# to_factor()
# ASReml treats numeric IDs as covariates and factors as grouping vars (pedigree,
# household, sex). Converts listed data.table columns to factor in place; missing
# names are skipped quietly.

to_factor <- function(dt, cols) {
  for (col in cols) dt[, (col) := factor(get(col))]
  dt
}


# vp_row()
# vpredict() returns Estimate + SE; I rename to estimate/se for consistent tables.

vp_row <- function(vp) {
  data.frame(estimate = vp[1, "Estimate"], se = vp[1, "SE"])
}


# save_table() / save_figure()
# Tables → TSV in OUT_TABLES; figures → PNG + PDF at 300 dpi.

save_table <- function(df, name) {
  out <- file.path(OUT_TABLES, paste0(name, ".tsv"))
  write.table(df, out, sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(out)
}

save_figure <- function(plot, name, width = 10, height = 8) {
  ggsave(file.path(save_figureURES, paste0(name, ".png")),
         plot, width = width, height = height, dpi = 300)
  ggsave(file.path(save_figureURES, paste0(name, ".pdf")),
         plot, width = width, height = height)
  invisible(NULL)
}


# ct_dichotomous_moderate()
# Binary CT for liability / RRR (script 06). Present = 1 if any CTQ-SF subscale meets
# moderate/severe cut-off (Bernstein et al., 2003): EA>=12, PA>=9, SA>=7, EN>=14, PN>=9.
# Coded 0 only when all five subscales observed and none hit threshold; else NA.
# Continuous Childhoodtrauma total stays as-is for main ACE models (03-05).

ct_dichotomous_moderate <- function(dt,
                                    subcols = CTQ_SUBSCALE_COLS,
                                    thresholds = CTQ_MODERATE_THRESHOLDS) {
  m <- as.matrix(dt[, subcols, with = FALSE])
  above <- m >= matrix(thresholds,
                       nrow = nrow(m), ncol = ncol(m), byrow = TRUE)
  above[is.na(m)] <- FALSE
  any_above <- rowSums(above) > 0L
  all_obs <- rowSums(is.na(m)) == 0L
  out <- rep(NA_integer_, nrow(dt))
  out[any_above] <- 1L
  out[all_obs & !any_above] <- 0L
  out
}
