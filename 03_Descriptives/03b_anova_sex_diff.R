#### =========================================================================
#### 03b_anova_sex_diff.R
#### -------------------------------------------------------------------------
#### Paper     : Reichelt F. et al. -- Familial Correlation between Personality
####             and Stress Exposure: Insights from the Multigenerational
####             Lifelines Cohort Study. European Journal of Personality
####             (manuscript EJP-26-4766, in revision).
#### Author    : Felix Reichelt  <f.reichelt@umcg.nl>
####
#### What this script does
#### ---------------------
#### Reproduces the sex-difference test on every phenotype reported in
#### Results 3.1 (and Suppl. Table S7c) using one-way ANOVA computed
#### directly from the published group means, SDs and Ns -- i.e. it
#### does not require access to the individual-level Lifelines data.
####
#### Why a "summary-stats ANOVA"?
#### ----------------------------
#### For two independent groups, the F statistic of a one-way ANOVA on
#### the raw data is *exactly* equal to the F you get from the formulas
####
####     SS_between = n1 * (m1 - m_bar)^2 + n2 * (m2 - m_bar)^2
####     SS_within  = (n1 - 1) * sd1^2 + (n2 - 1) * sd2^2
####     F          = (SS_between / 1) / (SS_within / (n1 + n2 - 2))
####
#### where m_bar is the grand mean weighted by n. Reviewers asked for
#### sex-difference p-values; this script lets anyone independently
#### check them without touching any Lifelines record.
####
#### Outputs   : output/tables/table_S3_sex_differences.tsv
####
#### Manuscript: Results 3.1; Suppl. Table S7c.
####
#### Acknowledgement
#### ---------------
#### Logic taken verbatim from Felix's earlier `ANOVA_sex_diff.R`.
#### Version history:
#### 2026-07-15  v1.1  Ani-aligned multiple testing (/16 uni, /75 bivar);
####                   per-component 1-df bivariate LRTs; INTERNAL tables
#### 2026-07-15  v1.0  Initial release
#### =========================================================================


## Required libraries -----------------------------------------------------
##
## No add-on libraries are required: everything below is base R
## (`data.frame`, `mapply`, `pf`, `save_table()` from `01_helpers.R`).


## Step 1 -- Embedded summary statistics by sex --------------------------
##
## Numbers below are the n/mean/SD reported in Suppl. Table S7a of the
## paper, separately for female and male respondents (1 = male,
## 2 = female in the raw Lifelines GENDER1M2F variable).

sex_summary <- data.frame(
  variable     = c("Age_years",
                    "Neuroticism", "Extraversion", "Conscientiousness",
                    "Loneliness", "Lack_of_social_support", "Childhoodtrauma",
                    "LTE_baseline_1A", "LTE_followup_1B", "LTE_second_1C",
                    "LTE_second_2A",   "LTE_latent_factor",
                    "LDI_baseline_1A", "LDI_followup_1B", "LDI_second_1C",
                    "LDI_second_2A",   "LDI_latent_factor"),
  n_female     = c(63386, 18656, 18650, 18655,
                    25941, 59996, 25897,
                    36072, 34440, 31414, 30916, 37446,
                    35947, 34503, 31392, 31044, 60790),
  mean_female  = c(45.30, 123.30, 40.00, 45.40,
                    11.20,  16.15, 32.30,
                     0.98,   0.69,  0.71,  0.70, -0.03,
                     2.00,   1.62,  1.67,  1.71, -0.21),
  sd_female    = c(13.20, 17.00,  6.00,  5.40,
                     3.40,  3.58,  6.30,
                     1.24,  1.00,  1.01,  1.01,  0.43,
                     2.10,  1.89,  2.00,  1.98,  1.50),
  n_male       = c(89336, 26440, 26433, 26438,
                    37611, 84546, 37605,
                    61053, 58225, 52751, 52212, 63280,
                    60822, 58328, 52755, 52381, 85648),
  mean_male    = c(44.17, 130.60, 39.60, 45.10,
                    11.02,  15.68, 32.90,
                     1.05,   0.80,  0.80,  0.81,  0.02,
                     2.45,   2.04,  2.03,  2.04,  0.15),
  sd_male      = c(13.10, 17.90,  5.80,  5.20,
                     3.44,  3.53,  8.00,
                     1.25,  1.06,  1.06,  1.07,  0.44,
                     2.41,  2.20,  2.26,  2.23,  1.80),
  stringsAsFactors = FALSE)


## Step 2 -- ANOVA from summary statistics (two-group case) -------------

anova_from_summary <- function(n1, m1, sd1, n2, m2, sd2) {
  N          <- n1 + n2                                              # total sample size across both groups
  grand_mean <- (n1 * m1 + n2 * m2) / N                              # sample-size-weighted overall mean
  ss_between <- n1 * (m1 - grand_mean)^2 + n2 * (m2 - grand_mean)^2  # between-group sum of squares (group means vs grand mean)
  ss_within  <- (n1 - 1) * sd1^2 + (n2 - 1) * sd2^2                  # within-group sum of squares (from the published SDs)
  df1        <- 1                                                    # k - 1 = 1 for two groups
  df2        <- N - 2                                                # N - k = N - 2 for two groups
  f_value    <- (ss_between / df1) / (ss_within / df2)               # ANOVA F = MS_between / MS_within
  list(F = f_value, df1 = df1, df2 = df2,
       p = pf(f_value, df1, df2, lower.tail = FALSE))                # upper-tail p-value
}


## Step 3 -- Apply to every variable and write the table ----------------

res <- mapply(anova_from_summary,
              sex_summary$n_female, sex_summary$mean_female, sex_summary$sd_female,
              sex_summary$n_male,   sex_summary$mean_male,   sex_summary$sd_male,
              SIMPLIFY = FALSE)

sex_summary$F_value         <- sapply(res, function(x) x$F)
sex_summary$df1             <- sapply(res, function(x) x$df1)
sex_summary$df2             <- sapply(res, function(x) x$df2)
sex_summary$p_value         <- sapply(res, function(x) x$p)
sex_summary$significant_001 <- sex_summary$p_value < 0.01

save_table(sex_summary, "table_S7c_sex_differences_anova")

message("Step 02b complete. Sex-difference ANOVA written to ",
        file.path(OUT_TABLES, "table_S3_sex_differences.tsv"))
