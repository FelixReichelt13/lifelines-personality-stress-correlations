#### =========================================================================
#### 03_descriptives.R
#### -------------------------------------------------------------------------
#### Paper     : Reichelt F. et al. -- Familial Correlation between Personality
####             and Stress Exposure: Insights from the Multigenerational
####             Lifelines Cohort Study. European Journal of Personality
####             (manuscript EJP-26-4766, in revision).
#### Author    : Felix Reichelt  <f.reichelt@umcg.nl>
####
#### What this script does
#### ---------------------
#### 1. Computes sample size and overall mean / SD / median / IQR for
####    every phenotype -- this is **Table 1** in the main paper.
#### 2. Computes the same numbers separately for men and women -- this
####    feeds **Supplementary Table S7** (descriptives by sex; age-group
####    descriptives are in `03c_age_brackets.R`).
#### 3. Tests sex differences with an age-adjusted linear model
####    (likelihood-ratio test between models with and without sex).
####
#### Inputs    : data/data_prepared.Rdata (from script 02).
#### Outputs   : three TSV tables in output/tables/:
####             - table_1_descriptives_overall.tsv         (main paper)
####             - table_S7a_descriptives_by_sex.tsv        (Suppl. S7)
####             - table_S7b_sex_differences_age_adjusted.tsv (Suppl. S7)
####
####             Childhood-trauma prevalence is reported in
####             Table S11 by `06_rrr_liability_heritability.R`.
####
#### Notes     : this script does not require ASReml-R; everything here is
####             plain base-R + dplyr, so it can be run by anyone reviewing
####             the paper without buying a license.
#### Version history:
#### 2026-07-15  v1.1  Ani-aligned multiple testing (/16 uni, /75 bivar);
####                   per-component 1-df bivariate LRTs; INTERNAL tables
#### 2026-07-15  v1.0  Initial release
#### =========================================================================


## Required libraries -----------------------------------------------------
##
## Loaded centrally by `01_setup.R`. The fixed-effects sex-difference
## test (Step 3) uses only base-R `lm()` and `anova()`, so this script
## has no ASReml dependency.

library(data.table)   # `data.table()`, `melt()` long format, `[i, j, by]`


load(file.path(DATA_DIR, "data_prepared.Rdata"))   # provides `baseline`, `baseline_int`, `family_data`

phenotypes <- intersect(c(STRESS_PHENOTYPES, PERSONALITY_PHENOTYPES),
                        names(baseline))


## Step 1 — Table 1: overall sample descriptives ------------------------
##
## Manuscript Table 1: N, mean (SD), median (IQR) per phenotype, plus
## age and sex composition of the cohort.

desc_overall <- data.table(
  Phenotype = c("AGE_1A1", phenotypes),
  N         = sapply(c("AGE_1A1", phenotypes),
                     function(p) sum(!is.na(baseline[[p]]))),
  Mean      = sapply(c("AGE_1A1", phenotypes),
                     function(p) mean(baseline[[p]], na.rm = TRUE)),
  SD        = sapply(c("AGE_1A1", phenotypes),
                     function(p) sd(baseline[[p]],   na.rm = TRUE)),
  Median    = sapply(c("AGE_1A1", phenotypes),
                     function(p) median(baseline[[p]], na.rm = TRUE)),
  IQR       = sapply(c("AGE_1A1", phenotypes),
                     function(p) IQR(baseline[[p]], na.rm = TRUE))
)

## Add sex composition as a separate row
sex_row <- data.table(
  Phenotype = "Female_pct",
  N         = sum(!is.na(baseline$GENDER)),
  Mean      = mean(baseline$GENDER == 2L, na.rm = TRUE) * 100,   # GENDER coding in Lifelines: 1 = male, 2 = female
  SD        = NA_real_, Median = NA_real_, IQR = NA_real_
)
desc_overall <- rbind(desc_overall, sex_row, fill = TRUE)
save_table(desc_overall, "table_1_descriptives_overall")


## Step 2 — Suppl. Table S7a: means and SDs by sex ---------------------

desc_by_sex <- baseline[, lapply(.SD, function(x) {
    c(N = sum(!is.na(x)),
      mean = mean(x, na.rm = TRUE),
      sd   = sd(x,   na.rm = TRUE))
  }), by = GENDER, .SDcols = phenotypes]

desc_long <- melt(desc_by_sex,
                  id.vars       = "GENDER",
                  variable.name = "Phenotype",
                  value.name    = "value")
save_table(desc_long, "table_S7a_descriptives_by_sex")


## Step 3 — Suppl. Table S7b: age-adjusted sex differences --------------
##
## We compare a model with sex against a model without sex (both adjust
## for AGE_c and AGE_c^2). The LRT-style ANOVA gives a single p-value
## per phenotype that is reported in Suppl. Table S7b.

sex_diff <- lapply(phenotypes, function(p) {
  m_full <- lm(reformulate(c("GENDER", "AGE_c", "AGE_c2"), response = p),  # phenotype ~ sex + age + age^2
               data = baseline)
  m_red  <- lm(reformulate(c("AGE_c", "AGE_c2"),           response = p),  # phenotype ~ age + age^2 (drop sex)
               data = baseline)
  lr     <- anova(m_red, m_full)                                            # F test for the sex term, conditional on age
  data.frame(Phenotype = p,
             coef_sex  = coef(m_full)[grep("GENDER", names(coef(m_full)))[1]],  # sex coefficient (in original phenotype units)
             p_value   = lr[2, "Pr(>F)"])
})
sex_diff <- do.call(rbind, sex_diff)
save_table(sex_diff, "table_S7b_sex_differences_age_adjusted")

message("Step 02 complete. Table 1 and S7a/S7b written to ", OUT_TABLES)
