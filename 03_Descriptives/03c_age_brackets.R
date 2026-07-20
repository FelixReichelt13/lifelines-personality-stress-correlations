#### =========================================================================
#### 03c_age_brackets.R
#### -------------------------------------------------------------------------
#### Paper     : Reichelt F. et al. -- Familial Correlation between Personality
####             and Stress Exposure: Insights from the Multigenerational
####             Lifelines Cohort Study. European Journal of Personality
####             (manuscript EJP-26-4766, in revision).
#### Author    : Felix Reichelt  <f.reichelt@umcg.nl>
####
#### What this script does
#### ---------------------
#### Computes prevalence/level of each stress exposure and personality
#### trait across seven age brackets (18-24, 25-34, 35-44, 45-54, 55-64,
#### 65-74, 75+). Writes one tidy file with the mean, SD, median, IQR,
#### min, max and N per bracket per phenotype. Used in:
####
####   - Supplementary Figure S3 (age-trend plots; rendered in
####     07b_figures_supplementary.R)
####   - Results section 3.1 ("phenotype values are largely stable
####     across the adult age range").
####
#### This script does *not* require ASReml-R.
####
#### Inputs    : data/data_prepared.Rdata (from script 02).
#### Outputs   : output/tables/table_S7_age_bracket_descriptives.tsv
####
#### Related outputs: Table S7d; Figure S2.
####
#### Pitfalls
#### --------
#### * Lifelines collects some phenotypes at the baseline assessment
####   (1A1, age column = AGE) and some at the second assessment (2B,
####   age column = AGE_2B). The age bracket used for childhood trauma
####   and loneliness should therefore be AGE_2B (when available),
####   because those questionnaires were only added at the second
####   assessment. The script uses AGE_2B when present and falls back
####   to AGE otherwise.
#### Version history:
#### 2026-07-15  v1.1  Ani-aligned multiple testing (/16 uni, /75 bivar);
####                   per-component 1-df bivariate LRTs; INTERNAL tables
#### 2026-07-15  v1.0  Initial release
#### =========================================================================


## Required libraries -----------------------------------------------------
##
## Loaded centrally by `01_setup.R`. This script has no ASReml
## dependency; it only summarises the prepared frame.

library(data.table)   # `data.table` syntax, `copy()`, `rbindlist()`


load(file.path(DATA_DIR, "data_prepared.Rdata"))   # provides `baseline` (the analytic data.table)

age_breaks <- c(18, 25, 35, 45, 55, 65, 75, Inf)                      # right-open intervals (see `right = FALSE` below)
age_labels <- c("18-24", "25-34", "35-44", "45-54", "55-64", "65-74", "75+")

dt <- copy(baseline)                                                  # `copy()` so we don't mutate the loaded object
dt[, age_for_bracket := if ("AGE_2B" %in% names(dt)) {
                          ifelse(!is.na(AGE_2B), AGE_2B, AGE)         # prefer age at 2B (CTQ, loneliness assessment)
                        } else AGE]                                   # fall back to baseline age (most other phenotypes)
dt[, age_bracket := cut(age_for_bracket,
                         breaks = age_breaks,
                         labels = age_labels,
                         right  = FALSE)]                             # `[18, 25)`, `[25, 35)`, ... -- left-closed, right-open


## Step 1 -- Long-format descriptives per phenotype x bracket ------------

phenos <- intersect(c(STRESS_PHENOTYPES, PERSONALITY_PHENOTYPES), names(dt))

summary_per_bracket <- rbindlist(lapply(phenos, function(p) {
  d <- dt[!is.na(get(p)) & !is.na(age_bracket),
          .(N      = .N,
            mean   = mean  (get(p), na.rm = TRUE),
            sd     = sd    (get(p), na.rm = TRUE),
            median = median(get(p), na.rm = TRUE),
            iqr    = IQR   (get(p), na.rm = TRUE),
            min    = min   (get(p), na.rm = TRUE),
            max    = max   (get(p), na.rm = TRUE)),
          by = age_bracket]
  d[, phenotype := p]
  d
}))

setcolorder(summary_per_bracket,
            c("phenotype", "age_bracket", "N",
              "mean", "sd", "median", "iqr", "min", "max"))

save_table(summary_per_bracket, "table_S7d_age_bracket_descriptives")

message("Step 09 complete. ", nrow(summary_per_bracket),
        " phenotype x age-bracket rows written.")
