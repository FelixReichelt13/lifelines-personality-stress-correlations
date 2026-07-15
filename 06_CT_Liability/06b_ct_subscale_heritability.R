#### =========================================================================
#### 06b_ct_subscale_heritability.R
#### -------------------------------------------------------------------------
#### Paper     : Reichelt F. et al. -- Familial Correlation between Personality
####             and Stress Exposure: Insights from the Multigenerational
####             Lifelines Cohort Study. European Journal of Personality
####             (manuscript EJP-26-4766, in revision).
#### Author    : Felix Reichelt  <f.reichelt@umcg.nl>
####
#### What this script does
#### ---------------------
#### **Standalone.** Reports the univariate heritability (ASReml) of each
#### of the five CTQ-SF subscales -- emotional abuse, physical abuse,
#### sexual abuse, emotional neglect and physical neglect -- on the
#### baseline (raw) scale and on the inverse-normal-transformed (INT)
#### scale. It also writes the *liability-scale* heritability of each
#### subscale at the moderate and at the low CTQ cut-offs (Bernstein
#### et al. 2003), obtained with the Wray & Gottesman (2012)
#### recurrence-risk-ratio method.
####
#### Two distinct sources are merged below:
####   * the ASReml continuous estimates are from script 04 (Suppl.
####     Table S8 baseline, Suppl. Table S10 INT).
####   * the liability h^2 numbers + 95 % CIs (CT total *and* the five
####     subscales) come from the cluster-robust `marginaleffects` +
####     `sandwich::vcovCL()` pipeline implemented in
####     `06_rrr_liability_heritability.R` (first-degree-relative
####     RRR -> Wray-Gottesman liability conversion). These are the
####     numbers reported in Suppl. Table S12. Emotional-neglect at
####     the moderate cut-off is degenerate (only n = 844 cases pass
####     threshold, CI spans negative values, see note below).
####   * the low-cut-off subscale rows use subscale-specific clinical
####     thresholds (Bernstein et al. 2003), not the total-score ≥ 30 rule
####     applied to each subscale sum.
####
#### Exploratory / supplementary — not part of the primary /16 or /75
#### Bonferroni families in scripts 04–04c.
####
#### Running this script does NOT require ASReml or Lifelines data --
#### it simply assembles the table into one tidy file for the supplement.
####
#### Manuscript: Suppl. Table S12 (CT liability h^2, total + subscales);
####             Discussion ("heritability ranging from 30 % to 62 %
####             for emotional neglect, emotional abuse, and physical
####             abuse" -- this script makes those numbers explicit).
####
#### Notes
#### -----
#### The liability-threshold pipeline (Wray & Gottesman, 2012) is
#### re-implemented in `06_rrr_liability_heritability.R`; the
#### numbers embedded below are this paper's final estimates for the
#### subscales (taken from the production run on the Lifelines HPC)
#### and the CT-total row.
####
#### CT coding  : Moderate/severe thresholds per Bernstein et al. (2003):
####              ea >= 12, pa >= 9, sa >= 7, en >= 14, pn >= 9 (each on
####              its own subscale sum). The CT-total row uses the same
####              any-subscale rule as `02_data_preparation.R`, not total >= 40.
#### Version history:
#### 2026-07-15  v1.2  Subscale-specific Bernstein thresholds documented
#### 2026-07-15  v1.1  Ani-aligned multiple testing (/16 uni, /75 bivar);
####                   per-component 1-df bivariate LRTs; INTERNAL tables
#### 2026-07-15  v1.0  Initial release
#### =========================================================================


## Required libraries -----------------------------------------------------
##
## No add-on libraries needed: this script just assembles a data.frame
## of pre-computed values and calls `save_table()` (defined in
## `01_helpers.R`). Useful as the absolute minimum reproduction script
## for the supplementary CT-subscale table.


## Step 1 -- Embedded subscale results --------------------------------------
##
## Columns
##   - subscale             : CTQ-SF subscale
##   - h2_asreml_baseline   : ASReml h^2 on raw continuous scale (script 04,
##                            Suppl. Table S8)
##   - h2_asreml_int        : ASReml h^2 on inverse-normal-transformed
##                            scale (script 04bc, Suppl. Table S10)
##   - n_above_mod          : cases above the **subscale-specific** moderate
##                            clinical threshold (Bernstein et al. 2003); NOT
##                            total CTQ ≥ 40 applied to each subscale sum
##   - h2_liability_mod     : Wray-Gottesman liability h^2 at moderate
##                            cut-off (first-degree-relative pipeline,
##                            Suppl. Table S12)
##   - h2_liab_mod_lo / _hi : 95 % CI on `h2_liability_mod`, taken from
##                            the cluster-robust `marginaleffects` +
##                            `sandwich::vcovCL()` pipeline reported in
##                            Suppl. Table S12
##   - n_above_low          : cases above the subscale-specific low threshold
##                            (Bernstein et al. 2003); NOT total CTQ ≥ 30 on
##                            each subscale unless noted for CT total row
##   - h2_liability_low     : Wray-Gottesman liability h^2 at low cut-off

ct_subscales <- data.frame(
  subscale           = c("Childhood_Trauma_total",
                          "Emotional_Abuse",
                          "Emotional_Neglect",
                          "Physical_Abuse",
                          "Sexual_Abuse",
                          "Physical_Neglect"),
  h2_asreml_baseline = c(0.507, 0.411, 0.307, 0.369, 0.334, 0.235),
  h2_asreml_int      = c(0.412, 0.224, 0.280, 0.313, 0.215, 0.069),
  n_above_mod        = c(3376, 3968,  844, 1964, 2735, 7519),
  h2_liability_mod   = c(0.572, 0.663, 0.023, 0.717, 0.439, 0.671),
  h2_liab_mod_lo     = c(0.487, 0.534, -0.327, 0.555, 0.328, 0.532),
  h2_liab_mod_hi     = c(0.658, 0.792,  0.372, 0.879, 0.549, 0.811),
  n_above_low        = c(34972, 10872, 27620, 3922, 8752, 9035),
  h2_liability_low   = c(0.258, 0.295, 0.313, 0.203, 0.246, 0.286),
  stringsAsFactors   = FALSE
)


## Step 2 -- Write the supplementary table ------------------------------------

save_table(ct_subscales, "table_S12b_ct_subscale_heritability")


## Step 3 -- A note on the degenerate emotional-neglect estimate ------------
##
## At the *moderate* cut-off, emotional neglect has very few cases above
## threshold (n = 844). The cluster-robust pipeline gives an h^2 point
## estimate of 2.3 % with a 95 % CI of [-32.7 %, 37.2 %] -- i.e. the CI
## spans zero and the estimate is uninformative. Applying the lower
## cut-off (>= 30) yields 31.3 %, very close to the ASReml continuous
## estimate (30.7 %), suggesting the small-sample behaviour at the
## moderate cut-off is the dominant cause. This is discussed in the
## manuscript supplement.

message("Step 07b complete. CT-subscale liability h^2 table written to ",
        file.path(OUT_TABLES, "table_S12b_ct_subscale_heritability.tsv"))
