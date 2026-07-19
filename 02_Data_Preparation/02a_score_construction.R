#### =========================================================================
#### 02a_score_construction.R
#### -------------------------------------------------------------------------
#### Paper     : Reichelt F. et al. -- Familial Correlation between Personality
####             and Stress Exposure: Insights from the Multigenerational
####             Lifelines Cohort Study. European Journal of Personality
####             (manuscript EJP-26-4766, in revision).
#### Author    : Felix Reichelt  <f.reichelt@umcg.nl>
####
#### Reads the four raw questionnaire files distributed by Lifelines
#### (CTQ, Loneliness, Social Support, Mental Health for LTE/LDI) and
#### writes a single tidy file `data/phenotypes_lifelines.Rdata` with one
#### row per PSEUDOIDEXT and one column per scale-level score:
####
####   Childhoodtrauma                   CTQ-SF total (25 clinical items)
####   ea_sum / pa_sum / sa_sum /
####     en_sum / pn_sum                 CTQ-SF subscales
####   Loneliness                        6-item De Jong Gierveld sum
####   Lack_of_social_support            9-item SPF-IL reversed sum
####   LTE_Baseline / _1B / _1C / _2A    LTE raw sum per wave
####   LDI_Baseline / _1B / _1C / _2A    LDI raw sum per wave
####
#### The cross-wave LTE / LDI latent factors used in the main analyses
#### are then constructed in `02_data_preparation.R`.
####
#### Scoring decisions:
####   - item recoding (reverse items) as in the questionnaire manuals
####   - mean-imputation of missing items when fewer than half are
####     missing per subscale (CTQ, LL, SS)
####   - complete-case sum scoring for LTE / LDI (imputing the occurrence
####     of random events would be theoretically implausible)
####
#### Public source references (Lifelines Wiki):
####   - CTQ-SF      https://wiki.lifelines.nl/doku.php?id=childhood_trauma_ctq
####   - Loneliness  https://wiki.lifelines.nl/doku.php?id=loneliness_ls
####   - SPF-IL      https://wiki.lifelines.nl/doku.php?id=social_functioning_spf-il
####   - LTE         https://wiki.lifelines.nl/doku.php?id=stress_lte
####   - LDI         https://wiki.lifelines.nl/doku.php?id=stress_ldi
####   - NEO         https://wiki.lifelines.nl/doku.php?id=personality_neo
####   - NEO SUMs    https://wiki.lifelines.nl/doku.php?id=personality_sum_scores_neo
####
#### Pitfalls:
####   - Lifelines distributes the questionnaires as labelled SPSS or as
####     CSV with "$7" placeholders for missing values. This script reads
####     the CSV; "$7" is recoded to NA before scoring.
####   - Some CTQ items are reverse-coded (q2, q5, q7, q13, q19, q26, q28);
####     the recoding (5 -> 1, 4 -> 2, 2 -> 4, 1 -> 5) is applied first.
####   - LTE / LDI use complete-case scoring per wave (no mean imputation).
####
#### CT coding  : Subscale sums (`ea_sum` … `pn_sum`) feed the dichotomous CT
####              flag in `02_data_preparation.R` (script 04bc liability h²).
####              Childhood trauma is present when at least one CTQ-SF subscale
####              meets the literature-based moderate/severe threshold used in
####              this project and Bernstein et al. (2003): emotional abuse
####              >= 12, physical abuse >= 9, sexual abuse >= 7, emotional
####              neglect >= 14 or physical neglect >= 9. A participant is
####              coded 0 only when all five subscales are observed and none
####              meets its threshold. The continuous total `Childhoodtrauma`
####              is unchanged and is used in primary ACE models (03–05).
####
#### Version history:
#### 2026-07-15  v1.2  Document subscale-specific Bernstein thresholds for CT
#### 2026-07-15  v1.1  Ani-aligned multiple testing (/16 uni, /75 bivar);
####                   per-component 1-df bivariate LRTs; INTERNAL tables
#### 2026-07-15  v1.0  Initial release
#### =========================================================================


## Required libraries -----------------------------------------------------
##
## Loaded centrally by `01_setup.R`; re-declared here so the script
## reads as self-contained.

library(data.table)   # `data.table()` constructor, `[i, j, by]`, in-place `:=`, `.SD` for column-wise ops
library(readr)        # `read_csv()` to ingest the four raw item-level CSV files


## Helpers used throughout this script ---------------------------------
##
##   reverse_5pt : Likert 1..5 -> 5..1 (CTQ reverse coding,
##                                       per Bernstein et al. 1998)
##   sum_with_mean_imputation : compute sum of selected items, but if
##     up to `max_missing` items are NA, replace them with the mean of
##     the observed items. Returns NA when too many items are missing.

reverse_5pt <- function(x) {
  x <- as.numeric(x)
  6 - x
}

sum_with_mean_imputation <- function(items, max_missing) {
  items   <- as.matrix(items)                                          # one row per participant, one column per item
  items   <- apply(items, 2, as.numeric)                               # CSV columns may read as character; coerce
  n_miss  <- rowSums(is.na(items))                                     # number of missing items per participant
  means   <- rowMeans(items, na.rm = TRUE)                             # per-person mean of observed items (used for imputation)
  items[is.na(items)] <- rep(means, times = sum(is.na(items)) / length(means))   # fill NAs with that person's own item mean
  s <- rowSums(items)                                                  # sum across items (after imputation)
  s[n_miss > max_missing] <- NA_real_                                  # too many missing? -> return NA (NOT 0)
  s
}


## Step 1 -- Childhood Trauma (CTQ-SF, Bernstein et al., 1998) ----------
##
## Lifelines administers the 28-item CTQ-SF only at assessment 2B
## (see https://wiki.lifelines.nl/doku.php?id=childhood_trauma_ctq).
## Item assignment to subscales is taken from the wiki:
##   emotional abuse    : q3, q8, q14, q18, q25
##   emotional neglect  : q5, q7, q13, q19, q28
##   physical abuse     : q9, q11, q12, q15, q17
##   physical neglect   : q1, q2, q4, q6, q26
##   sexual abuse       : q20, q21, q23, q24, q27
##   minimisation/denial: q10, q16, q22  (DISCARDED -- not summed)
##
## 5-point Likert: 1 Never true, 2 Rarely true, 3 Sometimes true,
## 4 Often true, 5 Very often true. The positively-worded items
## (q2, q5, q7, q13, q19, q26, q28) are reverse-coded so that a higher
## score consistently means *more* exposure to trauma.
##
## Moderate/severe subscale thresholds (Bernstein et al., 2003; used in
## `02_data_preparation.R` for `CT_dichotomous` and in `07b` for liability
## h²): ea >= 12, pa >= 9, sa >= 7, en >= 14, pn >= 9.

ctq <- as.data.table(read_csv(RAW_CTQ_FILE, show_col_types = FALSE))   # path defined in `01_setup.R` (RAW_CTQ_FILE)
setnames(ctq, "PROJECT_PSEUDO_ID", "ID")                                # harmonise the participant ID column name across files

rev_items <- c("q2", "q5", "q7", "q13", "q19", "q26", "q28")           # positively-worded CTQ items -> reverse before scoring
ctq[, (rev_items) := lapply(.SD, reverse_5pt), .SDcols = rev_items]    # `:=` updates in place; `.SD` = the columns listed in `.SDcols`

ea_items <- c("q3", "q8", "q14", "q18", "q25")             # emotional abuse
pa_items <- c("q9", "q11", "q12", "q15", "q17")            # physical abuse
sa_items <- c("q20", "q21", "q23", "q24", "q27")           # sexual abuse
en_items <- c("q5", "q7", "q13", "q19", "q28")             # emotional neglect
pn_items <- c("q1", "q2", "q4", "q6", "q26")               # physical neglect
all_25   <- c(ea_items, pa_items, sa_items, en_items, pn_items)

ctq[, ea_sum := sum_with_mean_imputation(.SD, max_missing = 2), .SDcols = ea_items]
ctq[, pa_sum := sum_with_mean_imputation(.SD, max_missing = 2), .SDcols = pa_items]
ctq[, sa_sum := sum_with_mean_imputation(.SD, max_missing = 2), .SDcols = sa_items]
ctq[, en_sum := sum_with_mean_imputation(.SD, max_missing = 2), .SDcols = en_items]
ctq[, pn_sum := sum_with_mean_imputation(.SD, max_missing = 2), .SDcols = pn_items]

ctq[, Childhoodtrauma := sum_with_mean_imputation(.SD, max_missing = 12),
     .SDcols = all_25]

ctq_out <- ctq[, .(ID, Childhoodtrauma, ea_sum, pa_sum, sa_sum, en_sum, pn_sum)]


## Step 2 -- Loneliness (De Jong Gierveld 6-item) ----------------------
##
## Lifelines variable names: loneliness_*_adu_q_01..06 (see
## https://wiki.lifelines.nl/doku.php?id=loneliness_ls). Items 1-3
## are negatively-worded (emptiness, missing people, abandoned) and
## items 4-6 are positively-worded (fall-back, trust, connection).
## Under the Lifelines coding convention (1 = "yes!" / strongly agree,
## 5 = "no!" / strongly disagree) we reverse items 1-3 so that high
## values consistently mean *more* loneliness. Mean-impute up to 3
## of the 6 items.

ll <- as.data.table(read_csv(RAW_LL_FILE, show_col_types = FALSE))
ll_items <- paste0("ll", 1:6)
setnames(ll,
         c("PROJECT_PSEUDO_ID",
           "loneliness_emptiness_adu_q_01",
           "loneliness_people_adu_q_02",
           "loneliness_abandoned_adu_q_03",
           "loneliness_fallback_adu_q_04",
           "loneliness_trust_adu_q_05",
           "loneliness_connection_adu_q_06"),
         c("ID", ll_items))

ll[, (ll_items[1:3]) := lapply(.SD, reverse_5pt), .SDcols = ll_items[1:3]]
ll[, Loneliness := sum_with_mean_imputation(.SD, max_missing = 3),
     .SDcols = ll_items]
ll_out <- ll[, .(ID, Loneliness)]


## Step 3 -- Social support (SPF-IL, 9 items) --------------------------
##
## Social Production Function Instrument for Level of Well-Being
## (Nieboer et al., 2005; see
## https://wiki.lifelines.nl/doku.php?id=social_functioning_spf-il).
## The FULL SPF-IL short form has 15 items (6 physical well-being +
## 9 social well-being). The general Lifelines assessments only
## include the **9 social well-being items** -- spfil1..9 (Affection,
## Confirmation, Status) -- which is the subset we use here.
## (The physical well-being items 12..17 are only collected at ACTQ
## for participants aged 60-80 and are not used in this paper.)
##
## Items rated on a 4-point scale: 1 not at all, 2 sometimes,
## 3 often, 4 always. Recoded 1..4 to 0..3, summed (range 0-27).
## We then flip the sign so that "lack of social support" runs in
## the same direction as the other stress phenotypes.

ss <- as.data.table(read_csv(RAW_SS_FILE, show_col_types = FALSE))
ss_items <- paste0("ss", 1:9)
setnames(ss,
         c("PROJECT_PSEUDO_ID",
           paste0("spfil_socsup_adu_q_0", 1:9)),
         c("ID", ss_items))
ss[, (ss_items) := lapply(.SD, function(x) as.numeric(x) - 1),
   .SDcols = ss_items]                                       # rescale Likert 1..4 -> 0..3 (so the total is 0-27 not 9-36)
ss[, ss_total := sum_with_mean_imputation(.SD, max_missing = 4),
   .SDcols = ss_items]                                       # allow up to 4 of 9 items missing (mean-imputed)
ss[, Lack_of_social_support := max(ss_total, na.rm = TRUE) - ss_total]  # flip the sign: high = LOW support (matches other stress phenotypes)
ss_out <- ss[, .(ID, Lack_of_social_support)]


## Step 4 -- LTE and LDI ------------------------------------------------
##
## LTE -- List of Threatening Experiences (Brugha & Cragg, 1990; see
## https://wiki.lifelines.nl/doku.php?id=stress_lte). 12 binary
## life-event categories + a 13th open-text "other" item which we
## DO NOT sum (per Lifelines convention). The 12 LTE items are
## administered at general assessments 1A, 1B, 1C and 2A (and at
## several additional questionnaires that are not used here).
##
## LDI -- Long-term Difficulties Inventory (Rosmalen et al., 2012;
## see https://wiki.lifelines.nl/doku.php?id=stress_ldi). A 3-point
## stress scale (0 not stressful, 1 slightly, 2 very stressful)
## across 12 life domains. **Important:** at the four waves used
## here (1A, 1B, 1C, 2A) only **items 1-9** are administered
## (housing, work, friends, partner, children, parents, relatives,
## free time, finances). Items 10-12 (health, school, religion)
## were added from assessment 3B onward and are *not* in any of the
## four waves used in this paper. The per-wave LDI sum is therefore
## over 9 items, range 0-18.
##
## Complete-case scoring per wave: any missing item -> NA. The
## cross-wave latent factor is fit in `02_data_preparation.R`.

mh <- as.data.table(read_csv(RAW_MENTAL_FILE, show_col_types = FALSE))
setnames(mh, "PROJECT_PSEUDO_ID", "ID")

score_wave <- function(dt, items, name) {
  ## Defensive: only keep columns that are actually present in the
  ## raw extraction (e.g. GLLM10..12 are not in any of the four
  ## waves we use, so they will be silently dropped if present in
  ## `items`). This keeps the script robust to small differences
  ## across Lifelines extractions.
  items_present <- intersect(items, names(dt))
  if (length(items_present) == 0L) {
    dt[[name]] <- NA_real_
    return(dt)
  }
  cols <- dt[, items_present, with = FALSE]
  cols[, (items_present) := lapply(.SD, as.numeric),
       .SDcols = items_present]
  if (grepl("LTE", name)) {
    ## BRUGHA: 1 = yes, 2 = no in the raw, recode 2 -> 0.
    cols[, (items_present) := lapply(.SD,
                                      function(x) ifelse(x == 2, 0, x)),
         .SDcols = items_present]
  } else {
    ## GLLM: 1 = not, 2 = somewhat, 3 = much in the raw,
    ## shift to 0..2.
    cols[, (items_present) := lapply(.SD, function(x) x - 1),
         .SDcols = items_present]
  }
  dt[[name]] <- ifelse(rowSums(is.na(cols)) > 0, NA_real_,
                        rowSums(cols))
  dt
}

lte_waves <- list(
  LTE_Baseline      = paste0("BRUGHA",   1:12),
  LTE_Followup_1B   = paste0("BRUGHA_b", 1:12),
  LTE_Followup_1C   = paste0("BRUGHA_c", 1:12),
  LTE_Second_Assess = paste0("BRUGHA_2", 1:12)
)

## LDI: items 1-9 only at waves 1A, 1B, 1C, 2A (see wiki note above).
ldi_waves <- list(
  LDI_Baseline      = paste0("GLLM",   1:9),
  LDI_Followup_1B   = paste0("GLLM_b", 1:9),
  LDI_Followup_1C   = paste0("GLLM_c", 1:9),
  LDI_Second_Assess = paste0("GLLM_2", 1:9)
)

for (nm in names(lte_waves)) mh <- score_wave(mh, lte_waves[[nm]], nm)
for (nm in names(ldi_waves)) mh <- score_wave(mh, ldi_waves[[nm]], nm)

mh_cols <- c("ID", names(lte_waves), names(ldi_waves))
mh_out  <- mh[, intersect(mh_cols, names(mh)), with = FALSE]


## Step 5 -- Merge everything on PROJECT_PSEUDO_ID --------------------

phenotypes_lifelines <- Reduce(
  function(x, y) merge(x, y, by = "ID", all = TRUE),
  list(ctq_out, ll_out, ss_out, mh_out)
)

save(phenotypes_lifelines,
     file = file.path(DATA_DIR, "phenotypes_lifelines.Rdata"))


## Step 6 -- Internal consistency (Suppl. Table S7e) --------------------
##
## Cronbach's alpha for every multi-item scale where item-level data are
## in scope. For the NEO scales we instead quote published alphas from
## the Dutch NEO-PI-R / NEO-FFI manual (Hoekstra, Ormel & de Fruyt, 2007),
## because the Lifelines personality release ships only the SUM* totals
## and not the items. The appropriate published reliabilities are:
##   - Neuroticism (NEO-PI-R 48-item domain)        : alpha ~ 0.92
##   - Extraversion (NEO-FFI 12-item domain; ceiling
##     for the 11-item version used here)           : alpha ~ 0.79
##   - Conscientiousness (NEO-FFI 12-item domain)   : alpha ~ 0.81

cronbach_alpha <- function(items) {
  m <- as.matrix(items)
  storage.mode(m) <- "double"
  if (ncol(m) < 2L) return(NA_real_)
  cc <- stats::cov(m, use = "pairwise.complete.obs")
  tot <- sum(cc, na.rm = TRUE)
  if (tot <= 0) return(NA_real_)
  k <- ncol(m)
  k / (k - 1) * (1 - sum(diag(cc), na.rm = TRUE) / tot)
}

## CTQ-SF: alpha on the 25 clinical items, after reverse-coding.
ctq_alpha_total <- cronbach_alpha(ctq[, ..all_25])
ctq_alpha_ea    <- cronbach_alpha(ctq[, ..ea_items])
ctq_alpha_pa    <- cronbach_alpha(ctq[, ..pa_items])
ctq_alpha_sa    <- cronbach_alpha(ctq[, ..sa_items])
ctq_alpha_en    <- cronbach_alpha(ctq[, ..en_items])
ctq_alpha_pn    <- cronbach_alpha(ctq[, ..pn_items])
ll_alpha        <- cronbach_alpha(ll[, ..ll_items])
ss_alpha        <- cronbach_alpha(ss[, ..ss_items])

## LTE: per-wave alpha across the 12 binary event items.
## LDI: per-wave alpha across the 9 domain items (items 1-9; the
## three later items are not in waves 1A-2A).
lte_alpha <- sapply(lte_waves, function(items) {
  cronbach_alpha(mh[, intersect(items, names(mh)), with = FALSE])
})
ldi_alpha <- sapply(ldi_waves, function(items) {
  cronbach_alpha(mh[, intersect(items, names(mh)), with = FALSE])
})

internal_consistency <- data.table(
  Scale = c("Childhoodtrauma_total",
            "CT_emotional_abuse",  "CT_physical_abuse",
            "CT_sexual_abuse",     "CT_emotional_neglect",
            "CT_physical_neglect",
            "Loneliness", "SPF_IL_social_support",
            names(lte_alpha), names(ldi_alpha),
            "Neuroticism_NEO_PI_R_48item",
            "Extraversion_NEO_FFI_11item",
            "Conscientiousness_NEO_FFI_12item"),
  n_items = c(25, 5, 5, 5, 5, 5, 6, 9,
              rep(12, length(lte_alpha)),   # LTE: 12 events/wave
              rep(9,  length(ldi_alpha)),   # LDI: 9 domains/wave (items 1-9 only at 1A-2A)
              48, 11, 12),                  # NEO: N is the 6-facet x 8-item
                                            # NEO-PI-R domain sum (= 48 items;
                                            # see 02_data_preparation.R header).
                                            # E is 11 NEO-FFI items (one
                                            # FFI item omitted by Lifelines);
                                            # C is 12 NEO-FFI items.
  Cronbach_alpha = c(ctq_alpha_total,
                     ctq_alpha_ea, ctq_alpha_pa, ctq_alpha_sa,
                     ctq_alpha_en, ctq_alpha_pn,
                     ll_alpha, ss_alpha,
                     lte_alpha, ldi_alpha,
                     ## Published Dutch NEO-PI-R / NEO-FFI alphas
                     ## (Hoekstra, Ormel & de Fruyt, 2007). The
                     ## Lifelines NEO release ships only the SUM*
                     ## totals, not the items, so we cannot compute
                     ## the Lifelines-specific alphas here. Replace
                     ## these with the Lifelines-computed values
                     ## once the items become available in the
                     ## extraction.
                     ##
                     ## Notes:
                     ##  (a) the Neuroticism phenotype used in this
                     ##      paper is the **full six-facet NEO-PI-R
                     ##      domain sum** (48 items), so the correct
                     ##      published reference is the NEO-PI-R
                     ##      Neuroticism domain alpha (~0.92,
                     ##      Hoekstra et al. 2007), NOT the 12-item
                     ##      NEO-FFI Neuroticism alpha (~0.85);
                     ##  (b) the Lifelines NEO-FFI Extraversion is
                     ##      11 items (one omitted by mistake), so
                     ##      the published 12-item alpha is a
                     ##      *ceiling* for the actual 11-item alpha.
                     0.92, 0.79, 0.81),
  Source = c(rep("Lifelines, this study", 8L),
             rep("Lifelines, this study (per wave)",
                 length(lte_alpha) + length(ldi_alpha)),
             "Hoekstra et al. 2007 (Dutch NEO-PI-R manual, 48-item domain alpha)",
             "Hoekstra et al. 2007 (Dutch NEO-FFI manual, 12-item ceiling for 11-item E)",
             "Hoekstra et al. 2007 (Dutch NEO-FFI manual)")
)

dir.create(OUT_TABLES, showWarnings = FALSE, recursive = TRUE)
save_table(internal_consistency, "table_S7e_internal_consistency")

message("Step 01a complete. N rows = ", nrow(phenotypes_lifelines),
        " | columns = ", ncol(phenotypes_lifelines),
        " | Suppl. Table S7e written.")
