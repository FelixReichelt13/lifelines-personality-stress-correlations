#### =========================================================================
#### 02_data_preparation.R
#### -------------------------------------------------------------------------
#### Paper     : Reichelt F. et al. -- Familial Correlation between Personality
####             and Stress Exposure: Insights from the Multigenerational
####             Lifelines Cohort Study. European Journal of Personality
####             (manuscript EJP-26-4766, in revision).
#### Author    : Felix Reichelt  <f.reichelt@umcg.nl>
####
#### Merges pedigree, household and phenotype files; restricts to adults
#### (>= 18 yr); builds the LTE / LDI cross-wave latent factors via a one-
#### factor SEM in `lavaan` with FIML; creates the dichotomous CT variable
#### for liability / RRR analyses (script 04bc) from CTQ-SF **subscale**
#### moderate/severe thresholds (Bernstein et al., 2003) 
####
#### CT coding  : Childhood trauma is present when at least one CTQ-SF subscale
####              meets the literature-based moderate/severe threshold used in
####              this project and Bernstein et al. (2003): emotional abuse
####              >= 12, physical abuse >= 9, sexual abuse >= 7, emotional
####              neglect >= 14 or physical neglect >= 9. A participant is
####              coded 0 only when all five subscales are observed and none
####              meets its threshold.
####
#### Inputs : the five files defined in 01_setup.R.
#### Output : data/data_prepared.Rdata, with three objects:
####            - baseline      : analysis dataset (baseline scale)
####            - baseline_int  : analysis dataset (INT scale, used in 06)
####            - family_data   : pedigree + household join (used in 03-05)
####
#### Personality scoring (consumed from the upstream Lifelines NEO release;
#### see the Lifelines NEO wiki and Schweren 2019 UGLI scoring script):
####   - Neuroticism      = NEO-PI-R 48-item six-facet domain sum
####                        (range 48-240); built as the sum of the six
####                        facet `_SUM3` scores with >= 5 / 8 valid items
####                        per facet.
####   - Extraversion / Conscientiousness = NEO-FFI domain SUM3
####                        (range 12-60); E uses 11 administered items
####                        rescaled to a 12-item-equivalent via the
####                        Lifelines rule SUM3 = round((SUM2/NVALID)*12);
####                        C uses 12 items administered.
#### All personality variables are z-standardised before variance
#### decomposition (scripts 03, 05), so the raw-scale asymmetry does not
#### affect estimates.
####
#### Pitfalls:
####   - Pedigree IDs must be character / factor, never integer.
####   - Unknown parents must be coded NA, not 0 or "".
####   - The LTE / LDI factor needs at least two waves; otherwise it
####     falls back to the single observed wave.
#### Version history:
#### 2026-07-15  v1.2  CT dichotomous via subscale-specific Bernstein thresholds
#### 2026-07-15  v1.1  Ani-aligned multiple testing (/16 uni, /75 bivar);
####                   per-component 1-df bivariate LRTs; INTERNAL tables
#### 2026-07-15  v1.0  Initial release
#### =========================================================================


## Required libraries -----------------------------------------------------
##
## All are loaded centrally by `01_setup.R`. The block below re-loads them
## explicitly so that this script is also self-contained when read on its
## own (calling `library()` on an already-loaded package is a no-op).

library(data.table)   # fast `data.table` syntax (`[i, j, by]`, `set*` helpers)
library(readr)        # `read_csv()` for the personality file
library(lavaan)       # one-factor CFA for LTE / LDI cross-wave latent score
                      # (`sem()`, FIML missing-data handling, `predict()`
                      # to obtain factor scores for the ASReml step)


## Step 1 — Pedigree and household ---------------------------------------
##
## The pedigree file contains one row per Lifelines participant with the
## participant id (PSEUDOIDEXT) and the ids of the father and mother. The
## household file maps each participant to a HouseID that we use as the
## "shared environment" effect (c) in the ACE model.

load(PEDIGREE_FILE)                                          # adds `Relatives_with_project_id_M_and_F` to the workspace
load(HOUSEHOLD_FILE)                                         # adds `HouseID_Final_with_project_id`

pedigree  <- as.data.table(Relatives_with_project_id_M_and_F) # one row per participant: ID, FatherID, MotherID
household <- as.data.table(HouseID_Final_with_project_id)     # one row per participant: ID, HouseID1, FamID, ...

## Rename the first column to "ID" so it lines up with the pedigree, and
## the fourth column to "FamID" so script 04bc can cluster cluster-robust
## standard errors by family.
setnames(household,
         c(names(household)[1], names(household)[4]),
         c("ID",                "FamID"))

family_data <- merge(pedigree, household, by = "ID")          # left join: every participant gets its FamID + HouseID1


## Step 2 — Phenotypes ---------------------------------------------------
##
## `data_for_baseline.Rdata` holds the raw, untransformed phenotypes that
## go into the main analyses. `data_for_baseline_transformed.Rdata` is the
## same table after pre-computed inverse normal transformation; we re-use
## it in script 04bc (sensitivity analysis) so that the transformation
## itself does not have to be repeated. Personality scores live in a
## separate CSV file because they come from a later Lifelines release.

load(DATA_FILE_BASELINE)
load(DATA_FILE_BASELINE_TRANSFORMED)

baseline     <- as.data.table(data_for_baseline)
baseline_int <- as.data.table(data_for_baseline_transformed)
personality  <- as.data.table(read_csv(PERSONALITY_FILE, show_col_types = FALSE))


## Step 3 — Harmonise IDs, merge family data, centre age ----------------
##
## The same routine is applied to both the baseline and the INT dataset
## so that downstream scripts can treat them interchangeably.

prepare <- function(dt) {
  dt[, PSEUDOIDEXT := as.character(PSEUDOIDEXT)]                 # ASReml requires factor/character IDs (integer IDs are silently treated as covariates)
  out <- merge(dt, family_data,
               by.x = "PSEUDOIDEXT", by.y = "ID", all.x = TRUE)  # left join keeps every phenotype row even if pedigree is missing
  mean_age <- mean(out$AGE, na.rm = TRUE)
  out[, AGE_c  := AGE - mean_age]     # centred age, used as fixed effect
  out[, AGE_c2 := AGE_c^2]            # quadratic age, also fixed (captures non-linear age trends; Suppl. Table S7d)
  to_factor(out, c("PSEUDOIDEXT", "FatherID", "MotherID", "GENDER", "HouseID1"))  # ASReml needs these as factors
}

baseline     <- prepare(baseline)
baseline_int <- prepare(baseline_int)


## Step 3b — Attach CTQ subscale sums (for dichotomous CT, script 04bc) ---
##
## Subscales are built in `02a_score_construction.R`. If they are not
## already in the baseline extract, merge them from phenotypes_lifelines.Rdata.

if (!all(CTQ_SUBSCALE_COLS %in% names(baseline))) {
  ph_path <- file.path(DATA_DIR, "phenotypes_lifelines.Rdata")
  if (file.exists(ph_path)) {
    load(ph_path)   # adds `phenotypes_lifelines`
    phenotypes_lifelines[, ID := as.character(ID)]
    merge_cols <- c("ID", intersect(CTQ_SUBSCALE_COLS, names(phenotypes_lifelines)))
    baseline <- merge(baseline,
                      phenotypes_lifelines[, ..merge_cols],
                      by.x = "PSEUDOIDEXT", by.y = "ID", all.x = TRUE)
    baseline_int <- merge(baseline_int,
                          phenotypes_lifelines[, ..merge_cols],
                          by.x = "PSEUDOIDEXT", by.y = "ID", all.x = TRUE)
  } else {
    warning("CTQ subscale columns not in baseline and ", ph_path,
            " not found — CT_dichotomous will be NA.")
  }
}


## Step 4 — Attach personality (Neuroticism, Extraversion, Conscientiousness)

personality[, PSEUDOIDEXT := as.character(PSEUDOIDEXT)]
trait_cols   <- intersect(PERSONALITY_PHENOTYPES, names(personality))

baseline     <- merge(baseline,
                      personality[, c("PSEUDOIDEXT", trait_cols), with = FALSE],
                      by = "PSEUDOIDEXT", all.x = TRUE)
baseline_int <- merge(baseline_int,
                      personality[, c("PSEUDOIDEXT", trait_cols), with = FALSE],
                      by = "PSEUDOIDEXT", all.x = TRUE)


## Step 5 — Adult sample (Methods 2.1) ----------------------------------

baseline     <- baseline[AGE >= 18]
baseline_int <- baseline_int[AGE >= 18]


## Step 6 — Latent factors for LTE and LDI ------------------------------
##
## Both questionnaires were administered at multiple waves. The
## Primary analysis uses a single congeneric factor across waves
## (Figure S1). When fewer than two waves are present in the file we
## fall back to the single available wave, otherwise we fit the SEM
## with FIML so that we keep participants with partial response.

lte_waves <- intersect(
  c("LTE_Baseline", "LTE_Followup_1B", "LTE_Followup_1C",
    "LTE_Second_Assessment", "LTE_1A", "LTE_1B", "LTE_1C", "LTE_2A"),
  names(baseline))

ldi_waves <- intersect(
  c("LDI_Baseline", "LDI_Followup_1B", "LDI_Followup_1C",
    "LDI_Second_Assessment", "LDI_1A", "LDI_1B", "LDI_1C", "LDI_2A"),
  names(baseline))

extract_factor <- function(dt, indicators, name) {
  if (length(indicators) < 2) {
    if (length(indicators) == 1) dt[[name]] <- dt[[indicators]]   # only one wave found -> use it directly
    return(dt)
  }
  model <- paste0("F =~ ", paste(indicators, collapse = " + "))   # `lavaan` one-factor formula: F =~ wave1 + wave2 + ...
  fit   <- sem(model, data = dt, missing = "fiml", estimator = "ML")  # FIML keeps participants with at least one wave observed
  dt[[name]] <- as.numeric(predict(fit))                          # factor scores, then standardised to SD = 1 by the ASReml step
  dt
}

baseline <- extract_factor(baseline, lte_waves, "Stressful_life_events")
baseline <- extract_factor(baseline, ldi_waves, "Life_difficulties")


## Step 7 — Dichotomous childhood trauma (used by script 04bc) ------------
##
## CT coding: present when at least one subscale meets its moderate/severe
## threshold (Bernstein et al., 2003; see header). Coded 0 only when all
## five subscales are observed and none meets its threshold. This differs
## from the continuous `Childhoodtrauma` total used in primary ACE models
## (03–05). Subscale-specific liability h² in `07b` uses the same per-
## subscale thresholds, not total >= 40 applied to each subscale sum.

baseline[, CT_dichotomous := ct_dichotomous_moderate(baseline)]


## Step 8 — Save the analysis dataset -----------------------------------

save(baseline, baseline_int, family_data,
     file = file.path(DATA_DIR, "data_prepared.Rdata"))

message("Step 01 complete. N(baseline) = ", nrow(baseline),
        " | N(INT) = ",                     nrow(baseline_int))
