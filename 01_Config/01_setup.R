#### =========================================================================
#### 01_setup.R
#### -------------------------------------------------------------------------
#### Paper     : Reichelt F. et al. -- Familial Correlation between Personality
####             and Stress Exposure: Insights from the Multigenerational
####             Lifelines Cohort Study. European Journal of Personality
####             (manuscript EJP-26-4766, in revision).
#### Author    : Felix Reichelt  <f.reichelt@umcg.nl>
####
#### Purpose   : Central configuration for paths, analysis variables and
####             multiple-testing thresholds used across the pipeline.
####             It loads the packages, fixes the random seed,
####             builds the output folders and defines two short vectors with
####             the names of the phenotypes used in the manuscript.
####
#### Dependencies & reproducibility: Individual-level Lifelines data are restricted 
####             and are not included in this repository. File names below
####             refer to the approved Lifelines extracts used for the manuscript.
####
#### Notes     : - ASReml-R is a commercial package
####               (https://vsni.co.uk/software/asreml-r/). Scripts 04-05
####               will not run without it. Folders 03, 06 and 07 only
####               need CRAN packages.
####             - All file paths below are written for the Lifelines HPC
####               environment. Edit them once and the rest of the pipeline
####               picks the new values up automatically.
####             - Tested with R >= 4.3 and ASReml-R 4.2 (Butler et al.,
####               2023). None of the analyses depend on a specific minor
####               version.
####
#### Version history:
#### 2026-07-15  v1.3  Script numbers aligned with folder numbers; 06_CT_Liability
#### 2026-07-15  v1.2  Numbered topic folders; script-by-script workflow
#### 2026-07-15  v1.1  Ani-aligned multiple testing (/16 uni, /75 bivar);
####                   per-component 1-df bivariate LRTs; INTERNAL tables
#### 2026-07-15  v1.0  Initial release
#### =========================================================================


## 1. Libraries -----------------------------------------------------------
##    Loaded in roughly the order in which they appear later in the pipeline.
##    ASReml-R is loaded first so that we fail fast if the license is missing.

library(asreml)          # commercial; see https://vsni.co.uk/software/asreml-r/
library(data.table)
library(dplyr)
library(tidyr)
library(readr)
library(haven)
library(lavaan)          # LTE / LDI latent factors with FIML
library(marginaleffects) # used in 07 for the RRR predictions
library(sandwich)        # cluster-robust SEs for the RRR
library(ggplot2)
library(patchwork)
library(scales)

## 2. Reproducibility -----------------------------------------------------
##    Used by inverse_normal() in 01_helpers.R. Any other random step should
##    set its own seed *inside* the function, not rely on this global one.

set.seed(12345)


## 3. Paths ---------------------------------------------------------------
##    Set working directory to the repository root, then:
##    source("01_Config/01_setup.R")

PROJECT_ROOT <- normalizePath(".", winslash = "/")
CONFIG_DIR   <- file.path(PROJECT_ROOT, "01_Config")

DATA_DIR     <- file.path(PROJECT_ROOT, "data")
OUTPUT_DIR   <- file.path(PROJECT_ROOT, "output")

PEDIGREE_FILE                  <- file.path(DATA_DIR, "Relatives_with_project_id_M_and_F")
HOUSEHOLD_FILE                 <- file.path(DATA_DIR, "HouseID_Final_with_project_id")
DATA_FILE_BASELINE             <- file.path(DATA_DIR, "data_for_baseline.Rdata")
DATA_FILE_BASELINE_TRANSFORMED <- file.path(DATA_DIR, "data_for_baseline_transformed.Rdata")
PERSONALITY_FILE               <- file.path(DATA_DIR, "adult_personality_OUTPUT.csv")

## Raw questionnaire files used by 02a_score_construction.R. These are
## the item-level CSVs distributed by Lifelines (CTQ-SF, Loneliness,
## SPF-IL Social Support, Mental Health for LTE/LDI). If you already
## have `data/phenotypes_lifelines.Rdata` from a previous run, you can
## skip 02a and start the pipeline at 02_data_preparation.R.
## Item-level CSVs for 02a_score_construction.R and Cronbach's alpha (S7e).
## If data/raw/ is absent, skip 02a or use pre-built phenotypes_lifelines.Rdata.
RAW_CTQ_FILE      <- file.path(DATA_DIR, "raw", "CTQ_items.csv")
RAW_LL_FILE       <- file.path(DATA_DIR, "raw", "Loneliness_items.csv")
RAW_SS_FILE       <- file.path(DATA_DIR, "raw", "SocialSupport_items.csv")
RAW_MENTAL_FILE   <- file.path(DATA_DIR, "raw", "Mental_Health_items.csv")


## 4. Output folders ------------------------------------------------------

OUT_TABLES  <- file.path(OUTPUT_DIR, "tables")
OUT_FIGURES <- file.path(OUTPUT_DIR, "figures")

for (d in c(OUTPUT_DIR, OUT_TABLES, OUT_FIGURES)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}


## 5. Phenotypes used in the manuscript ----------------------------------
##    Order here is the order used in every figure and table of the paper.
##    See Methods 2.2 and Supplementary Table S1 for the questionnaire-level
##    definitions of each phenotype.

STRESS_PHENOTYPES <- c(
  "Childhoodtrauma",         # CTQ-SF total score (Bernstein et al., 2003)
  "Loneliness",              # De Jong Gierveld Loneliness Scale
  "Lack_of_social_support",  # reversed SPF-IL
  "Stressful_life_events",   # LTE latent factor across waves
  "Life_difficulties"        # LDI latent factor across waves
)

PERSONALITY_PHENOTYPES <- c("Neuroticism", "Extraversion", "Conscientiousness")

## CTQ-SF subscale columns and moderate/severe cut-offs (Bernstein et al.,
## 2003) used for dichotomous CT in script 02 and liability h² in script 04bc.
## Primary ACE models (03–05) use the continuous CTQ-SF total (`Childhoodtrauma`).
CTQ_SUBSCALE_COLS <- c("ea_sum", "pa_sum", "sa_sum", "en_sum", "pn_sum")
CTQ_MODERATE_THRESHOLDS <- c(ea_sum = 12, pa_sum = 9, sa_sum = 7, en_sum = 14, pn_sum = 9)


## 6. Plotting palette ---------------------------------------------------
##    Kept here so every figure in script 04bc is colour-consistent.
##    Wes Anderson inspired palette for ACE (variance components) and correlation plots. 


ACE_COLOURS <- c(
  h2 = "#4C9FBC",  # additive genetic
  c2 = "#A8D5BA",  # shared environment
  e2 = "#F4A261"   # unique environment
)
COR_COLOURS <- c(
  rp = "#8E5E9B",  # phenotypic
  rg = "#4C9FBC",  # genetic
  rc = "#A8D5BA",  # shared environment
  re = "#F4A261"   # unique environment
)


## 7. Multiple-testing thresholds (Ani et al. 2026 protocol) ------------
##    Univariate: 8 phenotypes × 2 tested components (h², c²) = 16 tests.
##    Bivariate:  25 pairings × 3 correlations (rg, rc, re) = 75 tests.
##    e² is not tested univariately; rP is estimate + SE + CI only (no LRT).

UNIVARIATE_BONFERRONI_ALPHA <- 0.05 / 16
BIVARIATE_BONFERRONI_ALPHA  <- 0.05 / 75


## 8. Helper functions ---------------------------------------------------
source(file.path(CONFIG_DIR, "01_helpers.R"))

message("Setup loaded. Project root: ", PROJECT_ROOT)
