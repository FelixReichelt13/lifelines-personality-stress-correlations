# 01_setup.R -- paths, packages, phenotypes (Lifelines pipeline)
# Felix Reichelt, UMCG (f.reichelt@umcg.nl)

# libraries

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

# reproducibility (inverse_normal tie-breaking too)
set.seed(12345)


# paths

PROJECT_ROOT <- normalizePath(".", winslash = "/")
CONFIG_DIR   <- file.path(PROJECT_ROOT, "01_Config")

DATA_DIR     <- file.path(PROJECT_ROOT, "data")
OUTPUT_DIR   <- file.path(PROJECT_ROOT, "output")

PEDIGREE_FILE                  <- file.path(DATA_DIR, "Relatives_with_project_id_M_and_F")
HOUSEHOLD_FILE                 <- file.path(DATA_DIR, "HouseID_Final_with_project_id")
DATA_FILE_BASELINE             <- file.path(DATA_DIR, "data_for_baseline.Rdata")
DATA_FILE_BASELINE_TRANSFORMED <- file.path(DATA_DIR, "data_for_baseline_transformed.Rdata")
PERSONALITY_FILE               <- file.path(DATA_DIR, "adult_personality_OUTPUT.csv")

# Raw questionnaire files used by 02a_score_construction.R. These are
RAW_CTQ_FILE      <- file.path(DATA_DIR, "raw", "CTQ_items.csv")
RAW_LL_FILE       <- file.path(DATA_DIR, "raw", "Loneliness_items.csv")
RAW_SS_FILE       <- file.path(DATA_DIR, "raw", "SocialSupport_items.csv")
RAW_MENTAL_FILE   <- file.path(DATA_DIR, "raw", "Mental_Health_items.csv")


# output folders

OUT_TABLES  <- file.path(OUTPUT_DIR, "tables")
save_figureURES <- file.path(OUTPUT_DIR, "figures")

for (d in c(OUTPUT_DIR, OUT_TABLES, save_figureURES)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}


# phenotypes used in the manuscript

STRESS_PHENOTYPES <- c(
  "Childhoodtrauma",         # CTQ-SF total score (Bernstein et al., 2003)
  "Loneliness",              # De Jong Gierveld Loneliness Scale
  "Lack_of_social_support",  # reversed SPF-IL
  "Stressful_life_events",   # LTE latent factor across waves
  "Life_difficulties"        # LDI latent factor across waves
)

PERSONALITY_PHENOTYPES <- c("Neuroticism", "Extraversion", "Conscientiousness")

# CTQ-SF subscale columns and moderate/severe cut-offs (Bernstein et al.,
# 2003) used for dichotomous CT in script 02 and liability h² in script 04bc.
# Primary ACE models (03–05) use the continuous CTQ-SF total (`Childhoodtrauma`).
CTQ_SUBSCALE_COLS <- c("ea_sum", "pa_sum", "sa_sum", "en_sum", "pn_sum")
CTQ_MODERATE_THRESHOLDS <- c(
  ea_sum = 12, pa_sum = 9, sa_sum = 7, en_sum = 14, pn_sum = 9
)


# plotting palette


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


# multiple-testing thresholds (ani et al. 2026 protocol)

UNIVARIATE_BONFERRONI_ALPHA <- 0.05 / 16
BIVARIATE_BONFERRONI_ALPHA  <- 0.05 / 75


# helper functions
source(file.path(CONFIG_DIR, "01_helpers.R"))

message("Setup loaded. Project root: ", PROJECT_ROOT)
