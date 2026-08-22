#### =========================================================================
#### 04a_model_comparison_long_tables.R
#### -------------------------------------------------------------------------
#### Paper     : Reichelt F. et al. -- Familial Correlation between Personality
####             and Stress Exposure (EJP-26-4766, in revision).
#### Author    : Felix Reichelt  <f.reichelt@umcg.nl>
####
#### What this script does
#### ---------------------
#### Re-shapes the WIDE univariate model-comparison outputs from scripts 04
#### and 05 into the LONG "one row per model" layout in the manuscript:
####
####     Phenotype | Model | logLik | AIC | BIC | LR vs ACE | p
####       <name>    ACE       ..      ..    ..     —          —
####                 CE (tests A)  ..   ..    ..     LR_A       p_A
####                 AE (tests C)  ..   ..    ..     LR_C       p_C
####
#### Produces
####   Table S9  : eight main phenotypes (raw scale), from table_S8c + S8b.
####   Table S12 : INT phenotypes (mirrors S9), from table_S11c + S11b.
####   Table S11 : LTE/LDI per wave + latent factor (optional), from
####               table_S10c + S10b per-wave files when present.
####
#### Inputs    : table_S8b/S8c (script 04); table_S11b/S11c (script 05);
####             optional table_S10b/S10c per-wave (if per-wave ACE restored).
#### Outputs   : table_S9_univariate_model_comparison_long.tsv
####             table_S12_int_model_comparison_long.tsv
####             table_S11_perwave_model_comparison_long.tsv (optional)
####
#### This script does NOT call ASReml; it only reads TSVs. Run after 05 in
#### run_all.R so Table S12 can be built from the INT wide tables.
#### =========================================================================


## ---- Formatting helpers -------------------------------------------------

fmt_num <- function(x, d = 2)
  ifelse(is.na(x), "", formatC(as.numeric(x), format = "f", digits = d))

fmt_lr <- function(x)
  ifelse(is.na(x), "\u2014", formatC(as.numeric(x), format = "f", digits = 2))

fmt_p <- function(p) {
  p <- as.numeric(p)
  ifelse(is.na(p), "\u2014",
         ifelse(p < .001, "< .001",
                formatC(p, format = "f", digits = 3)))
}

read_tsv_req <- function(name) {
  f <- file.path(OUT_TABLES, paste0(name, ".tsv"))
  if (!file.exists(f)) stop("Required input not found: ", f)
  read.delim(f, stringsAsFactors = FALSE, check.names = FALSE)
}


## ---- Core reshape -------------------------------------------------------
##
## `ic_row`  must contain logLik_ACE/AE/CE, AIC_ACE/AE/CE, BIC_ACE/AE/CE.
## `lrt_row` must contain LR_ACE_vs_CE, p_ACE_vs_CE, LR_ACE_vs_AE, p_ACE_vs_AE.

three_model_rows <- function(ic_row, lrt_row, label) {
  data.frame(
    Phenotype = c(label, "", ""),
    Model     = c("ACE", "CE (tests A)", "AE (tests C)"),
    logLik    = fmt_num(c(ic_row$logLik_ACE, ic_row$logLik_CE, ic_row$logLik_AE)),
    AIC       = fmt_num(c(ic_row$AIC_ACE,    ic_row$AIC_CE,    ic_row$AIC_AE)),
    BIC       = fmt_num(c(ic_row$BIC_ACE,    ic_row$BIC_CE,    ic_row$BIC_AE)),
    `LR vs ACE` = c("\u2014",
                    fmt_lr(lrt_row$LR_ACE_vs_CE),
                    fmt_lr(lrt_row$LR_ACE_vs_AE)),
    p         = c("\u2014",
                  fmt_p(lrt_row$p_ACE_vs_CE),
                  fmt_p(lrt_row$p_ACE_vs_AE)),
    check.names = FALSE, stringsAsFactors = FALSE
  )
}


## =========================================================================
## Table S9 -- eight main phenotypes (raw scale)
## =========================================================================

s8c <- read_tsv_req("table_S8c_univariate_model_comparison")
s8b <- read_tsv_req("table_S8b_univariate_LRT")
main <- merge(s8c, s8b, by = "Phenotype", sort = FALSE)

main_labels <- c(
  Conscientiousness      = "Conscientiousness",
  Extraversion           = "Extraversion",
  Neuroticism            = "Neuroticism",
  Childhoodtrauma        = "Childhood Trauma",
  Life_difficulties      = "Latent life difficulties",
  Stressful_life_events  = "Latent stressful life events",
  Lack_of_social_support = "Lack of social support",
  Loneliness             = "Loneliness"
)
main_groups <- c(
  Conscientiousness      = "Personality traits",
  Extraversion           = "Personality traits",
  Neuroticism            = "Personality traits",
  Childhoodtrauma        = "Stress exposures",
  Life_difficulties      = "Stress exposures",
  Stressful_life_events  = "Stress exposures",
  Lack_of_social_support = "Stress exposures",
  Loneliness             = "Stress exposures"
)

s9_blocks <- lapply(names(main_labels), function(code) {
  r <- main[main$Phenotype == code, ]
  if (nrow(r) == 0) {
    message("Table S9: phenotype not found in inputs, skipping -> ", code)
    return(NULL)
  }
  out <- three_model_rows(r, r, main_labels[[code]])
  cbind(Group = c(main_groups[[code]], "", ""), out)
})
s9 <- do.call(rbind, Filter(Negate(is.null), s9_blocks))

save_table(s9, "table_S9_univariate_model_comparison_long")


## =========================================================================
## Table S12 -- INT phenotypes (mirrors S9)
## =========================================================================

s11c_path <- file.path(OUT_TABLES, "table_S11c_int_univariate_model_comparison.tsv")
s11b_path <- file.path(OUT_TABLES, "table_S11b_int_univariate_LRT.tsv")

if (file.exists(s11c_path) && file.exists(s11b_path)) {
  s11c <- read.delim(s11c_path, stringsAsFactors = FALSE, check.names = FALSE)
  s11b <- read.delim(s11b_path, stringsAsFactors = FALSE, check.names = FALSE)
  int_merged <- merge(s11c, s11b, by = "Phenotype", sort = FALSE)

  s12_blocks <- lapply(int_merged$Phenotype, function(code) {
    r <- int_merged[int_merged$Phenotype == code, ]
    lab <- paste0(gsub("_INT$", "", code), " (INT)")
    three_model_rows(r, r, lab)
  })
  s12 <- do.call(rbind, s12_blocks)
  save_table(s12, "table_S12_int_model_comparison_long")
} else {
  message("Table S12: table_S11b/S11c not found -- run 05_sensitivity_int.R first.")
}


## =========================================================================
## Table S11 -- LTE and LDI per wave + latent factor (optional)
## =========================================================================

perwave_ic  <- file.path(OUT_TABLES, "table_S10c_perwave_lte_ldi_model_comparison.tsv")
perwave_lrt <- file.path(OUT_TABLES, "table_S10b_perwave_lte_ldi_LRT.tsv")

if (file.exists(perwave_ic) && file.exists(perwave_lrt)) {
  wave_labels <- c(
    Baseline_1A      = "Baseline measurement (1A)",
    Followup_1B      = "Follow up measurement (1B)",
    Followup_1C      = "Second follow up measurement (1C)",
    Second_Assess_2A = "Second Assessment (2A)"
  )
  inv_labels <- c(
    LTE = "List of threatening experiences (LTE)",
    LDI = "Long-term difficulties (LDI)"
  )
  latent_code <- c(LTE = "Stressful_life_events", LDI = "Life_difficulties")

  s10c <- read.delim(perwave_ic,  stringsAsFactors = FALSE, check.names = FALSE)
  s10b <- read.delim(perwave_lrt, stringsAsFactors = FALSE, check.names = FALSE)
  pw   <- merge(s10c, s10b, by = c("Inventory", "Wave", "Column"), sort = FALSE)

  build_instrument <- function(inv) {
    blocks <- list()
    for (w in names(wave_labels)) {
      r <- pw[pw$Inventory == inv & pw$Wave == w, ]
      if (nrow(r) == 0) next
      lab <- paste0(inv_labels[[inv]], " \u2014 ", wave_labels[[w]])
      blocks[[length(blocks) + 1]] <- three_model_rows(r, r, lab)
    }
    lc <- latent_code[[inv]]
    rl <- main[main$Phenotype == lc, ]
    if (nrow(rl) > 0) {
      lab <- paste0(inv_labels[[inv]],
                    " \u2014 Latent factor score based on 1A+1B+1C+2A")
      blocks[[length(blocks) + 1]] <- three_model_rows(rl, rl, lab)
    }
    do.call(rbind, blocks)
  }

  s11_pw <- do.call(rbind, lapply(c("LTE", "LDI"), build_instrument))
  save_table(s11_pw, "table_S11_perwave_model_comparison_long")
} else {
  message("Table S11 (per-wave): S10b/S10c per-wave inputs not found -- skipped.")
}

message("04a_model_comparison_long_tables complete.")
