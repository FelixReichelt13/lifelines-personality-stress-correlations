#### =========================================================================
#### 04d_significance_audit.R
#### -------------------------------------------------------------------------
#### Paper     : Reichelt F. et al. -- Familial Correlation between Personality
####             and Stress Exposure: Insights from the Multigenerational
####             Lifelines Cohort Study. European Journal of Personality
####             (manuscript EJP-26-4766, in revision).
#### Author    : Felix Reichelt  <f.reichelt@umcg.nl>
####
#### What this script does
#### ---------------------
#### Reads the TSV files produced by scripts 04–05 and generates:
####
####   A. table_audit_univariate.tsv
####      Univariate ACE results with Bonferroni flags and LRT info.
####
####   B. table_audit_bivariate_full.tsv
####      All 25 bivariate pairings (10 stress×stress + 15 personality×
####      stress) with per-component 1-df LRT Bonferroni flags and
####      FDR-corrected q-values as sensitivity analysis.
####
####   C. table_audit_nonsig_components.tsv
####      Only the component-specific correlations that do NOT survive
####      Bonferroni correction (α = .05/75), for quick editorial reference.
####
####   D. table_audit_decomp_warnings.tsv
####      Bivariate pairings where the A/C/E decomposition is flagged
####      as unstable (small rp or out-of-range proportions).
####
####   E. audit_summary.txt
####      Human-readable summary: counts of sig / non-sig per component
####      and a plain-English interpretation guide.
####
#### Inputs (all under OUT_TABLES)
#### --------------------------------
####   table_S8_univariate_ace_full_INTERNAL.tsv   (script 04)
####   table_S8b_univariate_LRT.tsv                (script 04)
####   table_S13S15_bivariate_stress_stress_full_INTERNAL.tsv  (script 04b)
####   table_S14S16_bivariate_personality_stress_full_INTERNAL.tsv (script 04bc)
####
#### Multiple-testing thresholds used here
#### ----------------------------------------
####   Univariate: α_bonf = .05 / 16 = .003125 (8 phenotypes × 2 components
####               h²/c², parallel to Ani et al.; flags p_A and p_C < α)
####   Bivariate:  α_bonf = .05 / 75 = 6.67e-4 (25 pairings × 3 correlations;
####               each correlation tested by its own 1-df LRT, Ani et al.)
####   FDR (sensitivity, within bivariate component type):
####               Benjamini–Hochberg q < .05 computed separately for
####               p_rg, p_rc, p_re across all 25 pairings.
####
#### NOTE on philosophy
#### ------------------
#### The significance flags here are for INTERNAL REVIEW only.  Tables
#### that enter the manuscript supplement (S8, S13, S14, S15, S16) do NOT
#### include significance symbols.  The table notes state the Bonferroni
#### thresholds and direct the reader to this audit file.
#### Version history:
#### 2026-07-15  v1.1  Ani-aligned multiple testing (/16 uni, /75 bivar);
####                   per-component 1-df bivariate LRTs; INTERNAL tables
#### 2026-07-15  v1.0  Initial release
#### =========================================================================


suppressPackageStartupMessages({
  library(dplyr)
})


## ── helpers ────────────────────────────────────────────────────────────────

read_tsv_safe <- function(name) {
  f <- file.path(OUT_TABLES, paste0(name, ".tsv"))
  if (!file.exists(f)) {
    message("12_audit: file not found — ", f)
    return(NULL)
  }
  read.delim(f, stringsAsFactors = FALSE)
}

BONF_UNI  <- 0.05 / 16   # univariate: 8 phenotypes × 2 components (Ani et al.)
BONF_BIV  <- 0.05 / 75   # bivariate: 25 pairings × 3 per-component 1-df LRTs


## ── A. Univariate audit ────────────────────────────────────────────────────

uni <- read_tsv_safe("table_S8_univariate_ace_full_INTERNAL")

if (!is.null(uni)) {
  uni_audit <- uni %>%
    mutate(
      h2_pct      = round(h2 * 100, 1),
      c2_pct      = round(c2 * 100, 1),
      e2_pct      = round(e2 * 100, 1),
      h2_sig_bonf = p_A < BONF_UNI,      # additive genetic (ACE vs CE)
      c2_sig_bonf = p_C < BONF_UNI       # shared household (ACE vs AE)
    ) %>%
    select(Phenotype, N, h2_pct, h2_SE, c2_pct, c2_SE, e2_pct, e2_SE,
           p_A, p_C, h2_sig_bonf, c2_sig_bonf)

  write.table(uni_audit,
              file.path(OUT_TABLES, "table_audit_univariate.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  message("Univariate (alpha = ", signif(BONF_UNI, 3), ", 16 tests): h2 sig ",
          sum(uni_audit$h2_sig_bonf, na.rm = TRUE), " / ", nrow(uni_audit),
          "; c2 sig ", sum(uni_audit$c2_sig_bonf, na.rm = TRUE), " / ",
          nrow(uni_audit), ".")
}


## ── B. Bivariate full audit ────────────────────────────────────────────────

ss <- read_tsv_safe("table_S13S15_bivariate_stress_stress_full_INTERNAL")
ps <- read_tsv_safe("table_S14S16_bivariate_personality_stress_full_INTERNAL")

## Harmonise column names: stress×stress uses Phenotype1/Phenotype2,
## personality×stress uses Personality/Stress.
if (!is.null(ss)) {
  ss <- ss %>% rename(Trait1 = Phenotype1, Trait2 = Phenotype2)
  ss$Domain <- "Stress × Stress"
}
if (!is.null(ps)) {
  ps <- ps %>% rename(Trait1 = Personality, Trait2 = Stress)
  ps$Domain <- "Personality × Stress"
}

biv <- bind_rows(ss, ps)

if (nrow(biv) > 0) {

  ## Bonferroni flags (already computed in scripts 04b/04c, but recompute
  ## here for transparency in case the INTERNAL file was from an older run).
  ## Each correlation is tested by its own 1-df LRT (p_rg, p_rc, p_re) and
  ## corrected across 25 pairings × 3 correlations = 75 tests (Ani et al.).
  biv <- biv %>%
    mutate(
      rg_sig_bonf = p_rg < BONF_BIV,
      rc_sig_bonf = p_rc < BONF_BIV,
      re_sig_bonf = p_re < BONF_BIV
    )

  ## Benjamini–Hochberg FDR (sensitivity; within each component type
  ## across the 25 pairings; does NOT replace Bonferroni as primary).
  biv <- biv %>%
    mutate(
      rg_fdr = p.adjust(p_rg, method = "BH"),
      rc_fdr = p.adjust(p_rc, method = "BH"),
      re_fdr = p.adjust(p_re, method = "BH")
    )

  write.table(biv,
              file.path(OUT_TABLES, "table_audit_bivariate_full.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  message("Bivariate Bonferroni summary (α = ", signif(BONF_BIV, 3), ", 1-df LRT per component):")
  message("  rg : ", sum(biv$rg_sig_bonf, na.rm = TRUE), " / ", nrow(biv), " significant")
  message("  rc : ", sum(biv$rc_sig_bonf, na.rm = TRUE), " / ", nrow(biv), " significant")
  message("  re : ", sum(biv$re_sig_bonf, na.rm = TRUE), " / ", nrow(biv), " significant")


  ## ── C. Non-surviving components ──────────────────────────────────────────

  nonsig <- bind_rows(
    biv %>% filter(!rg_sig_bonf) %>%
      mutate(Component = "rg", Estimate = rg, SE = rg_SE, p = p_rg),
    biv %>% filter(!rc_sig_bonf) %>%
      mutate(Component = "rc", Estimate = rc, SE = rc_SE, p = p_rc),
    biv %>% filter(!re_sig_bonf) %>%
      mutate(Component = "re", Estimate = re, SE = re_SE, p = p_re)
  ) %>%
    select(Domain, Trait1, Trait2, Component, rp, Estimate, SE, p) %>%
    arrange(Domain, Component, Trait1, Trait2)

  write.table(nonsig,
              file.path(OUT_TABLES, "table_audit_nonsig_components.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  message("Non-surviving components written to table_audit_nonsig_components.tsv")
  message("  (", nrow(nonsig), " component-pairing combinations below α = ", signif(BONF_BIV, 3), ")")


  ## ── D. Decomposition warnings ────────────────────────────────────────────

  decomp_warn <- biv %>%
    filter(decomp_unstable) %>%
    select(Domain, Trait1, Trait2, rp, prop_A, prop_C, prop_E, decomp_unstable)

  write.table(decomp_warn,
              file.path(OUT_TABLES, "table_audit_decomp_warnings.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  message("Decomposition warnings: ", nrow(decomp_warn), " unstable pairings.")
  if (nrow(decomp_warn) > 0) {
    for (i in seq_len(nrow(decomp_warn))) {
      message("  ", decomp_warn$Trait1[i], " x ", decomp_warn$Trait2[i],
              "  rp=", round(decomp_warn$rp[i], 3),
              "  A=", round(decomp_warn$prop_A[i] * 100, 1), "%",
              "  C=", round(decomp_warn$prop_C[i] * 100, 1), "%",
              "  E=", round(decomp_warn$prop_E[i] * 100, 1), "%")
    }
  }


  ## ── E. Plain-English summary ──────────────────────────────────────────────

  n_ss_pairs <- sum(biv$Domain == "Stress × Stress")
  n_ps_pairs <- sum(biv$Domain == "Personality × Stress")

  ## Univariate counts (guarded: table_S8 INTERNAL may be absent).
  uni_h2_sig <- if (exists("uni_audit")) sum(uni_audit$h2_sig_bonf, na.rm = TRUE) else NA_integer_
  uni_c2_sig <- if (exists("uni_audit")) sum(uni_audit$c2_sig_bonf, na.rm = TRUE) else NA_integer_
  uni_n      <- if (exists("uni_audit")) nrow(uni_audit) else NA_integer_

  summary_lines <- c(
    "==============================================================",
    "SIGNIFICANCE AUDIT — EJP-26-4766 revision",
    sprintf("Generated: %s", Sys.time()),
    "==============================================================",
    "",
    "THRESHOLDS",
    sprintf("  Univariate: α = .05 / 16 = %.6f  (8 phenotypes × 2 components, Ani et al.)", BONF_UNI),
    sprintf("  Bivariate:  α = .05 / 75 = %.2e  (25 pairings × 3 correlations, Ani et al.)", BONF_BIV),
    "",
    "NOTE: following Ani et al. (2026, Diabetes Obes Metab), each of the",
    "three cross-trait correlations (rg, rc, re) is tested by its OWN 1-df",
    "likelihood-ratio test (full model vs a sub-model with that single",
    "covariance fixed to zero) and Bonferroni-corrected over all 75 tests.",
    "",
    "──────────────────────────────────────────────────────────────",
    "UNIVARIATE (Table S8)",
    "──────────────────────────────────────────────────────────────",
    sprintf("  h² sig: %s / %s ; c² sig: %s / %s  (after correction, α = %.6f).",
            uni_h2_sig, uni_n, uni_c2_sig, uni_n, BONF_UNI),
    "",
    "──────────────────────────────────────────────────────────────",
    sprintf("BIVARIATE: STRESS × STRESS (%d pairs, Table S13)", n_ss_pairs),
    "──────────────────────────────────────────────────────────────",
    "",
    sprintf("  Per-component 1-df LRTs (Bonferroni α = %.2e over 75 tests):", BONF_BIV),
    sprintf("  rg : %d / %d genetic correlations survive.",
            sum(biv$rg_sig_bonf[biv$Domain == "Stress × Stress"], na.rm = TRUE),
            n_ss_pairs),
    sprintf("  rc : %d / %d shared-env correlations survive.",
            sum(biv$rc_sig_bonf[biv$Domain == "Stress × Stress"], na.rm = TRUE),
            n_ss_pairs),
    sprintf("  re : %d / %d unique-env correlations survive.",
            sum(biv$re_sig_bonf[biv$Domain == "Stress × Stress"], na.rm = TRUE),
            n_ss_pairs),
    "",
    "  Non-surviving rc (stress × stress):",
    {
      nsc <- biv %>%
        filter(Domain == "Stress × Stress", !rc_sig_bonf) %>%
        mutate(msg = sprintf("    %s x %s  rc=%.3f (SE=%.3f, p=%.4g)",
                             Trait1, Trait2, rc, rc_SE, p_rc))
      if (nrow(nsc) == 0) "    (none)" else nsc$msg
    },
    "",
    "──────────────────────────────────────────────────────────────",
    sprintf("BIVARIATE: PERSONALITY × STRESS (%d pairs, Table S14)", n_ps_pairs),
    "──────────────────────────────────────────────────────────────",
    "",
    sprintf("  Per-component 1-df LRTs (Bonferroni α = %.2e over 75 tests):", BONF_BIV),
    sprintf("  rg : %d / %d genetic correlations survive.",
            sum(biv$rg_sig_bonf[biv$Domain == "Personality × Stress"], na.rm = TRUE),
            n_ps_pairs),
    sprintf("  rc : %d / %d shared-env correlations survive.",
            sum(biv$rc_sig_bonf[biv$Domain == "Personality × Stress"], na.rm = TRUE),
            n_ps_pairs),
    sprintf("  re : %d / %d unique-env correlations survive.",
            sum(biv$re_sig_bonf[biv$Domain == "Personality × Stress"], na.rm = TRUE),
            n_ps_pairs),
    "",
    "  Non-surviving rg (personality × stress):",
    {
      nsg <- biv %>%
        filter(Domain == "Personality × Stress", !rg_sig_bonf) %>%
        mutate(msg = sprintf("    %s x %s  rg=%.3f (SE=%.3f, p=%.4g)",
                             Trait1, Trait2, rg, rg_SE, p_rg))
      if (nrow(nsg) == 0) "    (none)" else nsg$msg
    },
    "",
    "  Non-surviving rc (personality × stress):",
    {
      nsc2 <- biv %>%
        filter(Domain == "Personality × Stress", !rc_sig_bonf) %>%
        mutate(msg = sprintf("    %s x %s  rc=%.3f (SE=%.3f, p=%.4g)",
                             Trait1, Trait2, rc, rc_SE, p_rc))
      if (nrow(nsc2) == 0) "    (none)" else nsc2$msg
    },
    "",
    "  Non-surviving re (personality × stress):",
    {
      nse <- biv %>%
        filter(Domain == "Personality × Stress", !re_sig_bonf) %>%
        mutate(msg = sprintf("    %s x %s  re=%.3f (SE=%.3f, p=%.4g)",
                             Trait1, Trait2, re, re_SE, p_re))
      if (nrow(nse) == 0) "    (none)" else nse$msg
    },
    "",
    "──────────────────────────────────────────────────────────────",
    "FDR SENSITIVITY (Benjamini-Hochberg q < .05)",
    "──────────────────────────────────────────────────────────────",
    sprintf("  rg (FDR): %d / %d survive",
            sum(biv$rg_fdr < .05, na.rm = TRUE), nrow(biv)),
    sprintf("  rc (FDR): %d / %d survive",
            sum(biv$rc_fdr < .05, na.rm = TRUE), nrow(biv)),
    sprintf("  re (FDR): %d / %d survive",
            sum(biv$re_fdr < .05, na.rm = TRUE), nrow(biv)),
    "",
    "  Note: FDR is reported as sensitivity only.  The primary correction",
    "  is Bonferroni over all 75 per-component 1-df LRTs (α = .05/75).",
    "",
    sprintf("DECOMPOSITION WARNINGS: %d pairings flagged as unstable.", nrow(decomp_warn)),
    "  (see table_audit_decomp_warnings.tsv)")

  writeLines(summary_lines,
             file.path(OUT_TABLES, "audit_summary.txt"))

  message("Audit summary written to audit_summary.txt")
}

message("Step 12 complete. Internal audit files written to ", OUT_TABLES)
