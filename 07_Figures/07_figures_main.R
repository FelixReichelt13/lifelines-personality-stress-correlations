#### =========================================================================
#### 07_figures_main.R
#### -------------------------------------------------------------------------
#### Paper     : Reichelt F. et al. -- Familial Correlation between Personality
####             and Stress Exposure: Insights from the Multigenerational
####             Lifelines Cohort Study. European Journal of Personality
####             (manuscript EJP-26-4766, in revision).
#### Author    : Felix Reichelt  <f.reichelt@umcg.nl>
####
#### What this script does
#### ---------------------
#### Builds the three main-text figures of the paper. **The script is
#### standalone**: all numbers used by the figures are encoded in the
#### data.frames near the top of the script. Source the file and you
#### get the same Figures 1-3 you see in the paper.
####
#### Figure mapping
#### --------------
####   figure_1_univariate_ace          Figure 1
####     stacked ACE bars (h2, c2, e2) per phenotype on the *baseline*
####     scale (5 stress phenotypes + 3 personality traits).
####
####   figure_2_stress_stress           Figure 2
####     5x5 panel grid for the stress phenotypes:
####       - lower triangle  = signed correlations (rp, rg, rc, re) +/- SE
####       - upper triangle  = unsigned ACE magnitude proportions (%) that
####                            decompose rp into its A, C and E parts
####       - diagonal        = empty
####
####   figure_3_personality_stress      Figure 3
####     5 (stress) x 3 (personality) grid, where each personality trait
####     gets *two* rows: a Correlation row (rp, rg, rc, re +/- SE) on a
####     [-1, +1] axis and a Magnitude row showing |A|, |C|, |E| as a
####     percentage of rp on a [0, 100] axis.
####
#### Outputs   : output/figures/figure_*.{png,pdf}
####
#### License   : This script does NOT require ASReml-R. CRAN packages
####             only: ggplot2, patchwork, dplyr, tidyr, scales,
####             grid, gtable.
#### Version history:
#### 2026-07-15  v1.1  Ani-aligned multiple testing (/16 uni, /75 bivar);
####                   per-component 1-df bivariate LRTs; INTERNAL tables
#### 2026-07-15  v1.0  Initial release
#### =========================================================================


## Required libraries -----------------------------------------------------
##
## `01_setup.R` already loads ggplot2, patchwork, scales, dplyr, tidyr.
## `grid` + `gtable` are figure-specific (used to add the right-hand
## "Correlation / Magnitude" strip column on Figure 3) and are loaded
## here so the script is self-contained when sourced directly.

library(ggplot2)      # all figure geoms / themes / scales
library(patchwork)    # `plot_baseline / plot_int` vertical composition + `plot_layout()`
library(scales)       # `percent_format()`, `number_format()` axis labellers
library(dplyr)        # `recode()`
library(tidyr)        # `pivot_longer()` (long -> wide for `geom_col` faceting)
library(grid)         # `textGrob()`, `gpar()`, `unit()` for the right-hand strip on Figure 3
library(gtable)       # `gtable_add_cols()`, `gtable_add_grob()` to modify the ggplot grob


## -----------------------------------------------------------------------
## 0. Embedded manuscript numbers
## -----------------------------------------------------------------------
## Univariate ACE point estimates with delta-method SEs (Table S8
## baseline, Table S10 INT). The stacked-bar figure shows only
## the proportions h^2 + c^2 + e^2 = 1; the SEs are listed alongside
## so the data.frame can be re-used to render error bars or to verify
## numbers against the tables
##
## Provenance: ASReml-R univariate ACE -- variance components extracted
## with `vpredict(model, h2 ~ V2 / (V1 + V2 + V3))` etc. (Reference
## Manual 4.2, Section 4.6). The SEs are the delta-method standard
## errors that `vpredict()` returns alongside each `Estimate`. The
## production run is `04_univariate_heritability.R` (Tables S8 / S10);
## numbers below are copied from that output.
##
## INT row covers only the five stress phenotypes (personality traits
## are not INT-transformed because they are already approximately
## normal -- Suppl. Fig. S17); the personality entries on the INT
## scale duplicate the baseline values so the data.frame stays
## rectangular.

univariate_baseline <- data.frame(
  Phenotype = c("Childhoodtrauma", "Life_difficulties", "Stressful_life_events",
                "Lack_of_social_support", "Loneliness",
                "Neuroticism", "Extraversion", "Conscientiousness"),
  Domain    = c(rep("Stress exposures", 5),
                rep("Personality traits", 3)),
  h2    = c(0.507, 0.223, 0.221, 0.176, 0.141,
            0.292, 0.279, 0.229),
  h2_SE = c(0.015, 0.008, 0.008, 0.008, 0.015,
            0.018, 0.018, 0.018),
  c2    = c(0.101, 0.324, 0.290, 0.207, 0.209,
            0.092, 0.085, 0.093),
  c2_SE = c(0.010, 0.005, 0.005, 0.005, 0.009,
            0.010, 0.010, 0.010),
  e2    = c(0.392, 0.453, 0.489, 0.617, 0.650,
            0.617, 0.636, 0.679),
  e2_SE = c(0.005, 0.005, 0.005, 0.003, 0.005,
            0.010, 0.010, 0.010),
  stringsAsFactors = FALSE)

univariate_int <- data.frame(
  Phenotype = c("Childhoodtrauma", "Life_difficulties", "Stressful_life_events",
                "Lack_of_social_support", "Loneliness",
                "Neuroticism", "Extraversion", "Conscientiousness"),
  Domain    = c(rep("Stress exposures", 5),
                rep("Personality traits", 3)),
  h2    = c(0.411, 0.237, 0.245, 0.175, 0.139,
            0.291, 0.275, 0.230),
  h2_SE = c(0.015, 0.008, 0.007, 0.008, 0.015,
            0.018, 0.018, 0.018),
  c2    = c(0.089, 0.320, 0.264, 0.204, 0.199,
            0.095, 0.085, 0.093),
  c2_SE = c(0.009, 0.005, 0.005, 0.005, 0.009,
            0.010, 0.010, 0.010),
  e2    = c(0.500, 0.443, 0.491, 0.621, 0.662,
            0.614, 0.640, 0.677),
  e2_SE = c(0.005, 0.005, 0.005, 0.003, 0.005,
            0.010, 0.010, 0.010),
  stringsAsFactors = FALSE)

## Bivariate stress x stress. Both directions of each pair are listed
## because the figure shows the 5x5 grid; SEs are from Table S13.

stress_stress <- read.table(text = "
Exp1                  Exp2                    rp     rp_se   rg     rg_se   rc     rc_se   re     re_se
Childhoodtrauma       Life_difficulties        0.33  0.0036   0.66  0.0220   0.61  0.0349   0.03  0.0142
Childhoodtrauma       Loneliness               0.26  0.0037   0.58  0.0428   0.62  0.0476   0.05  0.0216
Childhoodtrauma       Lack_of_social_support  -0.20  0.0039  -0.43  0.0273  -0.33  0.0392  -0.04  0.0141
Childhoodtrauma       Stressful_life_events    0.18  0.0039   0.40  0.0234   0.34  0.0319 -0.05  0.0171
Life_difficulties     Loneliness               0.30  0.0037   0.53  0.0435   0.52  0.0213   0.15  0.0112
Life_difficulties     Lack_of_social_support  -0.22  0.0026  -0.32  0.0267  -0.33  0.0145  -0.14  0.0066
Life_difficulties     Stressful_life_events    0.43  0.0022   0.51  0.0199   0.55  0.0095   0.31  0.0061
Loneliness            Lack_of_social_support  -0.36  0.0034  -0.85  0.0495  -0.61  0.0238  -0.16  0.0098
Loneliness            Stressful_life_events    0.12  0.0039   0.32  0.0455   0.22  0.0212   0.02  0.0130
Lack_of_social_support Stressful_life_events  -0.08  0.0027  -0.10  0.0278  -0.12  0.0144   0.05  0.0074
", header = TRUE, stringsAsFactors = FALSE)

## Mirror to get both directions for the 5x5 grid.
stress_stress_full <- rbind(
  stress_stress,
  setNames(stress_stress[, c("Exp2", "Exp1",
                              "rp", "rp_se", "rg", "rg_se",
                              "rc", "rc_se", "re", "re_se")],
           names(stress_stress)))

## Bivariate personality x stress (Figure 3 in the paper, Suppl. Table S14).

personality_stress <- read.table(text = "
Trait              Stress                  rp     rp_se   rg     rg_se   rc     rc_se   re     re_se
Conscientiousness  Childhoodtrauma         -0.14  0.0075  -0.18  0.0458  -0.44  0.1087  -0.06  0.0324
Conscientiousness  Life_difficulties       -0.19  0.0045  -0.23  0.0449  -0.59  0.0432  -0.07  0.0179
Conscientiousness  Loneliness              -0.23  0.0071  -0.21  0.0795  -0.42  0.0705  -0.20  0.0241
Conscientiousness  Lack_of_social_support  -0.37  0.0040  -0.51  0.0432  -0.63  0.0422  -0.28  0.0148
Conscientiousness  Stressful_life_events   -0.08  0.0047   0.01  0.0456  -0.17  0.0385  -0.09  0.0186
Extraversion       Childhoodtrauma         -0.15  0.0074  -0.37  0.0407  -0.13  0.1155   0.00  0.0336
Extraversion       Life_difficulties       -0.20  0.0045  -0.35  0.0399  -0.41  0.0430  -0.09  0.0185
Extraversion       Loneliness              -0.32  0.0066  -0.52  0.0702  -0.46  0.0706  -0.21  0.0241
Extraversion       Lack_of_social_support  -0.44  0.0037  -0.68  0.0361  -0.72  0.0414  -0.30  0.0147
Extraversion       Stressful_life_events   -0.04  0.0048   0.00  0.0407  -0.11  0.0407  -0.03  0.0194
Neuroticism        Childhoodtrauma          0.23  0.0072   0.41  0.0399   0.44  0.1033   0.06  0.0335
Neuroticism        Life_difficulties        0.41  0.0038   0.49  0.0343   0.67  0.0331   0.31  0.0165
Neuroticism        Loneliness               0.32  0.0066   0.61  0.0703   0.58  0.0658   0.18  0.0245
Neuroticism        Lack_of_social_support   0.36  0.0040   0.54  0.0394   0.56  0.0394   0.26  0.0156
Neuroticism        Stressful_life_events    0.21  0.0045   0.17  0.0392   0.42  0.0375   0.17  0.0189
", header = TRUE, stringsAsFactors = FALSE)


## -----------------------------------------------------------------------
## 1. Figure 1 -- Univariate ACE (baseline + INT, Stress | Personality)
## -----------------------------------------------------------------------
## Layout taken verbatim from NEWest_baseline_int_plus_perso_11_25.R:
##   - one panel for the five Stress exposures, one for the three
##     Personality traits (facet_grid(~ Domain, scales = "free_x"))
##   - bars are stacked bottom-to-top in the order Ve / Vc / h2
##   - in-bar percentage labels.

pheno_order <- c(
  ## stress
  "Childhoodtrauma", "Life_difficulties", "Stressful_life_events",
  "Lack_of_social_support", "Loneliness",
  ## personality
  "Neuroticism", "Extraversion", "Conscientiousness")
domain_order <- c("Stress exposures", "Personality traits")
stack_order  <- c("e2", "c2", "h2")  # bottom -> top

ace_long <- function(df, label) {
  d <- pivot_longer(df, c("h2", "c2", "e2"),
                     names_to  = "Variance_Component",
                     values_to = "Contribution")
  d <- as.data.frame(d)
  d$Variance_Component <- factor(d$Variance_Component, levels = stack_order)
  d$Phenotype <- factor(d$Phenotype, levels = pheno_order)
  d$Domain    <- factor(d$Domain,    levels = domain_order)
  d$Scale     <- label
  d
}

uni_long_raw <- ace_long(univariate_baseline, "Baseline (raw)")
uni_long_int <- ace_long(univariate_int,      "Inverse normal transformed")

ace_labels <- c(
  "h2" = expression(Genetic~(h^2)),
  "c2" = expression(Shared~(c^2)),
  "e2" = expression(Unique~(e^2)))

ace_bar_plot <- function(d, title) {
  ggplot(d, aes(Phenotype, Contribution, fill = Variance_Component)) +
    geom_col() +
    geom_text(aes(label = sprintf("%.1f%%", 100 * Contribution)),
              position = position_stack(vjust = 0.5),
              size = 4, fontface = "bold") +
    facet_grid(~ Domain, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = c("h2" = ACE_COLOURS[["h2"]],
                                  "c2" = ACE_COLOURS[["c2"]],
                                  "e2" = ACE_COLOURS[["e2"]]),
                       labels = ace_labels,
                       name   = "Variance components:") +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(x = NULL, y = "Contribution (%)", title = title,
          caption = "h2 = Va/Vp, c2 = Vc/Vp, e2 = Ve/Vp") +
    theme_minimal(base_size = 13) +
    theme(plot.title      = element_text(hjust = 0.5, face = "bold"),
          plot.caption    = element_text(hjust = 0.5, size = 8.5),
          axis.text.x     = element_text(angle = 50, hjust = 1, size = 10.5),
          strip.text      = element_text(face = "bold", size = 12.5),
          legend.position = "bottom")
}

plot_baseline <- ace_bar_plot(uni_long_raw, "Baseline data")
plot_int      <- ace_bar_plot(uni_long_int, "Inverse normal transformed (INT)")

fig1 <- plot_baseline / plot_int +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

save_figure(fig1, "figure_1_univariate_ace", width = 12, height = 10)


## -----------------------------------------------------------------------
## 2. Figure 2 -- Stress x Stress 5x5
## -----------------------------------------------------------------------
## Lower triangle  : signed correlations rp, rg, rc, re +/- SE.
## Upper triangle  : unsigned A/C/E magnitude proportions (decomposition
##                   of rp into the part coming from the genetic, the
##                   shared-environmental and the unique-environmental
##                   covariance). The decomposition uses the univariate
##                   h^2, c^2 and e^2 in `univariate_baseline`.
## Diagonal         : empty.

estimate_colors <- c("rp" = COR_COLOURS[["rp"]],
                      "rg" = COR_COLOURS[["rg"]],
                      "rc" = COR_COLOURS[["rc"]],
                      "re" = COR_COLOURS[["re"]])

stress_levels <- c("Childhoodtrauma", "Life_difficulties",
                    "Stressful_life_events", "Lack_of_social_support",
                    "Loneliness")

ss <- stress_stress_full
ss$Row <- factor(ss$Exp1, levels = stress_levels)
ss$Col <- factor(ss$Exp2, levels = stress_levels)

## Look-up vectors for the univariate proportions.
h2_s <- setNames(univariate_baseline$h2[univariate_baseline$Phenotype %in% stress_levels],
                  univariate_baseline$Phenotype[univariate_baseline$Phenotype %in% stress_levels])
c2_s <- setNames(univariate_baseline$c2[univariate_baseline$Phenotype %in% stress_levels],
                  univariate_baseline$Phenotype[univariate_baseline$Phenotype %in% stress_levels])
e2_s <- setNames(univariate_baseline$e2[univariate_baseline$Phenotype %in% stress_levels],
                  univariate_baseline$Phenotype[univariate_baseline$Phenotype %in% stress_levels])

ace_upper <- ss[as.integer(ss$Row) < as.integer(ss$Col), ]
ace_upper$h2x <- h2_s[as.character(ace_upper$Row)]
ace_upper$h2y <- h2_s[as.character(ace_upper$Col)]
ace_upper$c2x <- c2_s[as.character(ace_upper$Row)]
ace_upper$c2y <- c2_s[as.character(ace_upper$Col)]
ace_upper$e2x <- e2_s[as.character(ace_upper$Row)]
ace_upper$e2y <- e2_s[as.character(ace_upper$Col)]

ace_upper$A_part <- ace_upper$rg * sqrt(ace_upper$h2x * ace_upper$h2y)
ace_upper$C_part <- ace_upper$rc * sqrt(ace_upper$c2x * ace_upper$c2y)
ace_upper$E_part <- ace_upper$re * sqrt(ace_upper$e2x * ace_upper$e2y)
ace_upper$sum_abs <- abs(ace_upper$A_part) + abs(ace_upper$C_part) +
                      abs(ace_upper$E_part)
ace_upper$A_pct <- 100 * abs(ace_upper$A_part) / ace_upper$sum_abs
ace_upper$C_pct <- 100 * abs(ace_upper$C_part) / ace_upper$sum_abs
ace_upper$E_pct <- 100 - ace_upper$A_pct - ace_upper$C_pct
ace_upper$rp_lab <- sprintf("rp = %.2f", ace_upper$rp)

ace_upper_long <- pivot_longer(
  ace_upper[, c("Row", "Col", "A_pct", "C_pct", "E_pct", "rp_lab")],
  c("A_pct", "C_pct", "E_pct"),
  names_to = "Estimate_Type", values_to = "Value_pct")

ace_upper_long$Estimate_Type <- recode(ace_upper_long$Estimate_Type,
                                        "A_pct" = "rg",
                                        "C_pct" = "rc",
                                        "E_pct" = "re")
ace_upper_long$Estimate_Type <- factor(ace_upper_long$Estimate_Type,
                                        levels = c("rp", "rg", "rc", "re"))
ace_upper_long$Value_plot <- ace_upper_long$Value_pct / 100

## Add a dummy rp tick so rp shows on the x-axis of every upper panel.
ace_upper_full <- rbind(
  data.frame(unique(ace_upper_long[, c("Row", "Col", "rp_lab")]),
             Estimate_Type = factor("rp", levels = c("rp", "rg", "rc", "re")),
             Value_pct = NA_real_, Value_plot = NA_real_),
  ace_upper_long)

lower_df <- ss[as.integer(ss$Row) > as.integer(ss$Col), ]
lower_df <- pivot_longer(
  lower_df[, c("Row", "Col", "rp", "rp_se", "rg", "rg_se",
                "rc", "rc_se", "re", "re_se")],
  c("rp", "rg", "rc", "re"),
  names_to  = "Estimate_Type",
  values_to = "Value_signed")

lower_df$SE <- with(lower_df,
  ifelse(Estimate_Type == "rp", rp_se,
   ifelse(Estimate_Type == "rg", rg_se,
    ifelse(Estimate_Type == "rc", rc_se, re_se))))
lower_df$Estimate_Type <- factor(lower_df$Estimate_Type,
                                  levels = c("rp", "rg", "rc", "re"))
lower_df$upper_se <- pmax(lower_df$Value_signed + lower_df$SE,
                           lower_df$Value_signed - lower_df$SE)
lower_df$Label_text <- sprintf("%.2f", lower_df$Value_signed)
lower_df$Label_y    <- pmin(1, lower_df$upper_se + 0.10)

grid_all <- expand.grid(
  Row = factor(stress_levels, levels = stress_levels),
  Col = factor(stress_levels, levels = stress_levels))

fig2 <- ggplot() +
  geom_blank(data = grid_all,
              aes(x = factor("rp", levels = c("rp", "rg", "rc", "re")),
                  y = 0)) +
  geom_col(data = ace_upper_full,
           aes(x = Estimate_Type, y = Value_plot, fill = Estimate_Type),
           width = 0.85, position = position_dodge(width = 0.85),
           na.rm = TRUE) +
  geom_text(data = ace_upper_long,
            aes(x = Estimate_Type, y = Value_plot,
                label = paste0(round(Value_pct, 1), "%")),
            vjust = -0.3, size = 2.6, fontface = "bold") +
  geom_text(data = ace_upper[, c("Row", "Col", "rp_lab")],
            aes(x = 2.5, y = 0.90, label = rp_lab),
            inherit.aes = FALSE, size = 3.2, fontface = "bold") +
  geom_col(data = lower_df,
           aes(x = Estimate_Type, y = Value_signed, fill = Estimate_Type),
           width = 0.85, position = position_dodge(width = 0.85)) +
  geom_errorbar(data = lower_df,
                aes(x = Estimate_Type,
                    ymin = pmax(-0.2, Value_signed - SE),
                    ymax = pmin(1.0,  Value_signed + SE)),
                width = 0.2, linewidth = 0.4) +
  geom_hline(data = unique(lower_df[, c("Row", "Col")]),
              aes(yintercept = 0),
              linetype = "dotted", colour = "red", linewidth = 0.4) +
  geom_text(data = lower_df,
            aes(x = Estimate_Type, y = Label_y, label = Label_text),
            size = 2.6, fontface = "bold") +
  facet_grid(Row ~ Col, scales = "fixed", switch = "y", drop = FALSE) +
  scale_x_discrete(limits = c("rp", "rg", "rc", "re"), drop = FALSE) +
  scale_y_continuous(
    limits = c(-0.1, 1.0),
    breaks = c(0.0, 0.5, 1.0),
    labels = number_format(accuracy = 0.1),
    name   = "Correlation (lower triangle)",
    sec.axis = sec_axis(
      ~ . * 100,
      breaks = c(-20, 0, 50, 100),
      labels = number_format(accuracy = 1),
      name   = "ACE magnitude (%) (upper triangle)")) +
  scale_fill_manual(values = estimate_colors,
                     breaks = c("rp", "rg", "rc", "re"),
                     labels = c("Phenotypic (rp)", "Genetic (rg)",
                                "Shared env. (rc)", "Unique env. (re)")) +
  labs(title = "Stress x Stress -- correlations (lower) and ACE decomposition (upper)",
       x = NULL, y = NULL, fill = "Estimate type") +
  theme_minimal() +
  theme(legend.position  = "top",
        legend.title     = element_text(size = 9, face = "bold"),
        legend.text      = element_text(size = 8),
        axis.text.x      = element_text(size = 7, angle = 45, vjust = 1, hjust = 1),
        axis.text.y      = element_text(size = 7),
        panel.spacing    = unit(0.5, "lines"),
        panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        strip.text.x     = element_text(size = 8, face = "bold",
                                         margin = margin(t = 6, b = 6)),
        strip.text.y     = element_text(size = 8, face = "bold",
                                         margin = margin(l = 6, r = 6)),
        strip.placement  = "outside",
        strip.clip       = "off",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey85", linewidth = 0.3),
        plot.title       = element_text(size = 11, face = "bold", hjust = 0.5),
        plot.margin      = margin(t = 14, r = 18, b = 14, l = 18))

save_figure(fig2, "figure_2_stress_stress", width = 11, height = 10)


## -----------------------------------------------------------------------
## 3. Figure 3 -- Personality x Stress (paired Correlation + Magnitude rows)
## -----------------------------------------------------------------------

trait_levels  <- c("Conscientiousness", "Extraversion", "Neuroticism")

ps <- personality_stress
ps$Stress <- factor(ps$Stress, levels = stress_levels)
ps$Trait  <- factor(ps$Trait,  levels = trait_levels)

## Univariate ACE for the personality traits (Table 2 in the paper).
h2_p <- setNames(univariate_baseline$h2[univariate_baseline$Phenotype %in% trait_levels],
                  univariate_baseline$Phenotype[univariate_baseline$Phenotype %in% trait_levels])
c2_p <- setNames(univariate_baseline$c2[univariate_baseline$Phenotype %in% trait_levels],
                  univariate_baseline$Phenotype[univariate_baseline$Phenotype %in% trait_levels])
e2_p <- setNames(univariate_baseline$e2[univariate_baseline$Phenotype %in% trait_levels],
                  univariate_baseline$Phenotype[univariate_baseline$Phenotype %in% trait_levels])

ps$h2x <- h2_p[as.character(ps$Trait)]
ps$h2y <- h2_s[as.character(ps$Stress)]
ps$c2x <- c2_p[as.character(ps$Trait)]
ps$c2y <- c2_s[as.character(ps$Stress)]
ps$e2x <- e2_p[as.character(ps$Trait)]
ps$e2y <- e2_s[as.character(ps$Stress)]

ps$A_part   <- ps$rg * sqrt(ps$h2x * ps$h2y)
ps$C_part   <- ps$rc * sqrt(ps$c2x * ps$c2y)
ps$E_part   <- ps$re * sqrt(ps$e2x * ps$e2y)
ps$sum_abs  <- abs(ps$A_part) + abs(ps$C_part) + abs(ps$E_part)
ps$A_pct    <- 100 * abs(ps$A_part) / ps$sum_abs
ps$C_pct    <- 100 * abs(ps$C_part) / ps$sum_abs
ps$E_pct    <- 100 - ps$A_pct - ps$C_pct
ps$rp_lab   <- sprintf("rp = %.2f", ps$rp)

ace_df <- pivot_longer(
  ps[, c("Trait", "Stress", "A_pct", "C_pct", "E_pct", "rp_lab")],
  c("A_pct", "C_pct", "E_pct"),
  names_to = "Estimate_Type", values_to = "Value_pct")

ace_df$Estimate_Type <- recode(ace_df$Estimate_Type,
                                "A_pct" = "rg", "C_pct" = "rc", "E_pct" = "re")
ace_df$Estimate_Type <- factor(ace_df$Estimate_Type,
                                levels = c("rp", "rg", "rc", "re"))
ace_df$RowType <- factor("Magnitude",
                          levels = c("Correlation", "Magnitude"))
ace_df$Y <- ace_df$Value_pct

ace_dummy <- unique(ace_df[, c("Trait", "Stress", "RowType", "rp_lab")])
ace_dummy$Estimate_Type <- factor("rp", levels = c("rp", "rg", "rc", "re"))
ace_dummy$Value_pct <- NA_real_
ace_dummy$Y <- NA_real_
ace_plot <- rbind(ace_df, ace_dummy[, names(ace_df)])

corr_df <- pivot_longer(
  ps[, c("Trait", "Stress",
          "rp", "rp_se", "rg", "rg_se", "rc", "rc_se", "re", "re_se")],
  c("rp", "rg", "rc", "re"),
    names_to = "Estimate_Type", values_to = "Value_signed")
corr_df$SE <- with(corr_df,
  ifelse(Estimate_Type == "rp", rp_se,
   ifelse(Estimate_Type == "rg", rg_se,
    ifelse(Estimate_Type == "rc", rc_se, re_se))))
corr_df$Estimate_Type <- factor(corr_df$Estimate_Type,
                                 levels = c("rp", "rg", "rc", "re"))
corr_df$RowType <- factor("Correlation",
                           levels = c("Correlation", "Magnitude"))
corr_df$Y <- corr_df$Value_signed
corr_df$upper_se <- pmax(corr_df$Value_signed + corr_df$SE,
                          corr_df$Value_signed - corr_df$SE)
corr_df$lower_se <- pmin(corr_df$Value_signed + corr_df$SE,
                          corr_df$Value_signed - corr_df$SE)
corr_df$Label_text <- sprintf("%.2f", corr_df$Value_signed)
corr_df$Label_y    <- ifelse(corr_df$Value_signed >= 0,
                              pmin( 1.0, corr_df$upper_se + 0.20),
                              pmax(-1.0, corr_df$lower_se - 0.20))

row_levels <- c("Conscientiousness - Correlation",
                 "Conscientiousness - Magnitude",
                 "Extraversion - Correlation",
                 "Extraversion - Magnitude",
                 "Neuroticism - Correlation",
                 "Neuroticism - Magnitude")
corr_df$RowFacet <- factor(paste(corr_df$Trait, "-", corr_df$RowType),
                            levels = row_levels)
ace_plot$RowFacet <- factor(paste(ace_plot$Trait, "-", ace_plot$RowType),
                             levels = row_levels)

mag_rows  <- row_levels[grepl("Magnitude$",  row_levels)]
corr_rows <- row_levels[grepl("Correlation$", row_levels)]
bounds <- rbind(
  expand.grid(RowFacet = mag_rows,  Stress = stress_levels,
              Estimate_Type = c("rp", "rg", "rc", "re"),
              Y = c(0, 100), stringsAsFactors = FALSE),
  expand.grid(RowFacet = corr_rows, Stress = stress_levels,
              Estimate_Type = c("rp", "rg", "rc", "re"),
              Y = c(-1, 1),  stringsAsFactors = FALSE))

bounds$RowFacet      <- factor(bounds$RowFacet,      levels = row_levels)
bounds$Stress        <- factor(bounds$Stress,        levels = stress_levels)
bounds$Estimate_Type <- factor(bounds$Estimate_Type,
                                levels = c("rp", "rg", "rc", "re"))

nf <- function(z) {
  z <- ifelse(abs(z) < 1e-10, 0, z)
  number_format(accuracy = 1, trim = TRUE)(z)
}

p3 <- ggplot() +
  geom_point(data = bounds,
              aes(x = Estimate_Type, y = Y),
              alpha = 0, size = 0.001) +
  geom_col(data = ace_plot,
           aes(x = Estimate_Type, y = Y, fill = Estimate_Type),
           width = 0.85,
           position = position_dodge(width = 0.85), na.rm = TRUE) +
  geom_text(data = ace_plot[!is.na(ace_plot$Y), ],
            aes(x = Estimate_Type, y = Y,
                label = paste0(round(Value_pct, 1), "%")),
            vjust = -0.3, size = 2.3, fontface = "bold") +
  geom_text(data = ace_plot[ace_plot$Estimate_Type == "re", ],
            aes(x = 2.5, y = 90, label = rp_lab),
            inherit.aes = FALSE, size = 3.0, fontface = "bold") +
  geom_col(data = corr_df,
           aes(x = Estimate_Type, y = Y, fill = Estimate_Type),
           width = 0.85,
           position = position_dodge(width = 0.85)) +
  geom_errorbar(data = corr_df,
                aes(x = Estimate_Type,
                    ymin = pmax(-1.0, Y - SE),
                    ymax = pmin( 1.0, Y + SE)),
                width = 0.2, linewidth = 0.4) +
  geom_hline(data = unique(corr_df[, c("RowFacet", "Stress")]),
              aes(yintercept = 0),
              linetype = "dotted", colour = "red", linewidth = 0.4) +
  geom_text(data = corr_df,
            aes(x = Estimate_Type, y = Label_y, label = Label_text),
            size = 2.3, fontface = "bold") +
  facet_grid(
    RowFacet ~ Stress,
    scales = "free_y", switch = "y", drop = FALSE,
    labeller = labeller(
      RowFacet = function(x) sub(" - (Correlation|Magnitude)$", "", x),
      Stress   = label_value)) +
  scale_x_discrete(limits = c("rp", "rg", "rc", "re"), drop = FALSE,
                    expand = expansion(mult = c(0.02, 0.18))) +
  scale_y_continuous(
    breaks = function(lims) if (lims[1] < 0) c(-1, 0, 1) else c(0, 50, 100),
    labels = nf,
    expand = expansion(mult = c(0, 0.04)),
    name = NULL) +
  scale_fill_manual(values = estimate_colors,
                     breaks = c("rp", "rg", "rc", "re"),
                     labels = c("Phenotypic (rp)", "Genetic (rg)",
                                "Shared env. (rc)", "Unique env. (re)")) +
  labs(title = "Personality x Stress -- Correlation row + Magnitude row per trait",
       x = NULL, y = NULL, fill = "Estimate type") +
  theme_minimal() +
  theme(legend.position    = "top",
        legend.title       = element_text(size = 9, face = "bold"),
        legend.text        = element_text(size = 8),
        axis.text.x        = element_text(size = 7, angle = 45, vjust = 1, hjust = 1),
        axis.text.y        = element_text(size = 7),
        panel.spacing      = unit(0.5, "lines"),
        panel.border       = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        strip.text.x       = element_text(size = 8, face = "bold",
                                            margin = margin(t = 6, b = 6)),
        strip.text.y.left  = element_text(size = 7, face = "bold",
                                            margin = margin(l = 4, r = 4)),
        strip.placement    = "outside",
        strip.clip         = "off",
        panel.grid.minor   = element_blank(),
        panel.grid.major   = element_line(colour = "grey85", linewidth = 0.3),
        plot.title         = element_text(size = 11, face = "bold", hjust = 0.5),
        plot.margin        = margin(t = 14, r = 10, b = 14, l = 18)) +
  coord_cartesian(clip = "off")

## Add a right-hand strip column labelling "Correlation" / "Magnitude (%)"
gt <- ggplotGrob(p3)
right_col_width <- unit(9, "mm")
gt <- gtable_add_cols(gt, widths = right_col_width, pos = ncol(gt))

panel_layout <- gt$layout[grepl("^panel", gt$layout$name), ]
row_tops     <- sort(unique(panel_layout$t))
row_bottoms  <- sapply(row_tops, function(ti) max(panel_layout$b[panel_layout$t == ti]))

row_labels <- ifelse(grepl("Correlation$", row_levels), "Correlation", "Magnitude (%)")
for (i in seq_along(row_tops)) {
  gt <- gtable_add_grob(gt,
    grobs = textGrob(row_labels[i],
                       x = unit(0.5, "npc"), y = unit(0.5, "npc"),
                       just = c("centre", "centre"), rot = 270,
                       gp = gpar(fontface = "bold", fontsize = 7)),
    t = row_tops[i], b = row_bottoms[i],
    l = ncol(gt), r = ncol(gt),
    name = paste0("right-strip-", i),
    clip = "off")
}

ggsave(file.path(OUT_FIGURES, "figure_3_personality_stress.png"),
       gt, width = 12, height = 11, dpi = 300)
ggsave(file.path(OUT_FIGURES, "figure_3_personality_stress.pdf"),
       gt, width = 12, height = 11)

message("Step 08 complete. Figures 1-3 saved to ", OUT_FIGURES)
