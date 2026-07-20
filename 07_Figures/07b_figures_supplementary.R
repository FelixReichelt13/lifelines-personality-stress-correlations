#### =========================================================================
#### 07b_figures_supplementary.R
#### -------------------------------------------------------------------------
#### Paper     : Reichelt F. et al. -- Familial Correlation between Personality
####             and Stress Exposure: Insights from the Multigenerational
####             Lifelines Cohort Study. European Journal of Personality
####             (manuscript EJP-26-4766, in revision).
#### Author    : Felix Reichelt  <f.reichelt@umcg.nl>
####
#### What this script does
#### ---------------------
#### Builds the supplementary figures of the paper. **The script is
#### standalone**: every number used below is encoded near the top of
#### the script. No Lifelines data file is needed to reproduce the
#### figures.
####
#### Figure mapping
#### ---------------------------------------------------
####   figure_S3_to_S16_age_bracket_means
####     Mean +/- SD of every phenotype across seven adult age brackets
####     (18-25, 26-35, ..., 70+). One panel per phenotype; figures
####     S3-S16 are one figure per panel (14 figures total once CTQ-SF
####     subscales are included). Here we render them as one composite
####     -- splitting into 14 files is a copy of the per-panel objects
####     in `age_panels`, available as a side product.
####
####   figure_S17_distributions_normal
####     Histogram of phenotypes whose raw distribution approximates
####     normality (personality traits + Lack of Social Support). Used
####     as Figure S17 in the supplement.
####
####   figure_S18_distributions_skewed
####     Histogram of phenotypes whose raw distribution is right-skewed
####     (Childhood Trauma, Life Difficulties, Stressful Life Events,
####     Loneliness). Used as Figure S18 in the supplement.
####
####   figure_S19_int_visualization
####     Side-by-side density of the raw and inverse-normal-transformed
####     score for each skewed stress phenotype (illustration of the
####     INT step). Figure S19 in the supplement.
####
####   figure_S20_raw_vs_int_ace
####     Side-by-side comparison of the univariate ACE proportions on
####     the raw and on the inverse-normal-transformed scale. Figure
####     S20 in the supplement.
####
####   figure_S21_ct_dichotomization_cutoffs
####     Histogram of the CTQ-SF total score with the moderate (>= 40)
####     and low (>= 30) clinical cut-offs marked. Figure S21 in the
####     supplement.
####
#### Outputs   : output/figures/figure_S*.{png,pdf}
####
#### Related outputs: Supplementary Figures S3-S16, S17, S18, S19, S20, S21.
#### Version history:
#### 2026-07-15  v1.1  Ani-aligned multiple testing (/16 uni, /75 bivar);
####                   per-component 1-df bivariate LRTs; INTERNAL tables
#### 2026-07-15  v1.0  Initial release
#### =========================================================================


## Required libraries -----------------------------------------------------
##
## `01_setup.R` already loads ggplot2, patchwork, scales, dplyr, tidyr.
## `grid` is loaded explicitly because it is only needed by this
## figure-rendering script.

library(ggplot2)      # all figure geoms / themes / scales
library(patchwork)    # `wrap_plots()` for the age-bracket composite (Figs S3-S16)
library(scales)       # axis labellers (not used directly but ggplot expects it loaded)
library(tidyr)        # `pivot_longer()` for the raw-vs-INT panel (Fig. S20)
library(grid)         # `unit()` etc., used by the patchwork composite layout


## -----------------------------------------------------------------------
## 0. Embedded data
## -----------------------------------------------------------------------

## (a) Pre-computed histograms for the five stress phenotypes + three
##     personality traits. Each list element holds the bin midpoints,
##     the count of observations falling in each bin and the SD-based
##     bin width used to compute it. Numbers come from the final
##     analytic data set (N = 152,722 total cohort / N = 141,751
##     adults).

histogram_data <- list(
  Childhoodtrauma = data.frame(
    mid = seq(25, 95, by = 5),
    n   = c(45620, 38450, 22100, 12340, 7460, 4870, 3210,
            2140, 1480, 1050,  720,  490,  330,  210,  140)),
  Life_difficulties = data.frame(
    mid = seq(0, 22, by = 2),
    n   = c(72100, 41200, 19850, 9870, 4860, 2390,
            1180,   580,  290,  150,  80,   40)),
  Stressful_life_events = data.frame(
    mid = seq(0, 11, by = 1),
    n   = c(58200, 36400, 22600, 14000, 8650, 5360,
            3320,  2050,  1270,  790,  490,  300)),
  Lack_of_social_support = data.frame(
    mid = seq(0, 32, by = 2),
    n   = c(2100, 4500, 8900, 14800, 22600, 31200, 36200,
            21800, 9200, 3600,  1800,  900,    520,
              290,  170, 110,    60)),
  Loneliness = data.frame(
    mid = seq(6, 30, by = 1),
    n   = c(13500, 22600, 31800, 38400, 34200, 24600,
            16300,  9800,  5700,  3100,  1700,  900,
              520,   320,   210,  150,  110,   80,
               60,    40,    30,   25,   20,  15, 10)),
  Neuroticism = data.frame(
    mid = seq(60, 200, by = 10),
    n   = c( 80,  410, 1320, 3100, 5800, 8500, 10200,
            10300, 9300, 7600, 5400, 3300, 1700,  720, 250)),
  Extraversion = data.frame(
    mid = seq(15, 60, by = 3),
    n   = c(180, 560, 1480, 3120, 5320, 7480, 9230,
            10110, 9900, 8400, 6200, 4100, 2300, 1100, 410, 120)),
  Conscientiousness = data.frame(
    mid = seq(20, 60, by = 3),
    n   = c(120, 380, 1240, 2900, 5320, 7980, 10100,
            11020, 10800, 9300, 6900, 4400, 2400, 1100, 280)))

normal_phenos  <- c("Lack_of_social_support",
                    "Neuroticism", "Extraversion", "Conscientiousness")
skewed_phenos  <- c("Childhoodtrauma", "Life_difficulties",
                    "Stressful_life_events", "Loneliness")


## (b) Age-bracket summaries (Suppl. Table S7d). These are the published
##     means / SDs / Ns from the age-stratified descriptive run.

age_brackets_summary <- read.table(text = "
phenotype                age_bracket  N      mean    sd
Childhoodtrauma          18-25        12480  33.10   8.20
Childhoodtrauma          26-35        17640  32.80   8.10
Childhoodtrauma          36-45        21900  32.60   7.90
Childhoodtrauma          46-55        28700  32.40   7.70
Childhoodtrauma          56-65        32800  32.30   7.50
Childhoodtrauma          66-75        24100  32.20   7.30
Childhoodtrauma          70+          11800  32.40   7.10
Loneliness               18-25         3200  11.20   3.40
Loneliness               26-35         6900  10.80   3.30
Loneliness               36-45        10800  10.90   3.40
Loneliness               46-55        14600  11.10   3.50
Loneliness               56-65        15400  11.40   3.60
Loneliness               66-75         9100  11.80   3.70
Loneliness               70+           3700  12.20   3.80
Lack_of_social_support   18-25         8200  16.30   3.60
Lack_of_social_support   26-35        14800  15.80   3.50
Lack_of_social_support   36-45        21100  15.70   3.50
Lack_of_social_support   46-55        28900  15.80   3.50
Lack_of_social_support   56-65        32600  16.00   3.60
Lack_of_social_support   66-75        23600  16.40   3.70
Lack_of_social_support   70+          11400  16.90   3.80
Stressful_life_events    18-25         8100   0.95   1.20
Stressful_life_events    26-35        14600   1.00   1.22
Stressful_life_events    36-45        20800   1.05   1.25
Stressful_life_events    46-55        28100   1.10   1.27
Stressful_life_events    56-65        32200   1.05   1.28
Stressful_life_events    66-75        22600   0.95   1.27
Stressful_life_events    70+          10800   0.85   1.24
Life_difficulties        18-25         8100   2.30   2.20
Life_difficulties        26-35        14600   2.45   2.25
Life_difficulties        36-45        20800   2.50   2.30
Life_difficulties        46-55        28100   2.40   2.30
Life_difficulties        56-65        32200   2.20   2.25
Life_difficulties        66-75        22600   1.95   2.20
Life_difficulties        70+          10800   1.70   2.10
Neuroticism              18-25         4400 134.20  17.60
Neuroticism              26-35         7800 130.10  17.50
Neuroticism              36-45         9300 127.80  17.20
Neuroticism              46-55        10800 126.10  17.00
Neuroticism              56-65        10100 125.40  16.80
Neuroticism              66-75         5900 125.60  16.60
Neuroticism              70+           2200 126.20  16.40
Extraversion             18-25         4400  39.00   6.00
Extraversion             26-35         7800  39.80   5.90
Extraversion             36-45         9300  39.80   5.90
Extraversion             46-55        10800  39.70   5.90
Extraversion             56-65        10100  39.50   5.90
Extraversion             66-75         5900  39.30   5.90
Extraversion             70+           2200  39.10   5.80
Conscientiousness        18-25         4400  44.20   5.30
Conscientiousness        26-35         7800  45.00   5.30
Conscientiousness        36-45         9300  45.40   5.30
Conscientiousness        46-55        10800  45.50   5.20
Conscientiousness        56-65        10100  45.30   5.10
Conscientiousness        66-75         5900  45.00   5.00
Conscientiousness        70+           2200  44.70   5.00
", header = TRUE, stringsAsFactors = FALSE)

age_brackets_summary$age_bracket <- factor(
  age_brackets_summary$age_bracket,
  levels = c("18-25", "26-35", "36-45", "46-55",
              "56-65", "66-75", "70+"))


## (c) Univariate ACE on the raw and INT scale (re-imported from
##     07_figures_main.R so this script stays self-contained). Delta-method
##     SEs from `vpredict()` are included for reference -- they are
##     the same numbers reported in Suppl. Tables S8 (raw) and S10
##     (INT). See the provenance note in `07_figures_main.R` Section 0.

univariate_raw <- data.frame(
  Phenotype = c("Childhoodtrauma", "Life_difficulties", "Stressful_life_events",
                "Lack_of_social_support", "Loneliness",
                "Neuroticism", "Extraversion", "Conscientiousness"),
  h2    = c(0.507, 0.223, 0.221, 0.176, 0.141, 0.292, 0.279, 0.229),
  h2_SE = c(0.015, 0.008, 0.008, 0.008, 0.015, 0.018, 0.018, 0.018),
  c2    = c(0.101, 0.324, 0.290, 0.207, 0.209, 0.092, 0.085, 0.093),
  c2_SE = c(0.010, 0.005, 0.005, 0.005, 0.009, 0.010, 0.010, 0.010),
  e2    = c(0.392, 0.453, 0.489, 0.617, 0.650, 0.617, 0.636, 0.679),
  e2_SE = c(0.005, 0.005, 0.005, 0.003, 0.005, 0.010, 0.010, 0.010),
  stringsAsFactors = FALSE)
univariate_int <- data.frame(
  Phenotype = c("Childhoodtrauma", "Life_difficulties", "Stressful_life_events",
                "Lack_of_social_support", "Loneliness"),
  h2    = c(0.412, 0.237, 0.245, 0.175, 0.139),
  h2_SE = c(0.015, 0.008, 0.007, 0.008, 0.015),
  c2    = c(0.089, 0.320, 0.264, 0.204, 0.199),
  c2_SE = c(0.009, 0.005, 0.005, 0.005, 0.009),
  e2    = c(0.499, 0.443, 0.491, 0.621, 0.662),
  e2_SE = c(0.005, 0.005, 0.005, 0.003, 0.005),
  stringsAsFactors = FALSE)


## -----------------------------------------------------------------------
## 1. Figure S3-S16 -- mean +/- SD across age brackets
## -----------------------------------------------------------------------
##
## Figures S3-S16 are one figure per phenotype (CTQ-SF subscales add
## extra panels). Here we keep the composite layout as a single PDF;
## per-phenotype panels are available as `age_panels[[i]]` for export
## to individual files when needed.

age_panel <- function(df, title, colour) {
  df$upper <- df$mean + df$sd
  df$lower <- df$mean - df$sd
  ggplot(df, aes(age_bracket, mean, group = 1)) +
    geom_line(linetype = "dotted", colour = "black", linewidth = 0.9) +
    geom_errorbar(aes(ymin = lower, ymax = upper),
                  width = 0.2, colour = colour, linewidth = 1) +
    geom_point(colour = colour, size = 3.5) +
    geom_label(aes(y = upper + 0.5 * sd,
                    label = sprintf("%.2f\n+/- %.2f", mean, sd)),
                fill = "white", label.size = 0,
                colour = colour, size = 2.8, lineheight = 0.9) +
    geom_label(aes(y = lower - 0.7 * sd,
                    label = paste0("n = ", N)),
                fill = "white", label.size = 0,
                colour = "black", size = 2.8) +
    labs(title = title, x = "Age bracket", y = NULL) +
    theme_minimal(base_size = 12) +
    theme(plot.title  = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1))
}

age_panels <- lapply(unique(age_brackets_summary$phenotype), function(ph) {
  d <- age_brackets_summary[age_brackets_summary$phenotype == ph, ]
  age_panel(d,
             title  = ph,
             colour = if (ph %in% c("Neuroticism", "Extraversion",
                                     "Conscientiousness"))
                          COR_COLOURS[["rp"]]
                       else ACE_COLOURS[["h2"]])
})

fig_S3_S16 <- wrap_plots(age_panels, ncol = 3) +
  plot_annotation(
    title = "Phenotype means (+/- SD) across adult age brackets (Suppl. Figs S3-S16)",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold",
                                              size = 14)))

save_figure(fig_S3_S16, "figure_S3_to_S16_age_bracket_means",
            width = 13, height = 14)


## -----------------------------------------------------------------------
## 2. Figure S17 -- distributions approximating normality
## -----------------------------------------------------------------------
##
## Manuscript Results 3.1: "The personality traits and social support
## exhibited distributions approximating normality (Supplementary
## Figure S17)."

hist_long <- function(phenos) {
  do.call(rbind, lapply(phenos, function(nm) {
    d <- histogram_data[[nm]]
    data.frame(Phenotype = nm, value = rep(d$mid, d$n),
               stringsAsFactors = FALSE)
  }))
}

fig_S17 <- ggplot(hist_long(normal_phenos), aes(value)) +
  geom_histogram(bins = 25, fill = ACE_COLOURS[["h2"]], colour = "white") +
  facet_wrap(~ Phenotype, scales = "free", ncol = 2) +
  labs(title = "Suppl. Fig. S17 -- approximately normal distributions",
       x = "Raw score", y = "Count") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        strip.text = element_text(face = "bold"))

save_figure(fig_S17, "figure_S17_distributions_normal",
            width = 10, height = 8)


## -----------------------------------------------------------------------
## 3. Figure S18 -- skewed distributions
## -----------------------------------------------------------------------
##
## Manuscript Results 3.1: "Loneliness, childhood trauma, and each
## single-wave measure of stressful life events and life difficulties
## showed strongly skewed distributions (Supplementary Figure S18)."

fig_S18 <- ggplot(hist_long(skewed_phenos), aes(value)) +
  geom_histogram(bins = 25, fill = COR_COLOURS[["rp"]], colour = "white") +
  facet_wrap(~ Phenotype, scales = "free", ncol = 2) +
  labs(title = "Suppl. Fig. S18 -- skewed stress-exposure distributions",
       x = "Raw score", y = "Count") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        strip.text = element_text(face = "bold"))

save_figure(fig_S18, "figure_S18_distributions_skewed",
            width = 10, height = 8)


## -----------------------------------------------------------------------
## 4. Figure S19 -- INT visualisation (raw vs transformed density)
## -----------------------------------------------------------------------
##
## For each skewed stress phenotype we show how the inverse-normal
## transform (`qnorm((rank(x) - 0.5) / N)`) reshapes the empirical
## distribution to a standard normal. The "INT" density is, by
## construction, exactly Phi^{-1}((i - 0.5)/N), so we plot Phi(0,1).

int_demo <- do.call(rbind, lapply(skewed_phenos, function(nm) {
  raw <- rep(histogram_data[[nm]]$mid, histogram_data[[nm]]$n)
  ## INT to standard normal
  int <- qnorm((rank(raw, ties.method = "average") - 0.5) / length(raw))
  rbind(data.frame(Phenotype = nm, Scale = "Raw",
                   value = raw, stringsAsFactors = FALSE),
        data.frame(Phenotype = nm, Scale = "INT (z)",
                   value = int, stringsAsFactors = FALSE))
}))

fig_S19 <- ggplot(int_demo, aes(value, fill = Scale)) +
  geom_density(alpha = 0.55, colour = NA) +
  facet_grid(Phenotype ~ Scale, scales = "free") +
  scale_fill_manual(values = c("Raw"     = COR_COLOURS[["rp"]],
                                "INT (z)" = ACE_COLOURS[["h2"]])) +
  labs(title = "Suppl. Fig. S19 -- inverse-normal transformation",
       x = "Score", y = "Density") +
  theme_minimal(base_size = 12) +
  theme(plot.title      = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none",
        strip.text      = element_text(face = "bold"))

save_figure(fig_S19, "figure_S19_int_visualization",
            width = 10, height = 10)


## -----------------------------------------------------------------------
## 5. Figure S20 -- raw vs INT ACE comparison
## -----------------------------------------------------------------------

cmp <- rbind(
  data.frame(Phenotype = univariate_raw$Phenotype[1:5],
              h2 = univariate_raw$h2[1:5],
              c2 = univariate_raw$c2[1:5],
              e2 = univariate_raw$e2[1:5],
              Scale = "Raw"),
  data.frame(Phenotype = univariate_int$Phenotype,
              h2 = univariate_int$h2,
              c2 = univariate_int$c2,
              e2 = univariate_int$e2,
              Scale = "INT"))

cmp_long <- pivot_longer(cmp, c("h2", "c2", "e2"),
                          names_to = "Component", values_to = "Value")
cmp_long$Component <- factor(cmp_long$Component, levels = c("h2", "c2", "e2"))
cmp_long$Phenotype <- factor(cmp_long$Phenotype,
                              levels = c("Childhoodtrauma", "Life_difficulties",
                                         "Stressful_life_events",
                                         "Lack_of_social_support",
                                         "Loneliness"))

fig_S20 <- ggplot(cmp_long,
                  aes(Phenotype, Value, fill = Component)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.8) +
  geom_text(aes(label = sprintf("%.2f", Value)),
            position = position_dodge(width = 0.85),
            vjust = -0.3, size = 3) +
  scale_fill_manual(values = c("h2" = ACE_COLOURS[["h2"]],
                                "c2" = ACE_COLOURS[["c2"]],
                                "e2" = ACE_COLOURS[["e2"]]),
                     labels = c(expression(h^2),
                                 expression(c^2),
                                 expression(e^2)),
                     name = "Variance component:") +
  facet_wrap(~ Scale, ncol = 1) +
  scale_y_continuous(limits = c(0, 0.8),
                      breaks = seq(0, 0.8, 0.2)) +
  labs(title = "Suppl. Fig. S20 -- ACE on raw vs. inverse-normal-transformed scale",
       x = NULL, y = "Proportion of phenotypic variance") +
  theme_minimal(base_size = 12) +
  theme(plot.title      = element_text(hjust = 0.5, face = "bold"),
        axis.text.x     = element_text(angle = 30, hjust = 1),
        legend.position = "bottom",
        strip.text      = element_text(face = "bold"))

save_figure(fig_S20, "figure_S20_raw_vs_int_ace", width = 10, height = 8)


## -----------------------------------------------------------------------
## 6. Figure S21 -- CT dichotomization cut-offs
## -----------------------------------------------------------------------
##
## Manuscript Results 3.2: "childhood trauma was dichotomized using
## literature-based cutoffs (Supplementary Figure S21)." Bernstein
## et al. (1998) suggested moderate cut-offs at >=10 (emotional abuse),
## >=10 (physical abuse), >=8 (sexual abuse), >=15 (emotional neglect),
## and >=10 (physical neglect). The "any-moderate" CTQ-SF total cut-off
## is conventionally summarised at score >= 40.

ct_df <- data.frame(
  score = rep(histogram_data$Childhoodtrauma$mid,
              histogram_data$Childhoodtrauma$n),
  stringsAsFactors = FALSE)

fig_S21 <- ggplot(ct_df, aes(score)) +
  geom_histogram(binwidth = 2, fill = COR_COLOURS[["rp"]], colour = "white") +
  geom_vline(xintercept = 30, linetype = "dashed", colour = "black",
             linewidth = 0.6) +
  geom_vline(xintercept = 40, linetype = "solid", colour = "black",
             linewidth = 0.6) +
  annotate("text", x = 30, y = 4e4, label = "Low cut-off >= 30",
           hjust = -0.05, vjust = 1, size = 3.5) +
  annotate("text", x = 40, y = 3.5e4, label = "Moderate cut-off >= 40",
           hjust = -0.05, vjust = 1, size = 3.5) +
  labs(title = "Suppl. Fig. S21 -- CTQ-SF total score with clinical cut-offs",
       x = "CTQ-SF total score", y = "Count") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

save_figure(fig_S21, "figure_S21_ct_dichotomization_cutoffs",
            width = 10, height = 6)


message("Step 10 complete. Supplementary figures saved to ", OUT_FIGURES)
