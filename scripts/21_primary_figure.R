# ============================================
# 21_primary_figure.R  (primary analysis)
# PRIMARY figure: forest plot of sleep_disorder -> CV mortality across the
# adjustment ladder (minimal / confounder-only [PRIMARY] / +BMI / +clinical)
# PLUS the Fine-Gray subdistribution hazard (competing-risks confirmation).
#
# This is the clearer choice over a CIF curve because it shows the central
# robustness message of the reframe in one panel: the sleep->CV association is
# stable and significant across the over-adjustment ladder, and is confirmed by
# the competing-risks (Fine-Gray) estimator.
#
# Input : output/tables/table1_primary_sleep_cvd.csv
# Output: output/figures/fig_primary_sleep_cvd_forest.png  (300 dpi)
# ============================================

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(forcats)
})

cat("=====================================================\n")
cat(" PRIMARY figure: sleep -> CV forest (ladder + Fine-Gray)\n")
cat("=====================================================\n\n")

d <- read.csv("output/tables/table1_primary_sleep_cvd.csv", stringsAsFactors = FALSE)

# Keep only the CV-mortality rows (drop all-cause companion row)
d <- d %>% filter(outcome == "CV mortality")

# Build display labels in a sensible top-to-bottom order
lab_map <- c(
  "minimal"       = "Minimal (age + race)",
  "confounder"    = "Confounder-only (PRIMARY)",
  "plus_bmi"      = "+ BMI",
  "plus_clinical" = "+ BMI + clinical"
)
d$row_label <- NA_character_
cs <- d$estimator == "cause-specific Cox (survey-weighted)"
d$row_label[cs] <- lab_map[d$model[cs]]
fg <- d$estimator == "Fine-Gray subdistribution hazard"
d$row_label[fg] <- "Fine-Gray SHR (competing risks)"

# Order: ladder (cause-specific) then Fine-Gray at the bottom
ord <- c("Minimal (age + race)", "Confounder-only (PRIMARY)",
         "+ BMI", "+ BMI + clinical", "Fine-Gray SHR (competing risks)")
d <- d[match(ord, d$row_label), ]
d <- d[!is.na(d$row_label), ]
d$row_label <- factor(d$row_label, levels = rev(ord))

# Highlight the primary rows
d$is_prim <- d$is_primary == TRUE | d$is_primary == "TRUE"
d$col <- ifelse(d$is_prim, "#C0392B", "#34495E")
d$txt <- sprintf("%.2f (%.2f-%.2f)", d$HR, d$CI_low, d$CI_high)

# keep the HR = 1 reference line visible on the left, and leave room on the
# right for the in-plot HR (95% CI) text labels
xmin <- 0.9
xmax <- max(d$CI_high) * 1.35

p <- ggplot(d, aes(x = HR, y = row_label)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey55") +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high, color = is_prim),
                 height = 0.18, linewidth = 0.7) +
  geom_point(aes(color = is_prim, shape = estimator), size = 3.4) +
  geom_text(aes(x = xmax * 0.985, label = txt), hjust = 1, size = 3.2, color = "grey20") +
  scale_color_manual(values = c(`TRUE` = "#C0392B", `FALSE` = "#34495E"),
                     guide = "none") +
  scale_shape_manual(values = c("cause-specific Cox (survey-weighted)" = 16,
                                "Fine-Gray subdistribution hazard" = 17),
                     name = NULL,
                     labels = c("Cause-specific Cox", "Fine-Gray subdistribution")) +
  scale_x_log10(limits = c(xmin, xmax),
                breaks = c(0.9, 1, 1.5, 2, 2.5, 3),
                labels = function(x) sprintf("%.1f", x)) +
  labs(
    x = "Hazard ratio (95% CI), log scale  (← lower CV risk | higher CV risk →)",
    y = NULL,
    title = "Sleep disturbance and cardiovascular mortality",
    subtitle = "US women with a history of pregnancy, NHANES 2007-2018 (survey-weighted)",
    caption = paste0(
      "Primary model is confounder-only (age, race/ethnicity, education, PIR, smoking, alcohol, parity); ",
      "BMI and\nclinical factors (diabetes, hypertension, lipids) are mediators held out of the total-effect model. ",
      "Fine-Gray\ntreats non-CV death (n=642) as a competing risk. CV deaths n=249; EPV=20.8; E-value=2.70.")
) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank,
    panel.grid.major.y = element_blank,
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10, color = "grey35"),
    plot.caption = element_text(size = 7.5, color = "grey40", hjust = 0),
    legend.position = "top",
    axis.text.y = element_text(size = 10.5),
    plot.margin = margin(t = 12, r = 18, b = 8, l = 8)
)

if (!dir.exists("output/figures")) dir.create("output/figures", recursive = TRUE)
ggsave("output/figures/fig_primary_sleep_cvd_forest.png", plot = p,
       width = 8.8, height = 4.6, dpi = 300)

cat("[OK] output/figures/fig_primary_sleep_cvd_forest.png (300 dpi)\n")
cat(sprintf("Rows plotted: %d (ladder + Fine-Gray)\n", nrow(d)))
print(d[, c("row_label", "estimator", "HR", "CI_low", "CI_high", "P")])