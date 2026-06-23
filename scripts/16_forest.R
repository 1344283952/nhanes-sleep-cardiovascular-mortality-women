# ============================================
# 16_forest.R
# Forest plot for subgroup-stratified DII × all-cause mortality
# 包含 p_interaction 标注
# ============================================

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(forcats)
})

cat("========================================\n")
cat(" Forest plot 亚组\n")
cat("========================================\n\n")

strat <- read.csv("output/tables/table6_subgroup_strat.csv",
                  stringsAsFactors = FALSE)
inter <- read.csv("output/tables/table6b_subgroup_pinteraction.csv",
                  stringsAsFactors = FALSE)

if (!dir.exists("output/figures")) dir.create("output/figures", recursive = TRUE)

make_forest <- function(strat_df, inter_df, exposure_filter, outcome_filter,
                        title, fname) {
  d <- strat_df %>%
    filter(exposure == exposure_filter, outcome == outcome_filter) %>%
    mutate(label = paste0(subgroup_var, ": ", subgroup_level))

  if (nrow(d) == 0) {
    cat(sprintf("  (空) %s × %s\n", exposure_filter, outcome_filter))
    return(invisible(NULL))
  }

  p_int <- inter_df %>%
    filter(exposure == exposure_filter, outcome == outcome_filter)

  # 计算每个 subgroup 的 p_int 标注
  pint_labels <- p_int %>%
    mutate(p_label = sprintf("p_int = %.3f", p_interaction))

  d <- d %>%
    left_join(pint_labels %>% select(subgroup_var, p_label),
              by = "subgroup_var")

  # 美化 subgroup_var: 下划线 → 大写空格
  pretty_var <- function(v) {
    v <- gsub("_", " ", v)
    v <- tools::toTitleCase(v)
    v <- sub("Pir", "PIR", v)
    v <- sub("Bmi Cat", "BMI category", v)
    v <- sub("Smoke Status", "Smoking", v)
    v <- sub("Alcohol Any", "Alcohol use", v)
    v <- sub("Age Group", "Age group", v)
    v
  }
  d$label <- sprintf("%s: %s", pretty_var(d$subgroup_var), d$subgroup_level)

  # 排序：按 subgroup_var 然后 subgroup_level
  d$label <- factor(d$label, levels = unique(d$label))

  # Dynamic axis range to prevent CI whisker truncation, with soft caps to keep
  # ultra-sparse subgroups (CI span > 2 log units) from dominating axis scale.
  # fix: rows with CI span > 2 log units (i.e.
  # ratio conf.high / conf.low > 100) are FULL ROW omitted rather than retaining
  # an orphan label with clipped whisker. Full CIs available in Supp Table S9.
  log_span <- log10(d$conf.high / pmax(d$conf.low, 1e-10))
  ultra_sparse <- !is.finite(log_span) | log_span > 2
  d_omitted <- d[ultra_sparse, , drop = FALSE]
  d <- d[!ultra_sparse, , drop = FALSE]

  finite_lo <- d$conf.low[is.finite(d$conf.low) & d$conf.low > 0]
  finite_hi <- d$conf.high[is.finite(d$conf.high)]
  if (length(finite_lo) == 0) finite_lo <- 0.5
  if (length(finite_hi) == 0) finite_hi <- 3
  xmin <- max(0.1, min(0.5, min(finite_lo, na.rm = TRUE) * 0.9))
  xmax <- min(10,  max(3,   max(finite_hi, na.rm = TRUE) * 1.1))
  # After ultra-sparse row removal, remaining whisker clipping should be rare.
  n_clipped <- sum(d$conf.low < xmin | d$conf.high > xmax, na.rm = TRUE)
  n_omitted <- nrow(d_omitted)
  # Sensible log-scale breaks within the computed range
  candidate_breaks <- c(0.1, 0.2, 0.3, 0.5, 0.7, 1, 1.5, 2, 3, 5, 7, 10)
  brks <- candidate_breaks[candidate_breaks >= xmin & candidate_breaks <= xmax]
  if (length(brks) < 3) brks <- c(xmin, 1, xmax)

  p <- ggplot(d, aes(x = estimate, y = fct_rev(label))) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey60") +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                   height = 0.2, color = "#1A5276") +
    geom_point(color = "#E74C3C", size = 2.5) +
    scale_x_log10(limits = c(xmin, xmax), breaks = brks,
                  labels = function(x) sprintf("%.2g", x)) +
    labs(x = sprintf("HR (95%% CI) — %s per 1-unit, log scale (← favors low risk | favors high risk →)",
                     exposure_filter),
         y = NULL, title = title,
         caption = sprintf("Subgroup interaction P-values in Table 4%s",
                           if (n_omitted > 0)
                             sprintf("; %d ultra-sparse subgroup(s) (wide-CI artifact) omitted from plot, full CI in Supplementary Table S9",
                                     n_omitted)
                           else "")) +
    theme_minimal(base_size = 10) +
    theme(panel.grid.minor = element_blank,
          plot.title = element_text(size = 11, face = "bold"),
          plot.caption = element_text(size = 8, color = "grey40"),
          plot.margin = margin(t = 12, r = 16, b = 8, l = 8))

  ggsave(fname, plot = p, width = 8.5, height = max(6, 0.32 * nrow(d) + 1.4),
         dpi = 300, limitsize = FALSE)
  cat(sprintf("  OK %s (%d rows kept, %d ultra-sparse omitted, xrange %.2f-%.2f, %d clipped)\n",
              fname, nrow(d), n_omitted, xmin, xmax, n_clipped))
}

# 4 张 forest
make_forest(strat, inter, "DII", "mort_allcause",
            "DII × All-cause mortality by subgroup (M2 adjusted)",
            "output/figures/fig3a_forest_dii_allcause.png")

make_forest(strat, inter, "DII", "mort_cvd",
            "DII × CVD mortality by subgroup (M2 adjusted)",
            "output/figures/fig3b_forest_dii_cvd.png")

make_forest(strat, inter, "sleep_disorder", "mort_allcause",
            "Sleep disorder × All-cause mortality by subgroup (M2 adjusted)",
            "output/figures/fig3c_forest_sleep_allcause.png")

make_forest(strat, inter, "sleep_disorder", "mort_cvd",
            "Sleep disorder × CVD mortality by subgroup (M2 adjusted)",
            "output/figures/fig3d_forest_sleep_cvd.png")

cat("\n========================================\n")