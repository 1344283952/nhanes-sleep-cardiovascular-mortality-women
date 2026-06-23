# ============================================
# 09_rcs.R
# RCS 限制性立方样条非线性
# DII / sleep_score → 8 outcomes (2 mortality + 6 CVD)
# 用 rms::cph + 加权（同 a prior NHANES analysis 模板；rms 不支持 PSU/Strata，SE 偏乐观）
#
# 输出:
#   output/tables/rcs_pvalues.csv (p-overall + p-nonlinear)
#   output/figures/rcs_*.png (16 张 + overlay)
# ============================================

suppressPackageStartupMessages({
  library(survival)
  library(rms)
  library(dplyr)
  library(ggplot2)
})

cat("========================================\n")
cat(" RCS 非线性剂量-反应 \n")
cat("========================================\n\n")

load("data/processed/nhanes_design.RData")

if (!dir.exists("output/figures")) dir.create("output/figures", recursive = TRUE)
if (!dir.exists("output/tables"))  dir.create("output/tables",  recursive = TRUE)

# M2 协变量 (Ying et al. 2024 同款 lifestyle-only)
m2_covs <- c("age","race","education","pir","bmi",
             "smoke_status","alcohol_any","parity")

# ---- datadist ----
dd_cols <- unique(c("DII", "sleep_score", "WTDR_6YR",
                    m2_covs,
                    "permth", "mort_allcause", "mort_cvd",
                    "cvd_composite", "cvd_chf", "cvd_chd",
                    "cvd_angina", "cvd_mi", "cvd_stroke"))
nhanes_dd <- nhanes_final[, intersect(dd_cols, names(nhanes_final))]
ddist <- datadist(nhanes_dd)
options(datadist = "ddist")

# ---- outcomes & exposures ----
outcomes <- list(
  list(name="all_mortality", type="cox",   var="mort_allcause"),
  list(name="cvd_mortality", type="cox",   var="mort_cvd"),
  list(name="total_cvd",     type="logit", var="cvd_composite"),
  list(name="chf",           type="logit", var="cvd_chf"),
  list(name="chd",           type="logit", var="cvd_chd"),
  list(name="angina",        type="logit", var="cvd_angina"),
  list(name="mi",            type="logit", var="cvd_mi"),
  list(name="stroke",        type="logit", var="cvd_stroke")
)
exposures <- c("DII", "sleep_score")

# ---- 模型 ----
run_rcs_cox <- function(out_var, exp_name) {
  rhs <- paste(c(sprintf("rcs(%s, 4)", exp_name), m2_covs), collapse = " + ")
  f <- as.formula(sprintf("Surv(permth, %s) ~ %s", out_var, rhs))
  tryCatch(cph(f, data = nhanes_final, x = TRUE, y = TRUE,
               weights = nhanes_final$WTDR_6YR),
           error = function(e) { cat("cph ERR:", e$message, "\n"); NULL })
}
run_rcs_logit <- function(out_var, exp_name) {
  rhs <- paste(c(sprintf("rcs(%s, 4)", exp_name), m2_covs), collapse = " + ")
  f <- as.formula(sprintf("%s ~ %s", out_var, rhs))
  tryCatch(lrm(f, data = nhanes_final, x = TRUE, y = TRUE,
               weights = nhanes_final$WTDR_6YR),
           error = function(e) { cat("lrm ERR:", e$message, "\n"); NULL })
}

# ---- P-overall + P-nonlinear（按 rms::anova 输出格式提取）----
extract_p <- function(fit, exposure) {
  if (is.null(fit)) return(c(p_overall = NA, p_nonlin = NA))
  an <- tryCatch(anova(fit), error = function(e) NULL)
  if (is.null(an)) return(c(p_overall = NA, p_nonlin = NA))
  rn <- rownames(an)
  i_overall <- which(rn == exposure)
  i_nonlin  <- which(grepl("Nonlinear", rn))
  if (length(i_nonlin) > 1 && length(i_overall) > 0) {
    i_nonlin <- i_nonlin[i_nonlin > i_overall][1]
  } else if (length(i_nonlin) >= 1) {
    i_nonlin <- i_nonlin[1]
  } else {
    i_nonlin <- integer(0)
  }
  c(
    p_overall = if (length(i_overall) > 0) an[i_overall, "P"] else NA,
    p_nonlin  = if (length(i_nonlin)  > 0) an[i_nonlin,  "P"] else NA
)
}

# ---- 绘 RCS 图 ----
# fix:
#  (1) ref.zero = TRUE 让 Predict 把曲线锚定到 exposure 中位数 → exp(yhat)
#      产出真正的 HR/OR scale (中位数处 = 1.0, ref line 才有意义)
#  (2) scale_y_log10 让 multiplicative 倍数显示对称
#  (3) y label 根据 outcome_type 显示 "Hazard ratio" or "Odds ratio"
#  (4) annotate P_overall + P_nonlinear (从 ps 注入)
#  (5) ref line 现在真在 y = 1
plot_rcs <- function(fit, exposure, outcome_label, fname,
                     outcome_type = "cox", ps = c(p_overall = NA, p_nonlin = NA)) {
  if (is.null(fit)) return(invisible(NULL))
  pred <- tryCatch(Predict(fit, name = exposure, fun = exp, ref.zero = TRUE),
                   error = function(e) NULL)
  if (is.null(pred)) return(invisible(NULL))
  df <- as.data.frame(pred)
  names(df)[names(df) == exposure] <- "x"

  ylab_text <- if (outcome_type == "cox") "Hazard ratio (95% CI)" else "Odds ratio (95% CI)"

  # P-value 文本: 数值 0 → "< 0.001" (前面已说明 rms 加权 P 数值零)
  fmt_p <- function(p) {
    if (is.na(p)) "NA"
    else if (p < 0.001 || p == 0) "< 0.001"
    else sprintf("%.3f", p)
  }
  p_label <- sprintf("P_overall %s\nP_nonlinear %s",
                     fmt_p(ps["p_overall"]), fmt_p(ps["p_nonlin"]))

  # y range 留头: CI 可能压扁导致 ribbon 看不到; 强制 c(min,max) ± 头部
  y_lo <- min(c(df$lower, 1), na.rm = TRUE)
  y_hi <- max(c(df$upper, 1), na.rm = TRUE)
  y_lo <- max(0.05, y_lo * 0.9)
  y_hi <- y_hi * 1.1

  # annotate 放数据坐标里 (避免 -Inf/+Inf 与 log scale 互动诡异)
  x_lo_data <- min(df$x, na.rm = TRUE)
  y_top <- y_hi * 0.97

  p <- ggplot(df, aes(x = x, y = yhat)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey60") +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#5DADE2", alpha = 0.3) +
    geom_line(color = "#1A5276", linewidth = 1.1) +
    scale_y_log10(limits = c(y_lo, y_hi)) +
    labs(x = exposure, y = ylab_text,
         title = sprintf("RCS: %s -> %s", exposure, outcome_label)) +
    annotate("text", x = x_lo_data, y = y_top, hjust = 0, vjust = 1,
             label = p_label, size = 3.3) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank)
  ggsave(fname, plot = p, width = 5.5, height = 4, dpi = 200)
}

# ---- 主循环 ----
# fix: loop var was named `exp` which shadowed
# the built-in exp function, so Predict(..., fun = exp) silently used the
# loop's character vector instead of the exponential. Renamed to `exp_var`.
results <- list
for (out in outcomes) {
  for (exp_var in exposures) {
    cat(sprintf("--- %s × %s (%s) ---\n", exp_var, out$name, out$type))
    fit <- if (out$type == "cox") run_rcs_cox(out$var, exp_var)
           else                   run_rcs_logit(out$var, exp_var)
    ps  <- extract_p(fit, exp_var)
    results[[paste(exp_var, out$name, sep="_")]] <- data.frame(
      exposure       = exp_var,
      outcome        = out$name,
      outcome_type   = out$type,
      p_overall      = as.numeric(ps["p_overall"]),
      p_nonlinear    = as.numeric(ps["p_nonlin"]),
      stringsAsFactors = FALSE
)
    if (!is.null(fit)) {
      fname <- sprintf("output/figures/rcs_%s_%s.png", exp_var, out$name)
      plot_rcs(fit, exp_var, out$name, fname,
               outcome_type = out$type, ps = ps)
    }
  }
}

res_df <- bind_rows(results)
cat("\n--- RCS p-values ---\n")
print(res_df)

write.csv(res_df, "output/tables/rcs_pvalues.csv", row.names = FALSE)
cat("\n已保存:\n")
cat("  output/tables/rcs_pvalues.csv\n")
cat("  output/figures/rcs_DII_*.png (8 张)\n")
cat("  output/figures/rcs_sleep_score_*.png (8 张)\n")
cat("========================================\n")