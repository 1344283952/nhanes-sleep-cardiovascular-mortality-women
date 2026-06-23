# ============================================
# 11_mediation.R
# CMAverse 因果中介：DII → sleep_disorder → 死亡
#   - 主分析: 全因死亡 (events=894)
#   - 副分析: CVD 死亡 (events=251)
#   - 模型 "rb" (regression-based): mediator logistic + outcome coxph
#   - 未加权（CMAverse 限制；敏感性用 svyglm + difference method 交叉验证）
#
# 输出:
#   output/tables/table5_mediation_main.csv (全因)
#   output/tables/table5b_mediation_cvd.csv (CVD)
#   output/tables/table5_mediation.xlsx
#
# ⚠️ Limitation: CMAverse 不支持 svydesign 权重，主分析未加权。
#    Ying et al. 2024 + Shi B 2021 Epidemiology 同款处理；论文写明。
# ============================================

suppressPackageStartupMessages({
  library(CMAverse)
  library(survival)
  library(dplyr)
  library(openxlsx)
})

cat("========================================\n")
cat(" CMAverse 中介分析 \n")
cat("========================================\n\n")

load("data/processed/nhanes_design.RData")

# --------------------------------------------------
# 整理中介分析样本（完整 case，剔除 mediator/outcome/exposure/covar 缺失）
# --------------------------------------------------
m2_covs <- c("age","race","education","pir","bmi",
             "smoke_status","alcohol_any","parity")

nhanes_med <- nhanes_final %>%
  select(DII, sleep_disorder, mort_allcause, mort_cvd, permth,
         all_of(m2_covs)) %>%
  na.omit %>%
  as.data.frame
# 强制 factor 重新 droplevels（防 CMAverse 把空 level 当 NA）
factor_cols <- sapply(nhanes_med, is.factor)
nhanes_med[factor_cols] <- lapply(nhanes_med[factor_cols], droplevels)

cat(sprintf("中介分析样本 N: %d (从 %d 缩到 %d)\n",
            nrow(nhanes_med), nrow(nhanes_final),
            nrow(nhanes_med)))
cat(sprintf("  全因死亡 events: %d\n",
            sum(nhanes_med$mort_allcause)))
cat(sprintf("  CVD 死亡 events:  %d\n",
            sum(nhanes_med$mort_cvd)))
cat(sprintf("  sleep_disorder=1: %d (%.1f%%)\n",
            sum(nhanes_med$sleep_disorder),
            100 * mean(nhanes_med$sleep_disorder)))

# DII astar / a 设为 mean ± 0.5 SD（约 P25 到 P75 区间）
mu <- mean(nhanes_med$DII)
sigma <- sd(nhanes_med$DII)
astar <- mu - 0.5 * sigma
a_val <- mu + 0.5 * sigma
cat(sprintf("\nDII astar (低暴露): %.2f\n", astar))
cat(sprintf("DII a (高暴露):     %.2f\n", a_val))
cat(sprintf("DII 1 SD = %.2f\n", sigma))

# --------------------------------------------------
# 主分析: DII → sleep_disorder → 全因死亡
# --------------------------------------------------
cat("\n========================================\n")
cat("主分析: DII → sleep_disorder → 全因死亡\n")
cat("（200 bootstrap，估计 5-10 分钟）\n")
cat("========================================\n")

set.seed(20260513)

res_main <- tryCatch(
  cmest(data       = nhanes_med,
        model      = "rb",
        outcome    = "permth",
        event      = "mort_allcause",
        exposure   = "DII",
        mediator   = "sleep_disorder",
        basec      = m2_covs,
        EMint      = FALSE,
        mreg       = list("logistic"),
        yreg       = "coxph",
        astar      = astar,
        a          = a_val,
        mval       = list(0),
        estimation = "imputation",
        inference  = "bootstrap",
        nboot      = 1000,
        boot.ci.type = "per"),
  error = function(e) {
    cat("CMAverse ERR (全因):", e$message, "\n"); NULL
  }
)

if (!is.null(res_main)) {
  cat("\n--- 主分析 summary ---\n")
  print(summary(res_main))

  # 整理输出
  est <- res_main$effect.pe
  cil <- res_main$effect.ci.low
  ciu <- res_main$effect.ci.high
  pvl <- res_main$effect.pval

  out_main <- data.frame(
    effect    = names(est),
    estimate  = as.numeric(est),
    conf.low  = as.numeric(cil),
    conf.high = as.numeric(ciu),
    p.value   = as.numeric(pvl),
    outcome   = "all_cause",
    stringsAsFactors = FALSE
)
  cat("\n--- 主分析效应估计 ---\n")
  print(out_main)

  if (!dir.exists("output/tables")) dir.create("output/tables", recursive = TRUE)
  write.csv(out_main, "output/tables/table5_mediation_main.csv",
            row.names = FALSE)
} else {
  out_main <- data.frame
}

# --------------------------------------------------
# 副分析: DII → sleep_disorder → CVD 死亡
# --------------------------------------------------
cat("\n========================================\n")
cat("副分析: DII → sleep_disorder → CVD 死亡\n")
cat("========================================\n")

set.seed(20260513)

res_cvd <- tryCatch(
  cmest(data       = nhanes_med,
        model      = "rb",
        outcome    = "permth",
        event      = "mort_cvd",
        exposure   = "DII",
        mediator   = "sleep_disorder",
        basec      = m2_covs,
        EMint      = FALSE,
        mreg       = list("logistic"),
        yreg       = "coxph",
        astar      = astar,
        a          = a_val,
        mval       = list(0),
        estimation = "imputation",
        inference  = "bootstrap",
        nboot      = 1000,
        boot.ci.type = "per"),
  error = function(e) {
    cat("CMAverse ERR (CVD):", e$message, "\n"); NULL
  }
)

if (!is.null(res_cvd)) {
  cat("\n--- 副分析 summary ---\n")
  print(summary(res_cvd))

  est <- res_cvd$effect.pe
  cil <- res_cvd$effect.ci.low
  ciu <- res_cvd$effect.ci.high
  pvl <- res_cvd$effect.pval

  out_cvd <- data.frame(
    effect    = names(est),
    estimate  = as.numeric(est),
    conf.low  = as.numeric(cil),
    conf.high = as.numeric(ciu),
    p.value   = as.numeric(pvl),
    outcome   = "cvd",
    stringsAsFactors = FALSE
)
  cat("\n--- 副分析效应估计 ---\n")
  print(out_cvd)

  write.csv(out_cvd, "output/tables/table5b_mediation_cvd.csv",
            row.names = FALSE)
} else {
  out_cvd <- data.frame
}

# --------------------------------------------------
# Excel 合并
# --------------------------------------------------
wb <- createWorkbook
if (nrow(out_main) > 0) {
  addWorksheet(wb, "Mediation_AllCause")
  writeData(wb, "Mediation_AllCause", out_main)
}
if (nrow(out_cvd) > 0) {
  addWorksheet(wb, "Mediation_CVD")
  writeData(wb, "Mediation_CVD", out_cvd)
}
if (nrow(out_main) > 0 || nrow(out_cvd) > 0) {
  saveWorkbook(wb, "output/tables/table5_mediation.xlsx", overwrite = TRUE)
}

cat("\n========================================\n")
cat("已保存:\n")
cat("  output/tables/table5_mediation_main.csv\n")
cat("  output/tables/table5b_mediation_cvd.csv\n")
cat("  output/tables/table5_mediation.xlsx\n")
cat("========================================\n")