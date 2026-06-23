# ============================================
# run_all.R — full reproducible pipeline
# Sleep disturbance → cardiovascular mortality (NHANES 2007–2018)
# 用法: 在项目目录下执行
#   Rscript scripts/run_all.R
# 所有相对路径 (scripts/, data/, output/) 以当前 cwd 为基准
#
# 流程分两段:
#   01–17  数据准备 + 旧探索性分析 (部分已降级为补充材料)
#   20–24  PRIMARY 重构 (sleep → CV 死亡): 生成 manuscript 全部 headline 数字
#          + 出货主表 Table 1–3 + Figure 1
# ============================================

cat("╔══════════════════════════════════════════╗\n")
cat("║  Sleep → CV mortality (NHANES) — 全流程  ║\n")
cat("╚══════════════════════════════════════════╝\n\n")

start_time <- Sys.time

# Comment out 01 if raw data already in data/raw/ to skip re-download (~5-20 min).
# 01–17 准备数据并跑早期探索性分析 (其中部分输出已降级为补充材料, 不再是主结论);
# 20–24 是 sleep → CV 死亡的 PRIMARY 重构, 生成 manuscript 全部 headline 数字 + 出货主表/主图,
# 必须在 01–17 之后跑 (依赖 03/04 产出的 nhanes_final.RData + nhanes_design.RData)。
scripts <- c(
  "scripts/01_download_data.R",     # raw NHANES + mortality download
  "scripts/02_merge_data.R",        # 6 cycles + linked mortality
  "scripts/03_clean_data.R",        # variable construction + analytic sample
  "scripts/04_survey_design.R",     # svydesign with WTDR_6YR / SDMVPSU / SDMVSTRA
  "scripts/05_table1.R",            # (early) baseline characteristics
  "scripts/06_cox_main.R",          # (early) Cox + Schoenfeld PH tests
  "scripts/07_logistic_cvd.R",      # (early) logistic + BH FDR
  "scripts/08_subgroup.R",          # (early) subgroup × interactions
  "scripts/09_rcs.R",               # RCS dose-response (Supplementary Fig)
  "scripts/11_mediation.R",         # (early, exploratory) CMAverse mediation
  "scripts/12_sensitivity.R",       # sensitivity S1-S6 (incl. mice m=20, IPTW)
  "scripts/13_predictive.R",        # C-index + NRI/IDI + E-values
  "scripts/14_consort.R",           # sample-flow figure
  "scripts/15_dag.R",               # DAG figure
  "scripts/16_forest.R",            # subgroup forest plots
  "scripts/17_supp_sensitivity.R",  # Supplementary Fine-Gray + g-formula sensitivity
  # ---- PRIMARY (sleep -> CV mortality) reframe — produces EVERY shipped headline number ----
  "scripts/20_primary.R",       # Fine-Gray + cause-specific Cox primary + DII secondary
  "scripts/21_primary_figure.R",    # Figure 2 (primary forest: ladder + Fine-Gray)
  "scripts/22_baseline_table1.R",   # Table 1 baseline by sleep-disturbance status
  "scripts/23_demote_exploratory.R",# demote abandoned exploratory outputs to supplementary
  "scripts/24_display_tables.R",# shipped main Tables 1-3 + multiplicity display
  "scripts/25_primary_rcs_dag.R"    # Figure 3 (RCS dose-response) + Figure 4 (DAG)
)

n_ok <- 0
n_fail <- 0
fail_list <- character
for (s in scripts) {
  cat(paste0("\n>>> 正在执行: ", s, "\n"))
  cat(paste0(rep("-", 50), collapse = ""), "\n")

  tryCatch({
    source(s, echo = FALSE)
    cat(paste0(">>> ", s, " 执行成功 ✓\n"))
    n_ok <- n_ok + 1
  }, error = function(e) {
    cat(paste0(">>> ", s, " 执行失败 ✗\n"))
    cat(paste0("    错误: ", e$message, "\n"))
    n_fail <<- n_fail + 1
    fail_list <<- c(fail_list, s)
  })
}

end_time <- Sys.time
elapsed <- round(difftime(end_time, start_time, units = "mins"), 1)

cat(paste0("\n╔══════════════════════════════════════════╗\n"))
cat(sprintf("║  全部完成: 成功 %d / 失败 %d / 耗时 %s 分钟\n", n_ok, n_fail, elapsed))
cat(paste0("╚══════════════════════════════════════════╝\n"))

if (length(fail_list) > 0) {
  cat("\n失败脚本:\n")
  for (f in fail_list) cat(paste0("  ", f, "\n"))
}

cat("\n关键输出 — 出货主结论 (sleep -> CV mortality, 由 scripts 20-24 生成):\n")
cat("  output/tables/table1_baseline_by_sleep.xlsx     — Table 1 baseline by sleep-disturbance status\n")
cat("  output/tables/Table2_main_results.xlsx          — Table 2 sleep/DII vs CV/all-cause mortality (Fine-Gray + ladder + DII)\n")
cat("  output/tables/Table3_multiplicity.xlsx          — Table 3 multiple-comparison adjustment (4-test primary family)\n")
cat("  output/tables/table1_primary_sleep_cvd.csv      — primary Fine-Gray + cause-specific HRs (source CSV)\n")
cat("  output/tables/table2_secondary_dii.csv          — secondary DII associations (source CSV)\n")
cat("  output/tables/table3_multiplicity.csv           — family-wise correction (source CSV)\n")
cat("  output/figures/fig1_consort.png                 — Figure 1 sample-selection flow\n")
cat("  output/figures/fig_primary_sleep_cvd_forest.png — primary forest plot\n")
cat("\n  (01-17 的早期/探索性输出已降级为补充材料, 不再是 manuscript 主结论)\n")