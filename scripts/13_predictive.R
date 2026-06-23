# ============================================
# 13_predictive.R
# 预测层：C-index + ΔC + categorical NRI + IDI + E-value
#
# 比较 Base model (M2 协变量) vs Extended (M2 + DII + sleep_disorder)
# 输出 output/tables/table8_predictive.xlsx
# ============================================

suppressPackageStartupMessages({
  library(survival); library(survIDINRI); library(EValue)
  library(dplyr); library(openxlsx); library(broom)
})

cat("========================================\n")
cat(" 预测层 NRI / IDI / E-value \n")
cat("========================================\n\n")

load("data/processed/nhanes_design.RData")

# 未加权样本（survIDINRI 不支持 svy；与 CMAverse 同款 limitation）
nhanes_pred <- nhanes_final %>%
  select(permth, mort_allcause, mort_cvd,
         DII, sleep_disorder, gdm_binary,
         age, race, education, pir, bmi,
         smoke_status, alcohol_any, parity) %>%
  na.omit

cat(sprintf("预测分析样本 N = %d\n", nrow(nhanes_pred)))
cat(sprintf("  全因死亡 events: %d\n", sum(nhanes_pred$mort_allcause)))
cat(sprintf("  CVD 死亡 events:  %d\n\n", sum(nhanes_pred$mort_cvd)))

# ---- Base model: 仅传统协变量 ----
fit_base_all <- coxph(Surv(permth, mort_allcause) ~ age + race + education + pir +
                      bmi + smoke_status + alcohol_any + parity,
                      data = nhanes_pred, x = TRUE, y = TRUE)

# ---- Extended model: + DII + Sleep ----
fit_ext_all <- coxph(Surv(permth, mort_allcause) ~ age + race + education + pir +
                     bmi + smoke_status + alcohol_any + parity +
                     DII + sleep_disorder,
                     data = nhanes_pred, x = TRUE, y = TRUE)

# C-index
c_base_all <- summary(fit_base_all)$concordance[1]
c_ext_all  <- summary(fit_ext_all)$concordance[1]
delta_c_all <- c_ext_all - c_base_all

cat(sprintf("--- 全因死亡 ---\n"))
cat(sprintf("C-index Base:     %.4f\n", c_base_all))
cat(sprintf("C-index Extended: %.4f\n", c_ext_all))
cat(sprintf("ΔC-index:         %.4f\n\n", delta_c_all))

# survIDINRI: 时间依赖 NRI + IDI at t = 60 months (5 yr) 和 t = 120 months (10 yr)
indata_all <- as.matrix(nhanes_pred[, c("permth", "mort_allcause")])
storage.mode(indata_all) <- "numeric"

covs0_all <- model.matrix(~ age + race + education + pir + bmi +
                          smoke_status + alcohol_any + parity, data = nhanes_pred)
covs0_all <- covs0_all[, -1]   # drop intercept
covs1_all <- model.matrix(~ age + race + education + pir + bmi +
                          smoke_status + alcohol_any + parity +
                          DII + sleep_disorder, data = nhanes_pred)
covs1_all <- covs1_all[, -1]

nri_idi_results <- list
for (t0 in c(60, 120)) {
  cat(sprintf("--- 全因死亡 t = %d 个月 ---\n", t0))
  set.seed(20260513)
  res <- tryCatch(
    IDI.INF(indata = indata_all,
            covs0 = covs0_all, covs1 = covs1_all,
            t0 = t0, npert = 100),
    error = function(e) { cat("survIDINRI ERR:", e$message, "\n"); NULL })
  if (!is.null(res)) {
    out <- IDI.INF.OUT(res)
    cat(sprintf("IDI: %.4f (%.4f, %.4f)\n", out[1, 1], out[1, 2], out[1, 3]))
    cat(sprintf("Continuous NRI: %.4f (%.4f, %.4f)\n", out[2, 1], out[2, 2], out[2, 3]))
    cat(sprintf("Median improvement in risk score: %.4f\n\n", out[3, 1]))
    nri_idi_results[[paste0("t", t0)]] <- data.frame(
      outcome = "all_cause", time_months = t0,
      IDI    = out[1, 1], IDI_lower    = out[1, 2], IDI_upper    = out[1, 3],
      NRI    = out[2, 1], NRI_lower    = out[2, 2], NRI_upper    = out[2, 3]
)
  }
}
nri_idi_df <- bind_rows(nri_idi_results)

# ---- E-values: read principal HRs from table2_cox_main.csv (no hardcoding) ----
cat("--- E-values from current table2_cox_main.csv ---\n")
t2 <- read.csv("output/tables/table2_cox_main.csv", stringsAsFactors = FALSE)

get_hr <- function(out, exp, mod, term) {
  r <- t2[t2$outcome_label == out & t2$exposure_label == exp &
          t2$model == mod & t2$term == term, ]
  if (nrow(r) == 0) stop(sprintf("No row for %s %s %s %s", out, exp, mod, term))
  list(est = r$estimate[1], lo = r$conf.low[1], hi = r$conf.high[1])
}

hr_dii    <- get_hr("all_cause", "DII_cont",       "M2", "DII")
hr_diiQ4  <- get_hr("all_cause", "DII_Q",          "M2", "DII_QQ4")
hr_slpAC  <- get_hr("all_cause", "sleep_disorder", "M2", "sleep_disorder")
hr_slpCVD <- get_hr("cvd",       "sleep_disorder", "M2", "sleep_disorder")

ev      <- evalues.HR(est = hr_dii$est,    lo = hr_dii$lo,    hi = hr_dii$hi,    rare = FALSE)
ev_Q4   <- evalues.HR(est = hr_diiQ4$est,  lo = hr_diiQ4$lo,  hi = hr_diiQ4$hi,  rare = FALSE)
ev_slp_ac  <- evalues.HR(est = hr_slpAC$est,  lo = hr_slpAC$lo,  hi = hr_slpAC$hi,  rare = FALSE)
ev_sleep <- evalues.HR(est = hr_slpCVD$est, lo = hr_slpCVD$lo, hi = hr_slpCVD$hi, rare = FALSE)

cat(sprintf("DII per-unit all-cause M2 HR=%.3f: E=%.2f (lo %.2f)\n",
            hr_dii$est, ev[2, 1], ev[2, 2]))
cat(sprintf("DII Q4 vs Q1 all-cause M2 HR=%.3f: E=%.2f (lo %.2f)\n",
            hr_diiQ4$est, ev_Q4[2, 1], ev_Q4[2, 2]))
cat(sprintf("Sleep × all-cause M2 HR=%.3f: E=%.2f (lo %.2f)\n",
            hr_slpAC$est, ev_slp_ac[2, 1], ev_slp_ac[2, 2]))
cat(sprintf("Sleep × CVD M2 HR=%.3f: E=%.2f (lo %.2f)\n",
            hr_slpCVD$est, ev_sleep[2, 1], ev_sleep[2, 2]))

# ---- 汇总 + 保存 ----
pred_summary <- data.frame(
  Model    = c("Base (M2 covars only)", "Extended (+DII +Sleep_disorder)", "ΔC-index"),
  C_index  = c(c_base_all, c_ext_all, delta_c_all),
  outcome  = "all_cause"
)
print(pred_summary)

evalue_df <- data.frame(
  finding    = c(sprintf("DII × all-cause M2 (HR=%.2f)",   hr_dii$est),
                 sprintf("DII Q4 vs Q1 all-cause M2 (HR=%.2f)", hr_diiQ4$est),
                 sprintf("Sleep × all-cause M2 (HR=%.2f)",  hr_slpAC$est),
                 sprintf("Sleep × CVD M2 (HR=%.2f)",        hr_slpCVD$est)),
  HR         = c(hr_dii$est, hr_diiQ4$est, hr_slpAC$est, hr_slpCVD$est),
  HR_lo      = c(hr_dii$lo,  hr_diiQ4$lo,  hr_slpAC$lo,  hr_slpCVD$lo),
  HR_hi      = c(hr_dii$hi,  hr_diiQ4$hi,  hr_slpAC$hi,  hr_slpCVD$hi),
  E_value    = c(ev[2, 1], ev_Q4[2, 1], ev_slp_ac[2, 1], ev_sleep[2, 1]),
  E_value_CI = c(ev[2, 2], ev_Q4[2, 2], ev_slp_ac[2, 2], ev_sleep[2, 2])
)

if (!dir.exists("output/tables")) dir.create("output/tables", recursive = TRUE)

wb <- createWorkbook
addWorksheet(wb, "C_index");      writeData(wb, "C_index", pred_summary)
addWorksheet(wb, "NRI_IDI");      writeData(wb, "NRI_IDI", nri_idi_df)
addWorksheet(wb, "E_values");     writeData(wb, "E_values", evalue_df)
saveWorkbook(wb, "output/tables/table8_predictive.xlsx", overwrite = TRUE)

write.csv(pred_summary, "output/tables/table8a_cindex.csv", row.names = FALSE)
write.csv(nri_idi_df, "output/tables/table8b_nri_idi.csv", row.names = FALSE)
write.csv(evalue_df, "output/tables/table8c_evalue.csv", row.names = FALSE)

cat("\n已保存:\n")
cat("  output/tables/table8_predictive.xlsx (3 sheets)\n")
cat("========================================\n")