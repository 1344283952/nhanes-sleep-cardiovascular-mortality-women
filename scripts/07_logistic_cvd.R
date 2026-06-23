# ============================================
# 07_logistic_cvd.R
# Logistic 回归：DII / Sleep × 5 自报 CVD + 复合 CVD = 6 outcomes
# Crude / M1 / M2 (lifestyle 主) / M3 (full 敏感性)
# 用 svyglm + quasibinomial（survey-weighted）
#
# 输出: output/tables/table3_logistic_cvd.xlsx
# ============================================

suppressPackageStartupMessages({
  library(survey)
  library(dplyr)
  library(broom)
  library(openxlsx)
})

cat("========================================\n")
cat(" Logistic: DII / Sleep × 自报 CVD \n")
cat("========================================\n\n")

load("data/processed/nhanes_design.RData")

design <- update(design,
  DII_Q_numeric = as.numeric(DII_Q))

# ---- 协变量 ----
M1_covars <- c("age", "race")
M2_covars <- c("age", "race", "education", "pir", "bmi",
               "smoke_status", "alcohol_any", "parity")
M3_covars <- c(M2_covars, "diabetes", "hypertension", "ldl", "hdl", "tg")

outcomes <- c("cvd_composite", "cvd_chf", "cvd_chd",
              "cvd_angina", "cvd_mi", "cvd_stroke")
exposures <- list(
  DII_cont       = "DII",
  DII_Q          = "DII_Q",
  sleep_score    = "sleep_score",
  sleep_disorder = "sleep_disorder"
)
cov_str_of <- function(m) {
  switch(m, Crude = "",
            M1    = paste0(" + ", paste(M1_covars, collapse = " + ")),
            M2    = paste0(" + ", paste(M2_covars, collapse = " + ")),
            M3    = paste0(" + ", paste(M3_covars, collapse = " + ")))
}

# ---- Wrapper ----
fit_logit_one <- function(exp_var, out_var, mtype) {
  fs <- sprintf("%s ~ %s%s", out_var, exp_var, cov_str_of(mtype))
  fit <- tryCatch(
    svyglm(as.formula(fs), design = design, family = quasibinomial),
    error = function(e) {
      cat(sprintf("  [失败] %s ~ %s (%s): %s\n",
                  out_var, exp_var, mtype, e$message)); NULL
    }
)
  if (is.null(fit)) return(NULL)
  td <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)
  td$outcome  <- out_var
  td$exposure <- exp_var
  td$model    <- mtype
  td
}

# ---- 主循环：4 exposures × 6 outcomes × 4 models = 96 model fits ----
cat("--- 主效应（4 exposures × 6 outcomes × 4 models）---\n\n")
all_results <- list
for (out_var in outcomes) {
  for (exp_lab in names(exposures)) {
    for (m in c("Crude", "M1", "M2", "M3")) {
      res <- fit_logit_one(exposures[[exp_lab]], out_var, m)
      if (!is.null(res)) {
        res$exposure_label <- exp_lab
        all_results[[paste(out_var, exp_lab, m, sep="_")]] <- res
      }
    }
  }
}
results_df <- bind_rows(all_results)
results_main <- results_df %>%
  filter(grepl("^DII|^sleep_score|^sleep_disorder|^DII_Q", term))

cat("\n--- 主效应 OR (95% CI) ---\n")
print(results_main %>%
        select(outcome, exposure_label, model, term,
               estimate, conf.low, conf.high, p.value) %>%
        mutate(across(c(estimate, conf.low, conf.high),
                      ~ sprintf("%.3f", .)),
               p.value = sprintf("%.4f", p.value)),
      n = 200)

# ---- Benjamini-Hochberg FDR within 6-CVD endpoint family (per exposure × model) ----
cat("\n--- BH FDR (within 6 CVD endpoints per exposure × model) ---\n")
# For each (exposure_label × model) combine across 6 outcomes, take only the EXPOSURE term row
exp_terms <- c(DII = "^DII$", DII_Q = "^DII_QQ4$",
               sleep_score = "^sleep_score$", sleep_disorder = "^sleep_disorder$")

primary_terms <- results_main %>%
  group_by(outcome, exposure_label, model) %>%
  filter(case_when(
    exposure_label == "DII"            ~ term == "DII",
    exposure_label == "DII_Q"          ~ term == "DII_QQ4",
    exposure_label == "sleep_score"    ~ term == "sleep_score",
    exposure_label == "sleep_disorder" ~ term == "sleep_disorder",
    TRUE ~ FALSE)) %>%
  ungroup

# Apply BH within each (exposure × model) family across 6 outcomes
fdr_df <- primary_terms %>%
  group_by(exposure_label, model) %>%
  mutate(q.value = p.adjust(p.value, method = "BH")) %>%
  ungroup %>%
  arrange(exposure_label, model, outcome) %>%
  select(outcome, exposure_label, model, term,
         estimate, conf.low, conf.high, p.value, q.value)

cat("\nBH q-values (per exposure × model, across 6 CVD outcomes):\n")
print(fdr_df %>%
        mutate(across(c(estimate, conf.low, conf.high), ~ sprintf("%.3f", .)),
               across(c(p.value, q.value), ~ sprintf("%.4f", .))),
      n = 100)

# Merge q.value back into results_main (only where it exists)
results_main <- results_main %>%
  left_join(fdr_df %>% select(outcome, exposure_label, model, term, q.value),
            by = c("outcome", "exposure_label", "model", "term"))

write.csv(fdr_df, "output/tables/table3c_logistic_fdr.csv", row.names = FALSE)
cat("\n已保存 output/tables/table3c_logistic_fdr.csv\n\n")

# ---- P-trend for DII_Q ----
cat("\n--- DII P-trend (logistic) ---\n")
ptrend_results <- list
for (out_var in outcomes) {
  for (m in c("Crude", "M1", "M2", "M3")) {
    fs <- sprintf("%s ~ DII_Q_numeric%s", out_var, cov_str_of(m))
    fit <- tryCatch(svyglm(as.formula(fs), design = design,
                           family = quasibinomial),
                    error = function(e) NULL)
    if (!is.null(fit)) {
      td <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
        filter(term == "DII_Q_numeric")
      td$outcome <- out_var; td$model <- m
      ptrend_results[[paste(out_var, m, sep="_")]] <- td
    }
  }
}
ptrend_df <- bind_rows(ptrend_results)
print(ptrend_df %>%
        select(outcome, model, estimate, conf.low, conf.high, p.value))

# ---- 输出 ----
if (!dir.exists("output/tables")) dir.create("output/tables", recursive = TRUE)
write.csv(results_main, "output/tables/table3_logistic_cvd_main.csv", row.names = FALSE)
write.csv(ptrend_df,    "output/tables/table3b_logistic_dii_ptrend.csv", row.names = FALSE)

wb <- createWorkbook
addWorksheet(wb, "Logistic_Main")
writeData(wb, "Logistic_Main", results_main)
addWorksheet(wb, "DII_Ptrend")
writeData(wb, "DII_Ptrend", ptrend_df)
saveWorkbook(wb, "output/tables/table3_logistic_cvd.xlsx", overwrite = TRUE)

cat("\n========================================\n")
cat("已保存:\n")
cat("  output/tables/table3_logistic_cvd_main.csv\n")
cat("  output/tables/table3b_logistic_dii_ptrend.csv\n")
cat("  output/tables/table3_logistic_cvd.xlsx\n")
cat("========================================\n")