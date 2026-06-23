# ============================================
# 22_baseline_table1.R  (primary analysis)
# Survey-weighted baseline characteristics by SLEEP-DISTURBANCE status
# (PRIMARY exposure stratifier for the reframed paper).
#
# Input : data/processed/nhanes_design.RData
# Output: output/tables/table1_baseline_by_sleep.csv
#         output/tables/table1_baseline_by_sleep.xlsx
# ============================================

suppressPackageStartupMessages({
  library(survey)
  library(tableone)
  library(dplyr)
  library(openxlsx)
})

cat("=====================================================\n")
cat(" Baseline characteristics by sleep-disturbance status\n")
cat("=====================================================\n\n")

load("data/processed/nhanes_design.RData")

# Stratifier label
design <- update(design,
  sleep_label = factor(ifelse(sleep_disorder == 1,
                              "Sleep disturbance (score>=2)",
                              "No sleep disturbance (score<2)"),
                       levels = c("No sleep disturbance (score<2)",
                                  "Sleep disturbance (score>=2)")))

vars_cont <- c("age", "DII", "sleep_score", "parity",
               "bmi", "sbp", "dbp", "fbg", "hba1c", "homa_ir",
               "ldl", "hdl", "tc", "tg")
vars_cat  <- c("age_group", "race", "education", "marital",
               "pir_group", "bmi_cat", "smoke_status", "alcohol_any",
               "DII_Q", "gdm_history", "diabetes", "hypertension",
               "cvd_composite", "mort_allcause", "mort_cvd")
all_vars  <- c(vars_cont, vars_cat)
nonnormal <- c("DII", "sleep_score", "parity", "fbg", "hba1c", "homa_ir",
               "ldl", "hdl", "tc", "tg", "tg")

tab1 <- svyCreateTableOne(
  vars = all_vars, strata = "sleep_label", factorVars = vars_cat,
  data = design, test = TRUE)

print(tab1, showAllLevels = TRUE, nonnormal = nonnormal,
      formatOptions = list(big.mark = ","))

mat <- print(tab1, showAllLevels = TRUE, nonnormal = nonnormal,
             printToggle = FALSE, noSpaces = TRUE,
             formatOptions = list(big.mark = ","))

if (!dir.exists("output/tables")) dir.create("output/tables", recursive = TRUE)
write.csv(mat, "output/tables/table1_baseline_by_sleep.csv")

wb <- createWorkbook
addWorksheet(wb, "Baseline_by_sleep")
writeData(wb, "Baseline_by_sleep", mat, rowNames = TRUE)
saveWorkbook(wb, "output/tables/table1_baseline_by_sleep.xlsx", overwrite = TRUE)

cat("\n[OK] output/tables/table1_baseline_by_sleep.csv\n")
cat("[OK] output/tables/table1_baseline_by_sleep.xlsx\n")