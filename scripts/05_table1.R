# ============================================
# 05_table1.R (版)
# Table 1:  按 GDM 史分层加权基线
# Table 1b: 按全因死亡分层加权基线
#
# 输入:  data/processed/nhanes_design.RData
# 输出:  output/tables/table1_by_gdm.csv
#        output/tables/table1b_by_mortality.csv
#        output/tables/table1.xlsx (合并两 sheet)
# ============================================

suppressPackageStartupMessages({
  library(survey)
  library(tableone)
  library(dplyr)
  library(openxlsx)
})

cat("========================================\n")
cat("Table 1: 加权基线（）\n")
cat("========================================\n\n")

load("data/processed/nhanes_design.RData")

# 给 design 加 mort_label
design <- update(design,
  mort_label = factor(ifelse(mort_allcause == 1, "Death", "Alive"),
                      levels = c("Alive", "Death")))

# --------------------------------------------------
# 变量定义
# --------------------------------------------------
vars_cont <- c("age", "DII", "sleep_score",
               "bmi", "sbp", "dbp",
               "fbg", "hba1c", "homa_ir",
               "ldl", "hdl", "tc", "tg",
               "parity")

vars_cat  <- c("age_group", "race", "education",
               "marital", "pir_group", "bmi_cat",
               "smoke_status", "alcohol_any",
               "DII_Q", "sleep_disorder",
               "diabetes", "hypertension",
               "cvd_chf", "cvd_chd", "cvd_angina", "cvd_mi",
               "cvd_stroke", "cvd_composite",
               "mort_allcause", "mort_cvd")

all_vars <- c(vars_cont, vars_cat)
nonnormal_vars <- c("DII", "sleep_score", "fbg", "hba1c", "homa_ir",
                    "ldl", "hdl", "tc", "tg", "parity")

# --------------------------------------------------
# Table 1: 按 GDM 史
# --------------------------------------------------
cat("--- Table 1 (by GDM history) ---\n\n")
tab1_gdm <- svyCreateTableOne(
  vars       = all_vars,
  strata     = "gdm_history",
  factorVars = vars_cat,
  data       = design,
  test       = TRUE
)
print(tab1_gdm, showAllLevels = TRUE, nonnormal = nonnormal_vars,
      formatOptions = list(big.mark = ","))

tab1_mat <- print(tab1_gdm, showAllLevels = TRUE, nonnormal = nonnormal_vars,
                  printToggle = FALSE, noSpaces = TRUE,
                  formatOptions = list(big.mark = ","))

# --------------------------------------------------
# Table 1b: 按全因死亡
# --------------------------------------------------
cat("\n\n--- Table 1b (by all-cause mortality) ---\n\n")

# 死亡分层时不重复 mort_allcause / mort_cvd 列
vars_b      <- setdiff(all_vars, c("mort_allcause", "mort_cvd"))
vars_b_cat  <- setdiff(vars_cat, c("mort_allcause", "mort_cvd"))

tab1b_mort <- svyCreateTableOne(
  vars       = vars_b,
  strata     = "mort_label",
  factorVars = vars_b_cat,
  data       = design,
  test       = TRUE
)
print(tab1b_mort, showAllLevels = TRUE, nonnormal = nonnormal_vars,
      formatOptions = list(big.mark = ","))

tab1b_mat <- print(tab1b_mort, showAllLevels = TRUE, nonnormal = nonnormal_vars,
                   printToggle = FALSE, noSpaces = TRUE,
                   formatOptions = list(big.mark = ","))

# --------------------------------------------------
# 输出
# --------------------------------------------------
if (!dir.exists("output/tables")) dir.create("output/tables", recursive = TRUE)
write.csv(tab1_mat,  "output/tables/table1_by_gdm.csv")
write.csv(tab1b_mat, "output/tables/table1b_by_mortality.csv")

wb <- createWorkbook
addWorksheet(wb, "Table_1_by_GDM")
writeData(wb, "Table_1_by_GDM", tab1_mat, rowNames = TRUE)
addWorksheet(wb, "Table_1b_by_Mortality")
writeData(wb, "Table_1b_by_Mortality", tab1b_mat, rowNames = TRUE)
saveWorkbook(wb, "output/tables/table1.xlsx", overwrite = TRUE)

cat("\n已保存:\n")
cat("  output/tables/table1_by_gdm.csv\n")
cat("  output/tables/table1b_by_mortality.csv\n")
cat("  output/tables/table1.xlsx\n")
cat("========================================\n")