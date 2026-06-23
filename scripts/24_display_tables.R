# 24_display_tables.R — 生成投稿用 Table 2 / Table 3 展示表 (.xlsx)
# Table 1 (baseline) 已由 22_baseline_table1.R 生成; 本脚本补 Table 2 (主结果) + Table 3 (多重比较)
library(openxlsx)

fmt <- function(hr, lo, hi, p) {
  if (is.na(hr)) return(c(est="—", p=sprintf("%.3f", p)))
  c(est = sprintf("%.2f (%.2f–%.2f)", hr, lo, hi),
    p   = ifelse(p < 0.001, sprintf("%.1e", p), sprintf("%.3f", p)))
}

t1 <- read.csv("output/tables/table1_primary_sleep_cvd.csv", check.names = FALSE)
t2 <- read.csv("output/tables/table2_secondary_dii.csv", check.names = FALSE)
t3 <- read.csv("output/tables/table3_multiplicity.csv", check.names = FALSE)

# ---- Table 2: 主结果 (睡眠→CV 阶梯 + Fine-Gray + 全因; 次要 DII) ----
rows <- list
add <- function(outcome, exposure, model, estimator, hr, lo, hi, p, extra="") {
  f <- fmt(hr, lo, hi, p)
  rows[[length(rows)+1]] <<- data.frame(
    Outcome=outcome, Exposure=exposure, Model=model, Estimator=estimator,
    `HR/SHR (95% CI)`=f["est"], `P`=f["p"], Note=extra, check.names=FALSE)
}
# 睡眠 → CV: 阶梯 (cause-specific) + Fine-Gray
for (i in seq_len(nrow(t1))) {
  r <- t1[i,]
  if (r$outcome=="CV mortality" && grepl("Cox", r$estimator)) {
    lab <- c(minimal="Minimal (age, race)", confounder="Confounder-only (primary)",
             plus_bmi="+ BMI (mediator; sensitivity)", plus_clinical="+ clinical (mediators; sensitivity)")[r$model]
    add("Cardiovascular", "Sleep disturbance (score ≥2)", lab, "Cause-specific Cox",
        r$HR, r$CI_low, r$CI_high, r$P,
        ifelse(r$is_primary=="TRUE" | r$is_primary==TRUE, sprintf("E-value %.2f (CI %.2f); EPV %.0f; %d CV deaths", r$E_value_point, r$E_value_CIbound, r$EPV, r$N_events_CV), ""))
  }
  if (r$outcome=="CV mortality" && grepl("Fine-Gray", r$estimator))
    add("Cardiovascular", "Sleep disturbance (score ≥2)", "Confounder-only (PRIMARY estimator)", "Fine–Gray subdistribution", r$HR, r$CI_low, r$CI_high, r$P, "Competing-risks (non-CV death)")
  if (r$outcome=="all-cause mortality")
    add("All-cause", "Sleep disturbance (score ≥2)", "Confounder-only", "Cause-specific Cox", r$HR, r$CI_low, r$CI_high, r$P, "Companion")
}
# DII (次要)
for (i in seq_len(nrow(t2))) {
  r <- t2[i,]
  oc <- ifelse(grepl("CV", r$outcome), "Cardiovascular (null)", "All-cause")
  add(oc, paste0("DII (", r$contrast, ")"), "Confounder-only (secondary)", "Cause-specific Cox",
      r$HR, r$CI_low, r$CI_high, r$P, "")
}
T2 <- do.call(rbind, rows)

# ---- Table 3: 多重比较 ----
T3 <- data.frame(
  Test = t3$test,
  `HR (95% CI)` = mapply(function(h,l,u) ifelse(is.na(h),"—",sprintf("%.2f (%.2f–%.2f)",h,l,u)), t3$HR, t3$CI_low, t3$CI_high),
  `P (raw)` = sprintf("%.4f", t3$p_raw),
  `P (Bonferroni)` = sprintf("%.4f", t3$p_bonferroni),
  `P (BH-FDR)` = sprintf("%.4f", t3$p_BH_FDR),
  `Survives 0.05` = ifelse(t3$survives_bonferroni_0.05=="TRUE"|t3$survives_bonferroni_0.05==TRUE,"Yes","No (null)"),
  check.names = FALSE)

wb <- createWorkbook
addWorksheet(wb, "Table2_main_results"); writeData(wb, 1, T2); setColWidths(wb,1,1:ncol(T2),"auto")
saveWorkbook(wb, "output/tables/Table2_main_results.xlsx", overwrite=TRUE)
wb3 <- createWorkbook; addWorksheet(wb3,"Table3_multiplicity"); writeData(wb3,1,T3); setColWidths(wb3,1,1:ncol(T3),"auto")
saveWorkbook(wb3, "output/tables/Table3_multiplicity.xlsx", overwrite=TRUE)

cat("Table2_main_results.xlsx rows:", nrow(T2), "\n")
cat("Table3_multiplicity.xlsx rows:", nrow(T3), "\n")
print(T2[,c("Outcome","Model","HR/SHR (95% CI)","P")], row.names=FALSE)
cat("\n--- Table 3 ---\n"); print(T3, row.names=FALSE)