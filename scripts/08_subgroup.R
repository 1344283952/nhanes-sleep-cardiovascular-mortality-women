# ============================================
# 08_subgroup.R
# 亚组 + 交互：DII / sleep_disorder × 全因+CVD 死亡，9 项亚组
#
# 输出:
#   output/tables/table6_subgroup_strat.csv (各亚组 stratified HR)
#   output/tables/table6b_subgroup_pinteraction.csv (交互 p)
#   output/tables/table6_subgroup.xlsx (2 sheets)
# ============================================

suppressPackageStartupMessages({
  library(survey)
  library(survival)
  library(dplyr)
  library(broom)
  library(openxlsx)
})

cat("========================================\n")
cat(" 亚组 + 交互 \n")
cat("========================================\n\n")

load("data/processed/nhanes_design.RData")

# M2 协变量（在交互/分层中去掉 sg_var）
M2_covars <- c("age","race","education","pir","bmi",
               "smoke_status","alcohol_any","parity")

# 9 项亚组（gdm_history 已在  单独做）
subgroups <- list(
  age_group     = c("20-39", "40-59", ">=60"),
  race          = c("Non-Hispanic White", "Non-Hispanic Black",
                    "Mexican American", "Other Hispanic", "Other Race"),
  education     = c("Less than HS", "High school", "College or above"),
  pir_group     = c("<=1.3", "1.3-3.5", ">3.5"),
  bmi_cat       = c("<25", "25-29.9", ">=30"),
  smoke_status  = c("Never", "Ever"),
  alcohol_any   = c("No", "Yes"),
  diabetes      = c(0, 1),
  hypertension  = c(0, 1)
)

exposures <- c("DII", "sleep_disorder")
outcomes  <- c("mort_allcause", "mort_cvd")

# ---- Stratified Cox (M2 去 sg_var) ----
fit_strat <- function(exposure, outcome, sg_var, sg_level) {
  if (is.character(sg_level)) {
    expr_str <- sprintf("%s == '%s'", sg_var, sg_level)
  } else {
    expr_str <- sprintf("%s == %s", sg_var, sg_level)
  }
  sub_dsn <- tryCatch(
    eval(parse(text = sprintf("subset(design, %s)", expr_str))),
    error = function(e) NULL)
  if (is.null(sub_dsn)) return(NULL)

  cov_adj <- M2_covars[M2_covars != sg_var]
  cov_str <- if (length(cov_adj) > 0)
    paste0(" + ", paste(cov_adj, collapse = " + ")) else ""
  fs <- sprintf("Surv(permth, %s) ~ %s%s", outcome, exposure, cov_str)
  fit <- tryCatch(svycoxph(as.formula(fs), design = sub_dsn),
                  error = function(e) NULL,
                  warning = function(w) {
                    tryCatch(svycoxph(as.formula(fs), design = sub_dsn),
                             error = function(e) NULL)
                  })
  if (is.null(fit)) return(NULL)
  td <- tryCatch(
    broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term == exposure),
    error = function(e) NULL)
  if (is.null(td) || nrow(td) == 0) return(NULL)
  td$subgroup_var   <- sg_var
  td$subgroup_level <- as.character(sg_level)
  td$outcome        <- outcome
  td$exposure       <- exposure
  td
}

# ---- 整体交互项 p (全样本 Cox + DII × sg_var 交互) ----
fit_interaction_p <- function(exposure, outcome, sg_var) {
  cov_adj <- M2_covars[M2_covars != sg_var]
  cov_str <- if (length(cov_adj) > 0)
    paste0(" + ", paste(cov_adj, collapse = " + ")) else ""
  fs <- sprintf("Surv(permth, %s) ~ %s * %s%s",
                outcome, exposure, sg_var, cov_str)
  fit <- tryCatch(svycoxph(as.formula(fs), design = design),
                  error = function(e) NULL)
  if (is.null(fit)) return(NA)
  td <- tryCatch(broom::tidy(fit), error = function(e) NULL)
  if (is.null(td)) return(NA)
  inter_rows <- grepl(paste0("^", exposure, ":"), td$term)
  if (sum(inter_rows) == 0) return(NA)
  min(td$p.value[inter_rows], na.rm = TRUE)
}

# ---- 主循环 ----
cat("--- 跑 stratified Cox ---\n")
results <- list
for (exp in exposures) {
  for (out in outcomes) {
    for (sg_var in names(subgroups)) {
      cat(sprintf("  %s × %s × %s ...\n", exp, out, sg_var))
      for (level in subgroups[[sg_var]]) {
        res <- fit_strat(exp, out, sg_var, level)
        if (!is.null(res)) {
          results[[paste(exp, out, sg_var, level, sep="_")]] <- res
        }
      }
    }
  }
}
strat_df <- bind_rows(results)

cat("\n--- 跑全样本交互项 p ---\n")
inter_p <- list
for (exp in exposures) {
  for (out in outcomes) {
    for (sg_var in names(subgroups)) {
      p_int <- fit_interaction_p(exp, out, sg_var)
      inter_p[[paste(exp, out, sg_var, sep="_")]] <- data.frame(
        exposure       = exp,
        outcome        = out,
        subgroup_var   = sg_var,
        p_interaction  = p_int,
        stringsAsFactors = FALSE
)
    }
  }
}
inter_df <- bind_rows(inter_p)

# ---- 打印关键 ----
cat("\n--- Stratified HR (前 30) ---\n")
print(strat_df %>%
        select(exposure, outcome, subgroup_var, subgroup_level,
               estimate, conf.low, conf.high, p.value) %>%
        mutate(across(c(estimate, conf.low, conf.high),
                      ~ sprintf("%.3f", .)),
               p.value = sprintf("%.4f", p.value)) %>%
        head(30))

cat("\n--- 交互项 p (全样本) ---\n")
inter_print <- inter_df
inter_print$p_interaction[is.na(inter_print$p_interaction)] <- -1
tryCatch(
  print(inter_print %>% arrange(p_interaction)),
  error = function(e) cat("打印失败（不阻塞写 csv）:", e$message, "\n")
)

# ---- 输出 ----
if (!dir.exists("output/tables")) dir.create("output/tables", recursive = TRUE)
write.csv(strat_df, "output/tables/table6_subgroup_strat.csv", row.names = FALSE)
write.csv(inter_df, "output/tables/table6b_subgroup_pinteraction.csv",
          row.names = FALSE)

wb <- createWorkbook
addWorksheet(wb, "Stratified_HR");   writeData(wb, "Stratified_HR", strat_df)
addWorksheet(wb, "P_Interaction");   writeData(wb, "P_Interaction", inter_df)
saveWorkbook(wb, "output/tables/table6_subgroup.xlsx", overwrite = TRUE)

cat("\n========================================\n")
cat("已保存:\n")
cat("  output/tables/table6_subgroup_strat.csv\n")
cat("  output/tables/table6b_subgroup_pinteraction.csv\n")
cat("  output/tables/table6_subgroup.xlsx\n")
cat("========================================\n")