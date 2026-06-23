# ============================================
# 06_cox_main.R (-  主分析, M2 调整版)
# Cox 主模型：DII / Sleep × 全因死亡 + CVD 死亡
#   - Crude / M1 (age+race) / M2-lite (主, lifestyle) / M3-full (敏感性, +临床)
#   - 连续 + 四分位 + P-trend
#   - DII × GDM 交互项 (核心创新)
#   - GDM stratified Cox (M2-lite, 避免事件少拟合失败)
#
# M2 定义参考 Ying et al. 2024 Cardiovasc Diabetol 同款：
#   M2 = M1 + lifestyle，不含 hypertension/diabetes/血脂（这些是 DII 下游 mediator）
#   M3 = M2 + 临床调整（作为敏感性，对比 M2 看 dilute 程度）
#
# 输入:  data/processed/nhanes_design.RData
# 输出:  output/tables/table2_cox.xlsx (4 sheets)
# ============================================

suppressPackageStartupMessages({
  library(survey)
  library(survival)
  library(dplyr)
  library(broom)
  library(openxlsx)
})

cat("========================================\n")
cat(" Cox 主模型 (- M2 调整版)\n")
cat("========================================\n\n")

load("data/processed/nhanes_design.RData")

stopifnot(all(c("permth", "mort_allcause", "mort_cvd",
                "DII", "DII_Q", "sleep_score", "sleep_disorder",
                "gdm_history") %in% names(nhanes_final)))

design <- update(design,
  DII_Q_numeric  = as.numeric(DII_Q),
  sleep_Q_numeric = ntile(sleep_score, 4))

# --------------------------------------------------
# 协变量集（Ying et al. 2024 标准）
# --------------------------------------------------
M1_covars <- c("age", "race")
M2_covars <- c("age", "race", "education", "pir", "bmi",
               "smoke_status", "alcohol_any", "parity")
M3_covars <- c(M2_covars, "diabetes", "hypertension", "ldl", "hdl", "tg")

outcomes <- list(
  all_cause = "mort_allcause",
  cvd       = "mort_cvd"
)

exposures <- list(
  DII_cont       = "DII",
  DII_Q          = "DII_Q",
  sleep_score    = "sleep_score",
  sleep_disorder = "sleep_disorder"
)

cov_str_of <- function(mtype) {
  switch(mtype,
         Crude = "",
         M1    = paste0(" + ", paste(M1_covars, collapse = " + ")),
         M2    = paste0(" + ", paste(M2_covars, collapse = " + ")),
         M3    = paste0(" + ", paste(M3_covars, collapse = " + ")))
}

# --------------------------------------------------
# Cox wrapper
# --------------------------------------------------
fit_cox_one <- function(exp_var, out_var, mtype, dsn) {
  fs <- sprintf("Surv(permth, %s) ~ %s%s", out_var, exp_var, cov_str_of(mtype))
  fit <- tryCatch(
    svycoxph(as.formula(fs), design = dsn),
    error = function(e) {
      cat(sprintf("  [失败] %s ~ %s (%s): %s\n",
                  out_var, exp_var, mtype, e$message))
      NULL
    },
    warning = function(w) {
      # 拟合警告不阻塞，但记录
      tryCatch(svycoxph(as.formula(fs), design = dsn),
               error = function(e) NULL)
    }
)
  if (is.null(fit)) return(NULL)
  td <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)
  td$outcome  <- out_var
  td$exposure <- exp_var
  td$model    <- mtype
  td
}

# --------------------------------------------------
# Step 1: 主效应 — 4 暴露 × 2 outcome × 4 model = 32 主行
# --------------------------------------------------
cat("--- Step 1: 主效应 Cox (Crude/M1/M2-lite/M3-full) ---\n\n")

all_results <- list
for (out_lab in names(outcomes)) {
  for (exp_lab in names(exposures)) {
    for (m in c("Crude", "M1", "M2", "M3")) {
      res <- fit_cox_one(exposures[[exp_lab]], outcomes[[out_lab]], m, design)
      if (!is.null(res)) {
        res$exposure_label <- exp_lab
        res$outcome_label  <- out_lab
        all_results[[paste(out_lab, exp_lab, m, sep="_")]] <- res
      }
    }
  }
}
results_df <- bind_rows(all_results)
results_main <- results_df %>%
  filter(grepl("^DII|^sleep_score|^sleep_disorder|^DII_Q", term))

cat("\n--- 暴露主效应 HR (95% CI) ---\n")
print(results_main %>%
        select(outcome_label, exposure_label, model, term,
               estimate, conf.low, conf.high, p.value) %>%
        mutate(across(c(estimate, conf.low, conf.high),
                      ~ sprintf("%.3f", .)),
               p.value = sprintf("%.4f", p.value)),
      n = 200)

# --------------------------------------------------
# Step 1b: Schoenfeld 残差检验 (PH 假设, M2)
# --------------------------------------------------
cat("\n--- Step 1b: Schoenfeld 残差检验 (PH 假设, M2) ---\n")
ph_results <- list
for (out_lab in names(outcomes)) {
  for (exp_lab in names(exposures)) {
    fs <- sprintf("Surv(permth, %s) ~ %s%s",
                  outcomes[[out_lab]], exposures[[exp_lab]], cov_str_of("M2"))
    fit <- tryCatch(svycoxph(as.formula(fs), design = design),
                    error = function(e) NULL,
                    warning = function(w) tryCatch(svycoxph(as.formula(fs), design = design),
                                                   error = function(e) NULL))
    if (is.null(fit)) next
    zph <- tryCatch(cox.zph(fit), error = function(e) NULL)
    if (is.null(zph)) next
    z_tbl <- as.data.frame(zph$table)
    z_tbl$variable <- rownames(z_tbl)
    z_tbl$outcome  <- out_lab
    z_tbl$exposure <- exp_lab
    rownames(z_tbl) <- NULL
    ph_results[[paste(out_lab, exp_lab, sep = "_")]] <- z_tbl
  }
}
ph_df <- bind_rows(ph_results) %>%
  select(outcome, exposure, variable, chisq, df, p)
write.csv(ph_df, "output/tables/table_schoenfeld_ph.csv", row.names = FALSE)
cat(sprintf("Schoenfeld PH test rows: %d (per-variable + GLOBAL per fit)\n", nrow(ph_df)))
cat("PH 违反 (p < 0.05):\n")
print(ph_df %>% filter(p < 0.05) %>% arrange(p), n = 50)
cat(sprintf("\n违反数 / 总数: %d / %d (rate = %.1f%%)\n",
            sum(ph_df$p < 0.05, na.rm = TRUE), sum(!is.na(ph_df$p)),
            100 * mean(ph_df$p < 0.05, na.rm = TRUE)))

# --------------------------------------------------
# Step 2: DII P-trend
# --------------------------------------------------
cat("\n\n--- Step 2: DII P-trend ---\n")
ptrend_results <- list
for (out_lab in names(outcomes)) {
  for (m in c("Crude", "M1", "M2", "M3")) {
    fs <- sprintf("Surv(permth, %s) ~ DII_Q_numeric%s",
                  outcomes[[out_lab]], cov_str_of(m))
    fit <- tryCatch(svycoxph(as.formula(fs), design = design),
                    error = function(e) NULL)
    if (!is.null(fit)) {
      td <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
        filter(term == "DII_Q_numeric")
      td$outcome <- out_lab; td$model <- m
      ptrend_results[[paste(out_lab, m, sep="_")]] <- td
    }
  }
}
ptrend_df <- bind_rows(ptrend_results)
print(ptrend_df %>% select(outcome, model, estimate, conf.low, conf.high, p.value))

# --------------------------------------------------
# Step 3: DII × GDM 交互
# --------------------------------------------------
cat("\n\n--- Step 3: DII × GDM 交互（M2-lite + M3-full）---\n")
inter_results <- list
for (out_lab in names(outcomes)) {
  for (m in c("M2", "M3")) {
    fs <- sprintf("Surv(permth, %s) ~ DII * gdm_history%s",
                  outcomes[[out_lab]], cov_str_of(m))
    fit <- tryCatch(svycoxph(as.formula(fs), design = design),
                    error = function(e) NULL)
    if (!is.null(fit)) {
      td <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
        filter(grepl(":gdm_history|^DII$|^gdm_history", term))
      td$outcome <- out_lab; td$model <- m
      inter_results[[paste(out_lab, m, sep="_")]] <- td
    }
  }
}
inter_df <- bind_rows(inter_results)
print(inter_df %>% select(outcome, model, term, estimate, conf.low, conf.high, p.value))

# --------------------------------------------------
# Step 4: GDM stratified Cox (用 M2-lite 避免拟合失败，CVD 不分层)
# --------------------------------------------------
cat("\n\n--- Step 4: GDM-stratified Cox (M2-lite, 仅全因死亡) ---\n")
cat("（CVD 死亡 GDM=Yes 仅 ~10 events，拟合不可靠，不分层）\n")
strat_results <- list
for (g in c("Yes", "No")) {  # 去掉 Borderline (N=107)
  sub_dsn <- subset(design, gdm_history == g)
  fs <- sprintf("Surv(permth, mort_allcause) ~ DII%s", cov_str_of("M2"))
  fit <- tryCatch(svycoxph(as.formula(fs), design = sub_dsn),
                  error = function(e) NULL)
  if (!is.null(fit)) {
    td <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term == "DII")
    td$gdm <- g
    strat_results[[g]] <- td
  }
}
strat_df <- bind_rows(strat_results)
print(strat_df %>% select(gdm, estimate, conf.low, conf.high, p.value))

# --------------------------------------------------
# 保存
# --------------------------------------------------
if (!dir.exists("output/tables")) dir.create("output/tables", recursive = TRUE)

write.csv(results_main, "output/tables/table2_cox_main.csv", row.names = FALSE)
write.csv(ptrend_df,    "output/tables/table2b_dii_ptrend.csv", row.names = FALSE)
write.csv(inter_df,     "output/tables/table4_dii_gdm_interaction.csv", row.names = FALSE)
write.csv(strat_df,     "output/tables/table4b_dii_stratified_by_gdm.csv", row.names = FALSE)

wb <- createWorkbook
addWorksheet(wb, "Main_Cox");             writeData(wb, "Main_Cox", results_main)
addWorksheet(wb, "DII_Ptrend");           writeData(wb, "DII_Ptrend", ptrend_df)
addWorksheet(wb, "DII_GDM_Interaction");  writeData(wb, "DII_GDM_Interaction", inter_df)
addWorksheet(wb, "GDM_Stratified");       writeData(wb, "GDM_Stratified", strat_df)
saveWorkbook(wb, "output/tables/table2_cox.xlsx", overwrite = TRUE)

cat("\n========================================\n")
cat("已保存:\n")
cat("  output/tables/table2_cox_main.csv\n")
cat("  output/tables/table2b_dii_ptrend.csv\n")
cat("  output/tables/table4_dii_gdm_interaction.csv\n")
cat("  output/tables/table4b_dii_stratified_by_gdm.csv\n")
cat("  output/tables/table2_cox.xlsx (4 sheets)\n")
cat("========================================\n")