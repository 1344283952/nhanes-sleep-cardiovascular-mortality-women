# ============================================
# 12_sensitivity.R
# 5 套敏感性分析:
#   S1: 排除前 2 年内死亡（防反向因果）
#   S2: cycle 子集 (2007-2014 vs 2015-2018) 时间稳定性
#   S3: GDM 史子集（Yes / No / Borderline）
#   S5: mice 多重插补 (m=20)
#   S6: IPTW propensity score 加权
#
# (S4 DR1+DR2 双日均值需重 clean，作为 future work)
#
# 输出: output/tables/table7_sensitivity.xlsx
# ============================================

suppressPackageStartupMessages({
  library(survey); library(survival); library(dplyr); library(broom); library(openxlsx)
  library(WeightIt); library(cobalt); library(mice)
})

cat("========================================\n")
cat(" 敏感性分析 \n")
cat("========================================\n\n")

load("data/processed/nhanes_design.RData")

M2_covars <- c("age","race","education","pir","bmi","smoke_status","alcohol_any","parity")
cov_str <- paste0(" + ", paste(M2_covars, collapse = " + "))

fit_sens <- function(dsn, label, formula_str) {
  fit <- tryCatch(svycoxph(as.formula(formula_str), design = dsn),
                  error = function(e) NULL,
                  warning = function(w) tryCatch(svycoxph(as.formula(formula_str), design = dsn), error = function(e) NULL))
  if (is.null(fit)) return(data.frame(scenario = label, term = NA, estimate = NA,
                                       conf.low = NA, conf.high = NA, p.value = NA))
  broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(grepl("^DII$|^sleep_disorder$", term)) %>%
    mutate(scenario = label) %>%
    select(scenario, term, estimate, conf.low, conf.high, p.value)
}

results <- list

# ---- S1: 排除前 2 年内死亡 ----
cat("--- S1: 排除前 2 年内死亡 (permth > 24 m) ---\n")
design_S1 <- subset(design, !(mort_allcause == 1 & permth <= 24))
for (exp in c("DII", "sleep_disorder")) {
  for (out in c("mort_allcause", "mort_cvd")) {
    fs <- sprintf("Surv(permth, %s) ~ %s%s", out, exp, cov_str)
    results[[paste("S1", exp, out, sep="_")]] <- fit_sens(design_S1,
      sprintf("S1_%s_%s", out, exp), fs)
  }
}

# ---- S2a: 2007-2014 ----
cat("--- S2a: cycle 2007-2014 子集 ---\n")
design_S2a <- subset(design, cycle_year %in% c(2007, 2009, 2011, 2013))
for (exp in c("DII", "sleep_disorder")) {
  for (out in c("mort_allcause", "mort_cvd")) {
    fs <- sprintf("Surv(permth, %s) ~ %s%s", out, exp, cov_str)
    results[[paste("S2a", exp, out, sep="_")]] <- fit_sens(design_S2a,
      sprintf("S2a_2007_2014_%s_%s", out, exp), fs)
  }
}

# ---- S2b: 2015-2018 ----
cat("--- S2b: cycle 2015-2018 子集 ---\n")
design_S2b <- subset(design, cycle_year %in% c(2015, 2017))
for (exp in c("DII", "sleep_disorder")) {
  for (out in c("mort_allcause", "mort_cvd")) {
    fs <- sprintf("Surv(permth, %s) ~ %s%s", out, exp, cov_str)
    results[[paste("S2b", exp, out, sep="_")]] <- fit_sens(design_S2b,
      sprintf("S2b_2015_2018_%s_%s", out, exp), fs)
  }
}

# ---- S3: GDM 史子集 ----
for (g in c("Yes", "No", "Borderline")) {
  cat(sprintf("--- S3: GDM = %s ---\n", g))
  expr_str <- sprintf("gdm_history == '%s'", g)
  dsn_S3 <- tryCatch(eval(parse(text = sprintf("subset(design, %s)", expr_str))),
                     error = function(e) NULL)
  if (!is.null(dsn_S3)) {
    for (exp in c("DII", "sleep_disorder")) {
      for (out in c("mort_allcause")) {  # CVD events 太少
        fs <- sprintf("Surv(permth, %s) ~ %s%s", out, exp, cov_str)
        results[[paste("S3", g, exp, out, sep="_")]] <- fit_sens(dsn_S3,
          sprintf("S3_GDM%s_%s_%s", g, out, exp), fs)
      }
    }
  }
}

# ---- S5: mice m=20 多重插补 (M2 covariates + LDL/HDL/TG as auxiliary) ----
cat("\n--- S5: mice multiple imputation (m=20, maxit=10) ---\n")

mi_input <- nhanes_final %>%
  select(SEQN, SDMVPSU, SDMVSTRA, WTDR_6YR, permth, mort_allcause,
         DII, age, race, education, pir, bmi, smoke_status, alcohol_any, parity,
         ldl, hdl, tg)

cat("Missingness per column:\n")
miss_n <- colSums(is.na(mi_input))
print(miss_n)
n_with_miss <- sum(rowSums(is.na(mi_input)) > 0)
cat(sprintf("Total N = %d; rows with any missingness = %d (%.1f%%)\n\n",
            nrow(mi_input), n_with_miss, 100 * n_with_miss / nrow(mi_input)))

# Subset to rows with non-missing outcomes + survey design vars + DII
mi_complete_design <- mi_input %>%
  filter(!is.na(permth), !is.na(mort_allcause),
         !is.na(SDMVPSU), !is.na(SDMVSTRA), !is.na(WTDR_6YR),
         !is.na(DII))

cat(sprintf("After dropping rows missing outcome/PSU/strata/weight/DII: N = %d\n",
            nrow(mi_complete_design)))

# Build predictor matrix and method vector
set.seed(20260514)
init <- mice(mi_complete_design, maxit = 0, printFlag = FALSE)
pred_mat <- init$predictorMatrix
meth     <- init$method

# Don't impute design/outcome/ID
non_impute <- c("SEQN", "SDMVPSU", "SDMVSTRA", "WTDR_6YR",
                "permth", "mort_allcause")
pred_mat[, non_impute] <- 0
pred_mat[non_impute, ] <- 0
meth[non_impute]       <- ""

# Show what mice will do
cat("\nmice method per column:\n")
print(meth[meth != ""])

cat("\nRunning mice (m=20, maxit=10) ...\n")
t0 <- Sys.time
imp <- mice(mi_complete_design,
            m = 20, maxit = 10,
            method = meth, predictorMatrix = pred_mat,
            printFlag = FALSE)
cat(sprintf("mice complete in %.1f sec\n", as.numeric(Sys.time - t0, units = "secs")))

# Per-imputation svycoxph fit + Rubin's rules pooling
mi_per <- list
for (i in 1:imp$m) {
  d <- complete(imp, i)
  dsn_i <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA,
                     weights = ~WTDR_6YR, nest = TRUE, data = d)
  fit_i <- tryCatch(svycoxph(Surv(permth, mort_allcause) ~ DII + age + race +
                             education + pir + bmi + smoke_status +
                             alcohol_any + parity, design = dsn_i),
                    error = function(e) NULL,
                    warning = function(w) tryCatch(svycoxph(Surv(permth, mort_allcause) ~ DII + age + race +
                             education + pir + bmi + smoke_status +
                             alcohol_any + parity, design = dsn_i),
                             error = function(e) NULL))
  if (is.null(fit_i)) next
  tdy <- broom::tidy(fit_i) %>% filter(term == "DII")
  mi_per[[i]] <- data.frame(imp = i, log_hr = tdy$estimate, se = tdy$std.error)
}
mi_df <- bind_rows(mi_per)
cat(sprintf("Successful imputation fits: %d / %d\n", nrow(mi_df), imp$m))

if (nrow(mi_df) >= 2) {
  qbar  <- mean(mi_df$log_hr)
  ubar  <- mean(mi_df$se ^ 2)
  Bvar  <- var(mi_df$log_hr)
  M     <- nrow(mi_df)
  Tvar  <- ubar + (1 + 1 / M) * Bvar
  se_p  <- sqrt(Tvar)
  HR_p  <- exp(qbar)
  HR_lo <- exp(qbar - 1.96 * se_p)
  HR_hi <- exp(qbar + 1.96 * se_p)
  z_p   <- qbar / se_p
  p_pool <- 2 * pnorm(-abs(z_p))

  cat(sprintf("[S5 mice m=%d] Pooled DII × all-cause M2 HR = %.3f (%.3f - %.3f), P = %.4g\n",
              M, HR_p, HR_lo, HR_hi, p_pool))

  results[["S5_mice"]] <- data.frame(
    scenario  = sprintf("S5_mice_m%d_DII_mort_allcause_M2", M),
    term      = "DII",
    estimate  = HR_p,
    conf.low  = HR_lo,
    conf.high = HR_hi,
    p.value   = p_pool
)

  # Save per-imp results for transparency
  write.csv(mi_df, "output/tables/table7b_mice_per_imp.csv", row.names = FALSE)
}

# ---- S6: IPTW (Q4 vs Q1) ----
cat("--- S6: IPTW DII Q4 vs Q1 ---\n")
nhanes_iptw <- nhanes_final %>%
  filter(DII_Q %in% c("Q1", "Q4"),
         if_all(all_of(c("permth", "mort_allcause", M2_covars)), ~ !is.na(.))) %>%
  mutate(DII_high = as.integer(DII_Q == "Q4"))

cat(sprintf("IPTW 样本 N = %d (Q1: %d, Q4: %d)\n",
            nrow(nhanes_iptw),
            sum(nhanes_iptw$DII_high == 0),
            sum(nhanes_iptw$DII_high == 1)))

w_fit <- tryCatch(
  weightit(DII_high ~ age + race + education + pir + bmi +
                      smoke_status + alcohol_any + parity,
           data = nhanes_iptw, method = "ps", estimand = "ATE"),
  error = function(e) NULL)

if (!is.null(w_fit)) {
  bal <- bal.tab(w_fit, m.threshold = 0.1)
  cat("Balance summary:\n"); print(bal$Balance.Across.Pairs %||% bal$Balance %>% head(5))

  fit_S6_all <- coxph(Surv(permth, mort_allcause) ~ DII_high,
                      data = nhanes_iptw, weights = w_fit$weights)
  fit_S6_cvd <- coxph(Surv(permth, mort_cvd) ~ DII_high,
                      data = nhanes_iptw, weights = w_fit$weights)
  results[["S6_allcause"]] <- broom::tidy(fit_S6_all, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term == "DII_high") %>%
    mutate(scenario = "S6_IPTW_Q4vsQ1_allcause") %>%
    select(scenario, term, estimate, conf.low, conf.high, p.value)
  results[["S6_cvd"]] <- broom::tidy(fit_S6_cvd, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term == "DII_high") %>%
    mutate(scenario = "S6_IPTW_Q4vsQ1_cvd") %>%
    select(scenario, term, estimate, conf.low, conf.high, p.value)
}

# ---- 合并 + 输出 ----
all_sens <- bind_rows(results)
cat("\n--- 全部敏感性结果 ---\n")
print(all_sens, n = 50)

if (!dir.exists("output/tables")) dir.create("output/tables", recursive = TRUE)
write.csv(all_sens, "output/tables/table7_sensitivity.csv", row.names = FALSE)

wb <- createWorkbook
addWorksheet(wb, "Sensitivity")
writeData(wb, "Sensitivity", all_sens)
saveWorkbook(wb, "output/tables/table7_sensitivity.xlsx", overwrite = TRUE)

cat("\n已保存 output/tables/table7_sensitivity.xlsx\n")
cat("========================================\n")