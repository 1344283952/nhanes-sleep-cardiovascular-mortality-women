# ============================================
# _supp_finegray_gformula.R (— pre-submission sensitivity)
# Two reviewer-anticipated sensitivity analyses appended to Supplementary:
#   1. Fine-Gray subdistribution hazard (competing-risk) for CVD mortality
#      — addresses competing non-CVD deaths (n=643) biasing cause-specific Cox
#   2. CMAverse g-formula mediation (DII → sleep_disorder → all-cause mortality)
#      — addresses exposure-induced mediator-outcome confounding by BMI/diabetes
#
# Output:
#   output/tables/table_supp_S13_finegray.csv
#   output/tables/table_supp_S14_gformula.csv
# ============================================

suppressPackageStartupMessages({
  library(survival)
  library(cmprsk)
  library(CMAverse)
  library(dplyr)
  library(openxlsx)
})

cat("========================================\n")
cat("sensitivity analysis — Fine-Gray + g-formula sensitivity\n")
cat("========================================\n\n")

load("data/processed/nhanes_design.RData")

m2_covs <- c("age", "race", "education", "pir", "bmi",
             "smoke_status", "alcohol_any", "parity")

nhanes_sens <- nhanes_final %>%
  select(DII, sleep_disorder, mort_allcause, mort_cvd, permth,
         all_of(m2_covs)) %>%
  na.omit %>%
  as.data.frame

fac_cols <- sapply(nhanes_sens, is.factor)
nhanes_sens[fac_cols] <- lapply(nhanes_sens[fac_cols], droplevels)

cat(sprintf("Sensitivity-analysis sample N = %d\n", nrow(nhanes_sens)))
cat(sprintf("  all-cause deaths: %d\n", sum(nhanes_sens$mort_allcause)))
cat(sprintf("  CVD deaths:       %d\n", sum(nhanes_sens$mort_cvd)))

# -------------------------------------------
# 1. Fine-Gray competing-risk for CVD mortality
# -------------------------------------------
cat("\n--- 1. Fine-Gray subdistribution (CVD) ---\n")

# event coding: 1 = CVD death, 2 = non-CVD death, 0 = alive/censored
nhanes_sens$status_fg <- with(nhanes_sens,
  ifelse(mort_cvd == 1, 1L,
    ifelse(mort_allcause == 1 & mort_cvd == 0, 2L, 0L)))

cat("status_fg distribution:\n")
print(table(nhanes_sens$status_fg))

cov_form <- ~ DII + sleep_disorder + age + race + education + pir + bmi +
              smoke_status + alcohol_any + parity
X <- model.matrix(cov_form, data = nhanes_sens)[, -1]  # drop intercept

set.seed(20260515)
fg_fit <- crr(ftime   = nhanes_sens$permth,
              fstatus = nhanes_sens$status_fg,
              cov1    = X,
              failcode = 1L,
              cencode  = 0L)

fg_summary <- summary(fg_fit)
# crr summary stores conf.int with HR = exp(coef), and lower/upper 95%
ci_mat <- fg_summary$conf.int
coef_mat <- fg_summary$coef

# Find DII and sleep_disorder rows by row name
get_row <- function(mat, key) {
  matched <- grep(paste0("^", key, "$"), rownames(mat), value = TRUE)
  if (length(matched) == 0) {
    matched <- rownames(mat)[grepl(key, rownames(mat), ignore.case = TRUE)][1]
  }
  matched[1]
}

dii_row   <- get_row(ci_mat, "DII")
sleep_row <- get_row(ci_mat, "sleep_disorder")

fg_out <- data.frame(
  exposure = c("DII per 1-unit", "sleep_disorder (>=2)"),
  SHR      = c(ci_mat[dii_row, "exp(coef)"],     ci_mat[sleep_row, "exp(coef)"]),
  LCL      = c(ci_mat[dii_row, "2.5%"],          ci_mat[sleep_row, "2.5%"]),
  UCL      = c(ci_mat[dii_row, "97.5%"],         ci_mat[sleep_row, "97.5%"]),
  z        = c(coef_mat[dii_row, "z"],           coef_mat[sleep_row, "z"]),
  p_value  = c(coef_mat[dii_row, "p-value"],     coef_mat[sleep_row, "p-value"]),
  comparator_M2_cox = c("1.03 (0.96-1.12) cause-specific",
                        "1.62 (1.13-2.32) cause-specific"),
  method   = "Fine-Gray subdistribution hazard",
  n_total  = nrow(nhanes_sens),
  n_cvd_death = sum(nhanes_sens$status_fg == 1L),
  n_noncvd_death = sum(nhanes_sens$status_fg == 2L),
  stringsAsFactors = FALSE
)
cat("\nFine-Gray result:\n")
print(fg_out)
write.csv(fg_out, "output/tables/table_supp_S13_finegray.csv", row.names = FALSE)
cat("[OK] table_supp_S13_finegray.csv written\n")

# -------------------------------------------
# 2. g-formula mediation: DII -> sleep_disorder -> all-cause
# -------------------------------------------
cat("\n--- 2. g-formula causal mediation (bootstrap B=500) ---\n")

mu    <- mean(nhanes_sens$DII)
sigma <- sd(nhanes_sens$DII)
astar <- mu - 0.5 * sigma
a_val <- mu + 0.5 * sigma

set.seed(20260515)
res_gf <- tryCatch(
  cmest(data       = nhanes_sens,
        model      = "gformula",
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
        nboot      = 500,
        boot.ci.type = "per"),
  error = function(e) {
    cat("[ERR] g-formula failed:", conditionMessage(e), "\n")
    NULL
  }
)

if (!is.null(res_gf)) {
  s <- summary(res_gf)
  est <- s$summarydf

  cat("\nresult dataframe:\n")
  print(est)

  saveRDS(res_gf, "output/tables/supp_gformula_obj.rds")
  txt <- capture.output(print(s))
  writeLines(txt, "output/tables/table_supp_S14_gformula_print.txt")

  # tidy table for journal
  gf_tidy <- data.frame(
    parameter = rownames(est),
    estimate  = est[, 1],
    LCL       = est[, 2],
    UCL       = est[, 3],
    p_value   = est[, ncol(est)],
    method    = "g-formula (CMAverse), B = 500 bootstrap",
    comparator_rb = "regression-based: PM = 0.98% (95% CI -1.5% to 12.6%, P = 0.20)",
    n_total   = nrow(nhanes_sens),
    n_event   = sum(nhanes_sens$mort_allcause),
    contrast  = sprintf("DII a* = %.2f vs a = %.2f (1 SD shift)", astar, a_val),
    stringsAsFactors = FALSE
)
  write.csv(gf_tidy, "output/tables/table_supp_S14_gformula.csv", row.names = FALSE)
  cat("[OK] table_supp_S14_gformula.csv written\n")
}

# -------------------------------------------
# Combine into one XLSX
# -------------------------------------------
wb <- createWorkbook
addWorksheet(wb, "S13_FineGray_CVD")
writeData(wb, "S13_FineGray_CVD", fg_out)
if (exists("gf_tidy")) {
  addWorksheet(wb, "S14_gformula_allcause")
  writeData(wb, "S14_gformula_allcause", gf_tidy)
}
saveWorkbook(wb, "output/tables/table_supp_S13_S14_sensitivity.xlsx", overwrite = TRUE)
cat("\n[ALL DONE] Sensitivity tables saved\n")