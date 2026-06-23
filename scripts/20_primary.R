# ============================================
# 20_primary.R  (confounder-only primary analysis)
#
# NEW THESIS:
#   PRIMARY   = sleep disturbance (sleep_disorder, score>=2) -> CARDIOVASCULAR
#               mortality in US women with a history of pregnancy.
#   SECONDARY = pro-inflammatory diet (DII) -> ALL-CAUSE mortality.
#
# PRIMARY MODEL = CONFOUNDER-ONLY (total-effect) per the over-adjustment lesson:
#   BMI / diabetes / hypertension / lipids are MEDIATORS of sleep/diet -> CV
#   death and are HELD OUT of the primary model.
#     PRIMARY covars (C0) = age + race + education + pir + smoke_status
#                           + alcohol_any + parity                (NO BMI, NO clinical)
#   Sensitivity ladder (robustness):
#     minimal (M1)   = age + race
#     + BMI          = C0 + bmi                  (== old M2)
#     + clinical     = +BMI + diabetes + hypertension + ldl + hdl + tg   (== old M3)
#
# PRIMARY CV estimator = survey-weighted cause-specific Cox (svycoxph)
#                        + Fine-Gray subdistribution hazard (competing non-CV death).
#
# Inputs : data/processed/nhanes_design.RData
# Outputs:
#   output/tables/table1_primary_sleep_cvd.csv    (sleep->CV: ladder + Fine-Gray + E-value + EPV
#                                                   + sleep->all-cause companion)
#   output/tables/table2_secondary_dii.csv        (DII per-SD + quartiles -> all-cause; -> CV null)
#   output/tables/table3_multiplicity.csv         (Bonferroni + BH-FDR over 4-test PRIMARY family)
#   output/tables/_primary_verification.txt       (printed primary numbers + gates)
# ============================================

suppressPackageStartupMessages({
  library(survey)
  library(survival)
  library(cmprsk)
  library(EValue)
  library(dplyr)
  library(broom)
})

cat("=====================================================\n")
cat("  PRIMARY analysis (primary)\n")
cat(" PRIMARY: sleep_disorder -> CV mortality (confounder-only)\n")
cat(" SECONDARY: DII -> all-cause mortality (confounder-only)\n")
cat("=====================================================\n\n")

load("data/processed/nhanes_design.RData")

stopifnot(all(c("permth", "mort_allcause", "mort_cvd",
                "DII", "DII_Q", "sleep_disorder",
                "age", "race", "education", "pir",
                "smoke_status", "alcohol_any", "parity",
                "bmi") %in% names(nhanes_final)))

# DII per-SD exposure (standardised so HR is "per 1 SD")
dii_sd <- sd(nhanes_final$DII, na.rm = TRUE)
design <- update(design,
                 DII_z = DII / dii_sd,                 # HR per 1 SD of DII
                 DII_Q_numeric = as.numeric(DII_Q))    # for P-trend

# --------------------------------------------------
# Covariate sets  (PRIMARY = confounder-only, NO BMI)
# --------------------------------------------------
C0_covars  <- c("age", "race", "education", "pir",
                "smoke_status", "alcohol_any", "parity")   # PRIMARY confounder-only
M1_covars  <- c("age", "race")                              # minimal
BMI_covars <- c(C0_covars, "bmi")                           # + BMI (== old M2)
CLIN_covars<- c(BMI_covars, "diabetes", "hypertension",
                "ldl", "hdl", "tg")                         # + clinical (== old M3)

cov_str <- function(covs) if (length(covs)) paste0(" + ", paste(covs, collapse = " + ")) else ""

# --------------------------------------------------
# Survey-weighted Cox wrapper (cause-specific)
# --------------------------------------------------
fit_cox <- function(exp_var, out_var, covs, dsn = design, keep = exp_var) {
  fs <- sprintf("Surv(permth, %s) ~ %s%s", out_var, exp_var, cov_str(covs))
  fit <- tryCatch(svycoxph(as.formula(fs), design = dsn),
                  error = function(e) { cat("  [Cox fail]", fs, "::", e$message, "\n"); NULL })
  if (is.null(fit)) return(NULL)
  broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(grepl(keep, term))
}

# ==================================================
# PART A — PRIMARY: sleep_disorder -> CV mortality
# ==================================================
cat("--- PART A: sleep_disorder -> CV mortality ---\n\n")

# A1. Cause-specific Cox sensitivity ladder (CONFOUNDER-ONLY is PRIMARY)
ladder_specs <- list(
  minimal     = M1_covars,
  confounder  = C0_covars,     # <<< PRIMARY
  plus_bmi    = BMI_covars,
  plus_clinical = CLIN_covars
)

sleep_cvd_ladder <- bind_rows(lapply(names(ladder_specs), function(nm) {
  r <- fit_cox("sleep_disorder", "mort_cvd", ladder_specs[[nm]], keep = "^sleep_disorder$")
  if (is.null(r) || !nrow(r)) return(NULL)
  data.frame(
    outcome   = "CV mortality",
    exposure  = "sleep_disorder (score>=2)",
    model     = nm,
    estimator = "cause-specific Cox (survey-weighted)",
    HR        = r$estimate[1], CI_low = r$conf.low[1], CI_high = r$conf.high[1],
    P         = r$p.value[1],
    n_covars  = length(ladder_specs[[nm]]),
    is_primary = (nm == "confounder"),
    stringsAsFactors = FALSE
)
}))
cat("Sleep -> CV cause-specific Cox ladder:\n")
print(sleep_cvd_ladder[, c("model","HR","CI_low","CI_high","P","is_primary")], digits = 4)

# A2. Fine-Gray subdistribution hazard (CONFOUNDER-ONLY), competing non-CV death
cat("\n--- Fine-Gray subdistribution hazard (confounder-only, primary CV estimator) ---\n")
fg_dat <- nhanes_final %>%
  select(sleep_disorder, mort_cvd, mort_allcause, permth, all_of(C0_covars)) %>%
  na.omit
fac <- sapply(fg_dat, is.factor); fg_dat[fac] <- lapply(fg_dat[fac], droplevels)

# event coding: 1 = CV death, 2 = non-CV death, 0 = alive/censored
fg_dat$status_fg <- with(fg_dat,
  ifelse(mort_cvd == 1, 1L, ifelse(mort_allcause == 1 & mort_cvd == 0, 2L, 0L)))
cat("status_fg distribution (0=alive,1=CV death,2=non-CV death):\n")
print(table(fg_dat$status_fg))

Xfg <- model.matrix(~ sleep_disorder + age + race + education + pir +
                      smoke_status + alcohol_any + parity, data = fg_dat)[, -1]
set.seed(20260610)
fg_fit <- crr(ftime = fg_dat$permth, fstatus = fg_dat$status_fg,
              cov1 = Xfg, failcode = 1L, cencode = 0L)
fg_s <- summary(fg_fit)
fg_ci <- fg_s$conf.int; fg_co <- fg_s$coef
sl_row <- rownames(fg_ci)[grepl("sleep_disorder", rownames(fg_ci))][1]
fg_SHR <- fg_ci[sl_row, "exp(coef)"]
fg_lo  <- fg_ci[sl_row, "2.5%"]
fg_hi  <- fg_ci[sl_row, "97.5%"]
fg_p   <- fg_co[sl_row, "p-value"]
cat(sprintf("Fine-Gray SHR (sleep->CV, confounder-only): %.3f (%.3f-%.3f) P=%.4f\n",
            fg_SHR, fg_lo, fg_hi, fg_p))

# A3. E-value for the PRIMARY (confounder-only cause-specific Cox)
# CV death is a RARE outcome here (249/10923 = 2.3% prevalence), so HR ~ RR and
# we use rare=TRUE (VanderWeele: rare-outcome E-value treats HR directly as RR).
prim_row <- sleep_cvd_ladder[sleep_cvd_ladder$is_primary, ]
cv_prev  <- n_cv_events_tmp <- sum(nhanes_final$mort_cvd, na.rm = TRUE) /
            sum(!is.na(nhanes_final$mort_cvd))
ev <- evalues.HR(est = prim_row$HR, lo = prim_row$CI_low, hi = prim_row$CI_high, rare = TRUE)
ev_point <- ev[2, 1]; ev_ci <- ev[2, 2]
cat(sprintf("\nCV-death prevalence = %.4f (rare) -> E-value uses rare=TRUE\n", cv_prev))
cat(sprintf("E-value (point) = %.2f ; E-value (CI bound) = %.2f\n", ev_point, ev_ci))

# A4. EPV for the PRIMARY model
# parameters = number of estimated coefficients in confounder-only model
prim_form <- as.formula(sprintf("Surv(permth, mort_cvd) ~ sleep_disorder%s", cov_str(C0_covars)))
prim_fit  <- svycoxph(prim_form, design = design)
n_params  <- length(coef(prim_fit))
n_cv_events <- sum(fg_dat$status_fg == 1L)          # CV deaths in primary complete-case set
epv <- n_cv_events / n_params
cat(sprintf("EPV: CV events = %d / parameters = %d = %.1f\n", n_cv_events, n_params, epv))

# ==================================================
# PART B — sleep_disorder -> ALL-CAUSE (companion, confounder-only)
# ==================================================
cat("\n--- PART B: sleep_disorder -> all-cause mortality (confounder-only companion) ---\n")
sleep_allcause <- fit_cox("sleep_disorder", "mort_allcause", C0_covars, keep = "^sleep_disorder$")
sleep_ac_row <- data.frame(
  outcome = "all-cause mortality", exposure = "sleep_disorder (score>=2)",
  model = "confounder", estimator = "cause-specific Cox (survey-weighted)",
  HR = sleep_allcause$estimate[1], CI_low = sleep_allcause$conf.low[1],
  CI_high = sleep_allcause$conf.high[1], P = sleep_allcause$p.value[1],
  n_covars = length(C0_covars), is_primary = FALSE, stringsAsFactors = FALSE
)
cat(sprintf("Sleep -> all-cause (confounder-only): HR=%.3f (%.3f-%.3f) P=%.4f\n",
            sleep_ac_row$HR, sleep_ac_row$CI_low, sleep_ac_row$CI_high, sleep_ac_row$P))

# --------------------------------------------------
# Assemble Table 1 primary (ladder + Fine-Gray + companion + E-value + EPV)
# --------------------------------------------------
fg_block <- data.frame(
  outcome = "CV mortality", exposure = "sleep_disorder (score>=2)",
  model = "confounder", estimator = "Fine-Gray subdistribution hazard",
  HR = fg_SHR, CI_low = fg_lo, CI_high = fg_hi, P = fg_p,
  n_covars = length(C0_covars), is_primary = TRUE, stringsAsFactors = FALSE
)

table1_primary <- bind_rows(sleep_cvd_ladder, fg_block, sleep_ac_row)
# annotate E-value + EPV on the primary cause-specific row
table1_primary$E_value_point <- NA_real_
table1_primary$E_value_CIbound <- NA_real_
table1_primary$EPV <- NA_real_
table1_primary$N_events_CV <- NA_integer_
pidx <- which(table1_primary$is_primary &
              grepl("cause-specific", table1_primary$estimator))
table1_primary$E_value_point[pidx]   <- ev_point
table1_primary$E_value_CIbound[pidx] <- ev_ci
table1_primary$EPV[pidx]             <- epv
table1_primary$N_events_CV[pidx]     <- n_cv_events
# also annotate EPV/N on Fine-Gray primary row
fidx <- which(table1_primary$estimator == "Fine-Gray subdistribution hazard")
table1_primary$EPV[fidx]         <- epv
table1_primary$N_events_CV[fidx] <- n_cv_events

write.csv(table1_primary, "output/tables/table1_primary_sleep_cvd.csv", row.names = FALSE)
cat("\n[OK] output/tables/table1_primary_sleep_cvd.csv\n")

# ==================================================
# PART C — SECONDARY: DII -> all-cause (confounder-only)
# ==================================================
cat("\n--- PART C: DII -> all-cause (SECONDARY, confounder-only) ---\n")

# C1. per-SD
dii_ac_sd <- fit_cox("DII_z", "mort_allcause", C0_covars, keep = "^DII_z$")
cat(sprintf("DII per-SD -> all-cause: HR=%.3f (%.3f-%.3f) P=%.4f  [SD=%.3f]\n",
            dii_ac_sd$estimate[1], dii_ac_sd$conf.low[1], dii_ac_sd$conf.high[1],
            dii_ac_sd$p.value[1], dii_sd))

# C2. quartiles (Q2/Q3/Q4 vs Q1) + P-trend
dii_ac_q <- fit_cox("DII_Q", "mort_allcause", C0_covars, keep = "^DII_Q")
ptrend_fit <- svycoxph(as.formula(sprintf("Surv(permth, mort_allcause) ~ DII_Q_numeric%s",
                                          cov_str(C0_covars))), design = design)
ptrend_p <- broom::tidy(ptrend_fit) %>% filter(term == "DII_Q_numeric") %>% pull(p.value)

# C3. DII -> CV (expected NULL, report honestly) per-SD + quartiles
dii_cv_sd <- fit_cox("DII_z", "mort_cvd", C0_covars, keep = "^DII_z$")
cat(sprintf("DII per-SD -> CV (expected null): HR=%.3f (%.3f-%.3f) P=%.4f\n",
            dii_cv_sd$estimate[1], dii_cv_sd$conf.low[1], dii_cv_sd$conf.high[1],
            dii_cv_sd$p.value[1]))

table2_secondary <- bind_rows(
  data.frame(outcome = "all-cause", exposure = "DII per 1 SD", contrast = "per SD",
             HR = dii_ac_sd$estimate[1], CI_low = dii_ac_sd$conf.low[1],
             CI_high = dii_ac_sd$conf.high[1], P = dii_ac_sd$p.value[1],
             stringsAsFactors = FALSE),
  data.frame(outcome = "all-cause", exposure = "DII quartile",
             contrast = sub("DII_Q", "", dii_ac_q$term),
             HR = dii_ac_q$estimate, CI_low = dii_ac_q$conf.low,
             CI_high = dii_ac_q$conf.high, P = dii_ac_q$p.value,
             stringsAsFactors = FALSE),
  data.frame(outcome = "all-cause", exposure = "DII quartile", contrast = "P-trend",
             HR = NA, CI_low = NA, CI_high = NA, P = ptrend_p, stringsAsFactors = FALSE),
  data.frame(outcome = "CV (null)", exposure = "DII per 1 SD", contrast = "per SD",
             HR = dii_cv_sd$estimate[1], CI_low = dii_cv_sd$conf.low[1],
             CI_high = dii_cv_sd$conf.high[1], P = dii_cv_sd$p.value[1],
             stringsAsFactors = FALSE)
)
table2_secondary$model     <- "confounder-only"
table2_secondary$estimator <- "cause-specific Cox (survey-weighted)"
table2_secondary$DII_SD    <- dii_sd
write.csv(table2_secondary, "output/tables/table2_secondary_dii.csv", row.names = FALSE)
cat("[OK] output/tables/table2_secondary_dii.csv\n")
print(table2_secondary[, c("outcome","exposure","contrast","HR","CI_low","CI_high","P")], digits = 4)

# ==================================================
# PART D — MULTIPLICITY: 4-test PRIMARY family
#   1 sleep->CV (cause-specific, confounder-only)  <- headline
#   2 sleep->all-cause (confounder-only)
#   3 DII->all-cause (per SD, confounder-only)
#   4 DII->CV (per SD, confounder-only)
# ==================================================
cat("\n--- PART D: multiplicity over 4-test PRIMARY family ---\n")
fam_p <- c(
  `sleep->CV (cause-specific)`   = prim_row$P,
  `sleep->all-cause`             = sleep_ac_row$P,
  `DII->all-cause (per SD)`      = dii_ac_sd$p.value[1],
  `DII->CV (per SD)`             = dii_cv_sd$p.value[1]
)
m <- length(fam_p)
table3_mult <- data.frame(
  test       = names(fam_p),
  HR         = c(prim_row$HR, sleep_ac_row$HR, dii_ac_sd$estimate[1], dii_cv_sd$estimate[1]),
  CI_low     = c(prim_row$CI_low, sleep_ac_row$CI_low, dii_ac_sd$conf.low[1], dii_cv_sd$conf.low[1]),
  CI_high    = c(prim_row$CI_high, sleep_ac_row$CI_high, dii_ac_sd$conf.high[1], dii_cv_sd$conf.high[1]),
  p_raw      = as.numeric(fam_p),
  p_bonferroni = pmin(1, as.numeric(fam_p) * m),
  p_BH_FDR   = p.adjust(as.numeric(fam_p), method = "BH"),
  family_size = m,
  stringsAsFactors = FALSE
)
table3_mult$survives_bonferroni_0.05 <- table3_mult$p_bonferroni < 0.05
write.csv(table3_mult, "output/tables/table3_multiplicity.csv", row.names = FALSE)
cat("[OK] output/tables/table3_multiplicity.csv\n")
print(table3_mult, digits = 4)

# ==================================================
# PART E — VERIFICATION (print every primary number + gates)
# ==================================================
verify <- c(
  "================ PRIMARY VERIFICATION (primary) ================",
  sprintf("Analytic N (design)           : %d", nrow(nhanes_final)),
  sprintf("Fine-Gray complete-case N     : %d", nrow(fg_dat)),
  sprintf("CV deaths (primary set)       : %d", n_cv_events),
  sprintf("Non-CV deaths (competing)     : %d", sum(fg_dat$status_fg == 2L)),
  sprintf("All-cause deaths              : %d", sum(nhanes_final$mort_allcause, na.rm = TRUE)),
  "",
  "PRIMARY: sleep_disorder -> CV mortality (CONFOUNDER-ONLY)",
  sprintf("  Cause-specific Cox HR       : %.3f (%.3f-%.3f) P=%.4g",
          prim_row$HR, prim_row$CI_low, prim_row$CI_high, prim_row$P),
  sprintf("  Fine-Gray SHR               : %.3f (%.3f-%.3f) P=%.4g",
          fg_SHR, fg_lo, fg_hi, fg_p),
  sprintf("  E-value (point / CI bound)  : %.2f / %.2f", ev_point, ev_ci),
  sprintf("  EPV (CV events / params)    : %d / %d = %.1f", n_cv_events, n_params, epv),
  "  Sensitivity ladder (cause-specific Cox HR):",
  sprintf("    minimal (age+race)        : %.3f (%.3f-%.3f) P=%.4g",
          sleep_cvd_ladder$HR[sleep_cvd_ladder$model=="minimal"],
          sleep_cvd_ladder$CI_low[sleep_cvd_ladder$model=="minimal"],
          sleep_cvd_ladder$CI_high[sleep_cvd_ladder$model=="minimal"],
          sleep_cvd_ladder$P[sleep_cvd_ladder$model=="minimal"]),
  sprintf("    confounder-only (PRIMARY) : %.3f (%.3f-%.3f) P=%.4g",
          prim_row$HR, prim_row$CI_low, prim_row$CI_high, prim_row$P),
  sprintf("    + BMI                     : %.3f (%.3f-%.3f) P=%.4g",
          sleep_cvd_ladder$HR[sleep_cvd_ladder$model=="plus_bmi"],
          sleep_cvd_ladder$CI_low[sleep_cvd_ladder$model=="plus_bmi"],
          sleep_cvd_ladder$CI_high[sleep_cvd_ladder$model=="plus_bmi"],
          sleep_cvd_ladder$P[sleep_cvd_ladder$model=="plus_bmi"]),
  sprintf("    + clinical                : %.3f (%.3f-%.3f) P=%.4g",
          sleep_cvd_ladder$HR[sleep_cvd_ladder$model=="plus_clinical"],
          sleep_cvd_ladder$CI_low[sleep_cvd_ladder$model=="plus_clinical"],
          sleep_cvd_ladder$CI_high[sleep_cvd_ladder$model=="plus_clinical"],
          sleep_cvd_ladder$P[sleep_cvd_ladder$model=="plus_clinical"]),
  "",
  "COMPANION: sleep_disorder -> all-cause (confounder-only)",
  sprintf("  Cause-specific Cox HR       : %.3f (%.3f-%.3f) P=%.4g",
          sleep_ac_row$HR, sleep_ac_row$CI_low, sleep_ac_row$CI_high, sleep_ac_row$P),
  "",
  "SECONDARY: DII -> all-cause (confounder-only)",
  sprintf("  per 1 SD (SD=%.3f)          : HR=%.3f (%.3f-%.3f) P=%.4g",
          dii_sd, dii_ac_sd$estimate[1], dii_ac_sd$conf.low[1],
          dii_ac_sd$conf.high[1], dii_ac_sd$p.value[1]),
  sprintf("  Q4 vs Q1                    : HR=%.3f (%.3f-%.3f) P=%.4g",
          dii_ac_q$estimate[grepl("Q4", dii_ac_q$term)],
          dii_ac_q$conf.low[grepl("Q4", dii_ac_q$term)],
          dii_ac_q$conf.high[grepl("Q4", dii_ac_q$term)],
          dii_ac_q$p.value[grepl("Q4", dii_ac_q$term)]),
  sprintf("  P-trend                     : %.4g", ptrend_p),
  sprintf("  DII -> CV (null check)      : HR=%.3f (%.3f-%.3f) P=%.4g",
          dii_cv_sd$estimate[1], dii_cv_sd$conf.low[1],
          dii_cv_sd$conf.high[1], dii_cv_sd$p.value[1]),
  "",
  "MULTIPLICITY (4-test family, Bonferroni m=4):",
  sprintf("  sleep->CV        p_raw=%.4g  p_bonf=%.4g  survives=%s",
          table3_mult$p_raw[1], table3_mult$p_bonferroni[1], table3_mult$survives_bonferroni_0.05[1]),
  sprintf("  sleep->all-cause p_raw=%.4g  p_bonf=%.4g  survives=%s",
          table3_mult$p_raw[2], table3_mult$p_bonferroni[2], table3_mult$survives_bonferroni_0.05[2]),
  sprintf("  DII->all-cause   p_raw=%.4g  p_bonf=%.4g  survives=%s",
          table3_mult$p_raw[3], table3_mult$p_bonferroni[3], table3_mult$survives_bonferroni_0.05[3]),
  sprintf("  DII->CV          p_raw=%.4g  p_bonf=%.4g  survives=%s",
          table3_mult$p_raw[4], table3_mult$p_bonferroni[4], table3_mult$survives_bonferroni_0.05[4]),
  "",
  "GATES:",
  sprintf("  [GATE 1] sleep->CV survives Bonferroni(0.05)?  %s",
          ifelse(table3_mult$survives_bonferroni_0.05[1], "YES (PASS)", "NO (FAIL)")),
  sprintf("  [GATE 2] DII->CV null (p>0.05)?                %s",
          ifelse(dii_cv_sd$p.value[1] > 0.05, "YES (null, honest)", "NO")),
  sprintf("  [GATE 3] EPV >= 10 for primary?                %s (EPV=%.1f)",
          ifelse(epv >= 10, "YES (well-powered)", "NO"), epv),
  "=========================================================================="
)
writeLines(verify, "output/tables/_primary_verification.txt")
cat("\n")
cat(paste(verify, collapse = "\n"))
cat("\n\n[OK] output/tables/_primary_verification.txt\n")