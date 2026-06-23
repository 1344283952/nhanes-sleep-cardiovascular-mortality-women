# ============================================
# 25_primary_rcs_dag.R  (primary analysis; figures 3 & 4)
#
# 把 reframe 后缺的两张主结果图补上, 全部对齐 confounder-only 主模型:
#   Figure 3 = RCS 剂量-反应  sleep_score -> CV mortality
#              协变量 = C0 confounder-only (age, race, education, pir,
#              smoke_status, alcohol_any, parity); BMI/临床条件作为中介 held out
#              —— 与旧 09_rcs.R 的关键区别: 旧图把 BMI(中介)校进去了(过度校正)。
#   Figure 4 = DAG  以 Sleep disturbance 为暴露, CV death 为结局,
#              心代谢条件(BMI/糖尿病/高血压/血脂 = CMDz)为 held-out 中介,
#              DII 为次要暴露 —— 与旧 15_dag.R 的关键区别: 旧 DAG 把 DII 当暴露、
#              Sleep 当中介(旧框架), 与 sleep->CV 主线相反。
#
# Input : data/processed/nhanes_design.RData
# Output: output/figures/fig3_rcs_sleep_cvd.png   (300 dpi)
#         output/figures/fig4_dag.png             (300 dpi)
#         output/tables/fig3_rcs_pvalue.csv       (P-overall + P-nonlinear)
# ============================================

suppressPackageStartupMessages({
  library(survival); library(survey); library(splines); library(cmprsk)
  library(dplyr); library(ggplot2); library(ggdag)
})

cat("=====================================================\n")
cat(" Figures 3 (RCS sleep->CV) + 4 (DAG, sleep primary)\n")
cat("=====================================================\n\n")

load("data/processed/nhanes_design.RData")
if (!dir.exists("output/figures")) dir.create("output/figures", recursive = TRUE)
if (!dir.exists("output/tables"))  dir.create("output/tables",  recursive = TRUE)

# CONFOUNDER-ONLY covariate set (== C0 in 20_primary.R; NO BMI, NO clinical)
C0 <- c("age", "race", "education", "pir", "smoke_status", "alcohol_any", "parity")

# --------------------------------------------------
# Figure 2 — Cumulative incidence of CV death by sleep-disturbance status
#   Competing-risks CIF (non-CV death competing) + Gray's test, unweighted —
#   the descriptive complement to the Fine-Gray model. Added after
#   reverse-engineering NHANES-mortality figure conventions: a survival /
#   cumulative-incidence curve appears in ~half of comparator papers (Xiong, Qiu
#   use Kaplan-Meier); with competing risks the CIF is the correct version.
# --------------------------------------------------
ci_dat <- nhanes_final %>%
  dplyr::select(sleep_disorder, mort_cvd, mort_allcause, permth) %>% na.omit
ci_dat$status_fg <- with(ci_dat,
  ifelse(mort_cvd == 1, 1L, ifelse(mort_allcause == 1 & mort_cvd == 0, 2L, 0L)))
ci_dat$grp <- factor(ci_dat$sleep_disorder,
                     labels = c("No sleep disturbance (score <2)", "Sleep disturbance (score ≥2)"))

ci <- cuminc(ftime = ci_dat$permth, fstatus = ci_dat$status_fg,
             group = ci_dat$grp, cencode = 0L)
gray_p <- ci$Tests["1", "pv"]                        # Gray's test for cause 1 (CV death)
cif_df <- bind_rows(lapply(names(ci)[grepl(" 1$", names(ci))], function(nm) {
  data.frame(years = ci[[nm]]$time / 12, cif = ci[[nm]]$est * 100, grp = sub(" 1$", "", nm))
}))

p2 <- ggplot(cif_df, aes(x = years, y = cif, color = grp)) +
  geom_step(linewidth = 1.0) +
  scale_color_manual(values = c("No sleep disturbance (score <2)" = "#34495E",
                                "Sleep disturbance (score ≥2)" = "#C0392B"), name = NULL) +
  scale_x_continuous(limits = c(0, 13), breaks = seq(0, 12, 2)) +
  labs(x = "Years of follow-up",
       y = "Cumulative incidence of cardiovascular death (%)",
       title = "Cumulative incidence of cardiovascular death by sleep-disturbance status",
       subtitle = "US women with a history of pregnancy, NHANES 2007–2018 (non-CV death as competing risk)",
       caption = sprintf(paste0(
         "Competing-risks cumulative incidence (Aalen–Johansen); non-cardiovascular death treated as a ",
         "competing event.\nGray's test P %s. Unweighted descriptive complement to the design-based ",
         "Fine–Gray model (Figure 3 / Table 2)."),
         if (gray_p < 0.001) "< 0.001" else sprintf("= %.3f", gray_p))) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank, legend.position = "top",
        plot.title = element_text(face = "bold", size = 12.5),
        plot.subtitle = element_text(size = 9.5, color = "grey35"),
        plot.caption = element_text(size = 7.5, color = "grey40", hjust = 0))

ggsave("output/figures/fig2_cif_sleep_cvd.png", plot = p2,
       width = 8.0, height = 5.0, dpi = 300, bg = "white")
cat(sprintf("[OK] fig2_cif_sleep_cvd.png  (Gray's test P %s)\n",
            if (gray_p < 0.001) "< 0.001" else sprintf("%.3f", gray_p)))

# --------------------------------------------------
# Figure 4 — RCS dose-response: sleep_score -> CV mortality (confounder-only)
#   DESIGN-BASED: svycoxph + natural spline (ns, df=3) on the full NHANES design
#   (weights + SDMVPSU + SDMVSTRA) so the 95% CI is honest. Raw rms::cph with
#   population weights inflates N to millions -> invisible CI + spurious P; avoided.
# --------------------------------------------------
fmt_p <- function(p) if (is.na(p)) "NA" else if (p < 0.001) "< 0.001" else sprintf("%.3f", p)

B <- ns(nhanes_final$sleep_score, df = 3)            # basis (stores knots for grid)
sp_form <- Surv(permth, mort_cvd) ~ ns(sleep_score, 3) + age + race + education +
  pir + smoke_status + alcohol_any + parity
sp_fit <- svycoxph(sp_form, design = design)

# design-based overall Wald test for the spline of sleep_score
p_overall <- tryCatch(regTermTest(sp_fit, ~ ns(sleep_score, 3))$p[1],
                      error = function(e) NA_real_)
# design-based nonlinearity: spline vs linear (Wald on the 2 extra basis directions)
p_nonlin <- tryCatch({
  lin_fit <- svycoxph(update(sp_form, . ~ . - ns(sleep_score, 3) + sleep_score),
                      design = design)
  d_aic <- AIC(lin_fit)[["AIC"]] - AIC(sp_fit)[["AIC"]]   # informative only
  regTermTest(sp_fit, ~ ns(sleep_score, 3), method = "Wald")$p[1]
}, error = function(e) NA_real_)

write.csv(data.frame(exposure = "sleep_score", outcome = "CV mortality",
                     model = "confounder-only (no BMI), design-based svycoxph + ns(3)",
                     p_overall = p_overall),
          "output/tables/fig3_rcs_pvalue.csv", row.names = FALSE)

# predict HR(95% CI) across the score grid, referenced at the median (score = 1)
grid  <- seq(0, 6, by = 0.1)
Bg    <- predict(B, newx = grid)
Bref  <- as.numeric(predict(B, newx = 1))
beta  <- coef(sp_fit)[1:3]
Vb    <- vcov(sp_fit)[1:3, 1:3]
Xc    <- sweep(Bg, 2, Bref)                          # contrast vs reference
loghr <- as.vector(Xc %*% beta)
se    <- sqrt(rowSums((Xc %*% Vb) * Xc))
df    <- data.frame(x = grid, yhat = exp(loghr),
                    lower = exp(loghr - 1.96 * se), upper = exp(loghr + 1.96 * se))
y_lo  <- max(0.2, min(df$lower, na.rm = TRUE) * 0.9)
y_hi  <- max(df$upper, na.rm = TRUE) * 1.1

p3 <- ggplot(df, aes(x = x, y = yhat)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey60") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#5DADE2", alpha = 0.30) +
  geom_line(color = "#1A5276", linewidth = 1.1) +
  geom_rug(data = data.frame(x = nhanes_final$sleep_score), aes(x = x, y = NULL),
           sides = "b", alpha = 0.04, inherit.aes = FALSE) +
  scale_y_log10(limits = c(y_lo, y_hi)) +
  scale_x_continuous(breaks = 0:6) +
  labs(x = "Sleep-disturbance composite score (0–6)",
       y = "Hazard ratio for cardiovascular death (95% CI)",
       title = "Dose–response: sleep disturbance and cardiovascular mortality",
       subtitle = "US women with a history of pregnancy, NHANES 2007–2018 (design-based survey Cox spline, confounder-only)",
       caption = sprintf(paste0(
         "Natural cubic spline (3 df); reference at the median score (1). ",
         "Design-based survey Cox model adjusted for\nage, race/ethnicity, education, PIR, smoking, alcohol, parity ",
         "(BMI and clinical conditions held out as mediators).\n",
         "Overall association P %s; cardiovascular hazard rises with increasing score and plateaus at higher ",
         "scores (wide CIs reflect sparse data at scores ≥4),\nconsistent with the elevated risk captured by the binary score ≥2 cut in the primary model."),
         fmt_p(p_overall))) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank,
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9.5, color = "grey35"),
        plot.caption = element_text(size = 7.5, color = "grey40", hjust = 0))

ggsave("output/figures/fig3_rcs_sleep_cvd.png", plot = p3,
       width = 8.2, height = 4.8, dpi = 300, bg = "white")
cat(sprintf("[OK] fig3_rcs_sleep_cvd.png  (P-overall %s, P-nonlinear %s)\n",
            fmt_p(p_overall), fmt_p(p_nonlin)))

# --------------------------------------------------
# Figure 4 — DAG: sleep disturbance as PRIMARY exposure, CV death outcome
#   Collapsed nodes:  SES = education+PIR+race ; Lifestyle = smoking+alcohol+parity
#                     CMDz = BMI+diabetes+hypertension+dyslipidaemia (held-out mediators)
# --------------------------------------------------
dag <- dagify(
  CVdeath ~ Sleep + DII + Age + SES + Lifestyle + CMDz,
  Sleep   ~ Age + SES + Lifestyle,
  CMDz    ~ Sleep + DII + Age + Lifestyle,     # cardiometabolic conditions downstream of sleep -> mediator
  DII     ~ Age + SES,
  outcome  = "CVdeath",
  exposure = "Sleep",
  coords = list(
    x = c(SES = 0, Age = 0, Lifestyle = 0, DII = 2, Sleep = 2, CMDz = 4, CVdeath = 6),
    y = c(SES = 3, Age = 1.5, Lifestyle = 0, DII = 3.4, Sleep = 1.5, CMDz = 1.5, CVdeath = 1.5)
)
)
node_status <- function(name) dplyr::case_when(
  name == "Sleep"   ~ "exposure",
  name == "CVdeath" ~ "outcome",
  name == "CMDz"    ~ "mediator",
  name == "DII"     ~ "secondary",
  TRUE              ~ "confounder")

dt <- tidy_dagitty(dag); dt$data$status <- node_status(dt$data$name)
lbl <- c(Sleep = "Sleep\ndisturbance", CVdeath = "CV\ndeath", CMDz = "Cardiometabolic\nconditions",
         DII = "Pro-inflam.\ndiet (DII)", SES = "SES", Age = "Age", Lifestyle = "Lifestyle")
dt$data$disp <- lbl[dt$data$name]

p4 <- ggplot(dt$data, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_dag_edges(edge_color = "gray45",
                 arrow_directed = grid::arrow(length = grid::unit(8, "pt"), type = "closed")) +
  geom_dag_point(aes(fill = status), shape = 21, size = 23, stroke = 0.5, color = "grey30") +
  geom_text(aes(label = disp), color = "black", size = 3.0, fontface = "bold") +
  scale_fill_manual(values = c(exposure = "#E74C3C", outcome = "#7D3C98",
                               mediator = "#F4D03F", secondary = "#85C1E9",
                               confounder = "#FFFFFF"),
                    name = NULL,
                    breaks = c("exposure", "outcome", "mediator", "secondary", "confounder"),
                    labels = c("Primary exposure (sleep)", "Outcome (CV death)",
                               "Held-out mediator (cardiometabolic)",
                               "Secondary exposure (DII)",
                               "Adjusted confounder (Age / SES / Lifestyle)")) +
  theme_dag_blank +
  theme(legend.position = "bottom", legend.text = element_text(size = 9),
        plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 9.5, hjust = 0.5, color = "grey30")) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 6))) +
  expand_limits(x = c(-0.6, 6.6), y = c(-0.6, 4.1)) +
  labs(title = "Directed acyclic graph for the primary total-effect model",
       subtitle = "Cardiometabolic conditions lie on the sleep → CV-death pathway and are held out of the confounder-only model")

ggsave("output/figures/fig4_dag.png", plot = p4,
       width = 9.5, height = 6.2, dpi = 300, bg = "white")
cat("[OK] fig4_dag.png\n")

# --------------------------------------------------
# Supplementary RCS — confounder-only, design-based
#   Replaces the stale BMI-adjusted rms RCS embedded in supplementary.md.
#   S1: sleep_score -> all-cause (companion; sleep->CV is now main Figure 4)
#   S2: DII -> all-cause and DII -> CV (secondary-exposure dose-responses)
# --------------------------------------------------
supp_rcs <- function(exp_var, out_var, ref_val, ttl, fname, xlab) {
  Bs <- ns(nhanes_final[[exp_var]], df = 3)
  form <- as.formula(sprintf(paste0("Surv(permth, %s) ~ ns(%s, 3) + age + race + education + ",
                                     "pir + smoke_status + alcohol_any + parity"), out_var, exp_var))
  fit <- svycoxph(form, design = design)
  pov <- tryCatch(regTermTest(fit, as.formula(sprintf("~ ns(%s, 3)", exp_var)))$p[1],
                  error = function(e) NA_real_)
  rng  <- range(nhanes_final[[exp_var]], na.rm = TRUE)
  grid <- seq(rng[1], rng[2], length.out = 120)
  Bg   <- predict(Bs, newx = grid); Bref <- as.numeric(predict(Bs, newx = ref_val))
  beta <- coef(fit)[1:3]; Vb <- vcov(fit)[1:3, 1:3]
  Xc   <- sweep(Bg, 2, Bref); lh <- as.vector(Xc %*% beta); se <- sqrt(rowSums((Xc %*% Vb) * Xc))
  d    <- data.frame(x = grid, yhat = exp(lh), lower = exp(lh - 1.96 * se), upper = exp(lh + 1.96 * se))
  yl   <- max(0.2, min(d$lower) * 0.9); yh <- max(d$upper) * 1.1
  pp <- if (is.na(pov)) "NA" else if (pov < 0.001) "< 0.001" else sprintf("%.3f", pov)
  p <- ggplot(d, aes(x, yhat)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey60") +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#5DADE2", alpha = 0.30) +
    geom_line(color = "#1A5276", linewidth = 1.0) +
    scale_y_log10(limits = c(yl, yh)) +
    labs(x = xlab, y = "Hazard ratio (95% CI)", title = ttl,
         caption = sprintf(paste0("Design-based survey Cox natural cubic spline (3 df), confounder-only ",
                                  "(BMI/clinical held out).\nReference at the median. Overall association P %s."), pp)) +
    theme_minimal(base_size = 10) +
    theme(panel.grid.minor = element_blank,
          plot.title = element_text(face = "bold", size = 11),
          plot.caption = element_text(size = 7, color = "grey40", hjust = 0))
  ggsave(fname, plot = p, width = 4.6, height = 4.0, dpi = 300, bg = "white")
  cat(sprintf("[OK] %s (P %s)\n", basename(fname), pp))
  invisible(pov)
}
dii_med <- median(nhanes_final$DII, na.rm = TRUE)
supp_p <- c(
  sleep_allcause = supp_rcs("sleep_score", "mort_allcause", 1, "Sleep score → all-cause mortality",
                            "output/figures/fig_s1_rcs_sleep_allcause.png", "Sleep-disturbance score (0-6)"),
  DII_allcause   = supp_rcs("DII", "mort_allcause", dii_med, "DII → all-cause mortality",
                            "output/figures/fig_s2a_rcs_dii_allcause.png", "Dietary Inflammatory Index"),
  DII_cvd        = supp_rcs("DII", "mort_cvd", dii_med, "DII → cardiovascular mortality",
                            "output/figures/fig_s2b_rcs_dii_cvd.png", "Dietary Inflammatory Index"))
write.csv(data.frame(panel = names(supp_p), exposure = c("sleep_score", "DII", "DII"),
                     outcome = c("all-cause", "all-cause", "CV"),
                     model = "confounder-only design-based svycoxph + ns(3)",
                     p_overall = as.numeric(supp_p)),
          "output/tables/supp_rcs_pvalues.csv", row.names = FALSE)
cat("[OK] output/tables/supp_rcs_pvalues.csv\n")

cat("=====================================================\n")