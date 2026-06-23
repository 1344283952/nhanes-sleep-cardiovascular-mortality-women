# ============================================
# 03_clean_data.R (版)
# DII 27 项 + 睡眠紊乱评分 + GDM 史 + mortality outcome 编码
#
# 输入: data/processed/nhanes_raw_merged.RData
# 输出: data/processed/nhanes_final.RData + output/tables/flow_counts.csv
#
# 与 v1  差异:
#   - filter: 当前妊娠 (RIDEXPRG==1) → 曾妊娠 (RHQ160≥1) + GDM 报告 (RHQ162∈{1,2,3})
#   - outcome: PIH 复合 → mortality (全因 + CVD + 5 自报 CVD)
#   - 年龄: ≥18（v1 隐含）→ ≥20 (NHANES 标准成人)
#   - 不剔除 baseline CVD (因为是 mortality outcome；做敏感性 S2 排除)
#   - 周期: 5 → 6 (2007-2018)
#   - WTDR_5YR → WTDR_6YR
# ============================================

library(dplyr)
library(tidyr)
library(purrr)

cat("========================================\n")
cat("03_clean_data: DII + Sleep + GDM + mortality\n")
cat("========================================\n\n")

load("data/processed/nhanes_raw_merged.RData")
cat(sprintf("原始合并表: %d 行 × %d 列\n\n", nrow(nhanes_all), ncol(nhanes_all)))

# Flow log
flow <- list
log_flow <- function(label, n) {
  cat(sprintf("  [流程] %-55s n = %d\n", label, n))
  flow[[length(flow) + 1]] <<- data.frame(step = label, n = n,
                                          stringsAsFactors = FALSE)
}
log_flow("原始合并表（2007-2018 六周期）", nrow(nhanes_all))

# --------------------------------------------------
# Util
# --------------------------------------------------
coalesce_cols <- function(df, cols) {
  cols <- intersect(cols, names(df))
  if (length(cols) == 0) return(rep(NA_real_, nrow(df)))
  if (length(cols) == 1) return(df[[cols[1]]])
  Reduce(function(a, b) ifelse(is.na(a), b, a),
         lapply(cols, function(c) df[[c]]))
}
na_codes  <- function(x, codes) ifelse(x %in% codes, NA, x)
zero_to_na <- function(x) ifelse(!is.na(x) & x == 0, NA, x)

# --------------------------------------------------
# Step 0: 跨周期统一变量
# --------------------------------------------------
nhanes_all <- nhanes_all %>%
  mutate(
    sleep_hours = coalesce_cols(., c("SLD010H", "SLD012")),
    hdl_mgdl    = coalesce_cols(., c("LBDHDL", "LBXHDD", "LBDHDD")),
    vit_d_serum = coalesce_cols(., c("LBXVIDMS", "LBXVDMS",
                                     "LBXVD2MS", "LBXVD3MS"))
)

# --------------------------------------------------
# Step 1: 入选 — 女性 + 年龄 ≥ 20 + 曾妊娠 + GDM 报告完整
# --------------------------------------------------
df <- nhanes_all %>% filter(RIAGENDR == 2)
log_flow("女性 (RIAGENDR == 2)", nrow(df))

df <- df %>% filter(RIDAGEYR >= 20)
log_flow("年龄 ≥ 20", nrow(df))

df <- df %>%
  mutate(RHQ160 = na_codes(RHQ160, c(77, 99))) %>%
  filter(!is.na(RHQ160), RHQ160 >= 1)
log_flow("曾妊娠 (RHQ160 ≥ 1)", nrow(df))

df <- df %>%
  mutate(RHQ162 = na_codes(RHQ162, c(7, 9))) %>%
  filter(RHQ162 %in% c(1, 2, 3))
log_flow("GDM 史报告完整 (RHQ162 ∈ {1,2,3})", nrow(df))

# --------------------------------------------------
# Step 2: 膳食 + 睡眠数据可用性
# --------------------------------------------------
df <- df %>%
  filter(
    !is.na(DR1TKCAL),
    DR1TKCAL >= 500 & DR1TKCAL <= 5000,
    !is.na(sleep_hours) | !is.na(SLQ050) | !is.na(SLQ060)
)
log_flow("DR1TOT kcal[500,5000] + 至少 1 个 SLQ 可用", nrow(df))

# --------------------------------------------------
# Step 3: DII 27 项算法（Shivappa 2014, Public Health Nutr Table 4, v1 verbatim）
# --------------------------------------------------
dii_table <- data.frame(
  name = c("energy","protein","carbohydrate","total_fat","fiber","cholesterol",
           "sfa","mufa","pufa","n3","n6","niacin","thiamin","riboflavin",
           "vit_b6","vit_b12","folate","vit_a","beta_carotene","vit_c",
           "vit_d","vit_e","iron","magnesium","zinc","selenium","caffeine"),
  global_mean = c(2056, 79.4, 272.2, 71.4, 18.8, 279.4,
                  28.6, 27.0, 13.88, 1.06, 10.8,
                  25.90, 1.70, 1.74, 1.47, 5.15, 273,
                  983.9, 3718, 118.2, 6.26, 8.73,
                  13.35, 310.1, 9.84, 67, 8.05),
  global_sd = c(338, 13.9, 40.0, 19.4, 4.9, 51.2,
                8.0, 6.1, 3.76, 1.06, 7.5,
                11.77, 0.66, 0.79, 0.74, 2.7, 70.7,
                518.6, 1720, 43.46, 2.21, 1.49,
                3.71, 139.4, 2.19, 25.1, 6.67),
  inflam = c(0.180, 0.021, 0.097, 0.298, -0.663, 0.110,
             0.373, -0.009, -0.337, -0.436, -0.159,
             -0.246, -0.098, -0.068, -0.365, 0.106, -0.190,
             -0.401, -0.584, -0.424, -0.446, -0.419,
             0.032, -0.484, -0.313, -0.191, -0.110),
  stringsAsFactors = FALSE
)

extract_intakes <- function(df) {
  data.frame(
    energy        = df$DR1TKCAL,
    protein       = df$DR1TPROT,
    carbohydrate  = df$DR1TCARB,
    total_fat     = df$DR1TTFAT,
    fiber         = df$DR1TFIBE,
    cholesterol   = df$DR1TCHOL,
    sfa           = df$DR1TSFAT,
    mufa          = df$DR1TMFAT,
    pufa          = df$DR1TPFAT,
    n3            = rowSums(df[, c("DR1TP183","DR1TP205","DR1TP225","DR1TP226")]),
    n6            = rowSums(df[, c("DR1TP182","DR1TP204")]),
    niacin        = df$DR1TNIAC,
    thiamin       = df$DR1TVB1,
    riboflavin    = df$DR1TVB2,
    vit_b6        = df$DR1TVB6,
    vit_b12       = df$DR1TVB12,
    folate        = df$DR1TFOLA,
    vit_a         = df$DR1TVARA,
    beta_carotene = df$DR1TBCAR,
    vit_c         = df$DR1TVC,
    vit_d         = df$DR1TVD,
    vit_e         = df$DR1TATOC,
    iron          = df$DR1TIRON,
    magnesium     = df$DR1TMAGN,
    zinc          = df$DR1TZINC,
    selenium      = df$DR1TSELE,
    caffeine      = df$DR1TCAFF
)
}

compute_dii <- function(df) {
  intakes <- extract_intakes(df)
  contrib <- matrix(NA_real_, nrow = nrow(intakes), ncol = nrow(dii_table))
  colnames(contrib) <- dii_table$name
  for (i in seq_len(nrow(dii_table))) {
    nm <- dii_table$name[i]
    x  <- intakes[[nm]]
    z  <- (x - dii_table$global_mean[i]) / dii_table$global_sd[i]
    p  <- 2 * pnorm(z) - 1
    contrib[, i] <- p * dii_table$inflam[i]
  }
  rowSums(contrib)
}

df$DII <- compute_dii(df)
df <- df %>% filter(!is.na(DII))
log_flow("DII 27 项全部可算", nrow(df))

# DII 三/四分位
dii_T <- quantile(df$DII, probs = c(1/3, 2/3), na.rm = TRUE)
dii_Q <- quantile(df$DII, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
df <- df %>%
  mutate(
    DII_T = cut(DII, breaks = c(-Inf, dii_T[1], dii_T[2], Inf),
                labels = c("T1", "T2", "T3"), include.lowest = TRUE),
    DII_Q = cut(DII, breaks = c(-Inf, dii_Q[1], dii_Q[2], dii_Q[3], Inf),
                labels = c("Q1", "Q2", "Q3", "Q4"), include.lowest = TRUE)
)

# --------------------------------------------------
# Step 4: 睡眠紊乱评分（6 维 0-6 分，v1 verbatim）
# --------------------------------------------------
df <- df %>%
  mutate(
    SLQ050 = na_codes(SLQ050, c(7, 9)),
    SLQ060 = na_codes(SLQ060, c(7, 9)),
    SLQ030 = na_codes(SLQ030, c(7, 9)),
    SLQ040 = na_codes(SLQ040, c(7, 9)),
    SLQ120 = na_codes(SLQ120, c(7, 9))
) %>%
  mutate(
    sl_short      = as.integer(!is.na(sleep_hours) & sleep_hours < 6),
    sl_trouble    = as.integer(!is.na(SLQ050) & SLQ050 == 1),
    sl_disorderDx = as.integer(!is.na(SLQ060) & SLQ060 == 1),
    sl_snore      = as.integer(!is.na(SLQ030) & SLQ030 >= 2),
    sl_apnea      = as.integer(!is.na(SLQ040) & SLQ040 >= 2),
    sl_daytime    = as.integer(!is.na(SLQ120) & SLQ120 >= 3)
)

df$sl_dims_obs <- rowSums(!is.na(df[, c("sleep_hours","SLQ050","SLQ060",
                                        "SLQ030","SLQ040","SLQ120")]))
df <- df %>% filter(sl_dims_obs >= 3)
log_flow("≥3 个睡眠维度可评估", nrow(df))

df <- df %>%
  mutate(
    sleep_score    = sl_short + sl_trouble + sl_disorderDx +
                     sl_snore + sl_apnea + sl_daytime,
    sleep_disorder = as.integer(sleep_score >= 2),
    sleep_disorder_strict = as.integer(sleep_score >= 3)
)

# --------------------------------------------------
# Step 5: GDM 史变量
# --------------------------------------------------
df <- df %>%
  mutate(
    gdm_history = factor(
      case_when(
        RHQ162 == 1 ~ "Yes",
        RHQ162 == 2 ~ "No",
        RHQ162 == 3 ~ "Borderline"
),
      levels = c("No", "Yes", "Borderline")
),
    gdm_binary = as.integer(RHQ162 == 1)
)

# --------------------------------------------------
# Step 6: Mortality outcome 编码
# --------------------------------------------------
df <- df %>%
  mutate(
    mort_allcause = ifelse(!is.na(MORTSTAT), as.integer(MORTSTAT == 1), NA),
    mort_cvd      = ifelse(!is.na(MORTSTAT),
                           as.integer(MORTSTAT == 1 &
                                      UCOD_LEADING %in% c(1, 5)),
                           NA),
    permth = ifelse(!is.na(PERMTH_EXM), PERMTH_EXM, PERMTH_INT)
)

# --------------------------------------------------
# Step 7: 自报 CVD outcome（5 个 + 复合）
# --------------------------------------------------
df <- df %>%
  mutate(
    MCQ160B = na_codes(MCQ160B, c(7, 9)),
    MCQ160C = na_codes(MCQ160C, c(7, 9)),
    MCQ160D = na_codes(MCQ160D, c(7, 9)),
    MCQ160E = na_codes(MCQ160E, c(7, 9)),
    MCQ160F = na_codes(MCQ160F, c(7, 9))
) %>%
  mutate(
    cvd_chf       = as.integer(!is.na(MCQ160B) & MCQ160B == 1),
    cvd_chd       = as.integer(!is.na(MCQ160C) & MCQ160C == 1),
    cvd_angina    = as.integer(!is.na(MCQ160D) & MCQ160D == 1),
    cvd_mi        = as.integer(!is.na(MCQ160E) & MCQ160E == 1),
    cvd_stroke    = as.integer(!is.na(MCQ160F) & MCQ160F == 1),
    cvd_composite = as.integer(
      (cvd_chf == 1) | (cvd_chd == 1) | (cvd_angina == 1) |
      (cvd_mi == 1)  | (cvd_stroke == 1)
)
)

# --------------------------------------------------
# Step 8: 协变量（参考 Ying et al. 2024）
# --------------------------------------------------
df <- df %>%
  mutate(
    DMDEDUC2 = na_codes(DMDEDUC2, c(7, 9)),
    DMDMARTL = na_codes(DMDMARTL, c(77, 99)),
    SMQ020   = na_codes(SMQ020,   c(7, 9)),
    ALQ101   = na_codes(ALQ101,   c(7, 9)),
    ALQ111   = na_codes(ALQ111,   c(7, 9)),
    BPQ020   = na_codes(BPQ020,   c(7, 9)),
    DIQ010   = na_codes(DIQ010,   c(7, 9))
) %>%
  mutate(
    age           = RIDAGEYR,
    age_group     = factor(case_when(
                      age >= 20 & age < 40 ~ "20-39",
                      age >= 40 & age < 60 ~ "40-59",
                      age >= 60            ~ ">=60"
),
                    levels = c("20-39", "40-59", ">=60")),
    race          = factor(recode(as.character(RIDRETH1),
                                  "1" = "Mexican American",
                                  "2" = "Other Hispanic",
                                  "3" = "Non-Hispanic White",
                                  "4" = "Non-Hispanic Black",
                                  "5" = "Other Race"),
                           levels = c("Non-Hispanic White",
                                      "Non-Hispanic Black",
                                      "Mexican American",
                                      "Other Hispanic",
                                      "Other Race")),
    education     = factor(case_when(
                      DMDEDUC2 %in% c(1, 2) ~ "Less than HS",
                      DMDEDUC2 == 3          ~ "High school",
                      DMDEDUC2 %in% c(4, 5) ~ "College or above"
),
                    levels = c("Less than HS", "High school", "College or above")),
    marital       = factor(case_when(
                      DMDMARTL %in% c(1, 6) ~ "Married or partnered",
                      DMDMARTL %in% c(2, 3, 4, 5) ~ "Not married"
),
                    levels = c("Married or partnered", "Not married")),
    pir           = INDFMPIR,
    pir_group     = factor(case_when(
                      pir <= 1.3              ~ "<=1.3",
                      pir > 1.3 & pir <= 3.5  ~ "1.3-3.5",
                      pir > 3.5               ~ ">3.5"
),
                    levels = c("<=1.3", "1.3-3.5", ">3.5")),
    bmi           = BMXBMI,
    bmi_cat       = factor(case_when(
                      bmi < 25            ~ "<25",
                      bmi >= 25 & bmi < 30 ~ "25-29.9",
                      bmi >= 30           ~ ">=30"
), levels = c("<25", "25-29.9", ">=30")),
    smoke_status  = factor(case_when(
                      SMQ020 == 2 ~ "Never",
                      SMQ020 == 1 ~ "Ever"
), levels = c("Never", "Ever")),
    alcohol_any   = factor(case_when(
                      ALQ101 == 1 | ALQ111 == 1 ~ "Yes",
                      ALQ101 == 2 | ALQ111 == 2 ~ "No"
), levels = c("No", "Yes")),
    parity        = ifelse(!is.na(RHQ160), RHQ160, NA),

    # 临床协变量 / 中介
    fbg           = LBXGLU,
    hba1c         = LBXGH,
    fbs_ins       = LBXIN,
    homa_ir       = ifelse(!is.na(LBXGLU) & !is.na(LBXIN),
                           (LBXGLU * LBXIN) / 405, NA),
    ldl           = coalesce_cols(., c("LBDLDL", "LBDLDLM")),
    hdl           = hdl_mgdl,
    tc            = LBXTC,
    tg            = coalesce_cols(., c("LBXTR", "LBDTRSI")),
    cr            = LBXSCR,

    # 自报糖尿病 / 高血压
    diabetes      = as.integer(
                      (DIQ010 == 1) |
                      (!is.na(hba1c) & hba1c >= 6.5) |
                      (!is.na(fbg) & fbg >= 126)),
    hypertension  = as.integer(BPQ020 == 1)
)

# 血压平均（去 0 占位）
df <- df %>%
  mutate(
    BPXSY1 = zero_to_na(BPXSY1), BPXSY2 = zero_to_na(BPXSY2),
    BPXSY3 = zero_to_na(BPXSY3), BPXSY4 = zero_to_na(BPXSY4),
    BPXDI1 = zero_to_na(BPXDI1), BPXDI2 = zero_to_na(BPXDI2),
    BPXDI3 = zero_to_na(BPXDI3), BPXDI4 = zero_to_na(BPXDI4)
) %>%
  mutate(
    sbp = rowMeans(dplyr::select(., BPXSY2, BPXSY3, BPXSY4), na.rm = TRUE),
    dbp = rowMeans(dplyr::select(., BPXDI2, BPXDI3, BPXDI4), na.rm = TRUE),
    sbp = ifelse(is.nan(sbp), NA, sbp),
    dbp = ifelse(is.nan(dbp), NA, dbp)
)

# --------------------------------------------------
# Step 9: Mortality 合格性
# --------------------------------------------------
df <- df %>% filter(!is.na(ELIGSTAT), ELIGSTAT == 1)
log_flow("mortality ELIGSTAT == 1", nrow(df))

# --------------------------------------------------
# Step 10: 6 周期合并膳食权重
# --------------------------------------------------
df <- df %>% mutate(WTDR_6YR = WTDRD1 / 6)

# --------------------------------------------------
# Step 11: 核心协变量最后剔除
# --------------------------------------------------
core_covars <- c("age", "race", "education", "pir", "bmi",
                 "smoke_status", "DII", "sleep_score",
                 "gdm_history", "mort_allcause",
                 "WTDR_6YR", "SDMVPSU", "SDMVSTRA")

before_final <- nrow(df)
df <- df %>% filter(if_all(all_of(core_covars), ~ !is.na(.)))
cat(sprintf("\n核心协变量缺失剔除: -%d → %d\n", before_final - nrow(df), nrow(df)))
log_flow("核心协变量完整 (最终)", nrow(df))

# --------------------------------------------------
# 汇总
# --------------------------------------------------
nhanes_final <- df
cat(sprintf("\n========================================\n"))
cat(sprintf("最终分析样本: %d 行 × %d 列\n", nrow(nhanes_final), ncol(nhanes_final)))
cat(sprintf("  GDM 史阳性 (RHQ162==1): %d (%.1f%%)\n",
            sum(nhanes_final$gdm_binary), 100*mean(nhanes_final$gdm_binary)))
cat(sprintf("  全因死亡: %d (%.1f%%)\n",
            sum(nhanes_final$mort_allcause, na.rm = TRUE),
            100*mean(nhanes_final$mort_allcause, na.rm = TRUE)))
cat(sprintf("  CVD 死亡 (UCOD ∈ {1,5}): %d (%.1f%%)\n",
            sum(nhanes_final$mort_cvd, na.rm = TRUE),
            100*mean(nhanes_final$mort_cvd, na.rm = TRUE)))
cat(sprintf("  睡眠紊乱 (score≥2): %d (%.1f%%)\n",
            sum(nhanes_final$sleep_disorder),
            100*mean(nhanes_final$sleep_disorder)))
cat(sprintf("  DII 中位数: %.2f；范围: %.2f ~ %.2f\n",
            median(nhanes_final$DII),
            min(nhanes_final$DII), max(nhanes_final$DII)))
cat(sprintf("  自报复合 CVD: %d (%.1f%%)\n",
            sum(nhanes_final$cvd_composite, na.rm = TRUE),
            100*mean(nhanes_final$cvd_composite, na.rm = TRUE)))
cat(sprintf("========================================\n"))

# --------------------------------------------------
# 保存
# --------------------------------------------------
if (!dir.exists("data/processed")) dir.create("data/processed", recursive = TRUE)
save(nhanes_final, dii_table, file = "data/processed/nhanes_final.RData")
cat("已保存 data/processed/nhanes_final.RData\n")

if (!dir.exists("output/tables")) dir.create("output/tables", recursive = TRUE)
flow_df <- do.call(rbind, flow)
write.csv(flow_df, "output/tables/flow_counts.csv", row.names = FALSE)
cat("已保存 output/tables/flow_counts.csv\n")