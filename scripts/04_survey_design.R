# ============================================
# 04_survey_design.R (版)
# 复杂抽样设计（DII 主暴露 → 用单日膳食权重 WTDRD1/6）
# 输入:  data/processed/nhanes_final.RData
# 输出:  data/processed/nhanes_design.RData
# ============================================

library(survey)
library(dplyr)

cat("========================================\n")
cat("构建复杂抽样设计 \n")
cat("========================================\n\n")

load("data/processed/nhanes_final.RData")

# 重要：nest=TRUE 才能让 svyglm 在 1 个 PSU 的层正确处理
options(survey.lonely.psu = "adjust")

# --------------------------------------------------
# 主 design：DII 是单日 DR1 主暴露 → 用 WTDR_6YR 单日膳食权重
# --------------------------------------------------
design <- svydesign(
  ids     = ~SDMVPSU,
  strata  = ~SDMVSTRA,
  weights = ~WTDR_6YR,
  nest    = TRUE,
  data    = nhanes_final
)

cat(sprintf("svydesign 创建完成\n"))
cat(sprintf("  N: %d / PSU: %d / Strata: %d\n",
            nrow(nhanes_final),
            length(unique(nhanes_final$SDMVPSU)),
            length(unique(nhanes_final$SDMVSTRA))))
cat(sprintf("  加权总数: %.0f\n", sum(nhanes_final$WTDR_6YR)))

# --------------------------------------------------
# 子集（方便后续脚本直接 reload）
# --------------------------------------------------
# DII 四分位
design_dii_q1 <- subset(design, DII_Q == "Q1")
design_dii_q4 <- subset(design, DII_Q == "Q4")

# GDM 史（主修饰因子）
design_gdm_yes <- subset(design, gdm_history == "Yes")
design_gdm_no  <- subset(design, gdm_history == "No")

# 睡眠紊乱
design_sleep_dis <- subset(design, sleep_disorder == 1)
design_sleep_ok  <- subset(design, sleep_disorder == 0)

# Mortality 子集（已在清洗中限定 ELIGSTAT==1）
design_mort <- design

save(design, design_mort,
     design_dii_q1, design_dii_q4,
     design_gdm_yes, design_gdm_no,
     design_sleep_dis, design_sleep_ok,
     nhanes_final,
     file = "data/processed/nhanes_design.RData")

cat("\n已保存 data/processed/nhanes_design.RData\n")
cat("含: design / design_mort / design_dii_q[1,4] / design_gdm_[yes,no] / design_sleep_[dis,ok] / nhanes_final\n")
cat("========================================\n")