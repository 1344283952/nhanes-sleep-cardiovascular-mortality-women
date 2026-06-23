# ============================================
# 00_install_packages.R (版)
# 首次运行：装  全部依赖
# 需要 Rtools 4.6 才能编译 source-only 包（如 CMAverse）
# ============================================

cat("========================================\n")
cat(" 装包：基础 + Cox/中介/插补/IPTW/RCS/预测\n")
cat("========================================\n\n")

repos <- "https://cloud.r-project.org"

# --------------------------------------------------
# CRAN 包（按主题分组）
# --------------------------------------------------
packages <- c(
  # 数据基础（templates 通用）
  "tidyverse",    # dplyr, tidyr, ggplot2, purrr, stringr...
  "survey",       # 复杂抽样
  "haven",        # 读 .xpt
  "broom",        # tidy 模型输出
  "tableone",     # Table 1
  "openxlsx",     # Excel 输出
  "DiagrammeR",   # CONSORT 流程图
  "corrplot",     # 相关矩阵
  "remotes",      # github 装包

  # 生存分析 + Cox
  "survival",
  "survminer",    # ggsurvplot
  "rms",          # RCS, nomogram, val.surv

  # 中介分析
  "regmedint",    # rare-outcome closed-form (备用，CMAverse 主)
  "mediation",    # 老包，做加权敏感性

  # 多重插补（S5 敏感性）
  "mice",
  "VIM",          # missingness 可视化

  # IPTW / propensity score（S6 敏感性）
  "WeightIt",
  "cobalt",       # 平衡诊断
  "MatchIt",      # 备用

  # DAG / 因果识别
  "ggdag",
  "dagitty",

  # 预测层（Table 8）
  "pROC",         # ROC / AUC
  "nricens",      # NRI（含 binary + survival）
  "PredictABEL",  # IDI（备用，nricens 也有）

  # 多重比较 / E-value
  "EValue",       # VanderWeele 2017
  "Hmisc",        # rcs + describe

  # 绘图扩展
  "patchwork",
  "ggpubr",
  "RColorBrewer",

  # 报告辅助
  "rmarkdown",
  "knitr",
  "zip"
)

# --------------------------------------------------
# 装 CRAN 包（type="binary" 强制走预编译，避开 Rtools 慢编译）
# --------------------------------------------------
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("安装: %s ...\n", pkg))
    tryCatch(
      install.packages(pkg, repos = repos, type = "binary", quiet = TRUE),
      error = function(e) {
        cat(sprintf("  binary 失败，回退 source: %s\n", e$message))
        install.packages(pkg, repos = repos, quiet = TRUE)
      }
)
  } else {
    cat(sprintf("已装: %s\n", pkg))
  }
}

# --------------------------------------------------
# CMAverse 单独装（非 CRAN，GitHub source build）
# 需要 Rtools 4.6
# --------------------------------------------------
if (!requireNamespace("CMAverse", quietly = TRUE)) {
  cat("\n装 CMAverse from GitHub (BS1125/CMAverse) ...\n")
  cat("（需要 Rtools 4.6，首次约 5-10 分钟）\n")
  tryCatch(
    remotes::install_github("BS1125/CMAverse", upgrade = "never", quiet = FALSE),
    error = function(e) {
      cat(sprintf("CMAverse 装失败: %s\n", e$message))
      cat("→ 请确认 Rtools 4.6 已装：https://cran.r-project.org/bin/windows/Rtools/rtools46/\n")
    }
)
} else {
  cat("已装: CMAverse\n")
}

# --------------------------------------------------
# survIDINRI 单独装（biostat3 仓库，备用 NRI）
# --------------------------------------------------
if (!requireNamespace("survIDINRI", quietly = TRUE)) {
  cat("\n装 survIDINRI ...\n")
  tryCatch(
    install.packages("survIDINRI", repos = repos, type = "binary", quiet = TRUE),
    error = function(e) cat(sprintf("survIDINRI 装失败: %s\n", e$message))
)
}

cat("\n========================================\n")
cat("装包完成（如有失败请单独重装）\n")
cat("========================================\n")