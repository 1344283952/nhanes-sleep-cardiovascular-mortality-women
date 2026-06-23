# ============================================
# 02_merge_data.R  (版)
# 合并 NHANES 2007-2018 (6 周期) + NCHS Linked Mortality
# 输入:  data/raw/*.xpt + data/raw/mortality/*.dat
# 输出:  data/processed/nhanes_raw_merged.RData (含 nhanes_all + rx_all)
#
# 与 templates 原版 (a prior NHANES liver-disease analysis) 差异：
#   - 周期 1999-2018 → 2007-2018 (6 cycles)
#   - 模块清单：加 RHQ/SLQ/DR1TOT/DR2TOT/VID/PAQ；删 CBC（DII/Sleep/GDM 不用）
#   - 删 early_aliases (旧周期 LAB 系列)； 不需要 1999-2005
#   - VID_2_G 修订版单独 join（仅 2011-2012）
# ============================================

library(haven)
library(dplyr)
library(purrr)

cat("========================================\n")
cat("合并 NHANES 2007-2018  + 死亡链接\n")
cat("========================================\n\n")

raw_dir  <- "data/raw"
mort_dir <- "data/raw/mortality"

# --------------------------------------------------
# Step 1: 读所有 .xpt（glob，自动跳过坏文件）
# --------------------------------------------------
xpt_files <- list.files(raw_dir, pattern = "\\.xpt$", full.names = TRUE,
                        ignore.case = TRUE)
cat(sprintf("找到 %d 个 .xpt 文件\n\n", length(xpt_files)))

read_safe <- function(p) {
  tryCatch(read_xpt(p),
           error = function(e) {
             cat(sprintf("  [读失败] %s : %s\n", basename(p), e$message))
             NULL
           })
}

data_list <- map(xpt_files, read_safe)
names(data_list) <- gsub("\\.xpt$", "", basename(xpt_files), ignore.case = TRUE)
data_list <- data_list[!sapply(data_list, is.null)]
cat(sprintf("成功读入 %d 个数据帧\n\n", length(data_list)))

# --------------------------------------------------
# Step 2: 周期 ↔ 后缀 + 模块清单
# --------------------------------------------------
cycles <- data.frame(
  year   = c(2007, 2009, 2011, 2013, 2015, 2017),
  suffix = c("_E", "_F", "_G", "_H", "_I", "_J"),
  stringsAsFactors = FALSE
)

# 跨周期同名模块（每人 1 行）
modules_uniform <- c(
  "DEMO",
  "RHQ", "SLQ",                          # 生殖史 + 睡眠（ 关键）
  "DR1TOT", "DR2TOT",                    # 24h 膳食回顾（DII 主源）
  "VID",                                 # 维生素 D（DII 27 项之一）
  "BMX", "BPX",                          # 体测 + 血压
  "BPQ", "MCQ", "DIQ",                   # 慢病自报 + 自报 CVD
  "GHB", "GLU", "INS",                   # 血糖 / 胰岛素（HOMA-IR）
  "BIOPRO", "HDL", "TCHOL", "TRIGLY",    # 生化 + 血脂
  "SMQ", "ALQ", "PAQ"                    # 生活方式
)

# 长表（每人多行，单独处理）
modules_long <- c("RXQ_RX")

# --------------------------------------------------
# Step 3: 工具函数
# --------------------------------------------------
safe_left_join <- function(x, y, key = "SEQN") {
  if (is.null(y) || nrow(y) == 0) return(x)
  if (!key %in% names(y)) return(x)
  dup <- intersect(setdiff(names(x), key), names(y))
  if (length(dup) > 0) y <- y[, setdiff(names(y), dup), drop = FALSE]
  left_join(x, y, by = key)
}

merge_cycle <- function(year, suffix) {
  cat(sprintf("--- 周期 %d-%d (后缀 '%s') ---\n", year, year+1, suffix))
  base_name <- paste0("DEMO", suffix)
  if (!base_name %in% names(data_list)) {
    cat(sprintf("  [警告] 基表 %s 不存在，跳过\n\n", base_name))
    return(NULL)
  }
  result <- data_list[[base_name]]
  result$cycle_year <- year

  for (mod in modules_uniform[-1]) {  # skip DEMO
    key <- paste0(mod, suffix)
    if (key %in% names(data_list)) {
      result <- safe_left_join(result, data_list[[key]])
    }
  }

  # VID_2_G 修订版（仅 2011-2012），优先于 VID_G
  if (suffix == "_G" && "VID_2_G" %in% names(data_list)) {
    result <- safe_left_join(result, data_list[["VID_2_G"]])
  }

  cat(sprintf("  -> %d 行 × %d 列\n\n", nrow(result), ncol(result)))
  result
}

# --------------------------------------------------
# Step 4: 拼 6 个周期
# --------------------------------------------------
merged_list <- pmap(list(cycles$year, cycles$suffix), merge_cycle)
merged_list <- merged_list[!sapply(merged_list, is.null)]
nhanes_all <- bind_rows(merged_list)

cat(sprintf("========================================\n"))
cat(sprintf("6 周期合并完成: %d 行 × %d 列\n", nrow(nhanes_all), ncol(nhanes_all)))
cat(sprintf("========================================\n\n"))

# --------------------------------------------------
# Step 4.5: 长表 RXQ_RX 单独合并（用于 03 降压/降脂药识别）
# --------------------------------------------------
cat("--- Step 4.5: 长表 RXQ_RX 合并 ---\n")
rx_list <- list
for (suffix in cycles$suffix) {
  key <- paste0("RXQ_RX", suffix)
  if (key %in% names(data_list)) {
    rx_list[[key]] <- data_list[[key]]
  }
}
rx_all <- if (length(rx_list) > 0) bind_rows(rx_list) else data.frame(SEQN=integer(0))
cat(sprintf("rx_all: %d 行 × %d 列\n\n", nrow(rx_all), ncol(rx_all)))

# --------------------------------------------------
# Step 5: 读 NCHS Linked Mortality (固定宽度 .dat)
# 字段位置参考: https://www.cdc.gov/nchs/data-linkage/mortality-public.htm
# --------------------------------------------------
cat("--- Step 5: 接入死亡链接 ---\n")
mort_files <- list.files(mort_dir, pattern = "\\.dat$", full.names = TRUE)
cat(sprintf("找到 %d 个死亡链接文件\n", length(mort_files)))

read_mort <- function(path) {
  # NCHS Public-Use Linked Mortality File fixed-width layout (2019 release).
  # Bug fix DODQTR is 1 column (Q1/Q2/Q3/Q4 = '1'-'4'), not 2.
  # The previous DODQTR=2 shifted every subsequent field by one column,
  # truncating PERMTH_EXM to mod-100 (catastrophic for survival time).
  widths <- c(
    SEQN          = 6,
    PADDING1      = 8,
    ELIGSTAT      = 1,
    MORTSTAT      = 1,
    UCOD_LEADING  = 3,
    DIABETES      = 1,
    HYPERTEN      = 1,
    DODQTR        = 1,
    DODYEAR       = 4,
    WGT_NEW       = 8,
    SA_WGT_NEW    = 8,
    PERMTH_INT    = 3,
    PERMTH_EXM    = 3
)
  df <- tryCatch(
    read.fwf(path, widths = widths, header = FALSE,
             na.strings = c("", "."), stringsAsFactors = FALSE,
             col.names = names(widths)),
    error = function(e) {
      cat(sprintf("  [读失败] %s : %s\n", basename(path), e$message)); NULL
    }
)
  if (is.null(df)) return(NULL)
  df$PADDING1 <- NULL
  df$SEQN          <- as.integer(df$SEQN)
  df$ELIGSTAT      <- suppressWarnings(as.integer(df$ELIGSTAT))
  df$MORTSTAT      <- suppressWarnings(as.integer(df$MORTSTAT))
  df$UCOD_LEADING  <- suppressWarnings(as.integer(df$UCOD_LEADING))
  df$DIABETES      <- suppressWarnings(as.integer(df$DIABETES))
  df$HYPERTEN      <- suppressWarnings(as.integer(df$HYPERTEN))
  df$PERMTH_INT    <- suppressWarnings(as.integer(df$PERMTH_INT))
  df$PERMTH_EXM    <- suppressWarnings(as.integer(df$PERMTH_EXM))
  df
}

mort_all <- bind_rows(map(mort_files, read_mort))
cat(sprintf("死亡链接合并完成: %d 行\n", nrow(mort_all)))
cat(sprintf("  ELIGSTAT == 1 (合格): %d\n", sum(mort_all$ELIGSTAT == 1, na.rm=TRUE)))
cat(sprintf("  MORTSTAT == 1 (死亡): %d\n", sum(mort_all$MORTSTAT == 1, na.rm=TRUE)))

# --------------------------------------------------
# Step 6: 合并死亡数据进主表
# --------------------------------------------------
nhanes_all <- safe_left_join(nhanes_all, mort_all)
cat(sprintf("\n合并死亡后: %d 行 × %d 列\n", nrow(nhanes_all), ncol(nhanes_all)))

# --------------------------------------------------
# 保存
# --------------------------------------------------
if (!dir.exists("data/processed")) dir.create("data/processed", recursive = TRUE)
save(nhanes_all, rx_all, file = "data/processed/nhanes_raw_merged.RData")
cat("\n已保存 data/processed/nhanes_raw_merged.RData (含 rx_all)\n")
cat("========================================\n")