# ============================================
# 01_download_data.R （ 版）
# 下载 NHANES 2007-2018（6 周期）+ NCHS Linked Mortality File
# 论文：DII × Sleep × GDM 史 → 全因/CVD 死亡（套 Ying et al. 2024）
#
# 与 v1  差异：
#   - 加 _J 周期（2017-2018）
#   - 加 NCHS Linked Mortality（v1 不做死亡随访）
#
# 并行版（PSOCK cluster，Windows 兼容）:
#   - 默认 6 worker（CDC 服务器允许 6-8 并发；勿超过 10）
#   - 内容嗅探（首 6 字节必须是 "HEADER"，否则视为 HTML 假页删除）
#   - 进度日志同步写到 data/raw/_download.log
# ============================================

cat("========================================\n")
cat("下载 NHANES 2007-2018  + 死亡链接（并行）\n")
cat("========================================\n\n")

raw_dir  <- "data/raw"
mort_dir <- "data/raw/mortality"
log_path <- "data/raw/_download.log"
if (!dir.exists(raw_dir))  dir.create(raw_dir,  recursive = TRUE)
if (!dir.exists(mort_dir)) dir.create(mort_dir, recursive = TRUE)

writeLines(c("# NHANES download log", paste("# started:", Sys.time)),
           log_path)
log_line <- function(msg) {
  cat(msg, "\n", sep = "")
  cat(msg, "\n", sep = "", file = log_path, append = TRUE)
}

# --------------------------------------------------
# 周期映射（2007-2018，6 周期）
# 之所以从 2007 开始：DR1TVD（膳食维生素 D）从 2007-2008 才纳入 NHANES 24h 回顾
# 2005-2006 没有该字段，DII 27 项会缺一项 → 整 DII 算不出
# --------------------------------------------------
cycles <- data.frame(
  year   = c(2007, 2009, 2011, 2013, 2015, 2017),
  suffix = c("_E", "_F", "_G", "_H", "_I", "_J"),
  stringsAsFactors = FALSE
)

# --------------------------------------------------
# 模块清单（22 个跨周期同名模块）
# --------------------------------------------------
modules_uniform <- c(
  # 人口学与样本设计
  "DEMO",
  # 生殖史（暴露 + 修饰，关键  新增）
  "RHQ",
  # 24h 膳食回顾（DII 主源）
  "DR1TOT", "DR2TOT",
  # 睡眠
  "SLQ",
  # 体测 / 血压
  "BMX", "BPX",
  # 慢病自报 / 高血压问卷 / 肾脏
  "BPQ", "MCQ", "DIQ", "KIQ",
  # 实验室（HOMA-IR 中介 + 协变量）
  "GHB", "GLU", "INS", "BIOPRO",
  "HDL", "TCHOL", "TRIGLY",
  # 生活方式
  "SMQ", "ALQ", "PAQ",
  # 用药（保守降压/致脂药识别）
  "RXQ_RX"
)

# 维生素 D 文件命名特殊（DII 27 项之一）
# 2007-2010: VID_E/F；2011-2012: VID_G + VID_2_G（修订版）；2013-2018: VID_H/I/J
modules_vit_d <- list(
  list("VID_E",   2007), list("VID_F",   2009),
  list("VID_G",   2011), list("VID_2_G", 2011),
  list("VID_H",   2013), list("VID_I",   2015),
  list("VID_J",   2017)
)

# --------------------------------------------------
# 构造下载任务表
# --------------------------------------------------
tasks <- list
for (mod in modules_uniform) {
  for (i in seq_len(nrow(cycles))) {
    tasks[[length(tasks)+1]] <- list(
      filename = paste0(mod, cycles$suffix[i]),
      year = cycles$year[i],
      kind = "xpt"
)
  }
}
for (item in modules_vit_d) {
  tasks[[length(tasks)+1]] <- list(
    filename = item[[1]], year = as.integer(item[[2]]), kind = "xpt"
)
}

# 死亡链接任务（NCHS 2019 public file，6 周期各一个 .dat）
mort_base <- "https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality"
for (i in seq_len(nrow(cycles))) {
  yr <- cycles$year[i]
  fn <- paste0("NHANES_", yr, "_", yr+1, "_MORT_2019_PUBLIC.dat")
  tasks[[length(tasks)+1]] <- list(
    filename = fn, year = yr, kind = "mort",
    url = paste0(mort_base, "/", fn)
)
}

log_line(sprintf("总任务数: %d (xpt=%d, mort=%d)",
                 length(tasks),
                 sum(sapply(tasks, function(t) t$kind == "xpt")),
                 sum(sapply(tasks, function(t) t$kind == "mort"))))

# --------------------------------------------------
# Worker 函数（每个 worker 独立调用）
# --------------------------------------------------
download_task <- function(task, raw_dir, mort_dir) {
  is_real_xpt <- function(path) {
    if (!file.exists(path)) return(FALSE)
    if (file.size(path) < 1024) return(FALSE)
    con <- file(path, "rb"); on.exit(close(con))
    bytes <- readChar(con, 6, useBytes = TRUE)
    identical(bytes, "HEADER")
  }

  if (task$kind == "xpt") {
    destfile <- file.path(raw_dir, paste0(task$filename, ".xpt"))
    if (is_real_xpt(destfile)) {
      return(list(file = task$filename, status = "已存在"))
    }
    urls <- c(
      paste0("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/",
             task$year, "/DataFiles/", task$filename, ".xpt"),
      paste0("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/",
             task$year, "/DataFiles/", task$filename, ".XPT")
)
    for (url in urls) {
      ok <- tryCatch({
        download.file(url, destfile = destfile, mode = "wb",
                      quiet = TRUE, method = "libcurl")
        is_real_xpt(destfile)
      }, error = function(e) FALSE, warning = function(w) FALSE)
      if (isTRUE(ok)) return(list(file = task$filename, status = "ok"))
      if (file.exists(destfile)) file.remove(destfile)
    }
    return(list(file = task$filename, status = "FAIL"))
  } else {
    destfile <- file.path(mort_dir, task$filename)
    if (file.exists(destfile) && file.size(destfile) > 0) {
      return(list(file = task$filename, status = "已存在"))
    }
    ok <- tryCatch({
      download.file(task$url, destfile = destfile, mode = "wb",
                    quiet = TRUE, method = "libcurl")
      file.exists(destfile) && file.size(destfile) > 0
    }, error = function(e) FALSE, warning = function(w) FALSE)
    if (isTRUE(ok)) return(list(file = task$filename, status = "ok"))
    return(list(file = task$filename, status = "FAIL"))
  }
}

# --------------------------------------------------
# 并行执行（PSOCK cluster, Windows 兼容）
# --------------------------------------------------
library(parallel)

n_workers <- min(6, length(tasks))   # 6 并发
log_line(sprintf("启动 %d 个 worker 并行下载...", n_workers))

cl <- makeCluster(n_workers)
on.exit(stopCluster(cl), add = TRUE)

t0 <- Sys.time
results <- parLapplyLB(cl, tasks, download_task,
                       raw_dir = raw_dir, mort_dir = mort_dir)
elapsed <- round(as.numeric(difftime(Sys.time, t0, units = "secs")), 1)

stopCluster(cl)

# --------------------------------------------------
# 汇总
# --------------------------------------------------
status_tbl <- table(sapply(results, function(r) r$status))
log_line(sprintf("\n并行下载结束（%.1f 秒）：", elapsed))
for (s in names(status_tbl)) {
  log_line(sprintf("  %-10s : %d", s, status_tbl[[s]]))
}

fails <- sapply(results[sapply(results, function(r) r$status == "FAIL")],
                function(r) r$file)
if (length(fails) > 0) {
  log_line("\n失败文件（部分模块在某些周期不存在是正常的）:")
  for (f in fails) log_line(sprintf("  - %s", f))
}

n_xpt  <- length(list.files(raw_dir,  pattern = "\\.xpt$"))
n_mort <- length(list.files(mort_dir, pattern = "\\.dat$"))
log_line(sprintf("\ndata/raw/         %d 个 .xpt", n_xpt))
log_line(sprintf("data/raw/mortality/ %d 个 .dat", n_mort))
log_line("========================================")