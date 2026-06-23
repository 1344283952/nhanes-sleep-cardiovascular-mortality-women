# ============================================
# 14_consort.R
# CONSORT-like flowchart for sample selection
# 用 DiagrammeR (SVG) + ggplot2 兜底 (PNG)
# ============================================

suppressPackageStartupMessages({
  library(DiagrammeR); library(dplyr)
})

cat("========================================\n")
cat(" CONSORT 流程图\n")
cat("========================================\n\n")

flow <- read.csv("output/tables/flow_counts.csv")
print(flow)

# 给每行加 excluded n
flow$excluded <- c(0, diff(-flow$n))

# 用 DiagrammeR grViz
graph_str <- '
digraph CONSORT {
  graph [rankdir = TB, ranksep = 0.5, nodesep = 0.5, splines = ortho]
  node [shape = box, style = filled, fillcolor = "#E8F4FA", fontname = "Arial", fontsize = 10]

  total [label = "NHANES 2007-2018 (6 cycles)\\nN = 59,842"]
  female [label = "Women\\nN = 30,213"]
  adult [label = "Aged 20 years or older\\nN = 17,907"]
  pregnant [label = "One or more prior pregnancies\\nN = 12,952"]
  gdm [label = "Gestational-diabetes history reported\\nN = 12,923"]
  diet [label = "Plausible energy intake and complete diet data\\nN = 12,180"]
  sleep [label = "Three or more evaluable sleep dimensions\\nN = 12,165"]
  cov [label = "Eligible for mortality linkage\\nN = 12,136"]
  final [label = "Final analytic sample\\n(core covariates complete)\\nN = 10,931", fillcolor = "#FFE5A0", penwidth = 2]

  e1 [label = "Excluded: men\\n29,629", shape = plaintext]
  e2 [label = "Excluded: younger than 20 years\\n12,306", shape = plaintext]
  e3 [label = "Excluded: never pregnant\\n4,955", shape = plaintext]
  e4 [label = "Excluded: gestational-diabetes history not reported\\n29", shape = plaintext]
  e5 [label = "Excluded: implausible energy intake or incomplete diet\\n743", shape = plaintext]
  e6 [label = "Excluded: fewer than 3 sleep dimensions\\n15", shape = plaintext]
  e7 [label = "Excluded: ineligible for mortality linkage\\n29", shape = plaintext]
  e8 [label = "Excluded: missing covariates\\n1,205", shape = plaintext]

  total -> female; total -> e1
  female -> adult; female -> e2
  adult -> pregnant; adult -> e3
  pregnant -> gdm; pregnant -> e4
  gdm -> diet; gdm -> e5
  diet -> sleep; diet -> e6
  sleep -> cov; sleep -> e7
  cov -> final; cov -> e8

  outcomes [label = "Outcomes:\\n- All-cause death: 894 (8.2%)\\n- CVD death: 251 (2.3%)\\n- Composite CVD (self-rep): 1,098 (10.0%)", fillcolor = "#D4EFDF"]
  final -> outcomes
}
'

graph <- grViz(graph_str)

# Save SVG via export_svg if available
ok_svg <- tryCatch({
  library(DiagrammeRsvg)
  svg <- export_svg(graph)
  writeLines(svg, "output/figures/fig1_consort.svg")
  TRUE
}, error = function(e) FALSE)

# Save PNG via rsvg if available
ok_png <- tryCatch({
  library(rsvg)
  if (ok_svg) {
    rsvg::rsvg_png("output/figures/fig1_consort.svg",
                   "output/figures/fig1_consort.png",
                   width = 2400)  # ~300 DPI at 8 inches
    TRUE
  } else FALSE
}, error = function(e) FALSE)

# 兜底：保存 grViz 文本到 .dot 文件方便外部转换
writeLines(graph_str, "output/figures/fig1_consort.dot")

cat(sprintf("\nSVG: %s\nPNG: %s\nDOT: output/figures/fig1_consort.dot (兜底)\n",
            ifelse(ok_svg, "output/figures/fig1_consort.svg ✅", "❌ 缺 DiagrammeRsvg"),
            ifelse(ok_png, "output/figures/fig1_consort.png ✅", "❌ 缺 rsvg (用 .dot 在 graphviz 命令行转)")))
cat("========================================\n")