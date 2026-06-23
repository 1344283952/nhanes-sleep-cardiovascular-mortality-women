# ============================================
# 15_dag.R (- , simplified for clarity)
# DAG: simplified causal model with fewer collapsed nodes
# (fix: previous 17-node DAG had crossing edges + label overlap)
# ============================================

suppressPackageStartupMessages({
  library(ggdag); library(ggplot2); library(dplyr)
})

cat("========================================\n")
cat(" DAG (simplified)\n")
cat("========================================\n\n")

# Collapsed nodes for clarity:
#   SES        = Education + PIR + Race
#   Lifestyle  = BMI + Smoking + Alcohol + Parity
#   CMDz       = Diabetes + Hypertension + Dyslipidemia (downstream cardiometabolic)
dag <- dagify(
  Death     ~ DII + Sleep + Age + SES + Lifestyle + CMDz + GDMhx,
  Sleep     ~ DII + Age + Lifestyle + SES,
  CMDz      ~ DII + Age + Lifestyle + GDMhx,
  DII       ~ Age + SES,
  GDMhx     ~ Age + SES,
  outcome   = "Death",
  exposure  = "DII",
  coords = list(
    x = c(SES = 0, Age = 0, GDMhx = 0,
          DII = 2, Lifestyle = 2,
          Sleep = 4, CMDz = 4,
          Death = 6),
    y = c(SES = 3, Age = 1.5, GDMhx = 0,
          DII = 3, Lifestyle = 1,
          Sleep = 3, CMDz = 1,
          Death = 2)
)
)

node_status <- function(name) {
  case_when(
    name == "DII"   ~ "exposure",
    name == "Death" ~ "outcome",
    name == "Sleep" ~ "mediator",
    name == "CMDz"  ~ "mediator",
    name == "GDMhx" ~ "modifier",
    TRUE            ~ "confounder"
)
}

dag_tidy <- dag %>% tidy_dagitty
dag_tidy$data$status <- node_status(dag_tidy$data$name)

p <- ggplot(dag_tidy$data, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_dag_edges(edge_color = "gray45", edge_width = 0.6,
                 arrow_directed = grid::arrow(length = grid::unit(8, "pt"),
                                              type = "closed")) +
  geom_dag_point(aes(fill = status), shape = 21, size = 22,
                 stroke = 0.5, color = "grey30") +
  geom_dag_text(color = "black", size = 3.6, fontface = "bold") +
  scale_fill_manual(values = c(exposure   = "#E74C3C",
                               outcome    = "#E74C3C",
                               mediator   = "#F4D03F",
                               modifier   = "#85C1E9",
                               confounder = "#FFFFFF"),
                    name = NULL,
                    breaks = c("exposure", "outcome", "mediator",
                               "modifier", "confounder"),
                    labels = c("Exposure (DII)", "Outcome (Death)",
                               "Mediator (Sleep / CMDz)",
                               "Modifier (GDMhx)",
                               "Confounder (Age / SES / Lifestyle)")) +
  theme_dag_blank +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 9),
        legend.box = "horizontal") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE,
                             override.aes = list(size = 6))) +
  expand_limits(x = c(-0.5, 6.5), y = c(-0.7, 4.0)) +
  labs(title = "DAG: pro-inflammatory diet, sleep, and mortality",
       subtitle = "Nodes collapsed for clarity (see Figure 2 caption)") +
  theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey30"))

if (!dir.exists("output/figures")) dir.create("output/figures", recursive = TRUE)
ggsave("output/figures/fig2_dag.png", plot = p,
       width = 10, height = 7, dpi = 300, bg = "white")

cat("已保存 output/figures/fig2_dag.png\n")
cat("========================================\n")