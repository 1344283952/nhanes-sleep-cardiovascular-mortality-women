# ============================================
# 23_demote_exploratory.R  (primary analysis)
# Demote the ABANDONED framings to clearly-labeled EXPLORATORY supplementary
# outputs (prefix supp_exploratory_). These are NO LONGER headline but are kept
# for the honest "pre-registered hypotheses were null" disclosure.
#
# Abandoned framings (all null / implausible):
#   - GDM-history as effect modifier (DII x GDM interaction; GDM-stratified Cox)
#   - sleep mediation of DII -> mortality (CMAverse + g-formula)
#   - education / SES x exposure interaction (subgroup; education x sleep -> CV
#     fails Bonferroni in the 36-test subgroup family)
#
# This script COPIES (does not delete) the source CSVs into supp_exploratory_*.csv
# with a prepended provenance/demotion note row, and writes an index manifest.
#
# Input/Output: output/tables/*
# ============================================

suppressPackageStartupMessages({ library(dplyr) })

cat("=====================================================\n")
cat(" Demote abandoned framings -> supp_exploratory_*\n")
cat("=====================================================\n\n")

tdir <- "output/tables"

# source file -> (demoted name, framing label, status)
demote_map <- list(
  list(src = "table4_dii_gdm_interaction.csv",
       dst = "supp_exploratory_gdm_interaction.csv",
       framing = "GDM-history effect modification (DII x GDM)",
       status  = "EXPLORATORY - interaction null (all P>0.05 for DII:GDM Yes)"),
  list(src = "table4b_dii_stratified_by_gdm.csv",
       dst = "supp_exploratory_gdm_stratified.csv",
       framing = "GDM-history stratified Cox (DII -> all-cause)",
       status  = "EXPLORATORY - GDM=Yes stratum null, no modification"),
  list(src = "table5_mediation_main.csv",
       dst = "supp_exploratory_sleep_mediation_allcause.csv",
       framing = "Sleep mediation of DII -> all-cause (CMAverse)",
       status  = "EXPLORATORY - PM ~1%, mediation null (P=0.198)"),
  list(src = "table5b_mediation_cvd.csv",
       dst = "supp_exploratory_sleep_mediation_cvd.csv",
       framing = "Sleep mediation of DII -> CV (CMAverse)",
       status  = "EXPLORATORY - mediation null"),
  list(src = "table_supp_S14_gformula.csv",
       dst = "supp_exploratory_sleep_mediation_gformula.csv",
       framing = "Sleep mediation of DII -> all-cause (g-formula)",
       status  = "EXPLORATORY - PM ~1%, mediation null"),
  list(src = "table6b_subgroup_pinteraction.csv",
       dst = "supp_exploratory_subgroup_interaction_pvalues.csv",
       framing = "SES/education x exposure subgroup interactions",
       status  = "EXPLORATORY - education x sleep->CV P_int=0.023 raw FAILS Bonferroni (36-test family, threshold 0.00139)"),
  list(src = "table6_subgroup_strat.csv",
       dst = "supp_exploratory_subgroup_stratified.csv",
       framing = "Subgroup-stratified HRs (incl. education strata)",
       status  = "EXPLORATORY - hypothesis-generating only")
)

manifest <- list
for (m in demote_map) {
  src <- file.path(tdir, m$src)
  dst <- file.path(tdir, m$dst)
  if (!file.exists(src)) {
    cat(sprintf("  [skip] missing source: %s\n", m$src))
    manifest[[length(manifest)+1]] <- data.frame(
      demoted_file = m$dst, source_file = m$src, framing = m$framing,
      status = paste0("SOURCE MISSING - ", m$status), stringsAsFactors = FALSE)
    next
  }
  dat <- read.csv(src, stringsAsFactors = FALSE, check.names = FALSE)
  # write a provenance header comment + the data
  con <- file(dst, "w", encoding = "UTF-8")
  writeLines(sprintf("# EXPLORATORY (demoted, not headline) | framing: %s", m$framing), con)
  writeLines(sprintf("# status: %s", m$status), con)
  writeLines(sprintf("# source: %s", m$src), con)
  close(con)
  suppressWarnings(
    write.table(dat, dst, sep = ",", row.names = FALSE, col.names = TRUE,
                append = TRUE, qmethod = "double", fileEncoding = "UTF-8")
)
  cat(sprintf("  [OK] %-42s -> %s\n", m$src, m$dst))
  manifest[[length(manifest)+1]] <- data.frame(
    demoted_file = m$dst, source_file = m$src, framing = m$framing,
    status = m$status, stringsAsFactors = FALSE)
}

man_df <- bind_rows(manifest)
write.csv(man_df, file.path(tdir, "supp_exploratory_INDEX.csv"), row.names = FALSE)
cat("\n[OK] output/tables/supp_exploratory_INDEX.csv (manifest of demoted analyses)\n")
print(man_df[, c("demoted_file", "status")])
cat("\nNote: original source CSVs are RETAINED (not deleted) for full reproducibility.\n")