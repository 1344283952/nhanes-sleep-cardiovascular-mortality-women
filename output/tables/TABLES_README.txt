================================================================================
TABLES README — manuscript "Table N" / "Figure N" -> deliverable file mapping
Manuscript: "Sleep disturbance as a modifiable marker of cardiovascular
mortality risk in US women with a history of pregnancy: NHANES 2007-2018"
(primary endpoint = cardiovascular mortality)
================================================================================

MAIN DISPLAY TABLES (upload these to the submission system)
--------------------------------------------------------------------------------
Table 1  Baseline characteristics of the analytic sample (N = 10,931), by
         sleep-disturbance status (composite score <2 vs >=2). Survey-weighted.
           Display file (upload) : table1_display.xlsx
           Clean Markdown        : tables_display.md  (## Table 1 section)
           Underlying analysis   : table1_baseline_by_sleep.csv

Table 2  Survey-weighted associations of sleep disturbance and the Dietary
         Inflammatory Index with cardiovascular and all-cause mortality.
         Panel A primary Fine-Gray SHR; Panel B cause-specific Cox adjustment
         ladder; Panel C companion sleep->all-cause; Panels D-E secondary DII.
           Display file (upload) : table2_display.xlsx
           Clean Markdown        : tables_display.md  (## Table 2 section)
           Underlying analysis   : table1_primary_sleep_cvd.csv  (sleep->CV ladder,
                                   Fine-Gray, E-value, EPV, sleep->all-cause)
                                   table2_secondary_dii.csv      (DII per SD,
                                   quartiles, P-trend, DII->CV null)

Table 3  Multiple-comparison adjustment across the four-test primary family
         (sleep & DII x CV & all-cause): raw, Bonferroni, Benjamini-Hochberg.
           Display file (upload) : table3_display.xlsx
           Clean Markdown        : tables_display.md  (## Table 3 section)
           Underlying analysis   : table3_multiplicity.csv

MAIN FIGURE
--------------------------------------------------------------------------------
Figure 1 Sample-selection (CONSORT-style) flow diagram: pooled six-cycle NHANES
         2007-2018 -> final analytic sample N = 10,931.
           Deliverable           : ../figures/fig1_consort.png
           (vector source)       : ../figures/fig1_consort.svg / .dot

NOTE: The manuscript body and Figure Legends reference exactly ONE figure
(Figure 1 = fig1_consort.png) and three tables (Tables 1-3). The primary-
result forest plot ../figures/fig_primary_sleep_cvd_forest.png EXISTS but is
NOT cited in the manuscript; it is an optional visual, not a required
deliverable. No manuscript Table/Figure reference is missing its file.

KEY NUMBERS (cross-check vs manuscript; all verified equal)
--------------------------------------------------------------------------------
  Primary sleep->CV Fine-Gray SHR    1.58 (1.21-2.07), P = 0.0009
  Primary sleep->CV cause-specific   1.66 (1.16-2.35), P = 0.005
  E-value (point / CI bound)         2.70 / 1.60
  EPV / CV deaths                    20.8 / 249
  Adjustment ladder (cause-specific) 1.76 -> 1.66 -> 1.62 -> 2.13
  Sleep->all-cause                   1.34 (1.10-1.63), P = 0.004
  DII->all-cause per SD              1.16 (1.06-1.27), P = 0.002
  DII Q4 vs Q1                       1.52 (1.16-1.98), P-trend = 0.006
  DII->CV (null)                     1.08 (0.93-1.26), P = 0.32
  Bonferroni-adjusted (4-test)       sleep->CV 0.020; sleep->all 0.015;
                                     DII->all 0.008; DII->CV 1.00 (does not survive)

SUPPLEMENTARY / EXPLORATORY ANALYSIS FILES (see supplementary.md, not main tables)
--------------------------------------------------------------------------------
  Pre-specified NULL hypotheses and EXPLORATORY (not-a-claim) analyses:
    supp_exploratory_gdm_interaction.csv          GDM-modifier interaction (NULL)
    supp_exploratory_gdm_stratified.csv           GDM-stratified estimates
    supp_exploratory_sleep_mediation_allcause.csv sleep-mediation of DII (NULL)
    supp_exploratory_sleep_mediation_cvd.csv      sleep-mediation, CV outcome
    supp_exploratory_sleep_mediation_gformula.csv g-formula mediation sensitivity
    supp_exploratory_subgroup_interaction_pvalues.csv  subgroup interactions
    supp_exploratory_subgroup_stratified.csv      stratified estimates (incl. education, EXPLORATORY)
    supp_exploratory_INDEX.csv                    index of exploratory outputs
    table_supp_S13_finegray.csv                   Fine-Gray detail
    table_supp_S14_gformula.csv                   g-formula mediation detail
    flow_counts.csv                               CONSORT counts (Figure 1)

LEGACY FILES (from the abandoned DII-headline / SES-educational framings; NOT for
this submission): table1.xlsx, table1_by_gdm.csv, table1b_by_mortality.csv,
table2_cox*.{csv,xlsx}, table3_logistic_cvd*, table4_*, table5_*, table6_*,
table7_*, table8_*, rcs_pvalues.csv. Do NOT upload these as main tables.
================================================================================
