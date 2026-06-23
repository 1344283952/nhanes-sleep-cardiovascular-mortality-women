# Sleep disturbance and cardiovascular mortality in US women with a history of pregnancy (NHANES 2007–2018)

Reproducible analytic pipeline for the manuscript:

> **Sleep disturbance and cardiovascular mortality risk in US women with a history of pregnancy: a competing-risks analysis of NHANES 2007–2018**
> Submitted to *BMC Medicine*.

## What this repository contains

This is the full reproducible R pipeline for our prospective analysis of sleep disturbance, pro-inflammatory diet (Dietary Inflammatory Index, DII), and long-term mortality in **10,931 US women with a history of pregnancy** from six NHANES cycles (2007–2018), linked to the NCHS Public-Use Linked Mortality File.

- **Primary exposure → outcome**: sleep disturbance (6-item composite, score ≥2) → cardiovascular mortality
- **Secondary**: DII (Shivappa 2014, 27 NHANES-available components) → all-cause mortality
- **Primary estimator**: survey-weighted cause-specific Cox; unweighted Fine–Gray competing-risks subdistribution hazard as a concordant sensitivity analysis
- **Methods**: survey-weighted models (R `survey`), Fine–Gray competing risks, restricted cubic splines, E-values, multiple-comparison correction across the primary family, an adjustment-ladder sensitivity analysis

## Key findings

**Primary**:
- Sleep disturbance → cardiovascular mortality: design-primary survey-weighted cause-specific HR 1.65 (95% CI 1.16–2.35, *P* = 0.005); concordant unweighted Fine–Gray SHR 1.58 (1.21–2.07, *P* = 0.0009); E-value 2.70; robust across the adjustment ladder; survives multiple-comparison correction (249 CV deaths, well-powered)
- Sleep disturbance → all-cause mortality: HR 1.34 (1.10–1.63, *P* = 0.004)

**Secondary**:
- DII → all-cause mortality: HR 1.16 per SD (1.06–1.27, *P* = 0.002); Q4 vs Q1 1.52 (1.16–1.98); *P*-trend = 0.006. DII → cardiovascular mortality was null (HR 1.08, *P* = 0.32)

**Pre-specified hypotheses that did not hold (reported transparently)**:
- Prior gestational diabetes did not modify the associations (*P*~interaction~ ≥ 0.06; underpowered subset)
- Sleep did not measurably mediate the diet–mortality association (proportion mediated ≈ 1%, *P* = 0.20)
- An exploratory educational effect-modification analysis did not survive multiple-comparison correction and is reported as hypothesis-generating only

## Repository layout

```
scripts/                R analytic pipeline (run 00 → 16)
  00_install_packages.R    one-time install
  01_download_data.R       download NHANES + mortality
  02_merge_data.R          merge across 6 cycles
  03_clean_data.R          DII + sleep + GDM coding
  04_survey_design.R       svydesign setup
  05_table1.R              baseline descriptive table
  06_cox_main.R            Cox models (exploratory)
  07_logistic_cvd.R        logistic regression, self-reported CVD (exploratory)
  08_subgroup.R            stratified + interaction (exploratory)
  09_rcs.R                 restricted cubic splines (exploratory set)
  11_mediation.R           CMAverse causal mediation
  12_sensitivity.R         sensitivity analyses
  13_predictive.R          C-index + E-value
  14_consort.R             CONSORT flow chart (Figure 1)
  15_dag.R                 directed acyclic graph (legacy)
  16_forest.R              forest plots (exploratory)
  17_supp_sensitivity.R    Fine-Gray + g-formula supplementary sensitivity
  20_primary.R         PRIMARY: cause-specific Cox + Fine-Gray; DII secondary; multiplicity
  21_primary_figure.R      Figure 3 (primary forest: ladder + Fine-Gray)
  22_baseline_table1.R     Table 1 (baseline by sleep-disturbance status)
  23_demote_exploratory.R  demote exploratory analyses to supplementary
  24_display_tables.R  main display Tables 1-3
  25_primary_rcs_dag.R     Figure 2 (cumulative incidence), Figure 4 (RCS), Figure 5 (DAG), supplementary RCS
  run_all.R                end-to-end orchestrator (sources the scripts above)
data/processed/         intermediate .RData (re-created by scripts 01-04; not stored in this repo)
output/
  tables/               primary + supplementary tables (CSV + XLSX)
  figures/              CONSORT, cumulative-incidence curve, forest, RCS dose-response, DAG
task.md                 project specification
references.bib          36 BibTeX entries (cited in manuscript)
```

## How to reproduce

Software: **R 4.6** + **Rtools 4.6** (for compiling packages); **pandoc** (for documentation rendering).

```r
# 1. Install packages (one-time, ~30-60 min including CMAverse from GitHub)
Rscript scripts/00_install_packages.R

# 2. Download raw NHANES + mortality data (~15 min)
Rscript scripts/01_download_data.R

# 3. Run the full pipeline (~20-30 min)
Rscript scripts/run_all.R
```

The intermediate `.RData` files are not stored in this repository (they are `.gitignore`-d for size and re-identification hygiene). Scripts 01–04 re-create them from the public NHANES downloads; `scripts/run_all.R` then reproduces every table and figure in the submission.

## Reproducibility notes

- The full pipeline (data download → final figures) is deterministic; no randomized splits other than CMAverse bootstrap iterations.
- CMAverse uses B = 1000 bootstrap iterations with a fixed seed of `set.seed(20260513)` in `scripts/11_mediation.R` (reproducible to the 4th decimal).
- All survey-weighted analyses use the six-cycle dietary weight (WTDRD1 / 6) with `survey.lonely.psu = "adjust"`.
- The mediation analysis is unweighted (CMAverse limitation; documented in manuscript).

## Data

- **NHANES 2007–2018**: public domain. CDC/NCHS. https://wwwn.cdc.gov/nchs/nhanes/
- **NCHS Linked Mortality File**: public-use 2019 release. https://www.cdc.gov/nchs/data-linkage/mortality-public.htm

This repository does **not** include the raw `.XPT` NHANES files. `scripts/01_download_data.R` will download them into `data/raw/` (which is `.gitignore`-d).

## Reproducibility

This repository contains the full analytic pipeline as submitted to *BMC Medicine*. All scripts run end-to-end on a clean R 4.6 installation (after `00_install_packages.R`). The competing-risks and exploratory analyses use fixed seeds for reproducibility.

## Manuscript

See `manuscript.docx` (in the parent submission package).

## License

Code: MIT License (see `LICENSE`).
Data: NHANES is in the public domain (US government data).

## Citation

If you use this code, please cite:

> Li J, Sun X, Zhang J, Zhai L, Yu L. Sleep disturbance and cardiovascular mortality risk in US women with a history of pregnancy: a competing-risks analysis of NHANES 2007–2018. *BMC Medicine*. (Submitted; volume / pages / DOI to be assigned upon acceptance.)

## Contact

**Corresponding author**: Ling Yu (yulingyxb@jlu.edu.cn, ORCID 0000-0001-7362-3581), Department of Pharmacy, The Second Hospital of Jilin University, Changchun, Jilin Province, China.

**Co-authors**:
- Jie Li (first author) — Department of Obstetrics and Gynecology, The Second Hospital of Jilin University
- Xiubo Sun, Jing Zhang, Lijie Zhai — Department of Pharmacy, The Second Hospital of Jilin University

## Funding

Jilin Province Tianhua Health Public Welfare Foundation Project (Grant No. J2025JKJ012), awarded to Ling Yu.

## Declaration of Generative AI and AI-assisted technologies in the writing process

AI-assisted writing tools were used in two scope-limited ways:

(i) code implementation for the pre-submission sensitivity analyses (Fine-Gray competing-risk model and g-formula causal mediation) after the analytic plan was finalised by the authors;

(ii) sentence-level language polishing of the Methods and Results sections only.

The Background, Discussion, and Conclusions sections were authored without AI assistance. All study design choices, statistical model selection, scientific decisions, citation choices, numerical claims, and interpretations were independently made by the authors and verified against the underlying statistical outputs. Disclosure conforms to ICMJE Recommendations (Updated January 2026, Section V), consistent with COPE 2023 guidance on authorship and AI tools.