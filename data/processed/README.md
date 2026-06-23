# data/processed/ — analytic intermediates

Only small CSV summary files are checked in here.

Large binary intermediates (nhanes_design.RData, nhanes_final.RData,
nhanes_raw_merged.RData; ~80 MB combined) are excluded from this
repository (see `.gitignore`) for two reasons:
1. Repository size hygiene — keep `git clone` < 1 GB.
2. Re-identification gray zone — although NHANES raw fields are public
   domain, a re-assembled wide table with 1000+ joined fields is more
   conservatively re-built than redistributed in bulk.

## Rebuild instructions

From the repo root:
```
Rscript scripts/01_download_data.R   # ~700 MB to data/raw/
Rscript scripts/02_merge_data.R      # writes data/processed/nhanes_raw_merged.RData
Rscript scripts/03_clean_data.R      # writes data/processed/nhanes_final.RData + flow_counts.csv
Rscript scripts/04_survey_design.R   # writes data/processed/nhanes_design.RData
```
Then `Rscript scripts/run_all.R` re-runs the analytic pipeline.
