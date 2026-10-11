# Care-seeking for STIs and STI-related symptoms among Nigerian women (NDHS 2023-24)

Analysis code for Supplementary File 2.

- `STI_script.R`: the R analysis script (packages: `haven`, `dplyr`, `survey`)
- `index.html`: a readable, syntax-highlighted version of the same script, served via GitHub Pages

## Data

The Nigeria Demographic and Health Survey (NDHS) 2023-24 dataset is not distributed with this code. Request access through the DHS Program (https://dhsprogram.com/data/), then set `data_path` in `STI_script.R` to your local copy.

## Methods summary

Survey-weighted analyses account for the complex sampling design (clusters `V021`, strata `V022`, weights `V005/1,000,000`). Prevalence ratios are estimated with modified Poisson regression (`svyglm`, quasipoisson family, log link).
