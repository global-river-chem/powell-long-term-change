# Random forests for Si trend slopes

This workflow predicts site-level Si concentration and yield Sen slopes with
random forests. It compares static and average-driver predictors with the same
predictors plus 2002–2022 driver Sen slopes, evaluates predictions on testing
sites, and uses SHAP to explain testing-site predictions.

## Start here

Run `07_rf_trend_slopes.R` from the repository root.

Before running:

- change `data_root` to the local `long-term-change` folder
- change `master_data_root` to the local `master-datasets` folder
- export the current live site-reference table and set `site_reference_file` to
  that export
- use `stability_bootstraps <- 10L` for a quick trial or `100L` for a full run
- leave `make_plots <- TRUE` to create the three final plots
- set `refit_models <- TRUE` after changing data or model settings
- set `refit_models <- FALSE` to reuse matching saved results and recreate plots

The path comments near the top of the main script show which folders and live
export to change. Other input filenames are defined in
`rf_trend/specification.R`.

## First-time setup

Install the required packages once:

```r
install.packages(c(
  "dplyr", "tidyr", "ggplot2", "randomForest",
  "fastshap", "cowplot", "data.table"
))
```

The main script reports any missing package or input file before fitting a
model. A complete refit takes several minutes because tuning and stability
selection are repeated separately within each training set.

## Responses

The workflow fits two responses from Kathi Jo's supplied slope exports:

- Si concentration Sen slope
- Si yield Sen slope

The response is `estimates` for `chemical == "DSi"`. Each row is one site, so
the model predicts the direction and magnitude of long-term site-level change,
not an annual Si measurement.

All supplied Si slopes are retained regardless of statistical significance.
Significance is used only to distinguish points in the performance plot. NH4 is
not used, and sites are not removed solely because their Si slopes are extreme.

The response files contain 126 sites. Both models use the same 115 sites with
complete values for every selected predictor.

## Complete cases

Predictors are not imputed. A site is retained only when all 21 average
predictors and all eight trend predictors are finite. Using one complete cohort
for both predictor sets keeps their testing comparison paired.

Eleven sites are excluded:

- Ciudad Bolivar, Congo a Beach Brazzaville, SKEENA RIVER AT USK, and Saut
  Maripa lack P
- AMFREVILLE-SOUS-LES-MONTS and POSES 3 lack a specific-discharge trend
- CARRIERES-SOUS-POISSY and IVRY-SUR-SEINE lack an evapotranspiration trend
- Obidos and ST. LAWRENCE lack RCS
- Lower Atchafalaya has no valid drainage area in the live reference and
  therefore lacks specific discharge

Each saved model object records these sites and missing predictors under
`excluded_sites`. The workflow also stops if a missing or non-finite value
reaches a random forest.

## Predictor period

Annual predictors use 2002–2022. This is 21 annual time points spanning 20
elapsed years. There are no separate five-year or ten-year models.

Annual observations are first summarized to one average and one Sen slope per
site. They never become separate random-forest rows. Static watershed
characteristics are also summarized to one value per site.

The supplied response slopes were calculated previously and do not include
their record dates. This workflow standardizes the predictor period without
recalculating or replacing those response slopes.

## Predictors

The workflow compares two pre-specified predictor sets:

- `average_only`: 21 static or average-driver predictors
- `average_plus_trends`: the same 21 predictors plus eight 2002–2022 Sen slopes

The average predictors are median N, median P, NPP, evapotranspiration,
precipitation, temperature, snow cover, mean specific discharge, permafrost,
elevation, basin slope, RBI, RCS, four land-cover fractions, and four lithology
fractions.

The trend predictors are Sen slopes of N, P, NPP, evapotranspiration,
precipitation, temperature, snow cover, and specific discharge. Specific
discharge is discharge divided by drainage area.

At least five annual observations within 2002–2022 are required to calculate a
driver trend. This is an observation-coverage rule within the fixed period, not
a five-year model.

## Drainage areas

Specific discharge uses `drainSqKm` from active WRTDS rows in the latest live
site-reference export. Supplied yield slopes are rescaled by `old area / live
area`, which changes magnitude but not sign or significance. Lower Atchafalaya
is excluded because its area is undetermined.

## N and P

N combines NO3 and NOx. P combines the comparable dissolved P records. The
annual generalized WRTDS value is used first; an annual median from raw
chemistry fills a missing site-year. Raw N uses NO3 or NOx, and raw P uses SRP
or PO4. NH4 is excluded throughout.

This source combination provides N for all response sites. P remains missing
for four sites, which are excluded by the complete-case rule rather than
assigned estimated P values.

## Primary training and testing design

The primary analysis uses ten repeated outer site splits. Each split contains:

- 92 training sites, or 80% of the 115 complete sites
- 23 testing sites, or approximately 20%

The model is fitted with the 92 training sites and then predicts the 23 sites
that were not used for that fit. This process is repeated ten times with
different seeded site assignments.

There is no third site-level validation subset. The repeated testing results
provide out-of-sample validation, while tuning and stability selection remain
inside each training set.

These are independent repeated splits. A site can appear in testing data more
than once or not at all across the ten repetitions. Each testing metric is
calculated from all 23 testing sites in its corresponding split. The saved
result also reports the fraction of sites that appeared in testing at least
once.

Training and testing assignments are balanced within the same dominant
lithology groups used in the GRL workflow. Within each lithology group, sites
are also divided into three ranges of observed Si slopes. This helps each split
retain a mix of lithologies and low, middle, and high response values.

The site-average and average-plus-trend models receive identical assignments
for a response. Concentration and yield use separate seeded assignments.

## Training-only tuning and stability selection

All tuning and predictor selection occur within each outer training set:

- an initial random forest automatically selects the tree count and `mtry`
- the tree search evaluates 100 to 2,000 trees in increments of 100
- 100 bootstrap samples are drawn from the training sites only
- a predictor is selected in a bootstrap when its permutation importance is
  above the initial model's 75th-percentile importance threshold
- predictors selected in at least 95% of bootstraps are retained
- if no predictor reaches 95%, the five strongest initial predictors are used
- a second random forest retunes the tree count and `mtry` with the retained
  predictors

The second forest predicts only the corresponding testing sites. Testing sites
do not influence tuning, stability selection, or the retained predictor set.

Tree count and `mtry` are automated. Routine runs do not require manual tuning.
The preliminary plots used 10 bootstraps for a quick trial. The documented
workflow uses 100 for the full stability analysis. The stability controls are
defined in `default_rf_settings()` in
`rf_trend/specification.R`. Change them only for a planned methodological
sensitivity analysis, then set `refit_models <- TRUE`.

After validation, the workflow applies the same tuning and stability procedure
to all 115 sites and saves one final reusable model. Its out-of-bag result is
not substituted for the repeated-split testing performance.

## SHAP

SHAP values are calculated only for the average-plus-trends models. Within each
outer split, the training predictor data provide the SHAP background and SHAP
explains predictions for that split's testing sites only.

The Figure 2-style plot shows testing performance, mean absolute SHAP
importance, and SHAP distributions for concentration and yield. SHAP describes
how a fitted model used its predictors; it does not establish causation.
Interpret SHAP cautiously when testing performance is weak.

## Split sensitivity

The sensitivity analysis applies the complete tuning and stability pipeline to
four designs:

- 80% training and 20% testing, balanced by GRL dominant lithology
- 60% training and 40% testing, balanced by GRL dominant lithology
- 80% training and 20% testing, balanced by environmental cluster
- 60% training and 40% testing, balanced by environmental cluster

Each design is repeated ten times. The 60/40 design contains 69 training and 46
testing sites. Paired seeds generate one balanced site order for each
stratification method, so the 23 testing sites in an 80/20 repetition are
contained within the corresponding 46 testing sites in its 60/40 repetition.

Lithology and environmental clusters are alternative ways to balance the
splits. They are not random-forest predictors. The supplied cluster names are
retained, and missing membership is kept explicitly unassigned rather than
combined into an invented `Other` category.

## Performance metrics

Testing predictive R-squared is calculated separately for every repeated split:

```text
predictive R² = 1 - sum((observed - predicted)²) /
                    sum((observed - mean(observed))²)
```

A negative predictive R-squared means that predictions for those testing sites
were less accurate than using their mean observed slope. It is a valid testing
result and signals weak out-of-sample prediction. It should not be replaced by
squared Pearson correlation, which answers a different question. The saved
objects contain predictive R-squared, squared Pearson correlation, RMSE, and
MAE.

## Current results

The previous 100-bootstrap models used outdated drainage areas and must be
refit. In a 10-bootstrap check, the yield model with trends had median testing
R-squared 0.039 and mean -0.045; concentration remained weak. Do not report
final performance, SHAP, or sensitivity results until the full refit finishes.

## Inputs

Kathi Jo supplied the concentration and yield slope exports and the
environmental-cluster data through the [shared Google Drive results
folder](https://drive.google.com/drive/u/1/folders/10HMZLr9TAf2asrprZMyyuiYlIesm4Jug).
The workflow reads local copies and does not download directly from Google
Drive.

The results folder provides:

- `conc_slopes_export.csv`
- `yield_slopes_export.csv`
- `Si_sites_clusters_six_names.csv`

The [shared Google Drive data
folder](https://drive.google.com/drive/u/1/folders/1F16i4-dMvIvd_jTHKKnbdJALGa4PSPOJ)
needs these additional inputs:

- `all-data_si-extract_3_20260629.csv`
- `Full_Results_WRTDS_kalman_annual.csv`
- the latest export of the live `Site_Reference_Table`
- `20260105_masterdata_chem.csv`

The last file supplies raw NO3, NOx, SRP, and PO4 observations used when an
annual WRTDS nutrient value is missing.

Place local copies in this layout:

```text
SiSyn/
├── long-term-change/
│   └── data/
│       ├── conc_slopes_export.csv
│       ├── yield_slopes_export.csv
│       ├── all-data_si-extract_3_20260629.csv
│       └── Si_sites_clusters_six_names.csv
└── spatial-data-extractions/
    └── master-datasets/
        ├── Full_Results_WRTDS_kalman_annual.csv
        ├── 20260105_masterdata_chem.csv
        └── Site_Reference_Table_YYYYMMDD.csv
```

If the local files use different folders, change `data_root` and
`master_data_root` in `07_rf_trend_slopes.R`. Set `site_reference_file` there to
the latest live export. If another input filename changes, update it in
`rf_input_files()` in `rf_trend/specification.R`.

## Outputs

The workflow writes only three plots and compact RDS result objects to Box.
Nothing is written to the GitHub repository.

Each bootstrap setting has its own folder, so a full run cannot replace the
preliminary plots. Plot folders are:

- `rf_outputs/10_bootstraps/plots_only/`
- `rf_outputs/100_bootstraps/plots_only/`

Each contains:

- `Figure2_Si_slope_RF_SHAP.png`
- `Predictor_set_comparison.png`
- `Split_sensitivity_80-20_vs_60-40.png`

The matching `model_results/` folder allows plots to be recreated without
refitting. No intermediate CSV files or one-off analysis scripts are created.

## Safe edits

Common changes have one clear location:

- change local paths and run switches in `07_rf_trend_slopes.R`
- change `stability_bootstraps` there to select a 10- or 100-bootstrap run
- change years, seeds, repeated-split settings, tuning, or stability controls
  in `rf_trend/specification.R`
- change predictor lists in `predictor_specification()` in that same file
- change plot labels or appearance in `rf_trend/figure2.R`

Data, predictor, or model-setting changes require a refit. Plot-only changes do
not.

## Code map

- `07_rf_trend_slopes.R` runs the complete workflow
- `rf_trend/specification.R` defines inputs, responses, predictors, and settings
- `rf_trend/data-prep.R` creates the site-level complete-case table
- `rf_trend/resampling.R` creates splits, tunes models, selects stable
  predictors, and calculates testing SHAP
- `rf_trend/modeling.R` fits and stores the primary and final models
- `rf_trend/sensitivity.R` compares split proportions and stratification
- `rf_trend/figure2.R` creates the three final plots
