# Experiment 8: Ensemble Blending for Precipitation Forecasts

This directory contains an operational ensemble blending pipeline for monthly precipitation forecasts in Central Africa. The pipeline performs statistical post-processing of ensemble forecasts from the NMME project using an EMOS-LASSO (Ensemble Model Output Statistics with LASSO regularization) approach to generate calibrated predictive distributions.

## Pipeline Steps

The pipeline (WIP) consists of the following modular scripts:

- `config.R`: Shared configuration, paths, and helper functions (e.g., Gamma CRPS calculation).
- `00_baseline_stats.R`: Generates baseline ensemble statistics from raw NMME models.
- `01_load_and_align.R`: Loads hindcast (NMME) and observation (ERA5) data, filtering and aligning them to a common 1°x1° spatiotemporal grid.
- `02_compute_ensemble_stats.R`: Computes the ensemble mean and spread for each model.
- `03_training.R`: Fits per-cell, per-calendar-month, per-lead EMOS-LASSO models using leave-one-year-out cross-validation via `tidymodels` and `glmnet`. Outputs frozen model parameters.
- `04_blend_hindcasts.R`: Applies the frozen parameters to historical NMME data to generate the final, weighted EMOS ensembles.
- `06_evaluate_historical.qmd`: A Quarto document that generates visual maps and charts of the skill metrics (RMSE and CCC).

## Results

Results show an improvement over the raw NMME ensemble based on the Concordance Correlation Coefficient (CCC) across target months and lead times (each point is a grid cell of the AOI):

**January**

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-28-1.png)

**February**

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-30-1.png)

**March**

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-32-1.png)

**April**

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-34-1.png)

**May**

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-36-1.png)

**June**

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-38-1.png)

**July**

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-40-1.png)

**August**

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-42-1.png)

**September**

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-44-1.png)

**October**

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-46-1.png)

**November**

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-48-1.png)

**December**

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-50-1.png)
