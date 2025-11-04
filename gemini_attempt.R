# Load libraries
library(tidyverse)
library(stars)
library(furrr)
library(future)

# Source the general tools
source("functions/general_tools.R")

# Main parameters
models <- c(
  "canesm5",
  "cola-rsmas-ccsm4",
  "gfdl-spear",
  "nasa-geoss2s",
  "ncep-cfsv2"
)
years <- 2021:2024
target_month <- 7 # July
dest_dir <- "/mnt/pd-ssd-50/residuals/"

# Set up parallel processing
plan(multicore, workers = 8)

# Create destination directory if it doesn't exist
if (!dir.exists(dest_dir)) {
  dir.create(dest_dir, recursive = TRUE)
}

# Function to generate file URLs for a given model and year
generate_urls <- function(year, model) {
  # Observation for July of the given year
  obs_file <- sprintf(
    "gs://clim_data_reg_useast1/era5/monthly_means/total_precipitation/era5_total-precipitation_mon_%d-%02d-01.nc",
    year,
    target_month
  )

  # Observational climatology for July
  obs_clim_file <- sprintf(
    "gs://clim_data_reg_useast1/era5/monthly_means/climatologies/era5_total-precipitation_mon_gamma-params_1991-2020_%02d.nc",
    target_month
  )

  # Raw NMME forecast initialized in January of the given year
  raw_fcst_file <- sprintf(
    "gs://clim_data_reg_useast1/nmme/monthly/%s/precipitation/nmme_%s_precipitation_mon_ic-%d-01-01_leads-6.nc",
    model,
    model,
    year
  )

  # Bias-adjusted NMME forecast initialized in January of the given year
  bias_adj_fcst_file <- sprintf(
    "gs://clim_data_reg_useast1/nmme/monthly/%s/precipitation_biasadj/nmme_%s_precipitation_mon_ic-%d-01-01_leads-6_biasadj.nc",
    model,
    model,
    year
  )

  # NMME climatology for forecasts initialized in January
  nmme_clim_file <- sprintf(
    "gs://clim_data_reg_useast1/nmme/climatologies/%s/nmme_%s_precipitation_mon_gamma-params_1991-2020_01_leads-6.nc",
    model,
    model
  )

  list(
    obs = obs_file,
    obs_clim = obs_clim_file,
    raw_fcst = raw_fcst_file,
    bias_adj_fcst = bias_adj_fcst_file,
    nmme_clim = nmme_clim_file
  )
}

# --- Main Processing Loop ---
process_model_year <- function(model, year) {
  message(paste("Processing:", model, "-", year))

  # Generate file URLs
  urls <- generate_urls(year, model)

  # Download all necessary files
  unique_urls <- unique(unlist(urls))
  rt_gs_download_files(unique_urls, dest = dest_dir, quiet = TRUE)

  # Create a named list of local paths
  local_files <- purrr::map(urls, ~ file.path(dest_dir, fs::path_file(.)))

  # --- Read data ---
  obs_july <- read_ncdf(local_files$obs)
  obs_clim <- read_ncdf(local_files$obs_clim)
  raw_fcst <- read_ncdf(local_files$raw_fcst)
  bias_adj_fcst <- read_ncdf(local_files$bias_adj_fcst)
  nmme_clim <- read_ncdf(local_files$nmme_clim)

  # --- Pre-processing ---
  # Select July forecast from NMME data (lead time = 7 months from Jan)
  raw_fcst_july <- raw_fcst |> slice(L, target_month)
  bias_adj_fcst_july <- bias_adj_fcst |> slice(L, target_month)

  # --- Skill Comparison ---

  # 1. Mean Absolute Difference
  mae_raw <- abs(raw_fcst_july - obs_july)
  mae_bias_adj <- abs(bias_adj_fcst_july - obs_july)

  # 2. Anomaly Comparison
  # Calculate climatology means from gamma parameters (shape * scale)
  obs_clim_mean <- adrop((obs_clim[[1]] * obs_clim[[2]]))

  nmme_clim_july <- nmme_clim |> slice(L, target_month)
  nmme_clim_mean <- nmme_clim_july[[1]] * nmme_clim_july[[2]]

  # Calculate anomalies
  obs_anomaly <- obs_july - obs_clim_mean
  raw_fcst_anomaly <- raw_fcst_july - nmme_clim_mean
  bias_adj_fcst_anomaly <- bias_adj_fcst_july - nmme_clim_mean

  # Difference in anomalies
  anomaly_diff_raw <- abs(raw_fcst_anomaly - obs_anomaly)
  anomaly_diff_bias_adj <- abs(bias_adj_fcst_anomaly - obs_anomaly)

  # --- Save results ---
  # Combine results into a single stars object
  skill_results <- c(
    mae_raw = mae_raw,
    mae_bias_adj = mae_bias_adj,
    anomaly_diff_raw = anomaly_diff_raw,
    anomaly_diff_bias_adj = anomaly_diff_bias_adj,
    along = "variable"
  )

  output_filename <- file.path(
    dest_dir,
    paste0("skill_comparison_", model, "_", year, ".nc")
  )
  rt_write_nc(skill_results, filename = output_filename)

  return(output_filename)
}

# --- Run the analysis ---
# Create a grid of models and years
model_year_grid <- expand.grid(
  model = models,
  year = years,
  stringsAsFactors = FALSE
)

# Use future_map to run in parallel
future_map2(model_year_grid$model, model_year_grid$year, process_model_year)

message("Skill comparison script finished.")
