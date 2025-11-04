library(stars)
library(dplyr)
library(purrr)
library(stringr)
library(lubridate)

# Source the general tools
source("functions/general_tools.R")

# Define parameters
dest_dir <- "/mnt/pd-ssd-50/residuals/"
years <- 2021:2024
target_month <- 7 # July
models <- c("CFSv2", "CanCM4i", "GEM-NEMO") # Example models - adjust as needed

# Create destination directory if needed
if (!fs::dir_exists(dest_dir)) {
  fs::dir_create(dest_dir)
}

# Step 1: Download ERA5 observational data
# For July forecasts with 6-month leads, we need data from January (lead 6) to July (lead 0)
message("Downloading ERA5 observational data...")

era5_files <- map(years, \(year) {
  # For each year, download January through July
  months_to_download <- 1:7

  map_chr(months_to_download, \(month) {
    date_str <- sprintf("%04d-%02d-01", year, month)
    file_path <- str_glue(
      "gs://clim_data_reg_useast1/era5/monthly_means/total_precipitation/era5_total-precipitation_mon_{date_str}.nc"
    )
    file_path
  })
}) |>
  unlist()

era5_local <- rt_gs_download_files(era5_files, dest_dir)

# Step 2: Download ERA5 climatology for July
message("Downloading ERA5 climatology for July...")
era5_clim_file <- "gs://clim_data_reg_useast1/era5/monthly_means/climatologies/era5_total-precipitation_mon_gamma-params_1991-2020_07.nc"
era5_clim_local <- rt_gs_download_files(era5_clim_file, dest_dir)

# Step 3: Download NMME forecast data (raw and bias-adjusted)
message("Downloading NMME forecast data...")

nmme_files <- map(models, \(model) {
  # For each model and year, download forecasts initialized in January through July
  # that will forecast July
  map(years, \(year) {
    # January to July initial conditions
    init_months <- 1:7

    map(init_months, \(init_month) {
      init_date <- sprintf("%04d-%02d-01", year, init_month)

      # Raw forecast
      raw_file <- str_glue(
        "gs://clim_data_reg_useast1/nmme/monthly/{model}/precipitation/nmme_{model}_precipitation_mon_ic-{init_date}_leads-6.nc"
      )

      # Bias-adjusted forecast
      biasadj_file <- str_glue(
        "gs://clim_data_reg_useast1/nmme/monthly/{model}/precipitation_biasadj/nmme_{model}_precipitation_mon_ic-{init_date}_leads-6_biasadj.nc"
      )

      list(raw = raw_file, biasadj = biasadj_file)
    })
  })
}) |>
  set_names(models)

# Download raw files
nmme_raw_local <- map(models, \(model) {
  files <- nmme_files[[model]] |>
    unlist(recursive = FALSE) |>
    map_chr("raw")
  rt_gs_download_files(files, dest_dir)
}) |>
  set_names(models)

# Download bias-adjusted files
nmme_biasadj_local <- map(models, \(model) {
  files <- nmme_files[[model]] |>
    unlist(recursive = FALSE) |>
    map_chr("biasadj")
  rt_gs_download_files(files, dest_dir)
}) |>
  set_names(models)

# Step 4: Download NMME climatologies
message("Downloading NMME climatologies...")

nmme_clim_files <- map(models, \(model) {
  # Need climatologies for each initialization month (Jan-Jul)
  map_chr(1:7, \(init_month) {
    str_glue(
      "gs://clim_data_reg_useast1/nmme/climatologies/{model}/nmme_{model}_precipitation_mon_gamma-params_1991-2020_{sprintf('%02d', init_month)}_leads-6.nc"
    )
  })
}) |>
  set_names(models)

nmme_clim_local <- map(models, \(model) {
  rt_gs_download_files(nmme_clim_files[[model]], dest_dir)
}) |>
  set_names(models)

# Step 5: Load and process data
message("Processing data and calculating skill metrics...")

# Load ERA5 observations
era5_obs <- map(era5_local, read_stars)

# Load ERA5 climatology
era5_clim <- read_stars(era5_clim_local)

# Function to calculate lead time for July target
calc_lead_for_july <- function(init_month) {
  # Lead 0 = same month, lead 1 = next month, etc.
  lead <- (7 - init_month) %% 12
  if (init_month > 7) {
    lead <- lead + 12
  } # Handle previous year
  return(lead)
}

# Function to extract forecast for specific lead
extract_lead <- function(stars_obj, lead) {
  stars_obj |>
    filter(L == lead)
}

# Initialize results storage
skill_results <- list()

# Process each model
for (model in models) {
  message(str_glue("Processing model: {model}"))

  mae_raw <- list()
  mae_biasadj <- list()
  mae_anomaly_raw <- list()
  mae_anomaly_biasadj <- list()

  file_idx <- 1

  for (year in years) {
    for (init_month in 1:7) {
      # Calculate which lead month corresponds to July
      lead <- calc_lead_for_july(init_month)

      # Skip if lead is beyond 6 months (outside forecast range)
      if (lead > 6) {
        next
      }

      # Get observation for July of this year
      obs_idx <- which(str_detect(era5_local, sprintf("%04d-07-01", year)))
      if (length(obs_idx) == 0) {
        next
      }

      obs <- era5_obs[[obs_idx]]

      # Load forecasts
      raw_fc <- read_stars(nmme_raw_local[[model]][file_idx]) |>
        extract_lead(lead)

      biasadj_fc <- read_stars(nmme_biasadj_local[[model]][file_idx]) |>
        extract_lead(lead)

      # Calculate absolute errors
      mae_raw[[length(mae_raw) + 1]] <- abs(raw_fc - obs)
      mae_biasadj[[length(mae_biasadj) + 1]] <- abs(biasadj_fc - obs)

      # Load climatology for this initialization month
      clim_idx <- init_month
      nmme_clim <- read_stars(nmme_clim_local[[model]][clim_idx]) |>
        extract_lead(lead)

      # Calculate anomalies (forecast - climatology mean)
      # Assuming the climatology has shape and scale parameters
      # Mean of gamma distribution = shape * scale
      raw_anomaly <- raw_fc - (nmme_clim[[1]] * nmme_clim[[2]])
      biasadj_anomaly <- biasadj_fc - (nmme_clim[[1]] * nmme_clim[[2]])
      obs_anomaly <- obs - (era5_clim[[1]] * era5_clim[[2]])

      # Calculate absolute error of anomalies
      mae_anomaly_raw[[length(mae_anomaly_raw) + 1]] <- abs(
        raw_anomaly - obs_anomaly
      )
      mae_anomaly_biasadj[[length(mae_anomaly_biasadj) + 1]] <- abs(
        biasadj_anomaly - obs_anomaly
      )

      file_idx <- file_idx + 1
    }
  }

  # Calculate mean absolute error across all forecasts
  skill_results[[model]] <- list(
    mae_raw = Reduce(`+`, mae_raw) / length(mae_raw),
    mae_biasadj = Reduce(`+`, mae_biasadj) / length(mae_biasadj),
    mae_anomaly_raw = Reduce(`+`, mae_anomaly_raw) / length(mae_anomaly_raw),
    mae_anomaly_biasadj = Reduce(`+`, mae_anomaly_biasadj) /
      length(mae_anomaly_biasadj)
  )

  # Calculate skill improvement (negative values = bias adjustment improved)
  skill_results[[model]]$mae_improvement <-
    skill_results[[model]]$mae_biasadj - skill_results[[model]]$mae_raw

  skill_results[[model]]$mae_anomaly_improvement <-
    skill_results[[model]]$mae_anomaly_biasadj -
    skill_results[[model]]$mae_anomaly_raw
}

# Step 6: Save results
message("Saving results...")

for (model in models) {
  # Save MAE grids
  rt_write_nc(
    skill_results[[model]]$mae_raw,
    str_glue("{dest_dir}/mae_raw_{model}_july_2021-2024.nc")
  )

  rt_write_nc(
    skill_results[[model]]$mae_biasadj,
    str_glue("{dest_dir}/mae_biasadj_{model}_july_2021-2024.nc")
  )

  rt_write_nc(
    skill_results[[model]]$mae_improvement,
    str_glue("{dest_dir}/mae_improvement_{model}_july_2021-2024.nc")
  )

  # Save anomaly MAE grids
  rt_write_nc(
    skill_results[[model]]$mae_anomaly_raw,
    str_glue("{dest_dir}/mae_anomaly_raw_{model}_july_2021-2024.nc")
  )

  rt_write_nc(
    skill_results[[model]]$mae_anomaly_biasadj,
    str_glue("{dest_dir}/mae_anomaly_biasadj_{model}_july_2021-2024.nc")
  )

  rt_write_nc(
    skill_results[[model]]$mae_anomaly_improvement,
    str_glue("{dest_dir}/mae_anomaly_improvement_{model}_july_2021-2024.nc")
  )
}

message("Skill assessment complete!")

# Return results for inspection
skill_results
