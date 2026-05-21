#
# Script to generate NetCDFs that contain ensemble stats from all NMME models,
# one file per month. Based on all models and all their members, it calculates
# the mean, median, and quantiles 0.10, 0.25, 0.75, and 0.90. These files
# represent the baseline that the EMOS approach should beat.

source("ex_08_blend/config.R")
source("functions/general_tools.R")

fs::dir_create(str_glue("{OUTPUT_DIR}/baseline"))

# Obtain all file names
nmme_files <- fs::dir_ls(
  NMME_DIR,
  regexp = "\\.nc$",
  recurse = TRUE
)

# Create a dataframe of files
df_nmme <-
  tibble(file = nmme_files) |>
  mutate(
    model = str_extract(
      file,
      paste0("(?<=/)(", paste(MODELS, collapse = "|"), ")(?=/)")
    ),
    ic_date = as_date(str_extract(file, "\\d{4}-\\d{2}-\\d{2}"))
  ) |>
  filter(str_detect(file, "leads-7")) |>
  filter(!is.na(model), !is.na(ic_date)) |>
  # filter to baseline period 1991-01-01 to 2020-12-01
  filter(ic_date >= as_date("1991-01-01"), ic_date <= as_date("2020-12-01"))

dates_to_process <-
  unique(df_nmme$ic_date) |>
  sort()

message(str_glue(
  "Found {length(dates_to_process)} initialization dates to process."
))

# Loop over each month (initialization date)
for (d in dates_to_process) {
  d_str <- as.character(as_date(d))
  message(str_glue("Processing date: {d_str}"))

  out_file <- str_glue(
    "{OUTPUT_DIR}/baseline/nmme_ensemble_precipitation_mon_ic-{d_str}_leads-7_c-afr.nc"
  )

  files_for_date <- df_nmme |> filter(ic_date == d)

  if (nrow(files_for_date) == 0) {
    message("  No files found. Skipping.")
    next
  }

  # Load each model's data for this date
  res_list <-
    seq_len(nrow(files_for_date)) |>
    map(\(i) {
      f <- files_for_date$file[i]
      mod <- files_for_date$model[i]

      r <- read_mdim(f, proxy = TRUE) |>
        filter(
          X >= BBOX["xmin"],
          X <= BBOX["xmax"],
          Y >= BBOX["ymin"],
          Y <= BBOX["ymax"]
        ) |>
        st_as_stars() |>
        aperm(c("X", "Y", "M", "L"))

      # Handle geoss2s specific dimension (only 4 members in some years)
      # (not needed)
      # if (mod == "nasa-geoss2s") {
      #   r <- r |> slice(M, 1:4)
      # }

      return(r)
    })

  # Combine all models along the member dimension 'M'
  combined_stars <- do.call(c, c(res_list, along = "M"))

  # # Identify dimensions to apply margins
  # (not needed: aperm() above ensures dimensions have always the same order)
  # dim_names <- names(st_dimensions(combined_stars))
  # m_dim <- which(dim_names == "M")
  # margin_dims <- setdiff(seq_along(dim_names), m_dim)

  # Compute stats over the combined members
  baseline_stats <- st_apply(
    combined_stars,
    MARGIN = c("X", "Y", "L"),
    FUN = \(x) {
      x <- x[!is.na(x)]
      if (length(x) == 0) {
        return(c(
          pred_mean = NA_real_,
          pred_median = NA_real_,
          pred_q05 = NA_real_,
          pred_q20 = NA_real_,
          pred_q80 = NA_real_,
          pred_q95 = NA_real_
        ))
      } else {
        q <- quantile(x, probs = c(0.05, 0.2, 0.8, 0.95)) |> unname()
        return(c(
          pred_mean = mean(x),
          pred_median = median(x),
          pred_q05 = q[1],
          pred_q20 = q[2],
          pred_q80 = q[3],
          pred_q95 = q[4]
        ))
      }
    },
    .fname = "stats"
  )

  # Split the "stats" dimension into attributes
  baseline_stats <- baseline_stats |> split("stats")

  # Save to NetCDF
  rt_write_nc(baseline_stats, out_file)

  str_glue(
    "gcloud storage mv {out_file} gs://clim_data_reg_useast1/nmme/monthly/ensemble/central_africa/precipitation/"
  ) |>
    system(ignore.stdout = T, ignore.stderr = T)
}
