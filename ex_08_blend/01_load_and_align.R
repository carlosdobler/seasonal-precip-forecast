#' 01_load_and_align.R
#'
#' Purpose: Load hindcast (NMME models) and observation (ERA5) data into stars objects,
#' filter them to the Africa bounding box, load into memory, and align on a common
#' spatiotemporal grid (1x1 degree models' resolution) using st_warp(method="average").
#'
#' Inputs:
#' - NMME netcdf files in /mnt/pd-blend/nmme/
#' - ERA5 total-precipitation netcdf files in /mnt/pd-nmme-residuals/
#'
#' Outputs:
#' - ex_08_blend/data/nmme_aligned_{model}.rds
#' - ex_08_blend/data/era5_aligned.rds

source("ex_08_blend/config.R")

# Obtain all file names
nmme_files <- fs::dir_ls(
  NMME_DIR,
  regexp = "\\.nc$",
  recurse = T
)

# Create a dataframe of files
df_nmme <- tibble(file = nmme_files) |>
  mutate(
    model = str_extract(
      file,
      paste0("(?<=/)(", paste(MODELS, collapse = "|"), ")(?=/)")
    ),
    ic_date = as_date(str_extract(file, "\\d{4}-\\d{2}-\\d{2}"))
  ) |>
  filter(!is.na(model), !is.na(ic_date)) |>
  # filter to hindcast period 1991 to 2025
  # filter(year(ic_date) >= 1991, year(ic_date) <= 2025)
  filter(year(ic_date) >= 1991, year(ic_date) <= 2020) # 2021 - 2025 not ready


# target_grid <- NULL

# Load, crop, and save NMME data

plan(multicore, workers = 10)

for (mod in MODELS) {
  message(str_glue("Processing {mod}..."))
  files <- df_nmme |> filter(model == mod) |> pull(file)

  # Load each file, filter, and load into memory
  # Then bind them along the 'time' (ic_date) dimension
  res <- future_map(files, function(f) {
    r <-
      read_mdim(f, proxy = TRUE) |>
      filter(
        X >= BBOX["xmin"],
        X <= BBOX["xmax"],
        Y >= BBOX["ymin"],
        Y <= BBOX["ymax"]
      ) |>
      st_as_stars()

    if (mod == "nasa-geoss2s") {
      r <-
        r |>
        slice(M, 1:4)
    }

    return(r)
  })

  # Combine along time
  res <- do.call(c, c(res, along = "time")) |>
    st_set_dimensions(
      "time",
      values = df_nmme |> filter(model == mod) |> pull(ic_date)
    )

  write_rds(res, str_glue("{OUTPUT_DIR}/nmme_aligned_{mod}.rds"))
}


# Same for ERA5

era5_files <- fs::dir_ls(
  OBS_DIR,
  regexp = "era5_total-precipitation_mon_.*\\.nc$"
)

df_era5 <- tibble(file = era5_files) |>
  mutate(
    date = as_date(str_extract(file, "\\d{4}-\\d{2}-\\d{2}"))
  ) |>
  filter(!is.na(date)) |>
  # filter(year(date) >= 1991, year(date) <= 2025)
  filter(year(date) >= 1991, year(date) <= 2020)

res_era5 <- future_map(df_era5$file, function(f) {
  read_mdim(f, proxy = TRUE) |>
    filter(
      longitude >= BBOX["xmin"],
      longitude <= BBOX["xmax"],
      latitude >= BBOX["ymin"],
      latitude <= BBOX["ymax"]
    ) |>
    st_as_stars() |>
    adrop()
})

res_era5 <- do.call(c, c(res_era5, along = "time")) |>
  st_set_dimensions(
    "time",
    values = df_era5 |> pull(date)
  )

# Regrid ERA5 to the grid of the first NMME model
era5_regridded <- st_warp(
  res_era5,
  res |>
    slice(M, 1) |>
    slice(L, 1) |>
    slice(time, 1) |>
    adrop(),
  use_gdal = T,
  method = "average"
)

era5_regridded <-
  era5_regridded |>
  setNames("tp")

st_dimensions(era5_regridded) <- st_dimensions(
  res |>
    slice(M, 1) |>
    slice(L, 1)
)

# Convert ERA5 units
# Since it's m/day, multiply by 1000 to get mm/day
era5_regridded <- era5_regridded * 1000

write_rds(era5_regridded, str_glue("{OUTPUT_DIR}/era5_aligned.rds"))

plan(sequential)
message("Step 01 complete!")
