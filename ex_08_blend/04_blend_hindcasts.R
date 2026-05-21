#
# Script to apply frozen parameters to NMME hindcast (1991-2020) data
# Generates final, weighted ensembles with all stats

source("ex_08_blend/config.R")
source("functions/general_tools.R")
source("ex_08_blend/functions_blend.R")

dir_hindcasts <- "/mnt/pd-blend/hindcasts"
fs::dir_create(dir_hindcasts)

# Obtain all file names
nmme_files <- fs::dir_ls(
  NMME_DIR,
  regexp = "\\.nc$",
  recurse = TRUE
)

all_dates <-
  seq(
    as_date("1991-01-01"),
    as_date("2020-12-01"),
    by = "1 month"
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


walk(seq(12), \(init_month) {
  # load frozen model
  frozen_model <- read_rds(str_glue(
    "{OUTPUT_DIR}/frozen_model_initmonth_{str_pad(init_month, 2, 'left', '0')}.rds"
  ))

  # Get dates corresponding to this init_month
  dates_for_month <- all_dates[month(all_dates) == init_month]

  walk(dates_for_month, \(date) {
    #
    df_model_stats <-
      df_nmme |>
      filter(ic_date == date) |>
      pull(file) |>
      map_dfr(calculate_model_stats, bbox = BBOX)

    s_emos <-
      apply_emos(df_model_stats, date, frozen_model)

    out_file <- str_glue(
      "{dir_hindcasts}/nmme_ensemble_precipitation_mon_ic-{date}_leads-7_c-afr_emos-v01.nc"
    )
    rt_write_nc(s_emos, out_file)

    str_glue(
      "gcloud storage mv {out_file} gs://clim_data_reg_useast1/nmme/monthly/ensemble/central_africa/precipitation_emos_v01/"
    ) |>
      system(ignore.stdout = T, ignore.stderr = T)
  })
})

fs::dir_delete(dir_hindcasts)
