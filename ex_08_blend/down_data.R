# script to download NMME and ERA5 data from bucket to local dir

library(tidyverse)
rt <- new.env()
source("functions/general_tools.R", local = rt)

source("nmme_sources_df.R")

walk(df_sources$model, \(mod) {
  print(mod)
  ff <-
    rt$rt_gs_list_files(
      str_glue("gs://clim_data_reg_useast1/nmme/monthly/{mod}/precipitation/")
    ) |>
    str_subset(str_flatten(str_glue("ic-{seq(2021,2025)}"), "|"))

  dir_mod <- str_glue("/mnt/pd-blend/nmme/{mod}")
  fs::dir_create(dir_mod)

  ff |>
    rt$rt_gs_download_files(dir_mod) |>
    invisible()
})


ff <-
  rt$rt_gs_list_files(
    "gs://clim_data_reg_useast1/era5/monthly_means/total_precipitation"
  ) |>
  str_subset(str_flatten(str_glue("_{seq(2021,2025)}"), "|"))

dir_mod <- str_glue("/mnt/pd-nmme-residuals/")
# fs::dir_create(dir_mod)

ff |>
  rt$rt_gs_download_files(dir_mod) |>
  invisible()
