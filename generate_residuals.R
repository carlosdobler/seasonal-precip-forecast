# SCRIPT T0:
# (1) CALCULATE FORECAST DEVIATIONS FROM OBSERVED
# (2) SAVE DEVIATIONS FOR USE IN A MODEL TO IMPROVE SKILL

library(tidyverse)
library(stars)
# library(mirai)

sf_use_s2(F)

# daemons(6)

source("setup.R")

roi <-
  st_bbox(roi, crs = 4326)

source("functions/general_tools.R")

dir_data <- "/mnt/pd-ssd-50/residuals"
# fs::dir_create(dir_data)

ff_models <-
  rt_gs_list_files("gs://clim_data_reg_useast1/nmme/monthly") |>
  str_subset("ensemble", negate = T)

mod_names <- basename(ff_models)

dates <-
  seq(as_date("1991-01-01"), as_date("2024-12-01"), by = "1 month")


# ******

# Download and load ERA5 climatology

ff_era_clim <-
  "gs://clim_data_reg_useast1/era5/climatologies/*precipitation*" |>
  str_glue() |>
  rt_gs_list_files() |>
  str_subset("1991-2020") |>
  rt_gs_download_files(dir_data, quiet = T)

ss_era_clim <-
  ff_era_clim |>
  map(read_mdim) |>
  map(st_crop, roi, normalize = T) |>
  suppressWarnings() |>
  suppressMessages() |>
  map(\(s) mutate(s, across(everything(), \(x) if_else(x == -9999, NA, x)))) |>
  map(mutate, clim = alpha * beta)

# *******

# Loop over models

mod <- 1

# Download and load climatology

ff_nmme_clim <-
  "gs://clim_data_reg_useast1/nmme/climatologies/{mod_names[mod]}/*precipitation*" |>
  str_glue() |>
  rt_gs_list_files() |>
  rt_gs_download_files(dir_data, quiet = T)

ss_nmme_clim <-
  ff_nmme_clim |>
  map(read_mdim) |>
  map(st_warp, ss_era_clim[[1]]) |>
  suppressWarnings() |>
  suppressMessages() |>
  map(\(s) mutate(s, across(everything(), \(x) if_else(x == -9999, NA, x)))) |>
  map(mutate, clim = alpha * beta)


# *****

# Loop over dates

d <- dates[1]
d6 <- seq(as_date(d), as_date(d) + months(5), by = "1 month")


# *****

# Download and load ERA5 monthly data

era_available_dates <-
  fs::dir_ls(dir_data) |>
  str_subset("era5") |>
  str_subset("gamma-params", negate = T) |>
  str_sub(-13, -4)

# delete unnecessary files
needed_not <-
  era_available_dates[!era_available_dates %in% d6]

if (length(needed_not > 0)) {
  needed_not |>
    walk(\(f) {
      str_glue("{dir_data}/era5_total-precipitation_mon_{f}.nc") |>
        fs::file_delete()
    })
}


needed <-
  d6[!d6 %in% era_available_dates]

if (length(needed) > 0) {
  str_glue(
    "gs://clim_data_reg_useast1/era5/monthly_means/total_precipitation/era5_total-precipitation_mon_{needed}.nc"
  ) |>
    rt_gs_download_files(dir_data, parallel = F, quiet = T) |>
    invisible()
}


ss_era <-
  str_glue("{dir_data}/era5_total-precipitation_mon_{d6}.nc") |>
  set_names(d6) |>
  map(read_mdim) |>
  map(st_crop, roi, normalize = T) |>
  suppressWarnings() |>
  suppressMessages() |>
  map(adrop) |>
  map(units::drop_units)

# s_era <-
#   do.call(c, c(s_era, along = "L")) |>
#   st_set_dimensions("L", values = seq(6)) |>
#   units::drop_units()

# ******

# Download and load NMME monthly data

mon_nmme <-
  if_else(month(d) == 1, 12, month(d) - 1) |>
  str_pad(2, "left", "0")

f_nmme <-
  str_glue(
    "{ff_models[mod]}precipitation/nmme_{mod_names[mod]}_precipitation_mon_ic-{as_date(d) - months(1)}_leads-6.nc"
  ) |>
  rt_gs_download_files(dir_data, quiet = T)

f_nmme_ba <-
  str_glue(
    "{ff_models[mod]}precipitation_biasadj/nmme_{mod_names[mod]}_precipitation_mon_ic-{as_date(d) - months(1)}_leads-6_biasadj.nc"
  ) |>
  rt_gs_download_files(dir_data, quiet = T, update_only = T)

s_nmme <-
  f_nmme |>
  read_mdim() |>
  units::drop_units() |>
  mutate(prec = prec / 1000) |> # from mm/d to m/d
  st_warp(ss_era[[1]]) |>
  st_apply(c(1, 2, 3), mean, .fname = "prec")

s_nmme_ba <-
  f_nmme_ba |>
  read_mdim() |>
  units::drop_units() |>
  st_warp(ss_era[[1]]) |>
  st_apply(c(1, 2, 3), mean, .fname = "prec")


# *****

# Loop over lead times

lead_ <- 1
mon_i <- ((month(d) + lead_ - 2) %% 12) + 1

s <-
  c(
    # ERA5 monthly data
    pluck(ss_era, lead_) |>
      setNames("era5"),

    # NMME monthly data
    slice(s_nmme, L, lead_) |>
      setNames("nmme"),

    # NMME ba monthly data
    slice(s_nmme_ba, L, lead_) |>
      setNames("nmme_ba"),

    # clim era
    pluck(ss_era_clim, mon_i) |>
      setNames(c("era5_alpha", "era5_beta", "era5_clim")),

    #
    pluck(ss_nmme_clim, as.integer(mon_nmme)) |>
      slice(L, lead_) |>
      setNames(c("nmme_alpha", "nmme_beta", "nmme_clim"))
  )

s2 <-
  s |>
  merge() |>
  st_apply(
    c(1, 2),
    \(x) {
      if (any(is.na(x[c(1, 2, 4, 7)]))) {
        #
        era5_perc = NA
        nmme_perc = NA
        nmme_perc_ba = NA
        perc = NA
        perc_ba = NA
        #
      } else {
        #
        era5_perc = round(lmom::cdfgam(x[1], c(x[4], x[5])) * 100)
        nmme_perc = round(lmom::cdfgam(x[2], c(x[7], x[8])) * 100)
        nmme_perc_ba = round(lmom::cdfgam(x[3], c(x[4], x[5])) * 100)
        perc = era5_perc - nmme_perc
        perc_ba = era5_perc - nmme_perc_ba
        #
      }

      return(
        c(
          era5_perc = era5_perc,
          nmme_perc = nmme_perc,
          nmme_perc_ba = nmme_perc_ba,
          perc = perc,
          perc_ba = perc_ba
        )
      )
    },
    .fname = "p"
  ) |>
  split("p")

s <-
  s |>
  mutate(across(everything(), \(x) round(x * 1000)))


s1 <-
  s |>
  transmute(
    # MAE
    mae = era5 - nmme,
    mae_ba = era5 - nmme_ba,

    # ANOM
    anom_era5 = era5_clim - era5,
    anom_nmme = nmme_clim - nmme,
    anom_nmme_ba = era5_clim - nmme_ba,
    anom = anom_era5 - anom_nmme,
    anom_ba = anom_era5 - anom_nmme_ba
  )


c(s, s1, s2) |>
  as_tibble() |>
  select(
    era5,
    era5_alpha,
    era5_beta,
    era5_perc,
    nmme_ba,
    nmme_perc_ba,
    perc_ba
  ) |>
  slice_sample(n = 5)


# 3 ACC but with ba and both with era clim

# save

# clean up
fs::file_delete(c(f_nmme, f_nmme_ba))
