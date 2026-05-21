# Evaluation calibrated NMME


``` r
# load libraries
source(here::here("ex_08_blend/config.R"))
library(tidyverse)
library(patchwork)
library(fs)
library(colorspace)
library(yardstick)

# directories
dir_baseline <- file.path(OUTPUT_DIR, "eval_baseline")
dir_emos <- file.path(OUTPUT_DIR, "eval_emos")
dir_plots <- file.path(OUTPUT_DIR, "eval_plots")

# download NMME (baseline) and emos-corrected data (skip: already downloaded)
# fs::dir_create(c(dir_baseline, dir_emos, dir_plots))
# system(
#   str_glue(
#     "gcloud storage cp gs://clim_data_reg_useast1/nmme/monthly/ensemble/central_africa/precipitation/* {dir_baseline}/"
#   ),
#   ignore.stdout = TRUE,
#   ignore.stderr = TRUE
# )

# system(
#   str_glue(
#     "gcloud storage cp gs://clim_data_reg_useast1/nmme/monthly/ensemble/central_africa/precipitation_emos_v01/* {dir_emos}/"
#   ),
#   ignore.stdout = TRUE,
#   ignore.stderr = TRUE
# )

# load observational data (ERA5)
era5_df <-
  read_rds(file.path(OUTPUT_DIR, "era5_aligned.rds")) |>
  as.data.frame() |>
  rename(obs = tp, target_date = time) |>
  mutate(
    X = as.numeric(X),
    Y = as.numeric(Y)
  )

# calculate climatologies
era5_clim <-
  era5_df |>
  mutate(target_month = month(target_date), target_year = year(target_date)) |>
  filter(target_year >= 1991, target_year <= 2020) |>
  group_by(X, Y, target_month) |>
  summarise(clim_mean = mean(obs, na.rm = TRUE), .groups = "drop")

# load and prepare NMME (baseline)
base_files <- fs::dir_ls(dir_baseline, regexp = "\\.nc$")
df_base <- map_dfr(base_files, function(f) {
  ic_date <- as_date(str_extract(f, "\\d{4}-\\d{2}-\\d{2}"))
  # Only baseline period 1991-2020
  if (year(ic_date) < 1991 || year(ic_date) > 2020) {
    return(NULL)
  }
  r <- read_mdim(f, proxy = FALSE) |> slice(L, 2:7)
  as.data.frame(r) |>
    select(X, Y, L, pred_mean) |>
    mutate(
      ic_date = ic_date,
      X = as.numeric(X),
      Y = as.numeric(Y),
      L = as.numeric(L)
    )
}) |>
  rename(base_mean = pred_mean)

# load and prepare NMME+EMOS
emos_files <- fs::dir_ls(dir_emos, regexp = "\\.nc$")
df_emos <- map_dfr(emos_files, function(f) {
  ic_date <- as_date(str_extract(f, "\\d{4}-\\d{2}-\\d{2}"))
  if (year(ic_date) < 1991 || year(ic_date) > 2020) {
    return(NULL)
  }
  r <- read_mdim(f, proxy = FALSE) |> st_as_stars()
  as.data.frame(r) |>
    select(X, Y, L, pred_mean) |>
    # no data in NMME+EMOS means low precip (< 1 mm/month): flag w -9999
    mutate(pred_mean = if_else(is.na(pred_mean), -9999, pred_mean)) |>
    mutate(
      ic_date = ic_date,
      X = as.numeric(X),
      Y = as.numeric(Y),
      L = as.numeric(L)
    )
}) |>
  rename(emos_mean = pred_mean)

# join datasets
df_joined <-
  df_base |>
  inner_join(df_emos, by = c("X", "Y", "L", "ic_date")) |>
  mutate(
    target_date = ic_date %m+% months(L),
    target_month = month(target_date),
    target_year = year(target_date)
  ) |>
  filter(target_year >= 1991, target_year <= 2020)

df_joined <- df_joined |>
  inner_join(era5_df, by = c("X", "Y", "target_date")) |>
  inner_join(era5_clim, by = c("X", "Y", "target_month")) |>
  # if low precip, assign NA to that cell/target month
  mutate(
    across(
      c(obs, clim_mean, emos_mean, base_mean),
      ~ if_else(emos_mean == -9999, NA_real_, .x)
    )
  )


# if NMME+EMOS has no variability over the years it means
# the EMOS pipeline found no skill in any of the NMME models;
# fallback was to assign climatology, which means all years get
# the same value. Add some noise to be able to calculate metrics
# that require sd or variance
df_joined_metrics <-
  df_joined |>
  group_by(X, Y, target_month, L) |>
  mutate(
    emos_sd = sd(emos_mean, na.rm = TRUE),
    emos_mean = if_else(
      !is.na(emos_sd) & emos_sd == 0,
      emos_mean + rnorm(n(), mean = 0, sd = 1e-5),
      emos_mean
    )#,
    # base_sd = sd(base_mean, na.rm = TRUE),
    # base_mean = if_else(
    #   !is.na(base_sd) & base_sd == 0,
    #   base_mean + rnorm(n(), mean = 0, sd = 1e-5),
    #   base_mean
    # ),
    # clim_sd = sd(clim_mean, na.rm = TRUE),
    # clim_mean = if_else(
    #   !is.na(clim_sd) & clim_sd == 0,
    #   clim_mean + rnorm(n(), mean = 0, sd = 1e-5),
    #   clim_mean
    # )
  ) |>
  ungroup() |>
  select(-emos_sd)


# function to calculate skill metrics via yardstick
fn_compute_metrics <- function(
  df,
  metric_fn = yardstick::rmse,
  include_climatology = TRUE
) {
  # Group the dataframe once
  df_grouped <- df |> group_by(X, Y, target_month, L)

  # Yardstick calculates metrics for all groups natively and much faster
  res_base <- metric_fn(df_grouped, truth = obs, estimate = base_mean) |>
    rename(metric_base = .estimate) |>
    select(-.metric, -.estimator)

  res_emos <- metric_fn(df_grouped, truth = obs, estimate = emos_mean) |>
    rename(metric_emos = .estimate) |>
    select(-.metric, -.estimator)

  # Join the results back together
  res_out <- res_base |>
    left_join(res_emos, by = c("X", "Y", "target_month", "L"))

  if (include_climatology) {
    res_clim <- metric_fn(df_grouped, truth = obs, estimate = clim_mean) |>
      rename(metric_clim = .estimate) |>
      select(-.metric, -.estimator)

    res_out <- res_clim |>
      left_join(res_out, by = c("X", "Y", "target_month", "L"))
  }

  res_out
}

# function to plot skill maps
fn_plot_maps <- function(
  df,
  t_mon,
  lead,
  metric_name = "Metric",
  divergent_palette = FALSE
) {
  d_sub <- df |> filter(target_month == t_mon, L == lead)

  if (nrow(d_sub) == 0) {
    return(invisible(NULL))
  }

  has_clim <- "metric_clim" %in% names(d_sub)

  d_long <- d_sub |>
    select(X, Y, starts_with("metric_")) |>
    pivot_longer(
      cols = starts_with("metric_"),
      names_to = "model",
      values_to = "metric"
    ) |>
    mutate(
      model = case_match(
        model,
        "metric_clim" ~ "Climatology",
        "metric_base" ~ "NMME",
        "metric_emos" ~ "NMME+EMOS"
      )
    )

  if (has_clim) {
    d_long$model <- factor(
      d_long$model,
      levels = c("NMME", "NMME+EMOS", "Climatology")
    )
  } else {
    d_long$model <- factor(d_long$model, levels = c("NMME", "NMME+EMOS"))
  }

  p <- ggplot(d_long, aes(x = X, y = Y, fill = metric)) +
    geom_raster() +
    facet_wrap(~model) +
    labs(
      subtitle = str_glue(
        "Lead Time: {lead} months"
      ),
      fill = metric_name
    ) +
    theme(legend.position = "bottom",
          axis.title = element_blank()) +
    coord_fixed(expand = F)

  if (divergent_palette) {
    p <- p +
      scale_fill_continuous_diverging(
        palette = "Blue-Red 3",
        rev = TRUE,
        na.value = "grey90",
        limits = quantile(d_long$metric, c(0.02, 0.98), na.rm = T),
        oob = scales::squish
      )
  } else {
    p <- p +
      scale_fill_continuous_sequential(
        palette = "Viridis",
        rev = F,
        na.value = "grey90",
        limits = quantile(d_long$metric, c(0.02, 0.98), na.rm = T),
        trans = scales::modulus_trans(2),
        oob = scales::squish
      )
  }

  p + guides(fill = guide_colorbar(barwidth = 12, barheight = 0.8))
}

# function to plot skill over lead times
fn_plot_skill_leads <- function(df, t_mon, metric_name = "Metric") {
  d_sub <- df |> filter(target_month == t_mon)

  has_clim <- "metric_clim" %in% names(d_sub)

  d_long <- d_sub |>
    select(X, Y, L, starts_with("metric_")) |>
    pivot_longer(
      cols = starts_with("metric_"),
      names_to = "model",
      values_to = "metric"
    ) |>
    mutate(
      model = case_match(
        model,
        "metric_clim" ~ "Climatology",
        "metric_base" ~ "NMME",
        "metric_emos" ~ "NMME+EMOS"
      )
    )

  if (has_clim) {
    d_long$model <- factor(
      d_long$model,
      levels = c("NMME", "NMME+EMOS", "Climatology")
    )
  } else {
    d_long$model <- factor(d_long$model, levels = c("NMME", "NMME+EMOS"))
  }

  ggplot(mapping = aes(factor(L), y = metric, color = model)) +
    geom_jitter(
      data = d_long |> slice_sample(n = 5000),
      position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.7),
      alpha = 0.3,
      show.legend = F
    ) +
    stat_summary(
      data = d_long,
      aes(group = model, fill = model),
      fun = mean,
      geom = "point",
      shape = 23,
      size = 3,
      color = "black",
      stroke = 1,
      position = position_dodge(width = 0.7)
    ) +
    scale_color_discrete_qualitative(name = NULL, nmax = 3) +
    scale_fill_discrete_qualitative(name = NULL, nmax = 3) +
    labs(x = "Lead Time (months)", y = metric_name) +
    theme(legend.position = "bottom") +
    coord_cartesian(ylim = quantile(d_long$metric, c(0.03, 0.97), na.rm = T))
}
```

## RMSE (lower is better)

``` r
df_rmse <-
  fn_compute_metrics(
    df_joined_metrics,
    metric_fn = yardstick::rmse,
    include_climatology = T
  )
```

### Target month: January

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-2-1.png)

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-3-1.png)

(^ each point represents a grid cell)

### Target month: February

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-4-1.png)

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-5-1.png)

### Target month: March

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-6-1.png)

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-7-1.png)

### Target month: April

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-8-1.png)

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-9-1.png)

### Target month: May

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-10-1.png)

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-11-1.png)

### Target month: June

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-12-1.png)

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-13-1.png)

### Target month: July

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-14-1.png)

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-15-1.png)

### Target month: August

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-16-1.png)

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-17-1.png)

### Target month: September

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-18-1.png)

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-19-1.png)

### Target month: October

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-20-1.png)

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-21-1.png)

### Target month: November

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-22-1.png)

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-23-1.png)

### Target month: December

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-24-1.png)

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-25-1.png)

## CCC (Concordance Correlation Coefficient; higher is better)

``` r
df_ccc <-
  fn_compute_metrics(
    df_joined_metrics,
    metric_fn = yardstick::ccc,
    include_climatology = F
  )
```

### Target month: January

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-27-1.png)

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-28-1.png)

(^ each point represents a grid cell)

### Target month: February

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-29-1.png)

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-30-1.png)

### Target month: March

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-31-1.png)

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-32-1.png)

### Target month: April

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-33-1.png)

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-34-1.png)

### Target month: May

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-35-1.png)

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-36-1.png)

### Target month: June

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-37-1.png)

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-38-1.png)

### Target month: July

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-39-1.png)

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-40-1.png)

### Target month: August

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-41-1.png)

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-42-1.png)

### Target month: September

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-43-1.png)

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-44-1.png)

### Target month: October

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-45-1.png)

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-46-1.png)

### Target month: November

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-47-1.png)

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-48-1.png)

### Target month: December

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-49-1.png)

![](06_evaluate_historical_files/figure-commonmark/unnamed-chunk-50-1.png)
