calculate_model_stats <- function(fcst_file, bbox) {
  mod <-
    fcst_file |>
    fs::path_file() |>
    str_split_i("_", i = 2)

  r <- rlang::inject(
    read_mdim(fcst_file, proxy = TRUE) |>
      filter(
        X >= !!bbox[["xmin"]],
        X <= !!bbox[["xmax"]],
        Y >= !!bbox[["ymin"]],
        Y <= !!bbox[["ymax"]]
      ) |>
      st_as_stars()
  )

  # if (mod == "nasa-geoss2s") {
  #   r <- r |> slice(M, 1:4)
  # }

  dim_names <- names(st_dimensions(r))
  m_dim <- which(dim_names == "M")
  margin_dims <- setdiff(seq_along(dim_names), m_dim)

  mod_stats <- st_apply(
    r,
    MARGIN = margin_dims,
    FUN = \(x) c(mean = mean(x, na.rm = T), sd = sd(x, na.rm = T)),
    .fname = "stats"
  ) |>
    split("stats")

  df_mod <- as_tibble(mod_stats)

  # Remove NAs
  df_mod <- df_mod |> filter(!is.na(mean))

  df_mod |>
    mutate(model = str_replace_all(mod, "-", "_"))
}


# ************

apply_emos <- function(df, ic_date, frozen_model) {
  #
  all_models <- unique(df$model)

  df_fcst <-
    df |>
    pivot_wider(names_from = model, values_from = c("mean", "sd"))

  message(str_glue(
    "Processing forecast for ic date: {ic_date} | {length(all_models)} models found"
  ))

  df_fcst <-
    df_fcst |>
    filter(L > 0) |>
    mutate(
      target_date = as_date(ic_date) %m+% months(L)
    )

  # Convert arrays to data.frames and join
  # Note: frozen_model arrays have dimensions X, Y, L
  df_a <-
    as.data.frame.table(frozen_model$a_array, responseName = "a") |>
    as_tibble()

  df_c <-
    as.data.frame.table(frozen_model$c_array, responseName = "c") |>
    as_tibble()

  # b and d have an extra model dimension
  df_b <-
    as.data.frame.table(frozen_model$b_array, responseName = "b_val") |>
    as_tibble() |>
    mutate(model = str_replace_all(model, "-", "_")) |>
    pivot_wider(names_from = model, values_from = b_val, names_prefix = "b_")

  df_d <-
    as.data.frame.table(frozen_model$d_array, responseName = "d_val") |>
    as_tibble() |>
    mutate(model = str_replace_all(model, "-", "_")) |>
    pivot_wider(names_from = model, values_from = d_val, names_prefix = "d_")

  # Merge parameters
  df_params <-
    df_a |>
    inner_join(df_c, by = c("X", "Y", "L")) |>
    inner_join(df_b, by = c("X", "Y", "L")) |>
    inner_join(df_d, by = c("X", "Y", "L")) |>
    mutate(
      X = as.numeric(as.character(X)),
      Y = as.numeric(as.character(Y)),
      L = as.numeric(as.character(L))
    )

  df_result <-
    df_fcst |>
    inner_join(df_params, by = c("X", "Y", "L"))

  # Compute mu and phi (needed to get gamma params)
  # mu = softplus(a + (b1 * mean_model1) + (b2 * mean_model2) + ... + (b7 * mean_model7))
  # first, bx * mean_modelx and dx * sd_modelx
  for (mod in all_models) {
    df_result[[paste0("m_term_", mod)]] <-
      df_result[[paste0(
        "b_",
        mod
      )]] *
      df_result[[paste0(
        "mean_",
        mod
      )]]

    # phi
    df_result[[paste0("s_term_", mod)]] <-
      df_result[[paste0(
        "d_",
        mod
      )]] *
      df_result[[paste0(
        "sd_",
        mod
      )]]
  }

  # sum all terms
  m_terms <-
    df_result |>
    select(starts_with("m_term_")) |>
    rowSums(na.rm = TRUE)

  s_terms <-
    df_result |>
    select(starts_with("s_term_")) |>
    rowSums(na.rm = TRUE)

  df_result <- df_result |>
    mutate(
      # add intercept and apply softplus
      mu = softplus(a + m_terms),
      phi = softplus(c + s_terms),

      # get gamma params
      shape = 1 / phi,
      scale = mu * phi,

      # derive quantities
      pred_mean = mu,
      pred_median = qgamma(0.5, shape = shape, scale = scale),
      pred_q05 = qgamma(0.05, shape = shape, scale = scale),
      pred_q20 = qgamma(0.2, shape = shape, scale = scale),
      pred_q80 = qgamma(0.8, shape = shape, scale = scale),
      pred_q95 = qgamma(0.95, shape = shape, scale = scale)
    )

  # Convert back to stars object
  # rt_write_nc expects lon (x) and lat (y) to be the first two dimensions.
  # Our dimensions are X, Y, L.
  stars_out <-
    df_result |>
    select(
      X,
      Y,
      L,
      # mu,
      # phi,
      pred_mean,
      pred_median,
      pred_q05,
      pred_q20,
      pred_q80,
      pred_q95
    ) |>
    st_as_stars(dims = c("X", "Y", "L"))

  return(stars_out)
}
