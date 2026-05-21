#
# Script to fit EMOS-LASSO models for each grid cell, target calendar month, and lead time.
# Uses glmnet directly for LASSO fitting and rsample for leave-one-year-out (LOYO) CV.
# Generates frozen parameters (weights) arrays to be applied on new forecast runs

source("ex_08_blend/config.R")

plan(multicore, workers = parallelly::availableCores() - 1)

# Load and merge NMME stats

df_nmme <-
  map_dfr(MODELS, \(mod) {
    message(str_glue("  Loading {mod}..."))
    mod_stats <- read_rds(str_glue("{OUTPUT_DIR}/nmme_stats_{mod}.rds"))

    df_mod <- as_tibble(mod_stats) |>
      rename(ic_date = time) |>
      filter(!is.na(mean)) |>
      mutate(model = str_replace_all(mod, "-", "_"))

    return(df_mod)
  })

df_nmme <-
  pivot_wider(df_nmme, names_from = model, values_from = c(mean, sd))

# Add target_date to NMME based on lead time
df_nmme <-
  df_nmme |>
  filter(L > 0) |>
  mutate(target_date = ic_date %m+% months(L))

# Load ERA5 observations

era5 <-
  read_rds(str_glue("{OUTPUT_DIR}/era5_aligned.rds"))

# Merge NMME and ERA5
df_era5 <-
  as_tibble(era5) |>
  rename(obs = 4) |> # 4th column should be the variable name
  filter(!is.na(obs))

df_merged <-
  df_nmme |>
  inner_join(df_era5, by = join_by(X, Y, target_date == time)) |>
  mutate(
    target_month = month(target_date),
    init_month = month(ic_date),
    target_year = year(target_date)
  ) |>
  filter(target_year >= 1991, target_year <= 2020) # Training set

# Free up memory
rm(df_nmme, era5, df_era5)
gc()

# ***********

# Tuning stage

# tuning grid of 5 lambda values (lambda parameter: amount of regularization)
lambda_grid <-
  grid_regular(penalty(range = c(-4, 0)), levels = 5) |>
  pull(penalty)

mean_vars <- paste0("mean_", str_replace_all(MODELS, "-", "_"))
spread_vars <- paste0("sd_", str_replace_all(MODELS, "-", "_"))

# Helper: fit glmnet with Gamma, falling back to Gaussian on convergence failure
# LASSO: alpha = 1
fit_glmnet_safe <- function(X, y, lambda, alpha = 1, lower.limits = 0) {
  converged <- TRUE

  fit <- withCallingHandlers(
    tryCatch(
      glmnet::glmnet(
        X,
        y,
        family = Gamma(link = "identity"),
        lambda = lambda,
        alpha = alpha,
        lower.limits = lower.limits
      ),
      error = function(e) {
        converged <<- FALSE
        NULL
      }
    ),
    warning = function(w) {
      if (grepl("Convergence.*not reached", conditionMessage(w))) {
        converged <<- FALSE
        invokeRestart("muffleWarning")
      }
    }
  )

  # convergence failure
  if (!converged || is.null(fit)) {
    fit <- glmnet::glmnet(
      X,
      y,
      family = "gaussian",
      lambda = lambda,
      alpha = alpha,
      lower.limits = lower.limits
    )
  }

  fit
}

# Function to fit and evaluate one cell/month/lead
fit_cell_model <- function(df_subset) {
  if (nrow(df_subset) < 10) {
    return(NULL)
  }

  # less than 1 mm/month
  if (mean(df_subset$obs, na.rm = TRUE) < 0.033) {
    return(list(
      a = NA_real_,
      b = rep(NA_real_, length(mean_vars)),
      c = NA_real_,
      d = rep(NA_real_, length(spread_vars)),
      lambda = NA_real_
    ))
  }

  # LOYO CV
  # Tune lambda using CPRS as loss function
  splits <- loo_cv(df_subset)
  set.seed(9)
  splits <- splits |> slice_sample(n = 10) # only 10 samples

  crps_results <- numeric(length(lambda_grid))

  # loop through lambda values (grid)
  for (l_idx in seq_along(lambda_grid)) {
    lam <- lambda_grid[l_idx]
    crps_lam <- 0

    for (sp in splits$splits) {
      train_data <- training(sp)
      test_data <- testing(sp)

      if (nrow(test_data) == 0) {
        next
      }

      # 1. mean model
      X_train <- as.matrix(train_data[, mean_vars])
      y_train <- train_data$obs
      X_test <- as.matrix(test_data[, mean_vars])

      fit_m <- fit_glmnet_safe(X_train, y_train, lam)

      pred_m_train <- as.vector(predict(fit_m, newx = X_train, s = lam))
      pred_m_test <- as.vector(predict(fit_m, newx = X_test, s = lam))

      mu_train <- softplus(pred_m_train)
      mu_test <- softplus(pred_m_test)

      # 2. spread model
      y_var_train <- (y_train - mu_train)^2 # error of mean model
      X_spread_train <- as.matrix(train_data[, spread_vars])
      X_spread_test <- as.matrix(test_data[, spread_vars])

      # error is used as target for spread model
      fit_s <- fit_glmnet_safe(X_spread_train, y_var_train, lam)

      pred_s_test <- as.vector(predict(fit_s, newx = X_spread_test, s = lam))
      phi_test <- softplus(pred_s_test)

      # calculate CRPS
      crps_val <- crps_gamma(test_data$obs, mu_test, phi_test)
      crps_lam <- crps_lam + mean(crps_val, na.rm = TRUE)
    }
    # mean crps across folds
    crps_results[l_idx] <- crps_lam / length(splits$splits)
  }

  # identify best lamda value from grid
  min_indices <- which(crps_results == min(crps_results))
  best_lam <- max(lambda_grid[min_indices]) # max in case of ties

  # extract data to fit final mean model
  X_all <- as.matrix(df_subset[, mean_vars]) # features
  y_all <- df_subset$obs # target

  # fit final mean model
  fit_m_final <- fit_glmnet_safe(X_all, y_all, best_lam)

  # apply softplus to predictions
  mu_all <-
    predict(
      fit_m_final,
      newx = X_all,
      s = best_lam
    ) |>
    as.vector() |>
    softplus()

  # error = target for spread model
  y_var_all <- (y_all - mu_all)^2

  # features
  X_spread_all <- as.matrix(df_subset[, spread_vars])

  # fit final spread model
  fit_s_final <- fit_glmnet_safe(X_spread_all, y_var_all, best_lam)

  # get coefficients from path at best_lam value
  coef_m <- as.vector(coef(fit_m_final, s = best_lam))
  coef_s <- as.vector(coef(fit_s_final, s = best_lam))

  # final parameters
  list(
    a = coef_m[1], # intercept mean model
    b = coef_m[2:8], # coefficients/weights mean model
    c = coef_s[1], # intercept spread model
    d = coef_s[2:8], # coefficients/weights spread model
    lambda = best_lam
  )
}


# Fitting stage

for (init_month in seq(12)) {
  #seq(12)) {
  message(str_glue("   init month = {init_month}"))

  # subset data to init month
  nested_df <-
    df_merged |>
    filter(init_month == {{ init_month }}) |>
    group_by(X, Y, L) |>
    nest() |>
    ungroup()

  message(str_glue("   Total models to fit: {nrow(nested_df)}"))

  results <-
    nested_df |>
    mutate(model_params = future_map(data, fit_cell_model, .progress = TRUE))

  results <-
    results |>
    select(-data) |>
    filter(!map_lgl(model_params, is.null))

  message("   Done!")

  # Format frozen parameters

  lons <- sort(unique(results$X))
  lats <- sort(unique(results$Y))
  leads <- sort(unique(results$L))

  # Create empty arrays

  # Helper
  # (only for intercepts...edit to use for coefficients/weights)
  # (MOVE OUT OF FOR LOOP)
  create_empty_array <- function(dim_vec) {
    array(
      NA,
      dim = dim_vec,
      dimnames = list(
        X = lons,
        Y = lats,
        L = leads
      )
    )
  }

  # mean model intercepts
  a_array <- create_empty_array(c(length(lons), length(lats), length(leads)))

  # spread model intercepts
  c_array <- create_empty_array(c(length(lons), length(lats), length(leads)))

  # lambdas
  lambda_array <-
    create_empty_array(
      c(
        length(lons),
        length(lats),
        length(leads)
      )
    )

  # mean model weights
  b_array <-
    array(
      NA,
      dim = c(length(lons), length(lats), length(leads), length(MODELS)),
      dimnames = list(
        X = lons,
        Y = lats,
        L = leads,
        model = MODELS
      )
    )

  # spread model weights
  d_array <- array(
    NA,
    dim = c(length(lons), length(lats), length(leads), length(MODELS)),
    dimnames = list(
      X = lons,
      Y = lats,
      L = leads,
      model = MODELS
    )
  )

  # Populate empty arrays

  for (i in seq_len(nrow(results))) {
    r <- results[i, ]
    lon_idx <- as.character(r$X)
    lat_idx <- as.character(r$Y)
    l_idx <- as.character(r$L)

    p <- r$model_params[[1]]
    a_array[lon_idx, lat_idx, l_idx] <- p$a
    c_array[lon_idx, lat_idx, l_idx] <- p$c
    lambda_array[lon_idx, lat_idx, l_idx] <- p$lambda
    b_array[lon_idx, lat_idx, l_idx, ] <- p$b
    d_array[lon_idx, lat_idx, l_idx, ] <- p$d
  }

  frozen_model <-
    list(
      a_array = a_array,
      b_array = b_array,
      c_array = c_array,
      d_array = d_array,
      lambda_array = lambda_array,
      lons = lons,
      lats = lats,
      leads = leads,
      models = MODELS
    )

  write_rds(
    frozen_model,
    str_glue(
      "{OUTPUT_DIR}/frozen_model_initmonth_{str_pad(init_month, 2, 'left', '0')}.rds"
    )
  )
}
