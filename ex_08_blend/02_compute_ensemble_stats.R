#' 02_compute_ensemble_stats.R
#'
#' Purpose: Compute per-model ensemble mean and spread (standard deviation) from the aligned NMME data.
#'
#' Inputs:
#' - ex_08_blend/data/nmme_aligned_{model}.rds
#'
#' Outputs:
#' - ex_08_blend/data/nmme_stats_{model}.rds

source("ex_08_blend/config.R")

for (mod in MODELS) {
  message(str_glue("Processing {mod}"))

  # Load only the aligned data for this model
  stars_obj <- read_rds(
    str_glue(
      "{OUTPUT_DIR}/nmme_aligned_{mod}.rds"
    ),
    mod
  )

  # Identify the member dimension. Usually it's named "M".
  dim_names <- names(st_dimensions(stars_obj))
  m_dim <- which(dim_names == "M")

  margin_dims <- setdiff(seq_along(dim_names), m_dim)

  # Compute mean and sd
  mod_stats <- st_apply(
    stars_obj,
    MARGIN = margin_dims,
    FUN = \(x) c(mean = mean(x, na.rm = T), sd = sd(x, na.rm = T)),
    .fname = "stats"
  )

  mod_stats <-
    mod_stats |>
    split("stats")

  write_rds(mod_stats, str_glue("{OUTPUT_DIR}/nmme_stats_{mod}.rds"))
}

message("Step 02 complete!")
