#' Shared configuration and utilities for the blending pipeline
#'
#' This script contains shared paths, constants, and helper functions
#' used across the ensemble blending pipeline.

suppressPackageStartupMessages({
  library(tidyverse)
  library(stars)
  library(tidymodels)
  library(furrr)
  # library(glmnet)
  library(here)
})

# Source the general tools if needed
source(here::here("functions", "general_tools.R"))

# Source model names
source(here::here("nmme_sources_df.R"))

# Paths
NMME_DIR <- "/mnt/pd-blend/nmme/"
OBS_DIR <- "/mnt/pd-nmme-residuals/"
OUTPUT_DIR <- "/mnt/pd-blend/out/"

# Model names from nmme_sources_df
MODELS <- df_sources$model

# Bounding box for Africa domain
BBOX <- c(xmin = 7, ymin = -18, xmax = 50, ymax = 16)

#' Softplus activation function
#'
#' Enforces strict positivity.
#'
#' @param x Numeric vector
#' @return Numeric vector of softplus transformed values
#' @export
softplus <- function(x) {
  ifelse(x > 20, x, log1p(exp(x)))
}

#' CRPS for Gamma distribution
#'
#' Computes the closed-form Continuous Ranked Probability Score (CRPS)
#' for a Gamma distribution parameterized by mean and dispersion.
#' The variance is assumed to be phi * mu^2.
#'
#' @param y Numeric vector of observations
#' @param mu Numeric vector of predicted means (must be > 0)
#' @param phi Numeric vector of predicted dispersions (must be > 0)
#' @return Numeric vector of CRPS values
#' @export
crps_gamma <- function(y, mu, phi) {
  # Gamma parameterization:
  # shape alpha = 1 / phi
  # scale s = mu * phi
  shape <- 1 / phi
  scale <- mu * phi

  # Ensure strict positivity to avoid NaNs
  shape <- pmax(shape, .Machine$double.eps)
  scale <- pmax(scale, .Machine$double.eps)

  term1 <- y * (2 * pgamma(y, shape = shape, scale = scale) - 1)
  term2 <- shape * scale * (2 * pgamma(y, shape = shape + 1, scale = scale) - 1)
  term3 <- (scale / sqrt(pi)) * exp(lgamma(shape + 0.5) - lgamma(shape))

  term1 - term2 - term3
}
