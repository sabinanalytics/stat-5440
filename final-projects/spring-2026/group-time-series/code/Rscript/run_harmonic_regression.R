#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(tibble)
  library(broom)
})

args <- commandArgs(trailingOnly = TRUE)
script_dir <- local({
  args_all <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args_all, value = TRUE)
  if (length(file_arg) > 0L) {
    dirname(normalizePath(sub("^--file=", "", file_arg[1])))
  } else {
    normalizePath(getwd())
  }
})

source(file.path(script_dir, "_project_model_common.R"), local = FALSE)

code_dir <- find_code_dir()
config_path <- if (length(args) >= 1L) args[1] else file.path(code_dir, "config", "harmonic_config.json")
config_dir <- dirname(normalizePath(config_path, winslash = "/", mustWork = FALSE))
data_dir <- shared_noaa_data_dir_project(code_dir)
result_dir <- file.path(code_dir, "Result")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

cfg <- read_config_json(config_path)
model_name <- as.character(cfg$model_name %||% "Harmonic_Regression")
test_months <- as.integer(cfg$test_months %||% 24L)
start_date <- as.Date(cfg$start_date %||% "2000-01-01")
prediction_level <- as.numeric(cfg$prediction_level %||% 0.90)
if (!is.finite(prediction_level) || prediction_level <= 0 || prediction_level >= 1) {
  usage_error("`prediction_level` must be between 0 and 1.")
}

data_spec <- resolve_data_spec_project(
  data_mode = cfg$data %||% "default",
  data_dir = data_dir,
  base_dirs = c(".", config_dir, code_dir)
)

model_output_dir <- file.path(result_dir, model_name)
dir.create(model_output_dir, recursive = TRUE, showWarnings = FALSE)

result_subdir <- function(section) {
  path <- file.path(model_output_dir, section)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

result_path <- function(section, filename) {
  file.path(result_subdir(section), filename)
}

resolved_config <- list(
  model_name = model_name,
  data = cfg$data %||% "default",
  resolved_data = data_spec$description,
  test_months = test_months,
  start_date = format(start_date),
  prediction_level = prediction_level
)
save_json_pretty_project(resolved_config, file.path(model_output_dir, "resolved_config.json"))

cached_data <- NULL
if (tolower(cfg$data %||% "default") == "default") {
  cached_data <- load_prepared_data_cache_project(data_dir, start_date = start_date, test_months = test_months)
}

if (is.null(cached_data)) {
  dat4 <- prepare_panel_data_project(data_spec, start_date = start_date)
  split <- make_train_test_project(dat4, test_months = test_months)
  save_prepared_data_cache_project(data_dir, dat4, split, start_date = start_date, test_months = test_months)
  train <- split$train
  test <- split$test
} else {
  dat4 <- cached_data$aligned
  train <- cached_data$train
  test <- cached_data$test
}

fit_harm <- stats::lm(
  value ~ 0 + site + site:t + site:sin1 + site:cos1 + site:sin2 + site:cos2,
  data = train,
  weights = pmax(nvalue, 1)
)

pred_with_interval <- function(fit, newdat, level) {
  p <- predict(fit, newdata = newdat, interval = "prediction", level = level, se.fit = TRUE)
  sd_pred <- sqrt(p$se.fit^2 + p$residual.scale^2)
  tibble::tibble(
    pred_mean = p$fit[, "fit"],
    pred_lo = p$fit[, "lwr"],
    pred_hi = p$fit[, "upr"],
    pred_sd = sd_pred
  )
}

train_pred <- bind_cols(train, pred_with_interval(fit_harm, train, prediction_level))
test_pred <- bind_cols(test, pred_with_interval(fit_harm, test, prediction_level))

metrics <- tibble(
  model = model_name,
  train_rmse = rmse_project(train_pred$value, train_pred$pred_mean),
  train_mae = mae_project(train_pred$value, train_pred$pred_mean),
  train_coverage = coverage_project(train_pred$value, train_pred$pred_lo, train_pred$pred_hi),
  test_rmse = rmse_project(test_pred$value, test_pred$pred_mean),
  test_mae = mae_project(test_pred$value, test_pred$pred_mean),
  test_coverage = coverage_project(test_pred$value, test_pred$pred_lo, test_pred$pred_hi),
  test_mlpd = mlpd_gaussian_project(test_pred$value, test_pred$pred_mean, pmax(test_pred$pred_sd, 1e-8)),
  n_train = nrow(train_pred),
  n_test = nrow(test_pred)
)

coef_tbl <- broom::tidy(fit_harm)
readr::write_csv(coef_tbl, result_path("key_result", "coefficients.csv"))

cf <- stats::coef(fit_harm)
get_cf <- function(site, suffix) unname(cf[paste0("site", site, ":", suffix)])

seasonal_summary <- tibble(
  site = site_levels_project,
  beta_sin1 = vapply(site_levels_project, get_cf, numeric(1), suffix = "sin1"),
  beta_cos1 = vapply(site_levels_project, get_cf, numeric(1), suffix = "cos1"),
  beta_sin2 = vapply(site_levels_project, get_cf, numeric(1), suffix = "sin2"),
  beta_cos2 = vapply(site_levels_project, get_cf, numeric(1), suffix = "cos2")
) |>
  mutate(
    amp_h1 = harmonic_amplitude_project(beta_sin1, beta_cos1),
    amp_h2 = harmonic_amplitude_project(beta_sin2, beta_cos2),
    phase_h1 = atan2(beta_sin1, beta_cos1) * 180 / pi
  )

slopes <- tibble(
  site = site_levels_project,
  slope_ppm_per_year = vapply(
    site_levels_project,
    function(s) unname(cf[paste0("site", s, ":t")]),
    numeric(1)
  ) * 12
)

latent_trend <- bind_rows(train, test) |>
  select(site, date, t) |>
  mutate(
    trend = vapply(site, function(s) unname(cf[paste0("site", s)]), numeric(1)) +
      vapply(site, function(s) unname(cf[paste0("site", s, ":t")]), numeric(1)) * t
  )

predictions <- bind_rows(
  train_pred |>
    transmute(split = "train", site, date, actual = value, pred_mean, pred_lo, pred_hi, pred_sd),
  test_pred |>
    transmute(split = "test", site, date, actual = value, pred_mean, pred_lo, pred_hi, pred_sd)
)

acf_df <- predictions |>
  filter(split == "train") |>
  mutate(residual = actual - pred_mean) |>
  group_by(site) |>
  group_modify(~ {
    x <- .x$residual
    x <- x[is.finite(x)]
    if (length(x) < 3L) {
      return(tibble(lag = numeric(0), acf = numeric(0)))
    }
    acf_obj <- stats::acf(x, plot = FALSE, na.action = na.pass, lag.max = min(24L, length(x) - 1L))
    tibble(lag = as.numeric(acf_obj$lag)[-1], acf = as.numeric(acf_obj$acf)[-1])
  }) |>
  ungroup()

readr::write_csv(predictions, result_path("forecast", "predictions.csv"))
readr::write_csv(metrics, result_path("forecast", "metrics.csv"))
readr::write_csv(acf_df, result_path("time_series_diagnose", "residual_acf.csv"))
readr::write_csv(seasonal_summary, result_path("key_result", "seasonal_summary.csv"))
readr::write_csv(slopes, result_path("key_result", "slopes.csv"))
readr::write_csv(latent_trend, result_path("key_result", "latent_trend.csv"))

saveRDS(
  list(
    model = model_name,
    fit = fit_harm,
    train = train_pred,
    test = test_pred,
    metrics = metrics,
    coefficients = coef_tbl,
    seasonal_summary = seasonal_summary,
    slopes = slopes,
    latent_trend = latent_trend,
    residual_acf = acf_df
  ),
  result_path("key_result", "harmonic_regression_results.rds")
)

message("Finished. Results written to: ", normalize_maybe(model_output_dir))
