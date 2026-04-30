#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(tibble)
  library(forecast)
  library(purrr)
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
config_path <- if (length(args) >= 1L) args[1] else file.path(code_dir, "config", "sarima_config.json")
config_dir <- dirname(normalizePath(config_path, winslash = "/", mustWork = FALSE))
data_dir <- shared_noaa_data_dir_project(code_dir)
result_dir <- file.path(code_dir, "Result")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

cfg <- read_config_json(config_path)
model_name <- as.character(cfg$model_name %||% "SARIMA")
test_months <- as.integer(cfg$test_months %||% 24L)
start_date <- as.Date(cfg$start_date %||% "2000-01-01")
prediction_level <- as.numeric(cfg$prediction_level %||% 0.90)
seasonal_period <- as.integer(cfg$seasonal_period %||% 12L)
sarima_max_p <- as.integer(cfg$max_p %||% 3L)
sarima_max_q <- as.integer(cfg$max_q %||% 3L)
sarima_max_P <- as.integer(cfg$max_P %||% 3L)
sarima_max_Q <- as.integer(cfg$max_Q %||% 3L)
sarima_max_order <- as.integer(cfg$max_order %||% 8L)
sarima_verbose <- isTRUE(cfg$verbose %||% TRUE)

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
  prediction_level = prediction_level,
  seasonal_period = seasonal_period,
  max_p = sarima_max_p,
  max_q = sarima_max_q,
  max_P = sarima_max_P,
  max_Q = sarima_max_Q,
  max_order = sarima_max_order,
  verbose = sarima_verbose
)
save_json_pretty_project(resolved_config, file.path(model_output_dir, "resolved_config.json"))

sarima_log <- function(...) {
  if (sarima_verbose) cat(sprintf(...), "\n")
}

calc_aicc <- function(fit, n_obs) {
  k <- length(fit$coef) + 1L
  aic <- AIC(fit)
  if (n_obs <= k + 1L) {
    return(Inf)
  }
  aic + (2 * k * (k + 1)) / (n_obs - k - 1)
}

fit_sarima_grid <- function(y, site_code) {
  n_obs <- length(y)
  candidates <- expand.grid(
    p = 0:sarima_max_p,
    d = 0:2,
    q = 0:sarima_max_q,
    P = 0:sarima_max_P,
    D = 0:2,
    Q = 0:sarima_max_Q
  ) |>
    dplyr::filter((p + q + P + Q) <= sarima_max_order)

  sarima_log("[SARIMA][%s] Evaluating %d candidates", site_code, nrow(candidates))

  fits <- vector("list", nrow(candidates))
  n_converged <- 0L

  for (i in seq_len(nrow(candidates))) {
    spec <- candidates[i, ]
    fit_try <- tryCatch(
      stats::arima(
        y,
        order = c(spec$p, spec$d, spec$q),
        seasonal = list(order = c(spec$P, spec$D, spec$Q), period = seasonal_period),
        method = "ML"
      ),
      error = function(e) NULL
    )

    if (is.null(fit_try)) next

    n_converged <- n_converged + 1L
    fits[[i]] <- list(
      fit = fit_try,
      order = c(spec$p, spec$d, spec$q, spec$P, spec$D, spec$Q, seasonal_period),
      aicc = calc_aicc(fit_try, n_obs),
      bic = BIC(fit_try)
    )
  }

  fits <- purrr::compact(fits)
  if (!length(fits)) {
    usage_error("No SARIMA specification converged for site `", site_code, "`.")
  }

  best <- fits[[which.min(purrr::map_dbl(fits, "aicc"))]]
  sarima_log(
    "[SARIMA][%s] Selected (%d,%d,%d)(%d,%d,%d)[%d], AICc=%.3f, BIC=%.3f, converged=%d/%d",
    site_code,
    best$order[1], best$order[2], best$order[3],
    best$order[4], best$order[5], best$order[6], best$order[7],
    best$aicc, best$bic, n_converged, nrow(candidates)
  )
  best
}

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

z_alpha <- qnorm((1 + prediction_level) / 2)

fit_one_sarima <- function(df_train, df_test) {
  df_train <- df_train |> dplyr::arrange(date)
  df_test <- df_test |> dplyr::arrange(date)
  site_code <- as.character(df_train$site[1])

  y <- stats::ts(
    df_train$value,
    start = c(lubridate::year(min(df_train$date)), lubridate::month(min(df_train$date))),
    frequency = seasonal_period
  )

  best <- fit_sarima_grid(y, site_code)
  fit <- best$fit
  fc <- forecast::forecast(fit, h = nrow(df_test), level = prediction_level * 100)
  fitted_vals <- as.numeric(y - residuals(fit))
  sigma_hat <- sqrt(as.numeric(fit$sigma2))

  train_df <- tibble(
    site = site_code,
    date = df_train$date,
    actual = df_train$value,
    pred_mean = fitted_vals,
    pred_sd = sigma_hat,
    pred_lo = fitted_vals - z_alpha * sigma_hat,
    pred_hi = fitted_vals + z_alpha * sigma_hat,
    split = "train"
  )

  test_df <- tibble(
    site = site_code,
    date = df_test$date,
    actual = df_test$value,
    pred_mean = as.numeric(fc$mean),
    pred_sd = as.numeric(fc$upper[, 1] - fc$mean) / z_alpha,
    pred_lo = as.numeric(fc$lower[, 1]),
    pred_hi = as.numeric(fc$upper[, 1]),
    split = "test"
  )

  list(
    site = site_code,
    fit = fit,
    train = train_df,
    test = test_df,
    order = tibble(
      site = site_code,
      p = best$order[1],
      d = best$order[2],
      q = best$order[3],
      P = best$order[4],
      D = best$order[5],
      Q = best$order[6],
      period = best$order[7],
      aicc = best$aicc,
      bic = best$bic
    )
  )
}

sarima_runs <- purrr::map(site_levels_project, function(site_code) {
  fit_one_sarima(
    dplyr::filter(train, site == site_code),
    dplyr::filter(test, site == site_code)
  )
}) |> stats::setNames(site_levels_project)

train_pred <- purrr::map_dfr(sarima_runs, "train")
test_pred <- purrr::map_dfr(sarima_runs, "test")
predictions <- bind_rows(train_pred, test_pred)
order_summary <- purrr::map_dfr(sarima_runs, "order")

metrics <- tibble(
  model = model_name,
  train_rmse = rmse_project(train_pred$actual, train_pred$pred_mean),
  train_mae = mae_project(train_pred$actual, train_pred$pred_mean),
  train_coverage = coverage_project(train_pred$actual, train_pred$pred_lo, train_pred$pred_hi),
  test_rmse = rmse_project(test_pred$actual, test_pred$pred_mean),
  test_mae = mae_project(test_pred$actual, test_pred$pred_mean),
  test_coverage = coverage_project(test_pred$actual, test_pred$pred_lo, test_pred$pred_hi),
  test_mlpd = mlpd_gaussian_project(test_pred$actual, test_pred$pred_mean, pmax(test_pred$pred_sd, 1e-8)),
  n_train = nrow(train_pred),
  n_test = nrow(test_pred)
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
readr::write_csv(order_summary, result_path("key_result", "selected_orders.csv"))
readr::write_csv(acf_df, result_path("time_series_diagnose", "residual_acf.csv"))

saveRDS(
  list(
    model = model_name,
    train = train_pred,
    test = test_pred,
    predictions = predictions,
    metrics = metrics,
    selected_orders = order_summary,
    residual_acf = acf_df,
    runs = sarima_runs
  ),
  result_path("key_result", "sarima_results.rds")
)

message("Finished. Results written to: ", normalize_maybe(model_output_dir))
