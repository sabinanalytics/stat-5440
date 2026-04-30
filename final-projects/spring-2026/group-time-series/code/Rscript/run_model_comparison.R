#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(lubridate)
  library(stringr)
  library(scales)
  library(rstan)
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
config_path <- if (length(args) >= 1L) args[1] else file.path(code_dir, "config", "comparison_config.json")
result_root <- file.path(code_dir, "Result")
plot_root <- file.path(code_dir, "Plot")
dir.create(result_root, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_root, recursive = TRUE, showWarnings = FALSE)

site_levels <- site_levels_project
site_colors <- c(BRW = "#1f77b4", MLO = "#d62728", SMO = "#2ca02c", SPO = "#9467bd")

theme_project <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(face = "bold")
    )
}

usage_error <- function(...) {
  stop(..., call. = FALSE)
}

normalize_maybe <- function(path) {
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

require_file <- function(path, label) {
  if (!file.exists(path)) {
    usage_error(label, " not found: ", normalize_maybe(path))
  }
}

read_config <- function(path) {
  require_file(path, "Comparison config")
  cfg <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  if (!is.list(cfg) || is.null(cfg$models)) {
    usage_error("Comparison config must contain a `models` array.")
  }
  cfg
}

comparison_subdir <- function(root, compare_name, section) {
  path <- file.path(root, "comparison", compare_name, section)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

comparison_result_path <- function(compare_name, section, filename) {
  file.path(comparison_subdir(result_root, compare_name, section), filename)
}

comparison_plot_path <- function(compare_name, section, filename) {
  file.path(comparison_subdir(plot_root, compare_name, section), filename)
}

save_plot_file <- function(compare_name, section, filename, p, width = 12, height = 8, dpi = 150) {
  path <- comparison_plot_path(compare_name, section, filename)
  ggsave(path, p, width = width, height = height, dpi = dpi, bg = "white")
  invisible(path)
}

write_csv_file <- function(compare_name, section, filename, x) {
  path <- comparison_result_path(compare_name, section, filename)
  readr::write_csv(x, path)
  invisible(path)
}

result_file <- function(model_name, section, filename) {
  file.path(result_root, model_name, section, filename)
}

detect_model_type <- function(model_name) {
  if (file.exists(result_file(model_name, "forecast", "predictions.csv")) &&
      file.exists(result_file(model_name, "key_result", "coefficients.csv"))) {
    return("harmonic")
  }
  if (file.exists(result_file(model_name, "forecast", "predictions.csv")) &&
      file.exists(result_file(model_name, "key_result", "selected_orders.csv"))) {
    return("sarima")
  }
  if (file.exists(result_file(model_name, "forecast", "train_predictions.csv")) &&
      file.exists(result_file(model_name, "posterior", "posterior_summary.csv"))) {
    return("bayesian")
  }
  usage_error("Could not detect model type for result directory: ", model_name)
}

load_prediction_data <- function(model_name, model_type) {
  if (identical(model_type, "bayesian")) {
    train_path <- result_file(model_name, "forecast", "train_predictions.csv")
    test_path <- result_file(model_name, "forecast", "test_predictions.csv")
    require_file(train_path, "Bayesian train forecast file")
    require_file(test_path, "Bayesian test forecast file")

    train_df <- readr::read_csv(train_path, show_col_types = FALSE) |>
      transmute(
        split,
        site,
        date = as.Date(date),
        actual = observed,
        pred_mean,
        pred_lo = pred_lo90,
        pred_hi = pred_hi90,
        pred_sd = (pred_hi90 - pred_lo90) / (2 * qnorm(0.95))
      )
    test_df <- readr::read_csv(test_path, show_col_types = FALSE) |>
      transmute(
        split,
        site,
        date = as.Date(date),
        actual = observed,
        pred_mean,
        pred_lo = pred_lo90,
        pred_hi = pred_hi90,
        pred_sd = (pred_hi90 - pred_lo90) / (2 * qnorm(0.95))
      )
    return(bind_rows(train_df, test_df))
  }

  path <- result_file(model_name, "forecast", "predictions.csv")
  require_file(path, "Prediction file")
  pred_df <- readr::read_csv(path, show_col_types = FALSE) |>
    mutate(date = as.Date(date))

  if (!("pred_lo" %in% names(pred_df)) && "lo90" %in% names(pred_df)) {
    pred_df <- rename(pred_df, pred_lo = lo90)
  }
  if (!("pred_hi" %in% names(pred_df)) && "hi90" %in% names(pred_df)) {
    pred_df <- rename(pred_df, pred_hi = hi90)
  }

  pred_df
}

load_residual_acf <- function(model_name) {
  path <- result_file(model_name, "time_series_diagnose", "residual_acf.csv")
  require_file(path, "Residual ACF file")
  readr::read_csv(path, show_col_types = FALSE)
}

compute_site_metrics <- function(pred_df) {
  test_df <- pred_df |>
    filter(split == "test") |>
    mutate(site = factor(site, levels = site_levels))

  test_df |>
    group_by(site) |>
    summarise(
      mae = mean(abs(actual - pred_mean), na.rm = TRUE),
      rmse = sqrt(mean((actual - pred_mean)^2, na.rm = TRUE)),
      lpd = mean(dnorm(actual, mean = pred_mean, sd = pmax(pred_sd, 1e-8), log = TRUE), na.rm = TRUE),
      .groups = "drop"
    )
}

extract_trace_df <- function(fit_obj, variables) {
  arr <- rstan::extract(fit_obj, pars = variables, permuted = FALSE, inc_warmup = FALSE)
  var_names <- dimnames(arr)[[3]]

  bind_rows(lapply(seq_along(var_names), function(k) {
    mat <- arr[, , k, drop = TRUE]
    if (is.null(dim(mat))) {
      mat <- matrix(mat, ncol = 1)
    }
    out <- as.data.frame(mat)
    names(out) <- paste0("chain_", seq_len(ncol(out)))
    out$iteration <- seq_len(nrow(out))
    out |>
      pivot_longer(cols = starts_with("chain_"), names_to = "chain", values_to = "value") |>
      mutate(variable = var_names[k])
  }))
}

choose_trace_variables <- function(summary_df) {
  preferred <- c("alpha", "beta_global", "q_global", "level_1", "growth_1", "sigma_level", "sigma_growth")
  vars <- intersect(preferred, summary_df$variable)
  if (!length(vars)) {
    vars <- summary_df |>
      filter(!str_detect(variable, "mu\\[|y_rep\\[|log_lik")) |>
      slice_head(n = min(6L, n())) |>
      pull(variable)
  }
  vars
}

cfg <- read_config(config_path)
model_entries <- lapply(cfg$models, function(x) {
  if (is.null(x$model_name)) {
    usage_error("Each model in comparison config must include `model_name`.")
  }
  model_name <- as.character(x$model_name)
  model_label <- as.character(x$label %||% model_name)
  model_type <- as.character(x$model_type %||% detect_model_type(model_name))
  list(model_name = model_name, model_label = model_label, model_type = model_type)
})

default_compare_name <- paste(vapply(model_entries, `[[`, character(1), "model_name"), collapse = "_")
compare_name <- as.character(cfg$comparison_name %||% default_compare_name)

loaded_models <- lapply(model_entries, function(entry) {
  pred_df <- load_prediction_data(entry$model_name, entry$model_type) |>
    mutate(
      model_name = entry$model_name,
      model = factor(entry$model_label, levels = vapply(model_entries, `[[`, character(1), "model_label")),
      site = factor(site, levels = site_levels)
    )

  acf_df <- load_residual_acf(entry$model_name) |>
    mutate(
      model_name = entry$model_name,
      model = factor(entry$model_label, levels = vapply(model_entries, `[[`, character(1), "model_label")),
      site = factor(site, levels = site_levels)
    )

  list(
    model_name = entry$model_name,
    model_label = entry$model_label,
    model_type = entry$model_type,
    predictions = pred_df,
    acf = acf_df
  )
})

forecast_compare <- bind_rows(lapply(loaded_models, `[[`, "predictions"))
acf_compare <- bind_rows(lapply(loaded_models, `[[`, "acf"))

train_forecast_plot <- ggplot(filter(forecast_compare, split == "train"), aes(date)) +
  geom_ribbon(aes(ymin = pred_lo, ymax = pred_hi), fill = "grey80", alpha = 0.5) +
  geom_line(aes(y = actual), color = "black", linewidth = 0.35) +
  geom_line(aes(y = pred_mean, color = model), linewidth = 0.7) +
  facet_grid(model ~ site, scales = "free_y") +
  labs(
    title = "Training Forecast Comparison",
    subtitle = "Rows = models, columns = BRW/MLO/SMO/SPO",
    x = NULL,
    y = expression(CO[2] ~ "(ppm)"),
    color = NULL
  ) +
  theme_project()

test_forecast_plot <- ggplot(filter(forecast_compare, split == "test"), aes(date)) +
  geom_ribbon(aes(ymin = pred_lo, ymax = pred_hi), fill = "grey80", alpha = 0.5) +
  geom_line(aes(y = actual), color = "black", linewidth = 0.35) +
  geom_line(aes(y = pred_mean, color = model), linewidth = 0.7) +
  facet_grid(model ~ site, scales = "free_y") +
  labs(
    title = "Holdout Forecast Comparison",
    subtitle = "Rows = models, columns = BRW/MLO/SMO/SPO",
    x = NULL,
    y = expression(CO[2] ~ "(ppm)"),
    color = NULL
  ) +
  theme_project()

site_metrics <- bind_rows(lapply(loaded_models, function(x) {
  compute_site_metrics(x$predictions) |>
    mutate(
      model_name = x$model_name,
      model = factor(x$model_label, levels = vapply(model_entries, `[[`, character(1), "model_label"))
    )
}))

metric_long <- site_metrics |>
  pivot_longer(cols = c(mae, rmse, lpd), names_to = "metric", values_to = "value") |>
  mutate(
    metric = factor(metric, levels = c("mae", "rmse", "lpd"), labels = c("MAE", "RMSE", "LPD")),
    site = factor(site, levels = site_levels)
  )

metric_compare_plot <- ggplot(metric_long, aes(model, value, fill = model)) +
  geom_col(width = 0.7) +
  facet_grid(metric ~ site, scales = "free_y") +
  labs(
    title = "Forecast Metric Comparison",
    subtitle = "Rows = MAE / RMSE / LPD, columns = BRW / MLO / SMO / SPO",
    x = NULL,
    y = "Metric value",
    fill = NULL
  ) +
  theme_project() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

residual_acf_plot <- ggplot(acf_compare, aes(lag, acf)) +
  geom_hline(yintercept = 0, color = "grey55", linewidth = 0.35) +
  geom_col(fill = "#4c78a8", width = 0.06) +
  facet_grid(model ~ site, scales = "free_y") +
  labs(
    title = "Residual ACF Comparison",
    subtitle = "Rows = models, columns = BRW / MLO / SMO / SPO",
    x = "Lag",
    y = "ACF"
  ) +
  theme_project()

write_csv_file(compare_name, "forecast", "forecast_compare.csv", forecast_compare)
write_csv_file(compare_name, "forecast", "site_metric_compare.csv", site_metrics)
write_csv_file(compare_name, "time_series_diagnose", "residual_acf_compare.csv", acf_compare)

save_plot_file(compare_name, "forecast", "bundle_forecast_train_panels.png", train_forecast_plot, width = 14, height = max(6, 2.2 * length(model_entries)))
save_plot_file(compare_name, "forecast", "bundle_forecast_test_panels.png", test_forecast_plot, width = 14, height = max(6, 2.2 * length(model_entries)))
save_plot_file(compare_name, "forecast", "bundle_forecast_metric_compare.png", metric_compare_plot, width = 14, height = 9)
save_plot_file(compare_name, "time_series_diagnose", "bundle_residual_acf_panels.png", residual_acf_plot, width = 14, height = max(6, 2.2 * length(model_entries)))

bayesian_models <- Filter(function(x) identical(x$model_type, "bayesian"), model_entries)

if (length(bayesian_models) > 0L) {
  posterior_compare <- bind_rows(lapply(bayesian_models, function(entry) {
    path <- result_file(entry$model_name, "posterior", "posterior_summary.csv")
    require_file(path, "Posterior summary file")
    readr::read_csv(path, show_col_types = FALSE) |>
      transmute(
        model_name = entry$model_name,
        model = factor(entry$model_label, levels = vapply(bayesian_models, `[[`, character(1), "model_label")),
        variable,
        ess_bulk = if ("ess_bulk" %in% names(cur_data_all())) ess_bulk else .data[["n_eff"]],
        r_hat = if ("r_hat" %in% names(cur_data_all())) r_hat else .data[["Rhat"]]
      )
  }))

  rhat_ess_plot <- ggplot(posterior_compare, aes(ess_bulk, r_hat)) +
    geom_hline(yintercept = 1.01, linetype = "dashed", color = "#b22222", linewidth = 0.6) +
    geom_point(color = "#2f6c8f", alpha = 0.72, size = 1.8) +
    facet_wrap(~ model, nrow = 1) +
    scale_x_continuous(trans = "log10", breaks = c(10, 100, 400, 1000, 4000, 10000), labels = label_number()) +
    labs(
      title = "Bayesian Rhat vs ESS Comparison",
      subtitle = "One panel per Bayesian model",
      x = "ESS bulk",
      y = "R-hat"
    ) +
    theme_project()

  write_csv_file(compare_name, "posterior", "bundle_rhat_ess_compare.csv", posterior_compare)
  save_plot_file(compare_name, "posterior", "bundle_rhat_ess_panels.png", rhat_ess_plot, width = max(10, 4 * length(bayesian_models)), height = 5.5)

  trace_compare <- bind_rows(lapply(bayesian_models, function(entry) {
    summary_path <- result_file(entry$model_name, "posterior", "posterior_summary.csv")
    fit_path <- result_file(entry$model_name, "posterior", "stanfit.rds")
    require_file(summary_path, "Posterior summary file")
    require_file(fit_path, "Stan fit file")
    summary_df <- readr::read_csv(summary_path, show_col_types = FALSE)
    fit <- readRDS(fit_path)
    vars <- choose_trace_variables(summary_df)
    extract_trace_df(fit, vars) |>
      mutate(model = factor(entry$model_label, levels = vapply(bayesian_models, `[[`, character(1), "model_label")))
  }))

  write_csv_file(compare_name, "posterior", "bundle_trace_compare.csv", trace_compare)

  trace_plot <- ggplot(trace_compare, aes(iteration, value, color = chain)) +
    geom_line(alpha = 0.7, linewidth = 0.3) +
    facet_grid(model ~ variable, scales = "free_y") +
    labs(
      title = "Bayesian Trace Plot Comparison",
      subtitle = "Rows = Bayesian models; columns = selected parameters",
      x = "Iteration",
      y = "Value"
    ) +
    theme_project()

  save_plot_file(compare_name, "posterior", "bundle_trace_compare.png", trace_plot, width = 14, height = max(6, 2.6 * length(bayesian_models)))
}

message("Finished comparison. Outputs written under Result/comparison/", compare_name, " and Plot/comparison/", compare_name)
