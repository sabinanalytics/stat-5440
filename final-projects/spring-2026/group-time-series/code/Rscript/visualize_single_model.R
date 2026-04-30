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
result_root <- file.path(code_dir, "Result")
plot_root <- file.path(code_dir, "Plot")
dir.create(plot_root, recursive = TRUE, showWarnings = FALSE)

if (length(args) < 1L) {
  usage_error("Usage: Rscript visualize_single_model.R <model_name_or_result_path>")
}

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

resolve_model_result_dir <- function(input) {
  candidate <- trimws(input)
  direct <- resolve_existing_path(candidate, base_dirs = c(".", code_dir, result_root))
  if (!is.null(direct) && dir.exists(direct)) {
    return(direct)
  }

  by_name <- file.path(result_root, candidate)
  if (dir.exists(by_name)) {
    return(normalizePath(by_name, winslash = "/", mustWork = FALSE))
  }

  usage_error(
    "Could not resolve model result directory from input `", input, "`.\n",
    "Expected either a model name under ", normalize_maybe(result_root), " or a direct result path."
  )
}

require_result_files <- function(result_dir, files) {
  missing <- files[!file.exists(file.path(result_dir, files))]
  if (length(missing) > 0L) {
    usage_error(
      "Result directory ", normalize_maybe(result_dir), " is missing required files:\n",
      paste(missing, collapse = "\n")
    )
  }
}

result_file <- function(result_dir, section, filename) {
  file.path(result_dir, section, filename)
}

require_section_files <- function(result_dir, section, files) {
  missing <- files[!file.exists(result_file(result_dir, section, files))]
  if (length(missing) > 0L) {
    usage_error(
      "Result directory ", normalize_maybe(result_dir), " is missing required files under `", section, "`:\n",
      paste(missing, collapse = "\n")
    )
  }
}

save_plot_file <- function(plot_dir, filename, p, width = 11, height = 6.5, dpi = 150) {
  path <- file.path(plot_dir, filename)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path, p, width = width, height = height, dpi = dpi, bg = "white")
  invisible(path)
}

compute_acf_by_site <- function(pred_df) {
  pred_df |>
    filter(split == "train") |>
    mutate(residual = actual - pred_mean) |>
    group_by(site) |>
    group_modify(~ {
      x <- .x$residual
      x <- x[is.finite(x)]
      if (length(x) < 3L) {
        return(tibble(lag = numeric(0), acf = numeric(0), threshold = numeric(0)))
      }
      acf_obj <- stats::acf(x, plot = FALSE, na.action = na.pass, lag.max = min(24L, length(x) - 1L))
      threshold <- 1.96 / sqrt(length(x))
      tibble(
        lag = as.numeric(acf_obj$lag)[-1],
        acf = as.numeric(acf_obj$acf)[-1],
        threshold = threshold
      )
    }) |>
    ungroup()
}

plot_forecast_split <- function(df, split_name, model_name) {
  df_split <- df |>
    filter(split == split_name) |>
    mutate(site = factor(site, levels = site_levels))

  title_text <- if (split_name == "train") "Training Fit" else "Holdout Forecast"
  ggplot(df_split, aes(x = as.Date(date))) +
    geom_ribbon(aes(ymin = pred_lo, ymax = pred_hi), fill = "grey80", alpha = 0.6) +
    geom_line(aes(y = actual), color = "black", linewidth = 0.35) +
    geom_line(aes(y = pred_mean, color = site), linewidth = 0.8) +
    scale_color_manual(values = site_colors, guide = "none") +
    facet_wrap(~ site, ncol = 2, scales = "free_y") +
    labs(
      title = paste(model_name, "-", title_text),
      subtitle = "Black = observed, color = fitted/forecast mean, grey band = predictive interval",
      x = NULL,
      y = expression(CO[2] ~ "(ppm)")
    ) +
    theme_project()
}

plot_residual_acf <- function(acf_df, model_name) {
  threshold_df <- acf_df |>
    group_by(site) |>
    summarise(threshold = first(threshold), .groups = "drop")

  ggplot(acf_df, aes(lag, acf)) +
    geom_hline(yintercept = 0, color = "grey50") +
    geom_hline(data = threshold_df, aes(yintercept = threshold), color = "red", linetype = "dashed") +
    geom_hline(data = threshold_df, aes(yintercept = -threshold), color = "red", linetype = "dashed") +
    geom_col(fill = "#4c78a8", width = 0.06) +
    facet_wrap(~ site, ncol = 2, scales = "free_y") +
    labs(
      title = paste(model_name, "- Residual ACF"),
      subtitle = "Training residual autocorrelation by site; red dashed lines show approximate significance thresholds",
      x = "Lag",
      y = "ACF"
    ) +
    theme_project()
}

select_key_posterior <- function(summary_df) {
  variable_patterns <- c("^alpha$", "^beta_global$", "^q_global$", "^level_1$", "^growth_1$", "^sigma_")
  matched <- Reduce(`|`, lapply(variable_patterns, function(p) grepl(p, summary_df$variable)))
  out <- summary_df[matched, , drop = FALSE]
  if (nrow(out) == 0L) {
    out <- summary_df |>
      filter(!str_detect(variable, "log_lik|y_rep|mu\\[")) |>
      slice_head(n = min(15L, n()))
  }
  out |>
    mutate(variable = factor(variable, levels = rev(unique(variable))))
}

plot_key_results <- function(result_dir, model_name, model_type) {
  if (model_type == "harmonic") {
    require_section_files(result_dir, "key_result", c("slopes.csv", "seasonal_summary.csv"))
    slopes <- readr::read_csv(result_file(result_dir, "key_result", "slopes.csv"), show_col_types = FALSE)
    seasonal <- readr::read_csv(result_file(result_dir, "key_result", "seasonal_summary.csv"), show_col_types = FALSE)

    p1 <- ggplot(slopes, aes(site, slope_ppm_per_year, fill = site)) +
      geom_col(width = 0.65) +
      scale_fill_manual(values = site_colors, guide = "none") +
      labs(title = paste(model_name, "- Key Results"), subtitle = "Per-site linear trend slope", x = NULL, y = "ppm / year") +
      theme_project()

    seasonal_long <- seasonal |>
      select(site, amp_h1, amp_h2) |>
      pivot_longer(cols = c(amp_h1, amp_h2), names_to = "harmonic", values_to = "amplitude")

    p2 <- ggplot(seasonal_long, aes(site, amplitude, fill = harmonic)) +
      geom_col(position = "dodge", width = 0.7) +
      labs(title = paste(model_name, "- Seasonal Amplitudes"), x = NULL, y = "Amplitude") +
      theme_project()

    return(list(key_results_main = p1, key_results_secondary = p2))
  }

  if (model_type == "sarima") {
    require_section_files(result_dir, "key_result", "selected_orders.csv")
    orders <- readr::read_csv(result_file(result_dir, "key_result", "selected_orders.csv"), show_col_types = FALSE)
    orders <- mutate(orders, label = sprintf("(%d,%d,%d)(%d,%d,%d)[%d]", p, d, q, P, D, Q, period))

    p <- ggplot(orders, aes(site, aicc, fill = site)) +
      geom_col(width = 0.65) +
      geom_text(aes(label = label), vjust = -0.4, size = 3) +
      scale_fill_manual(values = site_colors, guide = "none") +
      labs(title = paste(model_name, "- Key Results"), subtitle = "Selected per-site SARIMA orders and AICc", x = NULL, y = "AICc") +
      theme_project()

    return(list(key_results_main = p))
  }

  require_section_files(result_dir, "posterior", "posterior_summary.csv")
  summary_df <- readr::read_csv(result_file(result_dir, "posterior", "posterior_summary.csv"), show_col_types = FALSE)
  key_df <- select_key_posterior(summary_df)

  q_low <- names(key_df)[str_detect(names(key_df), "^5(\\.|\\d)*5?%$|^2\\.5%$|^25%$|^10%$")]
  q_high <- names(key_df)[str_detect(names(key_df), "^94(\\.|\\d)*5?%$|^97\\.5%$|^75%$|^90%$")]
  if (!length(q_low) || !length(q_high)) {
    q_low <- "mean"
    q_high <- "mean"
  } else {
    q_low <- q_low[1]
    q_high <- q_high[1]
  }

  p <- ggplot(key_df, aes(x = mean, y = variable)) +
    geom_errorbarh(aes(xmin = .data[[q_low]], xmax = .data[[q_high]]), height = 0.2, color = "grey50") +
    geom_point(color = "#4c78a8", size = 2) +
    labs(title = paste(model_name, "- Key Posterior Results"), subtitle = "Posterior means with uncertainty intervals", x = "Estimate", y = NULL) +
    theme_project()

  list(key_results_main = p)
}

plot_rhat_ess <- function(result_dir, model_name) {
  require_section_files(result_dir, "posterior", "posterior_summary.csv")
  summary_df <- readr::read_csv(result_file(result_dir, "posterior", "posterior_summary.csv"), show_col_types = FALSE)
  ess_col <- names(summary_df)[str_detect(names(summary_df), "n_eff|ess_bulk")]
  rhat_col <- names(summary_df)[str_detect(names(summary_df), "Rhat|r_hat")]
  if (!length(ess_col) || !length(rhat_col)) {
    usage_error("posterior_summary.csv does not contain ESS/Rhat columns needed for convergence plots.")
  }
  ess_col <- ess_col[1]
  rhat_col <- rhat_col[1]

  plot_df <- summary_df |>
    filter(is.finite(.data[[ess_col]]), is.finite(.data[[rhat_col]]))

  ggplot(plot_df, aes(x = .data[[ess_col]], y = .data[[rhat_col]])) +
    geom_point(alpha = 0.7, color = "#e45756") +
    geom_hline(yintercept = 1.01, linetype = "dashed", color = "grey40") +
    scale_x_log10() +
    labs(title = paste(model_name, "- Rhat vs ESS"), subtitle = "Convergence diagnostics across posterior parameters", x = "Effective sample size", y = "Rhat") +
    theme_project()
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

plot_traceplot <- function(result_dir, model_name) {
  require_section_files(result_dir, "posterior", c("stanfit.rds", "posterior_summary.csv"))
  fit <- readRDS(result_file(result_dir, "posterior", "stanfit.rds"))
  summary_df <- readr::read_csv(result_file(result_dir, "posterior", "posterior_summary.csv"), show_col_types = FALSE)
  vars <- choose_trace_variables(summary_df)
  trace_df <- extract_trace_df(fit, vars)

  ggplot(trace_df, aes(iteration, value, color = chain)) +
    geom_line(alpha = 0.7, linewidth = 0.3) +
    facet_wrap(~ variable, scales = "free_y", ncol = 2) +
    labs(title = paste(model_name, "- Traceplot"), subtitle = "Selected posterior parameters by chain", x = "Iteration", y = "Value") +
    theme_project()
}

classify_model <- function(result_dir) {
  if (all(file.exists(c(
    result_file(result_dir, "forecast", "predictions.csv"),
    result_file(result_dir, "key_result", "coefficients.csv"),
    result_file(result_dir, "key_result", "slopes.csv"),
    result_file(result_dir, "key_result", "seasonal_summary.csv")
  )))) {
    return("harmonic")
  }
  if (all(file.exists(c(
    result_file(result_dir, "forecast", "predictions.csv"),
    result_file(result_dir, "key_result", "selected_orders.csv")
  )))) {
    return("sarima")
  }
  if (all(file.exists(c(
    result_file(result_dir, "forecast", "train_predictions.csv"),
    result_file(result_dir, "forecast", "test_predictions.csv"),
    result_file(result_dir, "posterior", "posterior_summary.csv")
  )))) {
    return("bayesian")
  }
  usage_error(
    "Could not classify result directory ", normalize_maybe(result_dir), ".\n",
    "Expected Harmonic, SARIMA, or Bayesian unified output files."
  )
}

load_prediction_frame <- function(result_dir, model_type) {
  if (model_type == "bayesian") {
    require_section_files(result_dir, "forecast", c("train_predictions.csv", "test_predictions.csv"))
    train_pred <- readr::read_csv(result_file(result_dir, "forecast", "train_predictions.csv"), show_col_types = FALSE) |>
      transmute(split, site, date = as.Date(date), actual = observed, pred_mean, pred_lo = pred_lo90, pred_hi = pred_hi90)
    test_pred <- readr::read_csv(result_file(result_dir, "forecast", "test_predictions.csv"), show_col_types = FALSE) |>
      transmute(split, site, date = as.Date(date), actual = observed, pred_mean, pred_lo = pred_lo90, pred_hi = pred_hi90)
    return(bind_rows(train_pred, test_pred))
  }

  require_section_files(result_dir, "forecast", "predictions.csv")
  pred_df <- readr::read_csv(result_file(result_dir, "forecast", "predictions.csv"), show_col_types = FALSE) |>
    mutate(date = as.Date(date))

  if (!("pred_lo" %in% names(pred_df)) && "lo90" %in% names(pred_df)) {
    pred_df <- dplyr::rename(pred_df, pred_lo = lo90)
  }
  if (!("pred_hi" %in% names(pred_df)) && "hi90" %in% names(pred_df)) {
    pred_df <- dplyr::rename(pred_df, pred_hi = hi90)
  }

  pred_df
}

result_dir <- resolve_model_result_dir(args[1])
model_name <- basename(result_dir)
plot_dir <- file.path(plot_root, model_name)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

model_type <- classify_model(result_dir)
pred_df <- load_prediction_frame(result_dir, model_type)

required_pred_cols <- c("split", "site", "date", "actual", "pred_mean", "pred_lo", "pred_hi")
missing_pred_cols <- setdiff(required_pred_cols, names(pred_df))
if (length(missing_pred_cols) > 0L) {
  usage_error(
    "Prediction data in ", normalize_maybe(result_dir), " is missing columns: ",
    paste(missing_pred_cols, collapse = ", ")
  )
}

acf_df <- compute_acf_by_site(pred_df)
if (nrow(acf_df) == 0L) {
  usage_error("Could not compute training residual ACF from result files in ", normalize_maybe(result_dir), ".")
}

save_plot_file(plot_dir, file.path("forecast", "forecast_train.png"), plot_forecast_split(pred_df, "train", model_name))
save_plot_file(plot_dir, file.path("forecast", "forecast_test.png"), plot_forecast_split(pred_df, "test", model_name))
save_plot_file(plot_dir, file.path("time_series_diagnose", "residual_acf.png"), plot_residual_acf(acf_df, model_name))

key_plots <- plot_key_results(result_dir, model_name, model_type)
for (plot_name in names(key_plots)) {
  save_plot_file(plot_dir, file.path("key_result", paste0(plot_name, ".png")), key_plots[[plot_name]], width = 10, height = 5.8)
}

if (model_type == "bayesian") {
  save_plot_file(plot_dir, file.path("posterior", "posterior_rhat_ess.png"), plot_rhat_ess(result_dir, model_name), width = 9, height = 6)
  save_plot_file(plot_dir, file.path("posterior", "posterior_traceplot.png"), plot_traceplot(result_dir, model_name), width = 11, height = 7)
}

message("Finished. Plots written to: ", normalize_maybe(plot_dir))
