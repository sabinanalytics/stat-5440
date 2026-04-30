#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(tibble)
  library(rstan)
  library(posterior)
})

args <- commandArgs(trailingOnly = TRUE)

find_script_dir <- function() {
  args_all <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args_all, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }

  frames <- sys.frames()
  for (fr in frames) {
    ofile <- tryCatch(fr$ofile, error = function(e) NULL)
    if (!is.null(ofile) && nzchar(ofile)) {
      return(dirname(normalizePath(ofile)))
    }
  }

  normalizePath(getwd())
}

code_dir <- dirname(find_script_dir())
data_dir <- file.path(code_dir, "Data", "NOAA_CO2_Monthly")
result_dir <- file.path(code_dir, "Result")
rstan_dir <- file.path(code_dir, "Rstan")
config_path <- if (length(args) >= 1L) args[1] else file.path(code_dir, "config", "bayesian_config.json")
config_dir <- dirname(normalizePath(config_path, winslash = "/", mustWork = FALSE))

dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

site_levels <- c("BRW", "MLO", "SMO", "SPO")
site_map <- c(BRW = "brw", MLO = "mlo", SMO = "smo", SPO = "spo")
start_date_default <- as.Date("2000-01-01")
noaa_base_url <- "https://gml.noaa.gov/aftp/data/trace_gases/co2/in-situ/surface/txt"

reserved_models <- c(
  BHRQ = file.path(rstan_dir, "BHRQ.stan"),
  BDLM = file.path(rstan_dir, "BDLM.stan"),
  BDLMQ = file.path(rstan_dir, "BDLMQ.stan"),
  BDLM_AR2 = file.path(rstan_dir, "BDLM_AR2.stan")
)

usage_error <- function(...) {
  stop(..., call. = FALSE)
}

require_file <- function(path, label) {
  if (!file.exists(path)) {
    usage_error(label, " not found: ", normalizePath(path, winslash = "/", mustWork = FALSE))
  }
}

normalize_maybe <- function(path) {
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

resolve_existing_path <- function(path_string, base_dirs = c(".", config_dir, code_dir)) {
  candidates <- unique(c(path_string, file.path(base_dirs, path_string)))
  for (candidate in candidates) {
    if (file.exists(candidate) || dir.exists(candidate)) {
      return(normalizePath(candidate, winslash = "/", mustWork = FALSE))
    }
  }
  NULL
}

read_config <- function(path) {
  require_file(path, "Config file")
  cfg <- jsonlite::fromJSON(path, simplifyVector = TRUE)
  if (!is.list(cfg)) {
    usage_error("Config JSON must decode to an object.")
  }
  cfg
}

standardize_control <- function(control) {
  if (is.null(control)) {
    return(list(adapt_delta = 0.98, max_treedepth = 13))
  }
  if (!is.list(control)) {
    usage_error("`control` in config must be a JSON object.")
  }
  control
}

resolve_model_spec <- function(stan_model) {
  if (is.null(stan_model) || !is.character(stan_model) || length(stan_model) != 1L) {
    usage_error("`stan_model` must be a single string.")
  }

  model_key <- trimws(stan_model)
  if (model_key %in% names(reserved_models)) {
    path <- reserved_models[[model_key]]
    require_file(path, paste0("Reserved Stan model `", model_key, "`"))
    return(list(
      model_name = model_key,
      stan_path = path,
      source_type = "reserved"
    ))
  }

  custom_path <- resolve_existing_path(model_key)
  if (!is.null(custom_path) && file.exists(custom_path) && !dir.exists(custom_path)) {
    return(list(
      model_name = tools::file_path_sans_ext(basename(custom_path)),
      stan_path = custom_path,
      source_type = "custom"
    ))
  }

  usage_error(
    "`stan_model` must be one of ",
    paste(sprintf("`%s`", names(reserved_models)), collapse = ", "),
    ", or an existing `.stan` file path."
  )
}

build_sampling_config <- function(cfg) {
  list(
    seed = as.integer(cfg$seed %||% 5440L),
    chains = as.integer(cfg$chains %||% 4L),
    cores = as.integer(cfg$cores %||% min(4L, getOption("mc.cores", 1L))),
    iter = as.integer(cfg$iter %||% 2000L),
    warmup = as.integer(cfg$warmup %||% 1000L),
    refresh = as.integer(cfg$refresh %||% 200L),
    control = standardize_control(cfg$control)
  )
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

site_filename <- function(site_code) {
  slug <- unname(site_map[[site_code]])
  sprintf("co2_%s_surface-insitu_1_ccgg_MonthlyData.txt", slug)
}

download_site_file <- function(site_code, overwrite = FALSE) {
  filename <- site_filename(site_code)
  target <- file.path(data_dir, filename)
  if (file.exists(target) && !overwrite) {
    return(target)
  }

  url <- sprintf("%s/%s", noaa_base_url, filename)
  download.file(url = url, destfile = target, mode = "wb", quiet = FALSE)
  target
}

validate_raw_noaa_file <- function(path) {
  cols <- c(
    "site_code", "year", "month", "day", "hour", "minute", "second",
    "datetime", "time_decimal", "midpoint_time",
    "value", "value_std_dev", "nvalue",
    "latitude", "longitude", "altitude", "elevation", "intake_height",
    "qcflag"
  )

  col_types_spec <- cols(
    site_code = col_character(),
    year = col_integer(),
    month = col_integer(),
    day = col_integer(),
    hour = col_integer(),
    minute = col_integer(),
    second = col_integer(),
    datetime = col_character(),
    time_decimal = col_double(),
    midpoint_time = col_character(),
    value = col_double(),
    value_std_dev = col_double(),
    nvalue = col_integer(),
    latitude = col_double(),
    longitude = col_double(),
    altitude = col_double(),
    elevation = col_double(),
    intake_height = col_double(),
    qcflag = col_character()
  )

  read_table(
    file = path,
    comment = "#",
    col_names = cols,
    col_types = col_types_spec,
    na = c("-999.99", "-999.999", "-99.99", "-9", "NA"),
    show_col_types = FALSE,
    progress = FALSE
  )
}

read_one_noaa_site <- function(site_code, path) {
  dat <- tryCatch(
    validate_raw_noaa_file(path),
    error = function(e) {
      usage_error(
        "Failed to read NOAA monthly file for site `", site_code, "` at ",
        normalize_maybe(path), ".\n",
        "Expected the original NOAA MonthlyData.txt fixed-width/table format.\n",
        "Suggestion: use the raw NOAA file directly, or provide a combined table with columns `site`, `date`, `value`."
      )
    }
  )

  dat |>
    mutate(
      site = site_code,
      date = make_date(year, month, 1),
      reject_flag = substr(qcflag, 1, 1),
      keep = is.na(reject_flag) | reject_flag == "."
    ) |>
    filter(keep, !is.na(value)) |>
    select(site, date, year, month, value, value_std_dev, nvalue, qcflag)
}

validate_combined_panel <- function(dat, label) {
  required_cols <- c("site", "date", "value")
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0L) {
    usage_error(
      "Custom data file ", normalize_maybe(label), " is missing required columns: ",
      paste(missing_cols, collapse = ", "), ".\n",
      "Suggestion: provide columns `site`, `date`, `value`, or `site`, `year`, `month`, `value`."
    )
  }

  dat <- as_tibble(dat)
  dat$site <- toupper(trimws(as.character(dat$site)))

  if (!inherits(dat$date, "Date")) {
    dat$date <- suppressWarnings(as.Date(dat$date))
  }

  if (anyNA(dat$date)) {
    usage_error(
      "Custom data file ", normalize_maybe(label), " has invalid `date` values.\n",
      "Suggestion: use ISO format like `2000-01-01`, or provide `year` and `month` columns instead."
    )
  }

  if (!all(dat$site %in% site_levels)) {
    usage_error(
      "Custom data file ", normalize_maybe(label), " contains unsupported site codes.\n",
      "Supported sites are: ", paste(site_levels, collapse = ", "), "."
    )
  }

  if (!is.numeric(dat$value)) {
    dat$value <- suppressWarnings(as.numeric(dat$value))
  }
  if (anyNA(dat$value)) {
    usage_error(
      "Custom data file ", normalize_maybe(label), " has non-numeric `value` entries.\n",
      "Suggestion: make `value` a numeric monthly CO2 concentration column."
    )
  }

  if (any(duplicated(dat[, c("site", "date")]))) {
    usage_error(
      "Custom data file ", normalize_maybe(label), " contains duplicate `site` + `date` rows.\n",
      "Suggestion: keep exactly one monthly observation per site and month."
    )
  }

  dat |>
    mutate(
      year = year(date),
      month = month(date),
      value_std_dev = if ("value_std_dev" %in% names(dat)) value_std_dev else NA_real_,
      nvalue = if ("nvalue" %in% names(dat)) nvalue else NA_integer_,
      qcflag = if ("qcflag" %in% names(dat)) qcflag else NA_character_
    ) |>
    select(site, date, year, month, value, value_std_dev, nvalue, qcflag) |>
    arrange(date, site)
}

read_combined_data_file <- function(path) {
  ext <- tolower(tools::file_ext(path))
  dat <- switch(
    ext,
    csv = readr::read_csv(path, show_col_types = FALSE, progress = FALSE),
    tsv = readr::read_tsv(path, show_col_types = FALSE, progress = FALSE),
    txt = readr::read_delim(path, delim = "\t", show_col_types = FALSE, progress = FALSE),
    rds = readRDS(path),
    usage_error(
      "Unsupported custom data file extension for ", normalize_maybe(path), ".\n",
      "Supported file types: `.csv`, `.tsv`, `.txt`, `.rds`, or a directory containing NOAA raw files."
    )
  )

  if (!("date" %in% names(dat)) && all(c("year", "month") %in% names(dat))) {
    dat <- mutate(dat, date = make_date(as.integer(year), as.integer(month), 1))
  }

  validate_combined_panel(dat, path)
}

resolve_data_files <- function(mode) {
  mode_clean <- trimws(as.character(mode))
  lower_mode <- tolower(mode_clean)

  if (lower_mode == "force download") {
    files <- vapply(site_levels, download_site_file, character(1), overwrite = TRUE)
    return(list(type = "noaa_files", data = unname(files), description = "force download"))
  }

  if (lower_mode == "default") {
    files <- file.path(data_dir, vapply(site_levels, site_filename, character(1)))
    missing <- !file.exists(files)
    if (any(missing)) {
      for (site_code in site_levels[missing]) {
        download_site_file(site_code, overwrite = FALSE)
      }
    }
    return(list(type = "noaa_files", data = files, description = "default"))
  }

  resolved_path <- resolve_existing_path(mode_clean)

  if (!is.null(resolved_path) && dir.exists(resolved_path)) {
    files <- file.path(resolved_path, vapply(site_levels, site_filename, character(1)))
    missing <- files[!file.exists(files)]
    if (length(missing) > 0L) {
      usage_error(
        "Custom data directory is missing required NOAA monthly files:\n",
        paste(normalize_maybe(missing), collapse = "\n")
      )
    }
    return(list(type = "noaa_files", data = files, description = normalize_maybe(resolved_path)))
  }

  if (!is.null(resolved_path) && file.exists(resolved_path)) {
    return(list(type = "combined_file", data = resolved_path, description = normalize_maybe(resolved_path)))
  }

  usage_error(
    "`data` must be `default`, `force download`, an existing data directory, or an existing data file path."
  )
}

prepare_panel_data <- function(data_spec, start_date = start_date_default) {
  raw_dat <- NULL

  if (identical(data_spec$type, "noaa_files")) {
    file_lookup <- setNames(as.list(data_spec$data), site_levels)
    raw_dat <- bind_rows(lapply(site_levels, function(site_code) {
      read_one_noaa_site(site_code, file_lookup[[site_code]])
    }))
  } else if (identical(data_spec$type, "combined_file")) {
    raw_dat <- read_combined_data_file(data_spec$data)
  } else {
    usage_error("Internal error: unsupported data specification type.")
  }

  dat4 <- raw_dat |>
    arrange(site, date)

  common_dates <- dat4 |>
    count(date) |>
    filter(n == length(site_levels)) |>
    pull(date)

  dat4 <- dat4 |>
    filter(date %in% common_dates, date >= start_date) |>
    arrange(date, site) |>
    mutate(site = factor(site, levels = site_levels))

  if (nrow(dat4) == 0L) {
    usage_error(
      "No usable aligned monthly panel remained after filtering.\n",
      "Suggestion: ensure all four sites have overlapping monthly data and dates on or after ",
      format(start_date), "."
    )
  }

  counts <- dat4 |>
    count(date)
  if (any(counts$n != length(site_levels))) {
    usage_error("Aligned panel validation failed: some dates do not contain all four sites.")
  }

  dat4
}

make_train_test <- function(dat, test_months = 24L) {
  dates_all <- sort(unique(dat$date))
  if (length(dates_all) <= test_months) {
    usage_error(
      "Not enough months in the aligned panel. Need more than `test_months` = ", test_months, "."
    )
  }
  test_dates <- tail(dates_all, test_months)
  train_dates <- dates_all[!dates_all %in% test_dates]
  list(
    train = dat |> filter(date %in% train_dates) |> arrange(date, site),
    test = dat |> filter(date %in% test_dates) |> arrange(date, site),
    train_dates = train_dates,
    test_dates = test_dates
  )
}

build_rows <- function(df, date_to_t) {
  df |>
    mutate(
      site_id = match(as.character(site), site_levels),
      t_id = unname(date_to_t[as.character(date)]),
      month_id = month(date)
    ) |>
    arrange(date, site)
}

build_stan_data <- function(model_name, tr, te, all_dates) {
  base <- list(
    N = nrow(tr),
    S = length(site_levels),
    T = length(all_dates),
    site = tr$site_id,
    t_id = tr$t_id,
    month_id = tr$month_id,
    y = tr$value,
    N_test = nrow(te),
    site_test = te$site_id,
    t_test = te$t_id,
    month_test = te$month_id
  )

  if (identical(model_name, "BHRQ")) {
    return(base)
  }

  base$y_test <- te$value
  base
}

extract_matrix <- function(fit, variable) {
  rstan::extract(fit, pars = variable, permuted = FALSE, inc_warmup = FALSE) |>
    posterior::as_draws_array() |>
    posterior::as_draws_matrix()
}

safe_extract_matrix <- function(fit, variable) {
  out <- tryCatch(extract_matrix(fit, variable), error = function(e) NULL)
  out
}

rmse <- function(y, y_hat) sqrt(mean((y - y_hat)^2, na.rm = TRUE))
mae <- function(y, y_hat) mean(abs(y - y_hat), na.rm = TRUE)
coverage <- function(y, y_lo, y_hi) mean(y >= y_lo & y <= y_hi, na.rm = TRUE)

.logsumexp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

mlpd_from_loglik_draws <- function(loglik_draws) {
  s <- nrow(loglik_draws)
  mean(apply(loglik_draws, 2, function(x) .logsumexp(x) - log(s)))
}

fit_summary_table <- function(fit) {
  smry <- rstan::summary(fit)$summary
  out <- as.data.frame(smry, check.names = FALSE)
  out$variable <- rownames(out)
  rownames(out) <- NULL
  tibble::as_tibble(out) |>
    relocate(variable)
}

sampler_diagnostics_table <- function(fit) {
  params <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  treedepth_limit <- as.numeric(sampling_cfg$control$max_treedepth %||% NA_real_)
  bind_rows(lapply(seq_along(params), function(idx) {
    p <- as.data.frame(params[[idx]])
    tibble(
      chain = idx,
      iterations = nrow(p),
      divergences = sum(p$divergent__, na.rm = TRUE),
      max_treedepth_hits = if (is.na(treedepth_limit)) NA_real_ else sum(p$treedepth__ >= treedepth_limit, na.rm = TRUE),
      mean_stepsize = mean(p$stepsize__, na.rm = TRUE),
      mean_accept_stat = mean(p$accept_stat__, na.rm = TRUE)
    )
  }))
}

posterior_prediction_table <- function(df, mu_draws, y_rep_draws, split_name) {
  tibble(
    split = split_name,
    site = as.character(df$site),
    date = as.character(df$date),
    observed = df$value,
    pred_mean = colMeans(mu_draws),
    pred_lo90 = apply(y_rep_draws, 2, quantile, probs = 0.05),
    pred_hi90 = apply(y_rep_draws, 2, quantile, probs = 0.95)
  )
}

save_json_pretty <- function(x, path) {
  writeLines(jsonlite::toJSON(x, auto_unbox = TRUE, pretty = TRUE, null = "null"), con = path)
}

shared_data_cache_paths <- function() {
  list(
    metadata = file.path(data_dir, "prepared_data_metadata.json"),
    aligned = file.path(data_dir, "aligned_panel_data.csv"),
    train = file.path(data_dir, "train_data.csv"),
    test = file.path(data_dir, "test_data.csv")
  )
}

load_prepared_data_cache <- function(start_date, test_months) {
  paths <- shared_data_cache_paths()
  required <- unlist(paths, use.names = FALSE)
  if (!all(file.exists(required))) {
    return(NULL)
  }

  meta <- jsonlite::fromJSON(paths$metadata, simplifyVector = TRUE)
  if (!identical(as.character(meta$start_date %||% ""), format(start_date))) {
    return(NULL)
  }
  if (!identical(as.integer(meta$test_months %||% NA_integer_), as.integer(test_months))) {
    return(NULL)
  }

  list(
    aligned = readr::read_csv(paths$aligned, show_col_types = FALSE) |>
      mutate(date = as.Date(date), site = factor(site, levels = site_levels)),
    train = readr::read_csv(paths$train, show_col_types = FALSE) |>
      mutate(date = as.Date(date), site = factor(site, levels = site_levels)),
    test = readr::read_csv(paths$test, show_col_types = FALSE) |>
      mutate(date = as.Date(date), site = factor(site, levels = site_levels))
  )
}

save_prepared_data_cache <- function(dat4, train, test, start_date, test_months) {
  paths <- shared_data_cache_paths()
  readr::write_csv(dat4, paths$aligned)
  readr::write_csv(train, paths$train)
  readr::write_csv(test, paths$test)
  save_json_pretty(
    list(
      start_date = format(start_date),
      test_months = as.integer(test_months),
      n_aligned = nrow(dat4),
      n_train = nrow(train),
      n_test = nrow(test)
    ),
    paths$metadata
  )
}

cfg <- read_config(config_path)
model_spec <- resolve_model_spec(cfg$stan_model)
sampling_cfg <- build_sampling_config(cfg)
test_months <- as.integer(cfg$test_months %||% 24L)
start_date <- as.Date(cfg$start_date %||% format(start_date_default))
data_mode <- cfg$data %||% "default"
data_spec <- resolve_data_files(data_mode)

model_output_dir <- file.path(result_dir, model_spec$model_name)
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
  stan_model = cfg$stan_model,
  resolved_model_name = model_spec$model_name,
  resolved_stan_path = normalize_maybe(model_spec$stan_path),
  data = data_mode,
  resolved_data = data_spec$description,
  seed = sampling_cfg$seed,
  chains = sampling_cfg$chains,
  cores = sampling_cfg$cores,
  iter = sampling_cfg$iter,
  warmup = sampling_cfg$warmup,
  refresh = sampling_cfg$refresh,
  control = sampling_cfg$control,
  test_months = test_months,
  start_date = format(start_date)
)
save_json_pretty(resolved_config, file.path(model_output_dir, "resolved_config.json"))

cached_data <- NULL
if (tolower(data_mode) == "default") {
  cached_data <- load_prepared_data_cache(start_date = start_date, test_months = test_months)
}

if (is.null(cached_data)) {
  panel_data <- prepare_panel_data(data_spec = data_spec, start_date = start_date)
  split <- make_train_test(panel_data, test_months = test_months)
  save_prepared_data_cache(panel_data, split$train, split$test, start_date = start_date, test_months = test_months)
  train_data <- split$train
  test_data <- split$test
} else {
  panel_data <- cached_data$aligned
  train_data <- cached_data$train
  test_data <- cached_data$test
}

all_dates <- sort(unique(panel_data$date))
date_to_t <- setNames(seq_along(all_dates), as.character(all_dates))
tr <- build_rows(train_data, date_to_t)
te <- build_rows(test_data, date_to_t)
stan_data <- build_stan_data(model_spec$model_name, tr, te, all_dates)

saveRDS(stan_data, result_path("posterior", "stan_data.rds"))

rstan::rstan_options(auto_write = TRUE)
options(mc.cores = max(1L, sampling_cfg$cores))

message("Compiling Stan model: ", model_spec$model_name)
sm <- rstan::stan_model(
  file = model_spec$stan_path,
  model_name = model_spec$model_name
)

message("Running sampling...")
fit <- rstan::sampling(
  object = sm,
  data = stan_data,
  seed = sampling_cfg$seed,
  chains = sampling_cfg$chains,
  cores = sampling_cfg$cores,
  iter = sampling_cfg$iter,
  warmup = sampling_cfg$warmup,
  refresh = sampling_cfg$refresh,
  control = sampling_cfg$control
)

saveRDS(fit, result_path("posterior", "stanfit.rds"))

summary_tbl <- fit_summary_table(fit)
readr::write_csv(summary_tbl, result_path("posterior", "posterior_summary.csv"))
readr::write_csv(
  sampler_diagnostics_table(fit),
  result_path("posterior", "sampler_diagnostics.csv")
)

mu_train_draws <- safe_extract_matrix(fit, "mu")
y_train_rep_draws <- safe_extract_matrix(fit, "y_rep")
mu_test_draws <- safe_extract_matrix(fit, "mu_test")
y_test_rep_draws <- safe_extract_matrix(fit, "y_test_rep")
log_lik_train_draws <- safe_extract_matrix(fit, "log_lik")
log_lik_test_draws <- safe_extract_matrix(fit, "log_lik_test")

if (is.null(mu_train_draws) || is.null(y_train_rep_draws) || is.null(mu_test_draws) || is.null(y_test_rep_draws)) {
  usage_error(
    "The Stan model did not expose expected generated quantities (`mu`, `y_rep`, `mu_test`, `y_test_rep`).\n",
    "Please update the Stan program or use one of the reserved models."
  )
}

train_pred <- posterior_prediction_table(tr, mu_train_draws, y_train_rep_draws, "train")
test_pred <- posterior_prediction_table(te, mu_test_draws, y_test_rep_draws, "test")
acf_df <- bind_rows(train_pred, test_pred) |>
  filter(split == "train") |>
  mutate(residual = observed - pred_mean) |>
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

readr::write_csv(train_pred, result_path("forecast", "train_predictions.csv"))
readr::write_csv(test_pred, result_path("forecast", "test_predictions.csv"))
readr::write_csv(acf_df, result_path("time_series_diagnose", "residual_acf.csv"))

metrics_tbl <- tibble(
  model = model_spec$model_name,
  train_rmse = rmse(tr$value, train_pred$pred_mean),
  train_mae = mae(tr$value, train_pred$pred_mean),
  train_coverage90 = coverage(tr$value, train_pred$pred_lo90, train_pred$pred_hi90),
  test_rmse = rmse(te$value, test_pred$pred_mean),
  test_mae = mae(te$value, test_pred$pred_mean),
  test_coverage90 = coverage(te$value, test_pred$pred_lo90, test_pred$pred_hi90),
  train_mlpd = if (is.null(log_lik_train_draws)) NA_real_ else mlpd_from_loglik_draws(log_lik_train_draws),
  test_mlpd = if (is.null(log_lik_test_draws)) NA_real_ else mlpd_from_loglik_draws(log_lik_test_draws),
  n_train = nrow(tr),
  n_test = nrow(te),
  months_total = length(all_dates)
)
readr::write_csv(metrics_tbl, result_path("forecast", "metrics.csv"))

message("Finished. Results written to: ", normalize_maybe(model_output_dir))
