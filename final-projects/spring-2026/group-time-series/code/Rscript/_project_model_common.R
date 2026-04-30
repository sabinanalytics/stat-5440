`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

find_code_dir <- function() {
  args_all <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args_all, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(dirname(normalizePath(sub("^--file=", "", file_arg[1])))))
  }

  frames <- sys.frames()
  for (fr in frames) {
    ofile <- tryCatch(fr$ofile, error = function(e) NULL)
    if (!is.null(ofile) && nzchar(ofile)) {
      return(dirname(dirname(normalizePath(ofile))))
    }
  }

  normalizePath(getwd())
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

read_config_json <- function(path) {
  require_file(path, "Config file")
  cfg <- jsonlite::fromJSON(path, simplifyVector = TRUE)
  if (!is.list(cfg)) {
    usage_error("Config JSON must decode to an object.")
  }
  cfg
}

site_levels_project <- c("BRW", "MLO", "SMO", "SPO")
site_map_project <- c(BRW = "brw", MLO = "mlo", SMO = "smo", SPO = "spo")
noaa_base_url_project <- "https://gml.noaa.gov/aftp/data/trace_gases/co2/in-situ/surface/txt"

shared_noaa_data_dir_project <- function(code_dir) {
  path <- file.path(code_dir, "Data", "NOAA_CO2_Monthly")
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

shared_data_cache_paths_project <- function(shared_data_dir) {
  list(
    metadata = file.path(shared_data_dir, "prepared_data_metadata.json"),
    aligned = file.path(shared_data_dir, "aligned_panel_data.csv"),
    train = file.path(shared_data_dir, "train_data.csv"),
    test = file.path(shared_data_dir, "test_data.csv")
  )
}

resolve_existing_path <- function(path_string, base_dirs) {
  candidates <- unique(c(path_string, file.path(base_dirs, path_string)))
  for (candidate in candidates) {
    if (file.exists(candidate) || dir.exists(candidate)) {
      return(normalizePath(candidate, winslash = "/", mustWork = FALSE))
    }
  }
  NULL
}

site_filename_project <- function(site_code) {
  slug <- unname(site_map_project[[site_code]])
  sprintf("co2_%s_surface-insitu_1_ccgg_MonthlyData.txt", slug)
}

download_site_file_project <- function(site_code, data_dir, overwrite = FALSE) {
  filename <- site_filename_project(site_code)
  target <- file.path(data_dir, filename)
  if (file.exists(target) && !overwrite) {
    return(target)
  }

  url <- sprintf("%s/%s", noaa_base_url_project, filename)
  download.file(url = url, destfile = target, mode = "wb", quiet = FALSE)
  target
}

validate_raw_noaa_file_project <- function(path) {
  cols <- c(
    "site_code", "year", "month", "day", "hour", "minute", "second",
    "datetime", "time_decimal", "midpoint_time",
    "value", "value_std_dev", "nvalue",
    "latitude", "longitude", "altitude", "elevation", "intake_height",
    "qcflag"
  )

  col_types_spec <- readr::cols(
    site_code = readr::col_character(),
    year = readr::col_integer(),
    month = readr::col_integer(),
    day = readr::col_integer(),
    hour = readr::col_integer(),
    minute = readr::col_integer(),
    second = readr::col_integer(),
    datetime = readr::col_character(),
    time_decimal = readr::col_double(),
    midpoint_time = readr::col_character(),
    value = readr::col_double(),
    value_std_dev = readr::col_double(),
    nvalue = readr::col_integer(),
    latitude = readr::col_double(),
    longitude = readr::col_double(),
    altitude = readr::col_double(),
    elevation = readr::col_double(),
    intake_height = readr::col_double(),
    qcflag = readr::col_character()
  )

  readr::read_table(
    file = path,
    comment = "#",
    col_names = cols,
    col_types = col_types_spec,
    na = c("-999.99", "-999.999", "-99.99", "-9", "NA"),
    show_col_types = FALSE,
    progress = FALSE
  )
}

read_one_noaa_site_project <- function(site_code, path) {
  dat <- tryCatch(
    validate_raw_noaa_file_project(path),
    error = function(e) {
      usage_error(
        "Failed to read NOAA monthly file for site `", site_code, "` at ",
        normalize_maybe(path), ".\n",
        "Expected the original NOAA MonthlyData.txt format.\n",
        "Suggestion: use the raw NOAA file directly, or provide a combined table with columns `site`, `date`, `value`."
      )
    }
  )

  dat |>
    dplyr::mutate(
      site = site_code,
      date = lubridate::make_date(year, month, 1),
      reject_flag = substr(qcflag, 1, 1),
      keep = is.na(reject_flag) | reject_flag == "."
    ) |>
    dplyr::filter(keep, !is.na(value)) |>
    dplyr::select(site, date, year, month, value, value_std_dev, nvalue, qcflag)
}

validate_combined_panel_project <- function(dat, label) {
  required_cols <- c("site", "date", "value")
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0L) {
    usage_error(
      "Custom data file ", normalize_maybe(label), " is missing required columns: ",
      paste(missing_cols, collapse = ", "), ".\n",
      "Suggestion: provide columns `site`, `date`, `value`, or `site`, `year`, `month`, `value`."
    )
  }

  dat <- tibble::as_tibble(dat)
  dat$site <- toupper(trimws(as.character(dat$site)))

  if (!inherits(dat$date, "Date")) {
    dat$date <- suppressWarnings(as.Date(dat$date))
  }
  if (anyNA(dat$date)) {
    usage_error(
      "Custom data file ", normalize_maybe(label), " has invalid `date` values.\n",
      "Suggestion: use ISO format like `2000-01-01`, or provide `year` and `month` columns."
    )
  }

  if (!all(dat$site %in% site_levels_project)) {
    usage_error(
      "Custom data file ", normalize_maybe(label), " contains unsupported site codes.\n",
      "Supported sites are: ", paste(site_levels_project, collapse = ", "), "."
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
    dplyr::mutate(
      year = lubridate::year(date),
      month = lubridate::month(date),
      value_std_dev = if ("value_std_dev" %in% names(dat)) value_std_dev else NA_real_,
      nvalue = if ("nvalue" %in% names(dat)) nvalue else NA_integer_,
      qcflag = if ("qcflag" %in% names(dat)) qcflag else NA_character_
    ) |>
    dplyr::select(site, date, year, month, value, value_std_dev, nvalue, qcflag) |>
    dplyr::arrange(date, site)
}

read_combined_data_file_project <- function(path) {
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
    dat <- dplyr::mutate(dat, date = lubridate::make_date(as.integer(year), as.integer(month), 1))
  }

  validate_combined_panel_project(dat, path)
}

resolve_data_spec_project <- function(data_mode, data_dir, base_dirs) {
  mode_clean <- trimws(as.character(data_mode))
  lower_mode <- tolower(mode_clean)

  if (lower_mode == "force download") {
    files <- vapply(site_levels_project, download_site_file_project, character(1), data_dir = data_dir, overwrite = TRUE)
    return(list(type = "noaa_files", data = unname(files), description = "force download"))
  }

  if (lower_mode == "default") {
    files <- file.path(data_dir, vapply(site_levels_project, site_filename_project, character(1)))
    missing <- !file.exists(files)
    if (any(missing)) {
      for (site_code in site_levels_project[missing]) {
        download_site_file_project(site_code, data_dir = data_dir, overwrite = FALSE)
      }
    }
    return(list(type = "noaa_files", data = files, description = "default"))
  }

  resolved_path <- resolve_existing_path(mode_clean, base_dirs = base_dirs)

  if (!is.null(resolved_path) && dir.exists(resolved_path)) {
    files <- file.path(resolved_path, vapply(site_levels_project, site_filename_project, character(1)))
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

  usage_error("`data` must be `default`, `force download`, an existing data directory, or an existing data file path.")
}

prepare_panel_data_project <- function(data_spec, start_date) {
  raw_dat <- NULL

  if (identical(data_spec$type, "noaa_files")) {
    file_lookup <- setNames(as.list(data_spec$data), site_levels_project)
    raw_dat <- dplyr::bind_rows(lapply(site_levels_project, function(site_code) {
      read_one_noaa_site_project(site_code, file_lookup[[site_code]])
    }))
  } else if (identical(data_spec$type, "combined_file")) {
    raw_dat <- read_combined_data_file_project(data_spec$data)
  } else {
    usage_error("Internal error: unsupported data specification type.")
  }

  dat4 <- raw_dat |>
    dplyr::arrange(site, date)

  common_dates <- dat4 |>
    dplyr::count(date) |>
    dplyr::filter(n == length(site_levels_project)) |>
    dplyr::pull(date)

  dat4 <- dat4 |>
    dplyr::filter(date %in% common_dates, date >= start_date) |>
    dplyr::arrange(date, site) |>
    dplyr::group_by(site) |>
    dplyr::mutate(
      t = dplyr::row_number(),
      sin1 = sin(2 * pi * t / 12),
      cos1 = cos(2 * pi * t / 12),
      sin2 = sin(4 * pi * t / 12),
      cos2 = cos(4 * pi * t / 12)
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(site = factor(site, levels = site_levels_project))

  if (nrow(dat4) == 0L) {
    usage_error(
      "No usable aligned monthly panel remained after filtering.\n",
      "Suggestion: ensure all four sites have overlapping monthly data and dates on or after ",
      format(start_date), "."
    )
  }

  dat4
}

make_train_test_project <- function(dat, test_months) {
  dates_all <- sort(unique(dat$date))
  if (length(dates_all) <= test_months) {
    usage_error("Not enough months in the aligned panel. Need more than `test_months` = ", test_months, ".")
  }

  test_dates <- tail(dates_all, test_months)
  train_dates <- dates_all[!dates_all %in% test_dates]
  list(
    train = dat |> dplyr::filter(date %in% train_dates) |> dplyr::arrange(date, site),
    test = dat |> dplyr::filter(date %in% test_dates) |> dplyr::arrange(date, site),
    train_dates = train_dates,
    test_dates = test_dates
  )
}

load_prepared_data_cache_project <- function(shared_data_dir, start_date, test_months) {
  paths <- shared_data_cache_paths_project(shared_data_dir)
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
    metadata = meta,
    aligned = readr::read_csv(paths$aligned, show_col_types = FALSE) |>
      dplyr::mutate(date = as.Date(date), site = factor(site, levels = site_levels_project)),
    train = readr::read_csv(paths$train, show_col_types = FALSE) |>
      dplyr::mutate(date = as.Date(date), site = factor(site, levels = site_levels_project)),
    test = readr::read_csv(paths$test, show_col_types = FALSE) |>
      dplyr::mutate(date = as.Date(date), site = factor(site, levels = site_levels_project))
  )
}

save_prepared_data_cache_project <- function(shared_data_dir, dat4, split, start_date, test_months) {
  paths <- shared_data_cache_paths_project(shared_data_dir)
  readr::write_csv(dat4, paths$aligned)
  readr::write_csv(split$train, paths$train)
  readr::write_csv(split$test, paths$test)
  save_json_pretty_project(
    list(
      start_date = format(start_date),
      test_months = as.integer(test_months),
      n_aligned = nrow(dat4),
      n_train = nrow(split$train),
      n_test = nrow(split$test)
    ),
    paths$metadata
  )
  invisible(paths)
}

rmse_project <- function(y, y_hat) sqrt(mean((y - y_hat)^2, na.rm = TRUE))
mae_project <- function(y, y_hat) mean(abs(y - y_hat), na.rm = TRUE)
coverage_project <- function(y, y_lo, y_hi) mean(y >= y_lo & y <= y_hi, na.rm = TRUE)
mlpd_gaussian_project <- function(y, mean_vec, sd_vec) mean(stats::dnorm(y, mean = mean_vec, sd = sd_vec, log = TRUE), na.rm = TRUE)
harmonic_amplitude_project <- function(beta_sin, beta_cos) sqrt(beta_sin^2 + beta_cos^2)

save_json_pretty_project <- function(x, path) {
  writeLines(jsonlite::toJSON(x, auto_unbox = TRUE, pretty = TRUE, null = "null"), con = path)
}
