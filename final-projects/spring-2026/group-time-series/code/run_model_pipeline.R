#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)

usage_error <- function(...) {
  stop(..., call. = FALSE)
}

normalize_maybe <- function(path) {
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

find_code_dir <- function() {
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

code_dir <- find_code_dir()
rscript_dir <- file.path(code_dir, "Rscript")
config_dir <- file.path(code_dir, "config")

resolve_existing_path <- function(path_string, base_dirs = c(".", code_dir, config_dir)) {
  candidates <- unique(c(path_string, file.path(base_dirs, path_string)))
  for (candidate in candidates) {
    if (file.exists(candidate) || dir.exists(candidate)) {
      return(normalizePath(candidate, winslash = "/", mustWork = FALSE))
    }
  }
  NULL
}

resolve_pipeline_target <- function(input) {
  key <- tolower(trimws(input))

  if (key %in% c("bayesian", "harmonic", "sarima")) {
    cfg_name <- switch(
      key,
      bayesian = "bayesian_config.json",
      harmonic = "harmonic_config.json",
      sarima = "sarima_config.json"
    )
    return(list(
      model_type = key,
      config_path = file.path(config_dir, cfg_name)
    ))
  }

  config_path <- resolve_existing_path(input)
  if (is.null(config_path) || !file.exists(config_path)) {
    usage_error(
      "Could not resolve pipeline target from `", input, "`.\n",
      "Use one of `bayesian`, `harmonic`, `sarima`, or provide a config file path."
    )
  }

  cfg_base <- tolower(basename(config_path))
  model_type <- if (grepl("bayesian", cfg_base)) {
    "bayesian"
  } else if (grepl("harmonic", cfg_base)) {
    "harmonic"
  } else if (grepl("sarima", cfg_base)) {
    "sarima"
  } else {
    usage_error(
      "Cannot infer model type from config file name `", basename(config_path), "`.\n",
      "Please use a config file whose name contains `bayesian`, `harmonic`, or `sarima`."
    )
  }

  list(
    model_type = model_type,
    config_path = config_path
  )
}

read_config <- function(path) {
  if (!file.exists(path)) {
    usage_error("Config file not found: ", normalize_maybe(path))
  }
  cfg <- jsonlite::fromJSON(path, simplifyVector = TRUE)
  if (!is.list(cfg)) {
    usage_error("Config JSON must decode to an object.")
  }
  cfg
}

resolve_visual_target <- function(model_type, cfg) {
  if (identical(model_type, "bayesian")) {
    stan_model <- as.character(cfg$stan_model %||% "")
    reserved <- c("BHRQ", "BDLM", "BDLMQ", "BDLM_AR2")
    if (stan_model %in% reserved) {
      return(stan_model)
    }
    if (!nzchar(stan_model)) {
      usage_error("Bayesian config must contain `stan_model`.")
    }
    return(tools::file_path_sans_ext(basename(stan_model)))
  }

  as.character(cfg$model_name %||% if (identical(model_type, "harmonic")) "Harmonic_Regression" else "SARIMA")
}

script_for_model <- function(model_type) {
  switch(
    model_type,
    bayesian = file.path(rscript_dir, "run_bayesian_model.R"),
    harmonic = file.path(rscript_dir, "run_harmonic_regression.R"),
    sarima = file.path(rscript_dir, "run_sarima_model.R"),
    usage_error("Unsupported model type: ", model_type)
  )
}

run_rscript <- function(script_path, script_args) {
  rscript_bin <- file.path(R.home("bin"), "Rscript")
  cmd <- c(script_path, script_args)
  status <- system2(rscript_bin, args = cmd)
  if (!identical(status, 0L)) {
    usage_error(
      "Command failed: ", normalize_maybe(rscript_bin), " ",
      paste(shQuote(cmd), collapse = " ")
    )
  }
}

if (length(args) < 1L) {
  usage_error(
    "Usage: Rscript run_model_pipeline.R <bayesian|harmonic|sarima|config_path>\n",
    "Example: Rscript run_model_pipeline.R bayesian"
  )
}

target <- resolve_pipeline_target(args[1])
cfg <- read_config(target$config_path)
visual_target <- resolve_visual_target(target$model_type, cfg)
train_script <- script_for_model(target$model_type)
visual_script <- file.path(rscript_dir, "visualize_single_model.R")

message("Step 0/1: data read / preprocessing / training / result saving")
run_rscript(train_script, script_args = c(target$config_path))

message("Step 2: plot saving")
run_rscript(visual_script, script_args = c(visual_target))

message("Finished pipeline for model: ", visual_target)
