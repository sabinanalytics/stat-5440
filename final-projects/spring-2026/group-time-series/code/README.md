# README

## 1. Model Name Mapping

The full names of the models and their abbreviations used in this project are:

- `HR`: Harmonic Regression
- `SARIMA`: Seasonal AutoRegressive Integrated Moving Average
- `BHRQ`: Bayesian Harmonic Regression with Quadratic Trend
- `BDLM`: Bayesian Dynamic Linear Model
- `BDLMQ`: Bayesian Dynamic Linear Model with Quadratic Observation Component
- `BDLM_AR2`: Bayesian Dynamic Linear Model with AR(2) Latent Residual Extension

Notes:

- `HR` and `SARIMA` are classical benchmark models.
- `BHRQ`, `BDLM`, `BDLMQ`, and `BDLM_AR2` are Bayesian models.
- In `Result/`, `Plot/`, and comparison outputs, these abbreviations or the configured `model_name` are used as the model identifiers.

## 2. `code` Folder Structure

The current structure of `group-time-series/code` is:

```text
group-time-series/code
|-- config/
|   |-- bayesian_config.json
|   |-- comparison_config.json
|   |-- harmonic_config.json
|   `-- sarima_config.json
|-- Data/
|   `-- NOAA_CO2_Monthly/
|-- Plot/
|-- Result/
|-- Rscript/
|   |-- _project_model_common.R
|   |-- run_bayesian_model.R
|   |-- run_harmonic_regression.R
|   |-- run_model_comparison.R
|   |-- run_sarima_model.R
|   `-- visualize_single_model.R
|-- Rstan/
|   |-- BDLM.stan
|   |-- BDLM_AR2.stan
|   |-- BDLMQ.stan
|   `-- BHRQ.stan
|-- run_model_pipeline.R
`-- README.md
```

The purpose of each part is:

- `config/`: stores all experiment configuration files in one place.
- `Data/NOAA_CO2_Monthly/`: the shared data directory for all models. On the first run, NOAA monthly CO2 data will be downloaded here automatically. Later runs will reuse the cached files and will not download them again by default. This folder also stores shared preprocessing outputs and train/test split caches.
- `Result/`: stores numerical outputs for each model.
- `Plot/`: stores single-model plots and model-comparison plots.
- `Rscript/`: stores the main R scripts for training, visualization, and comparison.
- `Rstan/`: stores the Stan programs used by the Bayesian models.
- `run_model_pipeline.R`: the one-stop pipeline controller for a single model run.

## 3. How To Run Each Script

### 3.1 `run_harmonic_regression.R`

Run:

```powershell
Rscript group-time-series/code/Rscript/run_harmonic_regression.R
```

By default it reads:

```text
group-time-series/code/config/harmonic_config.json
```

You can also pass a config path manually:

```powershell
Rscript group-time-series/code/Rscript/run_harmonic_regression.R group-time-series/code/config/harmonic_config.json
```

Main config fields:

- `model_name`: the output folder name for this experiment.
- `data`: the data source.
  - `"default"`: first checks the shared data under `group-time-series/code/Data/NOAA_CO2_Monthly/`; if the NOAA raw files or the required cached split do not exist, the script downloads or rebuilds them automatically.
  - `"force download"`: forces a fresh download of the NOAA raw files.
  - `"<path>"`: uses a user-provided data directory or data file.
- `test_months`: the number of final months used as the test set.
- `start_date`: the start date of the sample window.
- `prediction_level`: the prediction interval level, for example `0.9`.

Main functionality:

- Reads and preprocesses data.
- Fits the Harmonic Regression model.
- Saves forecasting outputs, metrics, residual ACF results, and key results.

To run different HR experiments with the same script, you only need to modify `harmonic_config.json`, for example:

- change `model_name` to separate different experiment outputs;
- change `start_date` to alter the sample window;
- change `test_months` to alter the train/test split;
- change `prediction_level` to alter interval width;
- change `data` to switch data sources.

### 3.2 `run_sarima_model.R`

Run:

```powershell
Rscript group-time-series/code/Rscript/run_sarima_model.R
```

By default it reads:

```text
group-time-series/code/config/sarima_config.json
```

Main config fields:

- `model_name`: the output folder name for this experiment.
- `data`: the data source, with the same logic described above.
- `test_months`: the number of final months used as the test set.
- `start_date`: the start date of the sample window.
- `prediction_level`: the prediction interval level.
- `seasonal_period`: the seasonal period, usually `12` for monthly data.
- `max_p`, `max_q`, `max_P`, `max_Q`: upper bounds for the SARIMA search.
- `max_order`: maximum total order constraint.
- `verbose`: whether to print search progress.

Main functionality:

- Reads and preprocesses data.
- Fits a SARIMA model separately for the 4 sites.
- Saves forecasting outputs, metrics, selected orders, residual ACF results, and key results.

To run different SARIMA experiments, the main way is to modify `sarima_config.json`, especially:

- `model_name`
- `start_date`
- `test_months`
- `seasonal_period`
- `max_p`, `max_q`, `max_P`, `max_Q`, `max_order`

### 3.3 `run_bayesian_model.R`

Run:

```powershell
Rscript group-time-series/code/Rscript/run_bayesian_model.R
```

By default it reads:

```text
group-time-series/code/config/bayesian_config.json
```

You can also pass a config path manually:

```powershell
Rscript group-time-series/code/Rscript/run_bayesian_model.R group-time-series/code/config/bayesian_config.json
```

Main config fields:

- `stan_model`: specifies the Bayesian model or a custom Stan file path.
  - `BHRQ`
  - `BDLM`
  - `BDLMQ`
  - `BDLM_AR2`
  - or a custom `.stan` file path
- `data`: the data source, with the same logic described above.
- `seed`: random seed.
- `chains`: number of MCMC chains.
- `cores`: number of parallel cores.
- `iter`: total number of iterations.
- `warmup`: number of warmup iterations.
- `refresh`: Stan log print frequency.
- `test_months`: the number of final months used as the test set.
- `start_date`: the start date of the sample window.
- `control`:
  - `adapt_delta`
  - `max_treedepth`

Main functionality:

- Reads and preprocesses data.
- Compiles and runs the selected Stan model.
- Saves forecasting outputs, posterior summaries, sampler diagnostics, residual ACF results, and Bayesian artifacts such as `stanfit.rds`.

To run different Bayesian experiments with the same script, you mainly modify `bayesian_config.json`, for example:

- change `stan_model` among `BHRQ`, `BDLM`, `BDLMQ`, and `BDLM_AR2`;
- change `iter`, `warmup`, `chains`, `cores`, and `control` to alter the sampling setup;
- change `start_date` and `test_months` to alter the sample window;
- change `data` to switch data sources.

Note: the Bayesian result folder name is usually determined automatically from `stan_model`.

### 3.4 `visualize_single_model.R`

Run:

```powershell
Rscript group-time-series/code/Rscript/visualize_single_model.R BDLM
```

or:

```powershell
Rscript group-time-series/code/Rscript/visualize_single_model.R group-time-series/code/Result/BDLM
```

Input logic:

- You can pass a model name directly, and the script will look under `Result/<model_name>/`.
- You can also pass a result-directory path directly.

Main functionality:

- Checks whether the result directory satisfies the plotting requirements.
- Generates, for all models:
  - train forecast plots;
  - test forecast plots;
  - residual ACF plots;
  - key-result plots.
- Generates, for Bayesian models only:
  - posterior Rhat vs ESS plots;
  - posterior traceplots.

### 3.5 `run_model_comparison.R`

Run:

```powershell
Rscript group-time-series/code/Rscript/run_model_comparison.R
```

By default it reads:

```text
group-time-series/code/config/comparison_config.json
```

Main config fields:

- `comparison_name`: the output folder name for the comparison.
- `models`: the list of models included in the comparison. Each entry should contain:
  - `model_name`: corresponding to `Result/<model_name>/`
  - `model_type`: `harmonic`, `sarima`, or `bayesian`
  - `label`: the display name used in the plots

Main functionality:

- Performs a unified comparison for the models specified in the config.
- For all models, generates:
  - forecast panels with dimensions model-count by 4 sites;
  - forecast metric comparison panels with dimensions 3 by 4, where the three rows correspond to `MAE`, `RMSE`, and `LPD`;
  - residual ACF panels with dimensions model-count by 4 sites.
- For Bayesian models included in the comparison, additionally generates:
  - `bundle_rhat_ess_panels`
  - traceplot comparison figures

Comparison plots are saved to:

```text
group-time-series/code/Plot/comparison/<comparison_name>/
```

### 3.6 `run_model_pipeline.R`

Run:

```powershell
Rscript group-time-series/code/run_model_pipeline.R harmonic
Rscript group-time-series/code/run_model_pipeline.R sarima
Rscript group-time-series/code/run_model_pipeline.R bayesian
```

You can also pass a config path directly:

```powershell
Rscript group-time-series/code/run_model_pipeline.R group-time-series/code/config/harmonic_config.json
Rscript group-time-series/code/run_model_pipeline.R group-time-series/code/config/sarima_config.json
Rscript group-time-series/code/run_model_pipeline.R group-time-series/code/config/bayesian_config.json
```

Main functionality:

1. If the shared data cache does not exist, it automatically performs data reading and preprocessing.
2. It trains the target model and saves its outputs under `Result/`.
3. It then calls `visualize_single_model.R` automatically and saves the corresponding plots under `Plot/`.

Important note: `run_model_pipeline.R` currently accepts `harmonic`, `sarima`, `bayesian`, or a config-file path.  
If you want to use it to run `BHRQ`, `BDLM`, `BDLMQ`, or `BDLM_AR2`, the recommended approach is to modify the `stan_model` field in `bayesian_config.json`, then rerun `run_model_pipeline.R`.

## 4. How To Use Configs For Different Experiments

The core design of this project is to minimize duplicated scripts and switch experiments mainly through the config files in `config/`.

The general rules are:

- modify `model_name` or `comparison_name` to avoid overwriting previous outputs;
- modify `data` to switch between default cached data, forced download, or custom input paths;
- modify `start_date` and `test_months` to change the sample window and the train/test split;
- modify model-specific parameters to change the actual experiment setup.

Therefore:

- the same `run_harmonic_regression.R` can be used for many HR experiments;
- the same `run_sarima_model.R` can be used for many SARIMA experiments;
- the same `run_bayesian_model.R` can be used for many Bayesian experiments;
- the same `run_model_comparison.R` can be used for many different model-comparison experiments by editing the comparison config.

## 5. Result And Plot Output Logic

### 5.1 Shared Data

The shared data directory for all models is:

```text
group-time-series/code/Data/NOAA_CO2_Monthly/
```

This folder stores:

- NOAA raw monthly CO2 files;
- aligned panel-data cache;
- train/test split cache;
- metadata describing the current cached preprocessing setup.

As a result:

- on the first run, data will be downloaded automatically into `group-time-series/code/Data/`;
- later runs will not download the same data again by default;
- model result folders do not duplicate `train_data.csv` or `test_data.csv`.

### 5.2 Model Result Output

The numerical outputs for each model are saved to:

```text
group-time-series/code/Result/<model_name>/
```

They are currently organized into the following thematic subfolders:

- `forecast/`
- `posterior/`
- `time_series_diagnose/`
- `key_result/`

In particular:

- `forecast/` stores prediction- and metric-related outputs;
- `posterior/` is mainly used by Bayesian models and stores posterior summaries, diagnostics, `stanfit.rds`, and related files;
- `time_series_diagnose/` stores time-series diagnostic outputs such as residual ACF results;
- `key_result/` stores model-specific key result tables.

### 5.3 Plot Output

Single-model plots are saved to:

```text
group-time-series/code/Plot/<model_name>/
```

Comparison plots are saved to:

```text
group-time-series/code/Plot/comparison/<comparison_name>/
```

They are also organized by theme, for example:

- `forecast/`
- `posterior/`
- `time_series_diagnose/`
- `key_result/`

## 6. Reproduction Note

It is important to note that:

- `group-time-series/code/` does not store the actual experiment outputs;
- instead, it presents a reasonably organized project structure for running and reproducing experiments.

If you want to reproduce the main contents shown in our slides, a reasonable pipeline is:

1. Use `group-time-series/code/run_model_pipeline.R` to run the following model experiments one by one:
   - `HR`
   - `SARIMA`
   - `BHRQ`
   - `BDLM`
   - `BDLMQ`
   - `BDLM_AR2`
2. Then run:

```powershell
Rscript group-time-series/code/Rscript/run_model_comparison.R
```

to compare the models.

More precisely:

- for `HR`, modify `harmonic_config.json` as needed and run `run_model_pipeline.R harmonic`;
- for `SARIMA`, modify `sarima_config.json` as needed and run `run_model_pipeline.R sarima`;
- for `BHRQ`, `BDLM`, `BDLMQ`, and `BDLM_AR2`, modify the `stan_model` field in `bayesian_config.json` and run `run_model_pipeline.R bayesian` separately for each case.

On the first run, data will be downloaded automatically into `group-time-series/code/Data/`. After that, later runs will not repeat the download by default.

## 7. Runtime Warning

The following runtimes are based on my local machine only and are given for reference. Different machines, parallel settings, and experiment configurations may lead to very different runtimes:

- `HR`, `BHRQ`: within 10 minutes
- `SARIMA`, `BDLM`, `BDLMQ`: around 2 to 3 hours
- `BDLM_AR2`: around 15 hours

Please therefore run some parts of the code cautiously, especially:

- `SARIMA`
- `BDLM`
- `BDLMQ`
- `BDLM_AR2`

## 8. Old Results

My previously saved local results are stored in:

```text
group-time-series/old_result
```

If you want to inspect historical outputs, you can check that folder directly.

## 9. Possible Bugs

These scripts were produced while I was reorganizing the project from a previously messy structure into the current organized version, so bugs may still exist.

If you find any bugs, please contact:

Tianyi Wang  
`wtyyy6@sas.upenn.edu`
