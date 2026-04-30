// co2_model_2.stan
// =============================================================
// Port of the Python notebook model (CO2/stan_models/co2_model_2.stan).
//
// Structure:
//   y[i, t] = level[t]
//           + site_intercept[i]
//           + site_slope[i] * x_time[t]
//           + seasonal[i, t]
//           + noise_i
//
//   level is a local-linear-trend random walk with a slowly drifting
//   growth[t] process; site_intercept / site_slope are demeaned soft
//   hierarchies; 2-harmonic seasonality is shared-scale-hierarchical;
//   per-site observation noise sigma_obs[i].
//
// Compared to the notebook original, this version additionally exports
// sigma_test and log_lik_test in generated quantities so R can compute
// exact holdout MLPD (matching the schema used by M4 / M4c).
// =============================================================

data {
  int<lower=1> N;
  int<lower=1> S;
  int<lower=1> T;

  array[N] int<lower=1, upper=S> site;
  array[N] int<lower=1, upper=T> t_id;
  array[N] int<lower=1, upper=12> month_id;
  vector[N] y;

  int<lower=0> N_test;
  array[N_test] int<lower=1, upper=S> site_test;
  array[N_test] int<lower=1, upper=T> t_test;
  array[N_test] int<lower=1, upper=12> month_test;
  vector[N_test] y_test;
}

transformed data {
  vector[12] sin1;
  vector[12] cos1;
  vector[12] sin2;
  vector[12] cos2;
  vector[T] x_time;

  for (m in 1:12) {
    real ang1 = 2 * pi() * (m - 1) / 12.0;
    real ang2 = 2 * ang1;
    sin1[m] = sin(ang1);
    cos1[m] = cos(ang1);
    sin2[m] = sin(ang2);
    cos2[m] = cos(ang2);
  }

  for (t in 1:T) {
    x_time[t] = (t - 1 - (T - 1) / 2.0) / 12.0;
  }
}

parameters {
  real level_1;
  real growth_1;
  vector[T - 1] z_level;
  vector[T - 1] z_growth;
  real<lower=0> sigma_level;
  real<lower=0> sigma_growth;

  vector[S] site_intercept_raw;
  real<lower=0> sigma_site_intercept;

  vector[S] site_slope_raw;
  real<lower=0> sigma_site_slope;

  vector[S] beta_sin1_raw;
  vector[S] beta_cos1_raw;
  vector[S] beta_sin2_raw;
  vector[S] beta_cos2_raw;
  real<lower=0> sigma_h1;
  real<lower=0> sigma_h2;

  vector<lower=0>[S] sigma_obs;
}

transformed parameters {
  vector[T] level;
  vector[T] growth;

  vector[S] site_intercept;
  vector[S] site_slope;

  vector[S] beta_sin1;
  vector[S] beta_cos1;
  vector[S] beta_sin2;
  vector[S] beta_cos2;

  level[1] = level_1;
  growth[1] = growth_1;

  for (t in 2:T) {
    growth[t] = growth[t - 1] + sigma_growth * z_growth[t - 1];
    level[t] = level[t - 1] + growth[t - 1] + sigma_level * z_level[t - 1];
  }

  site_intercept =
    sigma_site_intercept * (site_intercept_raw - mean(site_intercept_raw));

  site_slope =
    sigma_site_slope * (site_slope_raw - mean(site_slope_raw));

  beta_sin1 = sigma_h1 * beta_sin1_raw;
  beta_cos1 = sigma_h1 * beta_cos1_raw;
  beta_sin2 = sigma_h2 * beta_sin2_raw;
  beta_cos2 = sigma_h2 * beta_cos2_raw;
}

model {
  level_1 ~ normal(370, 20);
  growth_1 ~ normal(0.18, 0.05);

  z_level ~ std_normal();
  z_growth ~ std_normal();

  sigma_level ~ exponential(12);
  sigma_growth ~ exponential(200);

  site_intercept_raw ~ std_normal();
  sigma_site_intercept ~ exponential(0.3);

  site_slope_raw ~ std_normal();
  sigma_site_slope ~ exponential(2);

  beta_sin1_raw ~ std_normal();
  beta_cos1_raw ~ std_normal();
  beta_sin2_raw ~ std_normal();
  beta_cos2_raw ~ std_normal();

  sigma_h1 ~ exponential(0.4);
  sigma_h2 ~ exponential(1.0);

  sigma_obs ~ exponential(2.0);

  for (n in 1:N) {
    real seasonal_part =
      beta_sin1[site[n]] * sin1[month_id[n]] +
      beta_cos1[site[n]] * cos1[month_id[n]] +
      beta_sin2[site[n]] * sin2[month_id[n]] +
      beta_cos2[site[n]] * cos2[month_id[n]];

    real mu_n =
      level[t_id[n]] +
      site_intercept[site[n]] +
      site_slope[site[n]] * x_time[t_id[n]] +
      seasonal_part;

    y[n] ~ normal(mu_n, sigma_obs[site[n]]);
  }
}

generated quantities {
  vector[N] mu;
  vector[N] y_rep;
  vector[N] log_lik;

  vector[N_test] mu_test;
  vector[N_test] y_test_rep;
  vector<lower=0>[N_test] sigma_test;
  vector[N_test] log_lik_test;

  for (n in 1:N) {
    real seasonal_part =
      beta_sin1[site[n]] * sin1[month_id[n]] +
      beta_cos1[site[n]] * cos1[month_id[n]] +
      beta_sin2[site[n]] * sin2[month_id[n]] +
      beta_cos2[site[n]] * cos2[month_id[n]];

    mu[n] =
      level[t_id[n]] +
      site_intercept[site[n]] +
      site_slope[site[n]] * x_time[t_id[n]] +
      seasonal_part;

    y_rep[n] = normal_rng(mu[n], sigma_obs[site[n]]);
    log_lik[n] = normal_lpdf(y[n] | mu[n], sigma_obs[site[n]]);
  }

  for (n in 1:N_test) {
    real seasonal_part =
      beta_sin1[site_test[n]] * sin1[month_test[n]] +
      beta_cos1[site_test[n]] * cos1[month_test[n]] +
      beta_sin2[site_test[n]] * sin2[month_test[n]] +
      beta_cos2[site_test[n]] * cos2[month_test[n]];

    mu_test[n] =
      level[t_test[n]] +
      site_intercept[site_test[n]] +
      site_slope[site_test[n]] * x_time[t_test[n]] +
      seasonal_part;

    sigma_test[n] = sigma_obs[site_test[n]];
    y_test_rep[n] = normal_rng(mu_test[n], sigma_test[n]);
    log_lik_test[n] = normal_lpdf(y_test[n] | mu_test[n], sigma_test[n]);
  }
}
