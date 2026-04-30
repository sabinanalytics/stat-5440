// co2_model_2_upgrade_ar2.stan
// =============================================================
// Upgrade of co2_model_2.stan with a station-specific latent AR(2)
// residual process added on top of the original Model 2 backbone.
//
// Original backbone preserved:
//   y[s, t] =
//     level[t]
//     + site_intercept[s]
//     + site_slope[s] * x_time[t]
//     + harmonic_seasonality[s, month(t)]
//     + residual_state[s, t]
//     + observation_noise[s, t]
//
// New residual state:
//   r[s, t] = phi1[s] * r[s, t - 1] + phi2[s] * r[s, t - 2] + u[s, t]
//   u[s, t] ~ normal(0, sigma_ar[s])
//
// Stable parameterization:
//   station-specific PACF coefficients are mapped through tanh() and
//   then converted to AR(2) coefficients with the Durbin-Levinson map,
//   which guarantees stationarity.
// =============================================================

data {
  int<lower=1> N;
  int<lower=1> S;
  int<lower=2> T;

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

  real mu_pacf1;
  real mu_pacf2;
  real<lower=0> sigma_pacf1;
  real<lower=0> sigma_pacf2;
  vector[S] pacf1_raw;
  vector[S] pacf2_raw;
  vector<lower=0>[S] sigma_ar;
  vector[S] resid_init1_raw;
  vector[S] resid_init2_raw;
  matrix[S, T - 2] z_resid;

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

  vector[S] pacf1;
  vector[S] pacf2;
  vector[S] phi1;
  vector[S] phi2;
  matrix[S, T] ar2_residual;

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

  for (s in 1:S) {
    pacf1[s] = tanh(mu_pacf1 + sigma_pacf1 * pacf1_raw[s]);
    pacf2[s] = tanh(mu_pacf2 + sigma_pacf2 * pacf2_raw[s]);

    phi1[s] = pacf1[s] * (1 - pacf2[s]);
    phi2[s] = pacf2[s];

    ar2_residual[s, 1] = sigma_ar[s] * resid_init1_raw[s];
    ar2_residual[s, 2] = sigma_ar[s] * resid_init2_raw[s];

    for (t in 3:T) {
      ar2_residual[s, t] =
        phi1[s] * ar2_residual[s, t - 1]
        + phi2[s] * ar2_residual[s, t - 2]
        + sigma_ar[s] * z_resid[s, t - 2];
    }
  }
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

  mu_pacf1 ~ normal(0.4, 0.35);
  mu_pacf2 ~ normal(0, 0.25);
  sigma_pacf1 ~ normal(0, 0.25);
  sigma_pacf2 ~ normal(0, 0.20);
  pacf1_raw ~ std_normal();
  pacf2_raw ~ std_normal();
  sigma_ar ~ exponential(2.0);
  resid_init1_raw ~ std_normal();
  resid_init2_raw ~ std_normal();
  to_vector(z_resid) ~ std_normal();

  sigma_obs ~ exponential(2.0);

  for (n in 1:N) {
    int s = site[n];
    int t = t_id[n];
    real seasonal_part =
      beta_sin1[s] * sin1[month_id[n]] +
      beta_cos1[s] * cos1[month_id[n]] +
      beta_sin2[s] * sin2[month_id[n]] +
      beta_cos2[s] * cos2[month_id[n]];

    real mu_n =
      level[t] +
      site_intercept[s] +
      site_slope[s] * x_time[t] +
      seasonal_part +
      ar2_residual[s, t];

    y[n] ~ normal(mu_n, sigma_obs[s]);
  }
}

generated quantities {
  vector[N] mu;
  vector[N] trend_part;
  vector[N] seasonal_part;
  vector[N] residual_part;
  vector[N] y_rep;
  vector[N] log_lik;

  vector[N_test] mu_test;
  vector[N_test] trend_part_test;
  vector[N_test] seasonal_part_test;
  vector[N_test] residual_part_test;
  vector[N_test] y_test_rep;
  vector<lower=0>[N_test] sigma_test;
  vector[N_test] log_lik_test;

  for (n in 1:N) {
    int s = site[n];
    int t = t_id[n];

    trend_part[n] =
      level[t] + site_intercept[s] + site_slope[s] * x_time[t];
    seasonal_part[n] =
      beta_sin1[s] * sin1[month_id[n]] +
      beta_cos1[s] * cos1[month_id[n]] +
      beta_sin2[s] * sin2[month_id[n]] +
      beta_cos2[s] * cos2[month_id[n]];
    residual_part[n] = ar2_residual[s, t];
    mu[n] = trend_part[n] + seasonal_part[n] + residual_part[n];

    y_rep[n] = normal_rng(mu[n], sigma_obs[s]);
    log_lik[n] = normal_lpdf(y[n] | mu[n], sigma_obs[s]);
  }

  for (n in 1:N_test) {
    int s = site_test[n];
    int t = t_test[n];

    trend_part_test[n] =
      level[t] + site_intercept[s] + site_slope[s] * x_time[t];
    seasonal_part_test[n] =
      beta_sin1[s] * sin1[month_test[n]] +
      beta_cos1[s] * cos1[month_test[n]] +
      beta_sin2[s] * sin2[month_test[n]] +
      beta_cos2[s] * cos2[month_test[n]];
    residual_part_test[n] = ar2_residual[s, t];
    mu_test[n] =
      trend_part_test[n] + seasonal_part_test[n] + residual_part_test[n];

    sigma_test[n] = sigma_obs[s];
    y_test_rep[n] = normal_rng(mu_test[n], sigma_test[n]);
    log_lik_test[n] = normal_lpdf(y_test[n] | mu_test[n], sigma_test[n]);
  }
}
