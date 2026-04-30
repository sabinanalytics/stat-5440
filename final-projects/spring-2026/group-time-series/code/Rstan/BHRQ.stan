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
}

transformed data {
  vector[12] sin1;
  vector[12] cos1;
  vector[12] sin2;
  vector[12] cos2;

  vector[T] x_time;
  vector[T] x2_time;

  real t_center = (T + 1) / 2.0;

  for (m in 1:12) {
    real ang1 = 2 * pi() * (m - 1) / 12.0;
    real ang2 = 2 * ang1;

    sin1[m] = sin(ang1);
    cos1[m] = cos(ang1);
    sin2[m] = sin(ang2);
    cos2[m] = cos(ang2);
  }

  for (t in 1:T) {
    real x = (t - t_center) / 12.0;
    x_time[t] = x;
    x2_time[t] = square(x);
  }
}

parameters {
  real alpha;
  real beta_global;
  real q_global;

  vector[S] site_intercept_raw;
  vector[S] site_slope_raw;
  vector[S] site_curve_raw;

  real<lower=0> sigma_site_intercept;
  real<lower=0> sigma_site_slope;
  real<lower=0> sigma_site_curve;

  vector[S] beta_sin1_raw;
  vector[S] beta_cos1_raw;
  vector[S] beta_sin2_raw;
  vector[S] beta_cos2_raw;

  real<lower=0> sigma_h1;
  real<lower=0> sigma_h2;

  vector<lower=0>[S] sigma_obs;
}

transformed parameters {
  vector[S] site_intercept;
  vector[S] site_slope;
  vector[S] site_curve;

  vector[S] beta_sin1;
  vector[S] beta_cos1;
  vector[S] beta_sin2;
  vector[S] beta_cos2;

  site_intercept =
    sigma_site_intercept * (site_intercept_raw - mean(site_intercept_raw));

  site_slope =
    sigma_site_slope * (site_slope_raw - mean(site_slope_raw));

  site_curve =
    sigma_site_curve * (site_curve_raw - mean(site_curve_raw));

  beta_sin1 = sigma_h1 * beta_sin1_raw;
  beta_cos1 = sigma_h1 * beta_cos1_raw;
  beta_sin2 = sigma_h2 * beta_sin2_raw;
  beta_cos2 = sigma_h2 * beta_cos2_raw;
}

model {
  alpha ~ normal(400, 25);
  beta_global ~ normal(2.2, 0.8);
  q_global ~ normal(0, 0.08);

  site_intercept_raw ~ std_normal();
  site_slope_raw ~ std_normal();
  site_curve_raw ~ std_normal();

  sigma_site_intercept ~ normal(0, 1.5);
  sigma_site_slope ~ normal(0, 0.1);
  sigma_site_curve ~ normal(0, 0.04);

  beta_sin1_raw ~ std_normal();
  beta_cos1_raw ~ std_normal();
  beta_sin2_raw ~ std_normal();
  beta_cos2_raw ~ std_normal();

  sigma_h1 ~ normal(0, 3);
  sigma_h2 ~ normal(0, 1.5);

  sigma_obs ~ normal(0, 0.75);

  for (n in 1:N) {
    real seasonal_part =
      beta_sin1[site[n]] * sin1[month_id[n]] +
      beta_cos1[site[n]] * cos1[month_id[n]] +
      beta_sin2[site[n]] * sin2[month_id[n]] +
      beta_cos2[site[n]] * cos2[month_id[n]];

    real trend_part =
      alpha +
      site_intercept[site[n]] +
      (beta_global + site_slope[site[n]]) * x_time[t_id[n]] +
      (q_global + site_curve[site[n]]) * x2_time[t_id[n]];

    y[n] ~ normal(trend_part + seasonal_part, sigma_obs[site[n]]);
  }
}

generated quantities {
  vector[N] mu;
  vector[N] y_rep;
  vector[N] log_lik;

  vector[N_test] mu_test;
  vector[N_test] y_test_rep;

  for (n in 1:N) {
    real seasonal_part =
      beta_sin1[site[n]] * sin1[month_id[n]] +
      beta_cos1[site[n]] * cos1[month_id[n]] +
      beta_sin2[site[n]] * sin2[month_id[n]] +
      beta_cos2[site[n]] * cos2[month_id[n]];

    real trend_part =
      alpha +
      site_intercept[site[n]] +
      (beta_global + site_slope[site[n]]) * x_time[t_id[n]] +
      (q_global + site_curve[site[n]]) * x2_time[t_id[n]];

    mu[n] = trend_part + seasonal_part;
    y_rep[n] = normal_rng(mu[n], sigma_obs[site[n]]);
    log_lik[n] = normal_lpdf(y[n] | mu[n], sigma_obs[site[n]]);
  }

  for (n in 1:N_test) {
    real seasonal_part =
      beta_sin1[site_test[n]] * sin1[month_test[n]] +
      beta_cos1[site_test[n]] * cos1[month_test[n]] +
      beta_sin2[site_test[n]] * sin2[month_test[n]] +
      beta_cos2[site_test[n]] * cos2[month_test[n]];

    real trend_part =
      alpha +
      site_intercept[site_test[n]] +
      (beta_global + site_slope[site_test[n]]) * x_time[t_test[n]] +
      (q_global + site_curve[site_test[n]]) * x2_time[t_test[n]];

    mu_test[n] = trend_part + seasonal_part;
    y_test_rep[n] = normal_rng(mu_test[n], sigma_obs[site_test[n]]);
  }
}
