data {
  int<lower=1> N_train;
  int<lower=1> N_test;
  int<lower=1> D_loc;              // usually 2: latitude and longitude

  // Coordinate features for spatial RBF kernel
  array[N_train] vector[D_loc] X_loc_train;
  array[N_test] vector[D_loc] X_loc_test;

  // Standardized MedInc for linear income kernel
  vector[N_train] x_inc_train;
  vector[N_test] x_inc_test;

  vector[N_train] y_train;
}

transformed data { 
  real jitter = 1e-5;    // small diagonal term for numerical stability
}

parameters {
  // Spatial RBF kernel hyperparameters
  real<lower=1e-6> alpha_loc;
  real<lower=1e-6> rho_loc;

  // Linear MedInc kernel amplitude
  real<lower=1e-6> alpha_inc;

  // Observation noise
  real<lower=1e-6> sigma;
}

model {
  matrix[N_train, N_train] K_loc;
  matrix[N_train, N_train] K_inc;
  matrix[N_train, N_train] K;
  matrix[N_train, N_train] L_K;

  // Priors
  alpha_loc ~ lognormal(0, 0.5);
  rho_loc ~ lognormal(0, 0.5);
  alpha_inc ~ lognormal(0, 0.5);
  sigma ~ lognormal(-1, 0.5);

  // 1) Spatial RBF kernel on coordinates
  K_loc = gp_exp_quad_cov(X_loc_train, alpha_loc, rho_loc);

  // 2) Linear kernel on standardized MedInc:
  //    k_inc(m_i, m_j) = alpha_inc^2 * m_i * m_j
  K_inc = square(alpha_inc) * (x_inc_train * x_inc_train');

  // Additive kernel: K = K_loc + K_inc + noise
  K = K_loc + K_inc;

  for (i in 1:N_train) {
    K[i, i] = K[i, i] + square(sigma) + jitter;
  }

  L_K = cholesky_decompose(K);
  y_train ~ multi_normal_cholesky(rep_vector(0, N_train), L_K);
}

generated quantities {
  matrix[N_train, N_train] K_loc;
  matrix[N_train, N_train] K_inc;
  matrix[N_train, N_train] K;

  matrix[N_train, N_test] K_s_loc;
  matrix[N_train, N_test] K_s_inc;
  matrix[N_train, N_test] K_s;

  matrix[N_test, N_test] K_ss_loc;
  matrix[N_test, N_test] K_ss_inc;
  matrix[N_test, N_test] K_ss;
  matrix[N_test, N_test] K_pred;

  vector[N_train] alpha_vec;
  vector[N_test] f_test_mean;
  vector[N_test] y_test_rep;

  // Rebuild training covariance
  K_loc = gp_exp_quad_cov(X_loc_train, alpha_loc, rho_loc);
  K_inc = square(alpha_inc) * (x_inc_train * x_inc_train');
  K = K_loc + K_inc;

  for (i in 1:N_train) {
    K[i, i] = K[i, i] + square(sigma) + jitter;
  }

  // Cross covariance: train-test
  K_s_loc = gp_exp_quad_cov(X_loc_train, X_loc_test, alpha_loc, rho_loc);
  K_s_inc = square(alpha_inc) * (x_inc_train * x_inc_test');
  K_s = K_s_loc + K_s_inc;

  // Test-test covariance
  K_ss_loc = gp_exp_quad_cov(X_loc_test, alpha_loc, rho_loc);
  K_ss_inc = square(alpha_inc) * (x_inc_test * x_inc_test');
  K_ss = K_ss_loc + K_ss_inc;

  for (i in 1:N_test) {
    K_ss[i, i] = K_ss[i, i] + jitter;
  }

  // Posterior predictive mean and covariance
  alpha_vec = mdivide_left_spd(K, y_train);
  f_test_mean = K_s' * alpha_vec;
  K_pred = K_ss - K_s' * mdivide_left_spd(K, K_s);

  for (i in 1:N_test) {
    K_pred[i, i] = K_pred[i, i] + square(sigma) + jitter;
  }

  y_test_rep = multi_normal_rng(f_test_mean, K_pred);
}



