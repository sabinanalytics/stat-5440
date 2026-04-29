data {
  int<lower=1> N_train;
  int<lower=1> N_test;
  int<lower=1> D;

  // [Fix] Changed row_vector to vector. 
  // gp_exp_quad_cov strictly requires array[] vector.
  array[N_train] vector[D] X_train;
  array[N_test] vector[D] X_test;

  vector[N_train] y_train;
}

transformed data {
  real jitter = 1e-5;   // small diagonal term for numerical stability
}

parameters {
  real<lower=1e-6> alpha;   // kernel amplitude
  real<lower=1e-6> rho;     // length scale
  real<lower=1e-6> sigma;   // observation noise
}

model {
  matrix[N_train, N_train] K;
  matrix[N_train, N_train] L_K;

  // Weakly informative priors for positive GP hyperparameters
  alpha ~ lognormal(0, 0.5);
  rho ~ lognormal(0, 0.5);
  sigma ~ lognormal(-1, 0.5);

  // Use robust built-in function to avoid 0 * Inf gradient NaNs
  K = gp_exp_quad_cov(X_train, alpha, rho);

  for (i in 1:N_train) {
    K[i, i] = K[i, i] + square(sigma) + jitter;
  }

  L_K = cholesky_decompose(K);

  y_train ~ multi_normal_cholesky(rep_vector(0, N_train), L_K);
}

generated quantities {
  matrix[N_train, N_train] K;
  matrix[N_train, N_test] K_s;
  matrix[N_test, N_test] K_ss;
  matrix[N_test, N_test] K_pred;
  vector[N_train] alpha_vec;
  vector[N_test] f_test_mean;
  vector[N_test] y_test_rep;

  // Rebuild K using the built-in function
  K = gp_exp_quad_cov(X_train, alpha, rho);
  for (i in 1:N_train) {
    K[i, i] = K[i, i] + square(sigma) + jitter;
  }

  // Cross-covariance and test covariance built instantly
  K_s = gp_exp_quad_cov(X_train, X_test, alpha, rho);
  K_ss = gp_exp_quad_cov(X_test, alpha, rho);
  
  for (i in 1:N_test) {
    K_ss[i, i] = K_ss[i, i] + jitter;
  }

  // Posterior predictive mean and covariance
  alpha_vec = mdivide_left_spd(K, y_train);
  f_test_mean = K_s' * alpha_vec;
  K_pred = K_ss - K_s' * mdivide_left_spd(K, K_s);

  // Add observation noise and a small jitter term to stabilize predictive covariance
  for (i in 1:N_test) {
    K_pred[i, i] = K_pred[i, i] + square(sigma) + jitter;
  }

  // Draw posterior predictive samples for observed y values
  y_test_rep = multi_normal_rng(f_test_mean, K_pred);
}


