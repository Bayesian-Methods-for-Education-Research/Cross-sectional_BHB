#Stan files for multilevel study
#bdb.stan #stands for Bayesian Dynamic Borrowing stan codes
data {
  int<lower = 1> N;                         // number of students
  int<lower = 1> NN;                        // number of students in the current cycle
  int<lower = 1> M;                         // number of covariates
  int<lower = 1> MM;                        // number of covariates with random effects
  int<lower = 1> S;                         // number of schools
  int<lower = 1> C;                         // number of cycles
  matrix[N, M] x;                           // independent variables
  vector[N] y;                              // depenendent variables
  int<lower = 1, upper = S> sch[N];         // schools
  int<lower = 1, upper = C> cycle[N];       // cycles of students
  int<lower = 1, upper = C> cycle_sch[S];   // cycles of schools
  real<lower = 0> nu1;                      // parameter for tau_beta ~ inv_gamma(nu1, nu2)
  real<lower = 0> nu2;                      // parameter for tau_beta ~ inv_gamma(nu1, nu2)
  real<lower = 0> nu;                       // parameter for prec_u ~ wishart(nu, nu * sigma)
}

transformed data {
  vector[M] mu_x;
  vector[M] sd_x;
  matrix[N, M] x_std;
  
  real mu_y = mean(y);
  real sd_y = sd(y);
  vector[N] y_std = (y - mu_y) / sd_y;
  
  // x[, 1] is the intercept
  x_std[, 1] = x[, 1];
  for (m in 2:M) {
    mu_x[m] = mean(x[, m]);
    sd_x[m] = sd(x[, m]);
    x_std[, m] = (x[, m] - mu_x[m]) / sd_x[m];
  }
}

parameters {
  // fixed effects
  vector[M] beta_std[C];
  vector[M] mu_std;
  vector<lower = 0>[M] tau2_beta_std;
  
  // random effects
  vector[MM] u_std[S];
  cholesky_factor_corr[MM] L_omega_u;
  vector<lower = 0>[MM] tau_u_std;
  cov_matrix[MM] prec_u_std[C];
  
  real<lower = 0> sigma_e_std;
}

model {
  for (c in 1:C)
    beta_std[c] ~ normal(mu_std, sqrt(tau2_beta_std));
  mu_std ~ normal(0, 10);
  for (m in 1:M)
    tau2_beta_std[m] ~ inv_gamma(nu1, nu2);
  
  for (s in 1:S)
    u_std[s] ~ multi_normal_prec(rep_vector(0, MM), prec_u_std[cycle_sch[s]]);
  for (c in 1:C)
    prec_u_std[c] ~ wishart(nu, nu * inverse_spd(quad_form_diag(L_omega_u, tau_u_std)));
  L_omega_u ~ lkj_corr_cholesky(3);
  tau_u_std ~ cauchy(0, 1);
  
  for (n in 1:N)
    y_std[n] ~ normal(x_std[n, ] * beta_std[cycle[n]] + x_std[n, 1:MM] * u_std[sch[n]], sigma_e_std);
  sigma_e_std ~ cauchy(0, 1);
}

generated quantities {
  vector[M] beta;
  real<lower = 0> sigma_e = sigma_e_std * sd_y;
  vector[NN] log_lik;
  
  cov_matrix[MM] sigma_u[C];
  
  for (c in 1:C)
    sigma_u[c] = inverse_spd(prec_u_std[c]) * (sd_y ^ 2);
  
  beta[1] = sd_y * beta_std[1, 1] + mu_y;
  for (m in 2:M) {
    beta[m] = sd_y / sd_x[m] * beta_std[1, m];
    beta[1] -= beta[m] * mu_x[m];
  }
  
  for (n in 1:NN)
    log_lik[n] = normal_lpdf(y[n] | x[n, ] * beta + x[n, 1:MM] * (u_std[sch[n]] * sd_y), sigma_e);
}


#inf.stan # stands for informative prior stan codes

data {
  int<lower = 1> N;                         // number of students
  int<lower = 1> NN;                        // number of students in the current cycle
  int<lower = 1> M;                         // number of covariates
  int<lower = 1> MM;                        // number of covariates with random effects
  int<lower = 1> S;                         // number of schools
  int<lower = 1> C;                         // number of cycles
  matrix[N, M] x;                           // independent variables
  vector[N] y;                              // depenendent variables
  int<lower = 1, upper = S> sch[N];         // schools
  int<lower = 1, upper = C> cycle[N];       // cycles
  vector[M] mu;                             // prior mean of beta
  vector<lower = 0>[M] sigma;               // prior sd of beta
}

transformed data {
  vector[M] mu_x;
  vector[M] sd_x;
  matrix[N, M] x_std;
  
  real mu_y = mean(y);
  real sd_y = sd(y);
  vector[N] y_std = (y - mu_y) / sd_y;
  
  vector[M] mu_std;
  vector<lower = 0>[M] sigma_std;
  
  // x[, 1] is the intercept
  x_std[, 1] = x[, 1];
  for (m in 2:M) {
    mu_x[m] = mean(x[, m]);
    sd_x[m] = sd(x[, m]);
    x_std[, m] = (x[, m] - mu_x[m]) / sd_x[m];
  }
  
  mu_std[1] = mu[1] - mu_y;
  sigma_std[1] = sigma[1] ^ 2;
  for (m in 2:M) {
    real r = sd_x[m] / sd_y;
    mu_std[m] = r * mu[m];
    sigma_std[m] = r * sigma[m];
    mu_std[1] += mu_x[m] * mu[m];
    sigma_std[1] += (mu_x[m] * sigma[m]) ^ 2;
  }
  mu_std[1] /= sd_y;
  sigma_std[1] = sqrt(sigma_std[1]) / sd_y;
}

parameters {
  // fixed effects
  vector[M] beta_std;
  
  // random effects
  vector[MM] u_std[S];
  cholesky_factor_corr[MM] L_omega_u;
  vector<lower = 0>[MM] tau_u_std;
  
  real<lower = 0> sigma_e_std;
}

model {
  beta_std ~ normal(mu_std, sigma_std);
  
  for (s in 1:S)
    u_std[s] ~ multi_normal_cholesky(rep_vector(0, MM), diag_pre_multiply(tau_u_std, L_omega_u));
  L_omega_u ~ lkj_corr_cholesky(3);
  tau_u_std ~ cauchy(0, 1);
  
  for (n in 1:N)
    y_std[n] ~ normal(x_std[n, ] * beta_std + x_std[n, 1:MM] * u_std[sch[n]], sigma_e_std);
  sigma_e_std ~ cauchy(0, 1);
}

generated quantities {
  vector[M] beta;
  real<lower = 0> sigma_e = sigma_e_std * sd_y;
  vector[NN] log_lik;
  
  vector[MM] u[S];
  vector<lower = 0>[MM] tau_u = tau_u_std * sd_y;
  corr_matrix[MM] omega_u = multiply_lower_tri_self_transpose(L_omega_u);
  cov_matrix[MM] sigma_u = quad_form_diag(omega_u, tau_u);
  
  for (s in 1:S)
    u[s] = u_std[s] * sd_y;
  
  beta[1] = sd_y * beta_std[1] + mu_y;
  for (m in 2:M) {
    beta[m] = sd_y / sd_x[m] * beta_std[m];
    beta[1] -= beta[m] * mu_x[m];
  }
  
  for (n in 1:NN)
    log_lik[n] = normal_lpdf(y[n] | x[n, ] * beta + x[n, 1:MM] * u[sch[n]], sigma_e);
}

#pp.stan # stands for Power Priors stan codes
data {
  int<lower = 1> N;                         // number of students
  int<lower = 1> NN;                        // number of students in the current cycle
  int<lower = 1> M;                         // number of covariates
  int<lower = 1> MM;                        // number of covariates with random effects
  int<lower = 1> S;                         // number of schools
  int<lower = 1> C;                         // number of cycles
  matrix[N, M] x;                           // independent variables
  vector[N] y;                              // depenendent variables
  int<lower = 1, upper = S> sch[N];         // schools
  int<lower = 1, upper = C> cycle[N];       // cycles
  vector<lower = 0, upper = 1>[C] a;        // weight for each cycle
}

transformed data {
  vector[M] mu_x;
  vector[M] sd_x;
  matrix[N, M] x_std;
  
  real mu_y = mean(y);
  real sd_y = sd(y);
  vector[N] y_std = (y - mu_y) / sd_y;
  
  // x[, 1] is the intercept
  x_std[, 1] = x[, 1];
  for (m in 2:M) {
    mu_x[m] = mean(x[, m]);
    sd_x[m] = sd(x[, m]);
    x_std[, m] = (x[, m] - mu_x[m]) / sd_x[m];
  }
}

parameters {
  // fixed effects
  vector[M] beta_std;
  
  // random effects
  vector[MM] u_std[S];
  cholesky_factor_corr[MM] L_omega_u;
  vector<lower = 0>[MM] tau_u_std;
  
  real<lower = 0> sigma_e_std;
}

model {
  beta_std ~ normal(0, 10);
  
  for (s in 1:S)
    u_std[s] ~ multi_normal_cholesky(rep_vector(0, MM), diag_pre_multiply(tau_u_std, L_omega_u));
  L_omega_u ~ lkj_corr_cholesky(3);
  tau_u_std ~ cauchy(0, 1);
  
  for (n in 1:N)
    target += a[cycle[n]] * normal_lpdf(y_std[n] | x_std[n, ] * beta_std + x_std[n, 1:MM] * u_std[sch[n]],
                                        sigma_e_std);
    sigma_e_std ~ cauchy(0, 1);
}

generated quantities {
  vector[M] beta;
  real<lower = 0> sigma_e = sigma_e_std * sd_y;
  vector[NN] log_lik;
  
  vector<lower = 0>[MM] tau_u = tau_u_std * sd_y;
  corr_matrix[MM] omega_u = multiply_lower_tri_self_transpose(L_omega_u);
  cov_matrix[MM] sigma_u = quad_form_diag(omega_u, tau_u);
  
  beta[1] = sd_y * beta_std[1] + mu_y;
  for (m in 2:M) {
    beta[m] = sd_y / sd_x[m] * beta_std[m];
    beta[1] -= beta[m] * mu_x[m];
  }
  
  for (n in 1:NN)
    log_lik[n] = normal_lpdf(y[n] | x[n, ] * beta + x[n, 1:MM] * (u_std[sch[n]] * sd_y), sigma_e);
}

#reg.stan # noninformative and pooling analyses stan codes
data {
  int<lower = 1> N;                         // number of students
  int<lower = 1> NN;                        // number of students in the current cycle
  int<lower = 1> M;                         // number of covariates
  int<lower = 1> MM;                        // number of covariates with random effects
  int<lower = 1> S;                         // number of schools
  int<lower = 1> C;                         // number of cycles
  matrix[N, M] x;                           // independent variables
  vector[N] y;                              // depenendent variables
  int<lower = 1, upper = S> sch[N];         // schools
  int<lower = 1, upper = C> cycle[N];       // cycles
}

transformed data {
  vector[M] mu_x;
  vector[M] sd_x;
  matrix[N, M] x_std;
  
  real mu_y = mean(y);
  real sd_y = sd(y);
  vector[N] y_std = (y - mu_y) / sd_y;
  
  // x[, 1] is the intercept
  x_std[, 1] = x[, 1];
  for (m in 2:M) {
    mu_x[m] = mean(x[, m]);
    sd_x[m] = sd(x[, m]);
    x_std[, m] = (x[, m] - mu_x[m]) / sd_x[m];
  }
}

parameters {
  // fixed effects
  vector[M] beta_std;
  
  // random effects
  vector[MM] u_std[S];
  cholesky_factor_corr[MM] L_omega_u;
  vector<lower = 0>[MM] tau_u_std;
  
  real<lower = 0> sigma_e_std;
}

model {
  beta_std ~ normal(0, 10);
  
  for (s in 1:S)
    u_std[s] ~ multi_normal_cholesky(rep_vector(0, MM), diag_pre_multiply(tau_u_std, L_omega_u));
  L_omega_u ~ lkj_corr_cholesky(3);
  tau_u_std ~ cauchy(0, 1);
  
  for (n in 1:N)
    y_std[n] ~ normal(x_std[n, ] * beta_std + x_std[n, 1:MM] * u_std[sch[n]], sigma_e_std);
  sigma_e_std ~ cauchy(0, 1);
}

generated quantities {
  vector[M] beta;
  real<lower = 0> sigma_e = sigma_e_std * sd_y;
  vector[NN] log_lik;
  
  vector<lower = 0>[MM] tau_u = tau_u_std * sd_y;
  corr_matrix[MM] omega_u = multiply_lower_tri_self_transpose(L_omega_u);
  cov_matrix[MM] sigma_u = quad_form_diag(omega_u, tau_u);
  
  beta[1] = sd_y * beta_std[1] + mu_y;
  for (m in 2:M) {
    beta[m] = sd_y / sd_x[m] * beta_std[m];
    beta[1] -= beta[m] * mu_x[m];
  }
  
  for (n in 1:NN)
    log_lik[n] = normal_lpdf(y[n] | x[n, ] * beta + x[n, 1:MM] * (u_std[sch[n]] * sd_y), sigma_e);
}


