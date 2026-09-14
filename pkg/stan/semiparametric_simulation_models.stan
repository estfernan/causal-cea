/*
 * Semiparametric models for simulated insurance claims data.
 */

functions {
  vector gAR1_path(int J, vector eps, real nu, real rho, real tau) {
    // First-order autoregressive prior specification.
    vector[J] x;

    x[1] = nu + tau * eps[1];

    for (j in 2:J) {
      x[j] = nu * (1 - rho) + rho * x[j - 1] + tau * eps[j];
    }

    return x;
  }

  vector log_poisson_mean(vector loghaz, array[] int interval, vector duration, matrix S, vector psi) {
    // Log-mean for piecewise Poisson.
    return loghaz[interval] + S * psi + log(duration);
  }

  vector piecewise_gamma_mean(vector logmean, array[] int t, vector A, matrix H, real betaA, vector betaH) {
    // Mean specification for piecewise proportional-means Gamma cost model.
    return exp(logmean[t] + A * betaA + H * betaH);
  }
}

data {
  // core settings
  int<lower=1> n;   // number of patients
  int<lower=1> K;   // number of follow-up claims
  int<lower=1> J0;  // number of baseline time intervals
  int<lower=1> J;   // number of follow-up time intervals

  // baseline observations
  vector[n] U01;                       // continuous covariate
  vector[n] U02;                       // continuous covariate
  array[n] int<lower=0, upper=1> U03;  // binary covariate

  // baseline time-to-event observations
  int<lower=1> N0;                      // number of baseline partitions
  int<lower=1> P0;                      // number of baseline state covariates
  matrix[N0, P0] S0;                    // covariate matrix
  array[N0] int<lower=0, upper=1> d00;  // death-specific indicators
  array[N0] int<lower=0, upper=1> d01;  // visit-specific indicators
  vector<lower=0>[N0] D0;               // interval durations
  array[N0] int<lower=1, upper=J0> I0;   // interval IDs

  // follow-up time-to-event observations
  int<lower=1> N;                     // number of follow-up partitions
  int<lower=1> P;                     // number of follow-up state covariates
  matrix[N, P] S;                     // covariate matrix
  array[N] int<lower=0, upper=1> d0;  // death-specific indicators
  array[N] int<lower=0, upper=1> d1;  // visit-specific indicators
  vector<lower=0>[N] D;               // interval durations
  array[N] int<lower=1, upper=J> I;   // interval IDs

  // baseline covariate observations
  int<lower=1> R0;                     // number of baseline regression covariates
  matrix[n, R0] X0;                    // design matrix
  vector[n] L01;                       // continuous covariate
  array[n] int<lower=0, upper=1> L02;  // binary covariate

  // follow-up covariate observations
  int<lower=1> R;                     // number of follow-up regression covariates
  matrix[K, R] X;                     // design matrix
  vector[K] L1;                       // continuous covariate
  array[K] int<lower=0, upper=1> L2;  // binary covariate

  // baseline cost observations
  int<lower=1> C0;                    // number of baseline history covariates
  matrix[n, C0] H0;                   // history matrix
  array[n] int<lower=1, upper=J0> t0;  // interval indices
  vector[n] a0;                       // treatment decisions
  vector<lower=0>[n] y0;              // cost outcomes

  // follow-up cost observations
  int<lower=1> C;                    // number of follow-up history covariates
  matrix[K, C] H;                    // history matrix
  array[K] int<lower=1, upper=J> t;  // interval indices
  vector[K] a;                       // treatment decisions
  vector<lower=0>[K] y;              // cost outcomes
}

parameters {
  // baseline covariates
  real mu01;
  real<lower=0> sigma01;
  real mu02;
  real<lower=0> sigma02;
  real<lower=0, upper=1> p03;

  // log-gAR(1) baseline hazards: death and visit
  vector[P0] psi00;
  vector[J0] eps_h00;
  real nu_h00;
  real<lower=-1, upper=1> rho_h00;
  real<lower=0> tau_h00;
  vector[P0] psi01;
  vector[J0] eps_h01;
  real nu_h01;
  real<lower=-1, upper=1> rho_h01;
  real<lower=0> tau_h01;

  // log-gAR(1) follow-up hazards: death and visit
  vector[P] psi0;
  vector[J] eps_h0;
  real nu_h0;
  real<lower=-1, upper=1> rho_h0;
  real<lower=0> tau_h0;
  vector[P] psi1;
  vector[J] eps_h1;
  real nu_h1;
  real<lower=-1, upper=1> rho_h1;
  real<lower=0> tau_h1;

  // baseline covariate GLMs
  vector[R0] omegaL01;
  vector[R0] omegaL02;
  real<lower=0> sigmaL01;

  // follow-up covariate GLMs
  vector[R] omegaL1;
  vector[R] omegaL2;
  real<lower=0> sigmaL1;

  // gAR(1) baseline means
  vector[C0] betaH0;
  real betaA0;
  real<lower=0> gamma0;
  vector[J0] eps_m0;
  real nu_m0;
  real<lower=-1, upper=1> rho_m0;
  real<lower=0> tau_m0;

  // gAR(1) follow-up means
  vector[C] betaH;
  real betaA;
  real<lower=0> gamma;
  vector[J] eps_m;
  real nu_m;
  real<lower=-1, upper=1> rho_m;
  real<lower=0> tau_m;
}

transformed parameters {
  // log-gAR(1) hazards: death and visit
  vector[J0] log_h00 = gAR1_path(J0, eps_h00, nu_h00, rho_h00, tau_h00);
  vector[J0] log_h01 = gAR1_path(J0, eps_h01, nu_h01, rho_h01, tau_h01);
  vector[J] log_h0 = gAR1_path(J, eps_h0, nu_h0, rho_h0, tau_h0);
  vector[J] log_h1 = gAR1_path(J, eps_h1, nu_h1, rho_h1, tau_h1);

  // gAR(1) means
  vector[J0] log_m0 = gAR1_path(J0, eps_m0, nu_m0, rho_m0, tau_m0);
  vector[J] log_m = gAR1_path(J, eps_m, nu_m, rho_m, tau_m);
}

model {
  // weakly informative priors on baseline covariate parameters
  mu01 ~ normal(0, 1);
  sigma01 ~ normal(0, 1);
  mu02 ~ normal(0, 1);
  sigma02 ~ normal(0, 1);
  p03 ~ uniform(0, 1);

  // baseline models
  U01 ~ normal(mu01, sigma01);
  U02 ~ normal(mu02, sigma02);
  U03 ~ bernoulli(p03);

  // weakly informative priors on baseline time-to-event parameters
  psi00 ~ normal(0, 1);
  eps_h00 ~ std_normal();
  nu_h00 ~ normal(0, 1);
  rho_h00 ~ uniform(-1, 1);
  tau_h00 ~ normal(0, 1);
  psi01 ~ normal(0, 1);
  eps_h01 ~ std_normal();
  nu_h01 ~ normal(0, 1);
  rho_h01 ~ uniform(-1, 1);
  tau_h01 ~ normal(0, 1);

  // baseline time-to-event Poisson-based hazards
  vector[N0] eta00 = log_poisson_mean(log_h00, I0, D0, S0, psi00);
  vector[N0] eta01 = log_poisson_mean(log_h01, I0, D0, S0, psi01);

  // observed data likelihood: baseline time-to-events
  target += poisson_log_lpmf(d00 | eta00);
  target += poisson_log_lpmf(d01 | eta01);

  // weakly informative priors on follow-up time-to-event parameters
  psi0 ~ normal(0, 1);
  eps_h0 ~ std_normal();
  nu_h0 ~ normal(0, 1);
  rho_h0 ~ uniform(-1, 1);
  tau_h0 ~ normal(0, 1);
  psi1 ~ normal(0, 1);
  eps_h1 ~ std_normal();
  nu_h1 ~ normal(0, 1);
  rho_h1 ~ uniform(-1, 1);
  tau_h1 ~ normal(0, 1);

  // follow-up time-to-event Poisson-based hazards
  vector[N] eta0 = log_poisson_mean(log_h0, I, D, S, psi0);
  vector[N] eta1 = log_poisson_mean(log_h1, I, D, S, psi1);

  // observed data likelihood: follow-up time-to-events
  target += poisson_log_lpmf(d0 | eta0);
  target += poisson_log_lpmf(d1 | eta1);

  // weakly informative priors on covariate parameters
  omegaL01 ~ normal(0, 1);
  omegaL02 ~ normal(0, 1);
  sigmaL01 ~ normal(0, 1);
  omegaL1 ~ normal(0, 1);
  omegaL2 ~ normal(0, 1);
  sigmaL1 ~ normal(0, 1);

  // covariate GLMs
  L01 ~ normal_id_glm(X0, 0.0, omegaL01, sigmaL01);
  L02 ~ bernoulli_logit_glm(X0, 0.0, omegaL02);
  L1 ~ normal_id_glm(X, 0.0, omegaL1, sigmaL1);
  L2 ~ bernoulli_logit_glm(X, 0.0, omegaL2);

  // weakly informative priors on baseline cost parameters
  betaH0 ~ normal(0, 1);
  betaA0 ~ normal(0, 1);
  gamma0 ~ normal(0, 1);
  eps_m0 ~ std_normal();
  nu_m0 ~ normal(0, 1);
  rho_m0 ~ uniform(-1, 1);
  tau_m0 ~ normal(0, 1);

  // weakly informative priors on follow-up cost parameters
  betaH ~ normal(0, 1);
  betaA ~ normal(0, 1);
  gamma ~ normal(0, 1);
  eps_m ~ std_normal();
  nu_m ~ normal(0, 1);
  rho_m ~ uniform(-1, 1);
  tau_m ~ normal(0, 1);

  // piecewise cost means
  vector[n] mu0 = piecewise_gamma_mean(log_m0, t0, a0, H0, betaA0, betaH0);
  vector[K] mu = piecewise_gamma_mean(log_m, t, a, H, betaA, betaH);

  // observed data likelihood: cost outcomes
  y0 ~ gamma((mu0 .* mu0) / gamma0,  mu0 / gamma0);
  y ~ gamma((mu .* mu) / gamma,  mu / gamma);
}

generated quantities {
  // inverse log-transform unconstrained variables
  vector<lower=0>[J0] h00 = exp(log_h00);
  vector<lower=0>[J0] h01 = exp(log_h01);
  vector<lower=0>[J] h0 = exp(log_h0);
  vector<lower=0>[J] h1 = exp(log_h1);
  vector<lower=0>[J0] m0 = exp(log_m0);
  vector<lower=0>[J] m = exp(log_m);
}
