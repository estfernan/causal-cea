/*
 * Semiparametric models for simulated insurance claims data.
 */

functions {
  vector neg_cumhaz(vector t, matrix S, real lambda, vector psi) {
    // Negative cumulative hazard for the Exponential model.
    return -lambda * exp(S * psi) .* t;
  }

  vector log_haz(vector t, matrix S, real lambda, vector psi) {
    // Log-hazard for the Exponential model.
    return log(lambda) + S * psi;
  }

  vector weibull_gamma_mean(vector t, vector A, matrix H, real eta, real alpha, real betaA, vector betaH) {
    // Mean specification for Weibull-based proportional-means Gamma cost model.
    return eta * alpha * pow(t, alpha - 1) .* exp(A * betaA + H * betaH);
  }
}

data {
  // core settings
  int<lower=1> n;  // number of patients
  int<lower=1> K;  // number of follow-up claims

  // baseline observations
  vector[n] U01;                       // continuous covariate
  vector[n] U02;                       // continuous covariate
  array[n] int<lower=0, upper=1> U03;  // binary covariate

  // baseline time-to-event observations
  int<lower=1> P0;                      // number of baseline state covariates
  matrix[n, P0] S0;                    // covariate matrix
  array[n] int<lower=0, upper=1> d00;  // death-specific indicators
  array[n] int<lower=0, upper=1> d01;  // visit-specific indicators
  vector<lower=0>[n] t0;               // event times

  // follow-up time-to-event observations
  int<lower=1> P;                     // number of follow-up state covariates
  matrix[K, P] S;                     // covariate matrix
  array[K] int<lower=0, upper=1> d0;  // death-specific indicators
  array[K] int<lower=0, upper=1> d1;  // visit-specific indicators
  vector<lower=0>[K] t;               // event times

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
  int<lower=1> C0;        // number of baseline history covariates
  matrix[n, C0] H0;       // history matrix
  vector[n] a0;           // treatment decisions
  vector<lower=0>[n] y0;  // cost outcomes

  // follow-up cost observations
  int<lower=1> C;        // number of follow-up history covariates
  matrix[K, C] H;        // history matrix
  vector[K] a;           // treatment decisions
  vector<lower=0>[K] y;  // cost outcomes
}

parameters {
  // baseline covariates
  real mu01;
  real<lower=0> sigma01;
  real mu02;
  real<lower=0> sigma02;
  real<lower=0, upper=1> p03;

  // baseline hazards: death and visit
  vector[P0] psi00;
  real<lower=0> lambda00;
  vector[P0] psi01;
  real<lower=0> lambda01;

  // follow-up hazards: death and visit
  vector[P] psi0;
  real<lower=0> lambda0;
  vector[P] psi1;
  real<lower=0> lambda1;

  // baseline covariate GLMs
  vector[R0] omegaL01;
  vector[R0] omegaL02;
  real<lower=0> sigmaL01;

  // follow-up covariate GLMs
  vector[R] omegaL1;
  vector[R] omegaL2;
  real<lower=0> sigmaL1;

  // baseline means
  vector[C0] betaH0;
  real betaA0;
  real<lower=0> gamma0;
  real<lower=0> eta0;
  real<lower=0> alpha0;

  // follow-up means
  vector[C] betaH;
  real betaA;
  real<lower=0> gamma;
  real<lower=0> eta;
  real<lower=0> alpha;
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
  lambda00 ~ normal(0, 1);
  psi01 ~ normal(0, 1);
  lambda01 ~ normal(0, 1);

  // baseline time-to-event computations
  vector[n] neg_H00 = neg_cumhaz(t0, S0, lambda00, psi00);
  vector[n] log_h00 = log_haz(t0, S0, lambda00, psi00);
  vector[n] neg_H01 = neg_cumhaz(t0, S0, lambda01, psi01);
  vector[n] log_h01 = log_haz(t0, S0, lambda01, psi01);

  // observed data likelihood: baseline time-to-events
  for (i in 1:n) {
    if (d00[i] == 1) {
      target += log_h00[i];
    }

    if (d01[i] == 1) {
      target += log_h01[i];
    }

    target += neg_H00[i] + neg_H01[i];
  }

  // weakly informative priors on follow-up time-to-event parameters
  psi0 ~ normal(0, 1);
  lambda0 ~ normal(0, 1);
  psi1 ~ normal(0, 1);
  lambda1 ~ normal(0, 1);

  // follow-up time-to-event computations
  vector[K] neg_H0 = neg_cumhaz(t, S, lambda0, psi0);
  vector[K] log_h0 = log_haz(t, S, lambda0, psi0);
  vector[K] neg_H1 = neg_cumhaz(t, S, lambda1, psi1);
  vector[K] log_h1 = log_haz(t, S, lambda1, psi1);

  // observed data likelihood: follow-up time-to-events
  for (i in 1:K) {
    if (d0[i] == 1) {
      target += log_h0[i];
    }

    if (d1[i] == 1) {
      target += log_h1[i];
    }

    target += neg_H0[i] + neg_H1[i];
  }

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
  eta0 ~ normal(0, 1);
  alpha0 ~ normal(0, 1);

  // weakly informative priors on follow-up cost parameters
  betaH ~ normal(0, 1);
  betaA ~ normal(0, 1);
  gamma ~ normal(0, 1);
  eta ~ normal(0, 1);
  alpha ~ normal(0, 1);

  // piecewise cost means
  vector[n] mu0 = weibull_gamma_mean(t0, a0, H0, eta0, alpha0, betaA0, betaH0);
  vector[K] mu = weibull_gamma_mean(t, a, H, eta, alpha, betaA, betaH);

  // observed data likelihood: cost outcomes
  y0 ~ gamma((mu0 .* mu0) / gamma0,  mu0 / gamma0);
  y ~ gamma((mu .* mu) / gamma,  mu / gamma);
}
