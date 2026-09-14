/*
 * Discrete-time models for simulated insurance claims data.
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
}

data {
  // core settings
  int<lower=1> n;                       // number of patients
  int<lower=1> K0;                      // number of initial observations
  int<lower=1> K;                       // number of follow-up observations
  int<lower=1> J;                       // number of time intervals
  array[K] int<lower=1, upper=J> idxs;  // follow-up interval IDs

  // baseline observations
  vector[n] U01;                       // continuous covariate
  vector[n] U02;                       // continuous covariate
  array[n] int<lower=0, upper=1> U03;  // binary covariate

  // initial discrete time-to-event observations
  int<lower=1> P0;                     // number of baseline state covariates
  matrix[K0, P0] S0;                   // covariate matrix
  array[K0] int<lower=0, upper=1> d0;  // event indicators

  // follow-up discrete time-to-event observations
  int<lower=1> P;                    // number of follow-up state covariates
  matrix[K, P] S;                    // covariate matrix
  array[K] int<lower=0, upper=1> d;  // event indicators

  // baseline covariate observations
  int<lower=1> R0;                      // number of baseline regression covariates
  matrix[K0, R0] X0;                    // design matrix
  vector[K0] L01;                       // continuous covariate
  array[K0] int<lower=0, upper=1> L02;  // binary covariate

  // follow-up covariate observations
  int<lower=1> R;                     // number of follow-up regression covariates
  matrix[K, R] X;                     // design matrix
  vector[K] L1;                       // continuous covariate
  array[K] int<lower=0, upper=1> L2;  // binary covariate

  // baseline cost observations
  int<lower=1> C0;                          // number of baseline history covariates
  matrix[K0, C0] H0;                        // history matrix
  vector[K0] a0;                            // treatment decisions
  vector<lower=0>[K0] y0;                   // cost outcomes
  array[K0] int<lower=0, upper=1> is_pos0;  // non-zero cost indicator

  // follow-up cost observations
  int<lower=1> C;                         // number of follow-up history covariates
  matrix[K, C] H;                         // history matrix
  vector[K] a;                            // treatment decisions
  vector<lower=0>[K] y;                   // cost outcomes
  array[K] int<lower=0, upper=1> is_pos;  // non-zero cost indicator
}

parameters {
  // baseline covariates
  real mu01;
  real<lower=0> sigma01;
  real mu02;
  real<lower=0> sigma02;
  real<lower=0, upper=1> p03;

  // gAR(1) means logit-hazards
  vector[P0] psi0;
  vector[P] psi;
  vector[J] eps_h;
  real nu_h;
  real<lower=-1, upper=1> rho_h;
  real<lower=0> tau_h;

  // baseline covariate GLMs
  vector[R0] omegaL01;
  vector[R0] omegaL02;
  real<lower=0> sigmaL01;

  // follow-up covariate GLMs
  vector[R] omegaL1;
  vector[R] omegaL2;
  real<lower=0> sigmaL1;

  // hurdle component
  vector[C0] phiH0;
  real phiA0;
  vector[C] phiH;
  real phiA;
  vector[J] theta;

  // gAR(1) means
  vector[C0] betaH0;
  real betaA0;
  real<lower=0> gamma0;
  vector[C] betaH;
  real betaA;
  real<lower=0> gamma;
  vector[J] eps_m;
  real nu_m;
  real<lower=-1, upper=1> rho_m;
  real<lower=0> tau_m;
}

transformed parameters {
  // gAR(1) discrete intercepts
  vector[J] alpha_h = gAR1_path(J, eps_h, nu_h, rho_h, tau_h);
  vector[J] alpha_m = gAR1_path(J, eps_m, nu_m, rho_m, tau_m);
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

  // weakly informative priors on discrete time-to-event parameters
  psi0 ~ normal(0, 1);
  psi ~ normal(0, 1);
  eps_h ~ std_normal();
  nu_h ~ normal(0, 1);
  rho_h ~ uniform(-1, 1);
  tau_h ~ normal(0, 1);

  // observed data likelihood: discrete time-to-events
  d0 ~ bernoulli_logit_glm(S0, alpha_h[1], psi0);
  d ~ bernoulli_logit_glm(S, alpha_h[idxs], psi);

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

  // weakly informative priors on cost parameters for hurdle portion
  phiH0 ~ normal(0, 1);
  phiA0 ~ normal(0, 1);
  phiH ~ normal(0, 1);
  phiA ~ normal(0, 1);
  theta ~ normal(0, 1);

  // non-zero logit-scale probabilities
  vector[K0] eta0 = theta[1] + a0 * phiA0 + H0 * phiH0;
  vector[K] eta = theta[idxs]  + a * phiA + H * phiH;

  // observed data likelihood: cost outcomes - hurdle model component
  is_pos0 ~ bernoulli_logit(eta0);
  is_pos ~ bernoulli_logit(eta);

  // weakly informative priors on cost parameters
  betaH0 ~ normal(0, 1);
  betaA0 ~ normal(0, 1);
  gamma0 ~ normal(0, 1);
  betaH ~ normal(0, 1);
  betaA ~ normal(0, 1);
  gamma ~ normal(0, 1);
  eps_m ~ std_normal();
  nu_m ~ normal(0, 1);
  rho_m ~ uniform(-1, 1);
  tau_m ~ normal(0, 1);

  // cost means
  vector[K0] mu0 = exp(alpha_m[1] + a0 * betaA0 + H0 * betaH0);
  vector[K] mu = exp(alpha_m[idxs]  + a * betaA + H * betaH);

  // observed data likelihood: initial non-zero cost outcomes
  for (i in 1:K0) {
    if (is_pos0[i] == 1) {
      y0[i] ~ gamma((mu0[i] .* mu0[i]) / gamma0,  mu0[i] / gamma0);
    }
  }

  // observed data likelihood: follow-up non-zero cost outcomes
  for (i in 1:K) {
    if (is_pos[i] == 1) {
      y[i] ~ gamma((mu[i] .* mu[i]) / gamma,  mu[i] / gamma);
    }
  }
}
