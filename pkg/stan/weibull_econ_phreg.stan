/*
 * Weibull Proportional-Hazards Modeling
 */

data
{
  int<lower=1> n;   // sample size
  int<lower=1> p;   // number of predictors
  vector[n] t;      // failure times
  vector[n] delta;  // event indicator (1 - observed, 0 - censored)
  matrix[n, p] X;   // design matrix
}

parameters
{
  real<lower=0> lambda;  // alternative scale parameter
  real<lower=0> nu;      // shape parameter
  vector[p] beta;      // scaled regression coefficients
}

model
{
  // prior distributions for model parameters
  nu ~ lognormal(0, 1);
  lambda ~ lognormal(0, 1);
  beta ~ normal(0, 1);

  // temporary variables
  real obs;
  real b;
  real log_pdf;
  real log_surv;
  vector[n] sigma = pow(lambda * exp(X * beta), -1 / nu);

  // increment log-likelihood
  for (i in 1:n)
  {
    obs = delta[i];
    b = sigma[i];

    log_pdf = weibull_lpdf(t[i] | nu, b);
    log_surv = weibull_lccdf(t[i] | nu, b);

    target += obs * log_pdf + (1 - obs) * log_surv;
  }
}
