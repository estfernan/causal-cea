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
  real<lower=0> lambda;  // inverse scale (or rate) parameter
  vector[p] beta;        // regression coefficients
}

model
{
  // prior distributions for model parameters
  lambda ~ normal(0, 1);
  beta ~ normal(0, 1);

  // temporary variables
  real obs;
  real b;
  real log_pdf;
  real log_surv;
  vector[n] Lambda = lambda * exp(X * beta);

  // increment log-likelihood
  for (i in 1:n)
  {
    obs = delta[i];
    b = Lambda[i];

    log_pdf = exponential_lpdf(t[i] | b);
    log_surv = exponential_lccdf(t[i] | b);

    target += obs * log_pdf + (1 - obs) * log_surv;
  }
}
