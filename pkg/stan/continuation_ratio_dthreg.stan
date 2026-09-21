/*
 * Continuation Ratio Regression Modeling
 */

data
{
  int<lower=1> n;                    // sample size
  int<lower=1> p;                    // number of predictors
  int<lower=1> K;                     // number of indices
  matrix[n, p] X;                    // design matrix
  array[n] int<lower=0, upper=1> y;  // response vector
  array[n] int<lower=1,upper=K> idx;  // indices for observations
}

parameters
{
  real gamma_0;              // location parameter for logit baseline hazards
  vector<lower=0>[K] sigma;  // standard deviation for logit baseline hazards
  vector[K] eps;             // error terms for logit baseline hazards
  vector[p] beta;            // regression coefficients
}

transformed parameters
{
  vector[K] gamma = gamma_0 + sigma .* eps;  // logit baseline hazards
}

model
{
  // prior distributions
  gamma_0 ~ normal(0, 1);
  sigma ~ gamma(1, 1);
  eps ~ normal(0, 1);
  beta ~ normal(0, 3);

  // temporary variables
  vector[n] alpha = gamma[idx] + X * beta;

  // observed data likelihood
  y ~ bernoulli_logit(alpha);
}
