/*
 * Gamma-Gompertz Proportional-Means Modeling
 */

functions
{
  // baseline mean (intensity) function
  vector mu_0(vector t, real eta, real xi)
  {
    return eta * exp(xi * t);
  }
}

data
{
  int<lower=1> n;        // sample size
  int<lower=1> p;        // number of predictors
  vector<lower=0>[n] y;  // response vector
  vector<lower=0>[n] t;  // observed times
  matrix[n, p] X;        // design matrix
}

parameters
{
  real<lower=0> alpha;  // Gamma distribution: shape parameter
  real<lower=0> eta;    // mean function: shape parameter
  real<lower=0> xi;     // mean function: scale parameter
  vector[p] beta;       // regression coefficients
}

model
{
  // prior distributions for model parameters
  alpha ~ lognormal(0, 1);
  eta ~ lognormal(0, 1);
  xi ~ lognormal(0, 1);
  beta ~ normal(0, 1);

  // temporary variables
  vector[n] mu = mu_0(t, eta, xi) .* exp(X * beta);

  // observed data likelihood
  y ~ gamma(mu^2 ./ alpha, mu ./ alpha);
}
