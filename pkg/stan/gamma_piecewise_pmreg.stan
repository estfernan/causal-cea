/*
 * Gamma-Piecewise Proportional-Means Modeling
 */

data
{
  int<lower=1> n;                     // sample size
  int<lower=1> p;                     // number of predictors
  int<lower=1> K;                     // number of cuts in the time-axis
  vector<lower=0>[n] y;               // response vector
  vector<lower=0>[n] t;               // observed times
  matrix[n, p] X;                     // design matrix
  array[n] int<lower=1,upper=K> idx;  // indices for assigned time-axis interval
}

parameters
{
  real mu_0;                 // location parameter for log-intensity rates
  vector<lower=0>[K] sigma;  // standard deviation for log-intensity rates
  vector[K] eps;             // error terms for log-intensity rates
  vector[p] beta;            // regression coefficients
  real<lower=0> tau;         // variance parameter
}

transformed parameters
{
  vector[K] log_mu = mu_0 + sigma .* eps;  // log-intensity rates
  vector<lower=0>[K] mu = exp(log_mu);    // intensity rates
}

model
{
  // prior distributions
  mu_0 ~ normal(0, 1);
  sigma ~ gamma(1, 1);
  eps ~ normal(0, 1);
  beta ~ normal(0, 1);
  tau ~ gamma(1, 1);

  // temporary variables
  vector[n] lambda = mu[idx] .* exp(X * beta);

  // observed data likelihood
  y ~ gamma(lambda^2 ./ tau, lambda ./ tau);
}
