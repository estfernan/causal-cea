/*
 * Gaussian Modeling
 */

data
{
  int<lower=1> n;  // sample size
  vector[n] y;     // response vector
}

parameters
{
  real mu;              // location parameter
  real<lower=0> sigma;  // standard deviation parameter
}

model
{
  // prior distributions
  mu ~ normal(0, 3);
  sigma ~ gamma(1, 1);

  // observed data likelihood
  y ~ normal(mu, sigma);
}
