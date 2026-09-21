/*
 * Bernoulli Modeling
 */

data
{
  int<lower=1> n;                    // sample size
  array[n] int<lower=0, upper=1> y;  // response vector
}

parameters
{
  real theta;  // log-odds
}

transformed parameters
{
  real<lower=0, upper=1> p = inv_logit(theta);  // success probability
}

model
{
  theta ~ normal(0, 1);
  y ~ bernoulli(p);
}
