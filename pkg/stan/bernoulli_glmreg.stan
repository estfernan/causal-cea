/*
 * Bernoulli Generalized Linear Modeling
 */

data
{
  int<lower=1> n;                    // sample size
  int<lower=1> p;                    // number of predictors
  matrix[n, p] X;                    // design matrix
  array[n] int<lower=0, upper=1> y;  // response vector
}

parameters
{
  vector[p] beta;  // regression coefficients
}

model
{
  beta ~ normal(0, 3);
  y ~ bernoulli_logit_glm(X, 0, beta);
}
