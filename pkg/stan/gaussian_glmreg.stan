/*
 * Gaussian Generalized Linear Modeling
 */

data
{
  int<lower=1> n;  // sample size
  int<lower=1> p;  // number of predictors
  vector[n] y;     // response vector
  matrix[n, p] X;  // design matrix
}

parameters
{
  vector[p] beta;       // regression coefficients
  real<lower=0> sigma;  // standard deviation parameter
}

transformed parameters
{
  real<lower=0> tau = 1 / sigma^2;  // precision parameter
}

model
{
  // prior distributions
  beta ~ normal(0, 3);
  sigma ~ gamma(1, 1);

  // observed data likelihood
  y ~ normal_id_glm(X, 0, beta, sigma);
}
