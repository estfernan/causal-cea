/*
 * Categorical Generalized Linear Modeling
 */

data
{
  int<lower=1> n;                    // sample size
  int<lower=1> p;                    // number of predictors
  int<lower=1> K;                    // number of categories
  matrix[n, p] X;                    // design matrix
  array[n] int<lower=1, upper=K> y;  // response vector
}

transformed data
{
  vector[K] alpha = rep_vector(0, K);  // zero-vector for intercepts
  vector[p] empty = rep_vector(0, p);  // zero-vector for reference category
}

parameters
{
  matrix[p, K - 1] beta;  // regression coefficients
}

transformed parameters
{
  matrix[p, K] beta_all = append_col(beta, empty);
}

model
{
  to_vector(beta) ~ normal(0, 3);                 // prior for category-specific coefficients
  y ~ categorical_logit_glm(X, alpha, beta_all);  // observed data likelihood
}
