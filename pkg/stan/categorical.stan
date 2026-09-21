/*
 * Categorical Modeling
 */

data
{
  int<lower=1> n;                    // sample size
  int<lower=1> K;                    // number of categories
  array[n] int<lower=1, upper=K> y;  // response vector
}

parameters
{
  vector[K] beta;  // log-odds for each category
}

transformed parameters
{
  simplex[K] theta = softmax(beta);  // event probabilities
}

model
{
  beta ~ normal(0, 1);     // prior for log-odds
  y ~ categorical(theta);  // observed data likelihood
}
