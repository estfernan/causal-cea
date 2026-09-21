/*
 * Piecewise Exponential Modeling
 */

data
{
  int<lower=0> n;                       // sample size
  int<lower=0> p;                       // number of predictors
  int<lower=0> K;                       // number of time-axis intervals
  vector[n] H;                          // duration of each interval
  matrix[n, p] X;                       // design matrix
  array[n] int<lower=0,upper=1> delta;  // interval-specific failure indicator
  array[n] int<lower=1,upper=K> idx;    // interval assignment indicator
}

parameters
{
  vector[K] eps;                // error terms for log-hazard
  real mu;                      // location parameter for log-hazard
  vector<lower=0>[K] sigma;     // standard deviation for error terms
  real<lower=-1, upper=1> rho;  // correlation between sequential hazards
  vector[p] beta;               // regression coefficients
}

transformed parameters
{
  vector[K] log_lambda;
  log_lambda[1] = mu + sigma[1] * eps[1];

  // gAR(1) prior over the log-hazard rates
  for (k in 2:K)
  {
    log_lambda[k] = mu * (1 - rho) + rho * log_lambda[k - 1] + sigma[k] * eps[k];
  }

  vector[n] log_hazard = log_lambda[idx] + X * beta + log(H);
}

model
{
  // prior distributions
  eps ~ normal(0, 1);
  mu ~ normal(0, 1);
  sigma ~ lognormal(0, 1);
  rho ~ uniform(-1, 1);
  beta ~ normal(0, 1);

  // observed data likelihood
  delta ~ poisson_log(log_hazard);
}

generated quantities
{
  vector<lower=0>[K] lambda = exp(log_lambda);
}
