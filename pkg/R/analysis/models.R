#
# Statistical models using Stan.
#

#'
#' Bayesian Continuation Ratio Modeling
#'
bayesian_continuation_ratio <- function(formula, data, index_weights, method = "sample", chains = 1, parallel_chains = 1, iter = 1000, warmup = 1000, refresh = 100, max_treedepth = NULL, format = "draws_matrix", stan_name = CONTINUATION_RATIO, stan_path = STAN_PATH) {
  stan_model <- .init_cmdstanr(stan_name, stan_path)

  f_call <- match.call(expand.dots = FALSE)
  env <- parent.frame()

  mf <- .extract_model_frame(f_call, env)

  # run Stan's MCMC algorithm with extracted components
  fit <- .run_stan_algorithms(
    data = .extract_discrete_tte_objects(mf, index_weights),
    stan_model = stan_model,
    method = method,
    chains = chains,
    parallel_chains = parallel_chains,
    refresh = refresh,
    warmup = warmup, iter = iter,
    max_treedepth = max_treedepth
  )

  draws <- fit$draws(variables = c("gamma", "beta"), format = format)

  if (method == "optimize") {
    draws <- draws[1, , drop = TRUE]
  }

  return(draws)
}

# """
# Generalized linear models
# """

#'
#' Bayesian Bernoulli Generalized Linear Modeling
#'
bayesian_bernoulli_glm <- function(formula, data, method = "sample", chains = 1, parallel_chains = 1, iter = 1000, warmup = 1000, refresh = 100, format = "draws_matrix", stan_name = BERNOULLI_GLM, stan_path = STAN_PATH) {
  stan_model <- .init_cmdstanr(stan_name, stan_path)

  f_call <- match.call(expand.dots = FALSE)
  env <- parent.frame()

  mf <- .extract_model_frame(f_call, env)

  # run stan's MCMC algorithm with extracted components
  fit <- .run_stan_algorithms(
    data = .extract_regression_objects(mf),
    stan_model = stan_model,
    method = method,
    chains = chains,
    parallel_chains = parallel_chains,
    refresh = refresh,
    warmup = warmup, iter = iter
  )

  draws <- fit$draws(variables = "beta", format = format)

  if (method == "optimize") {
    draws <- draws[1, , drop = TRUE]
  }

  return(draws)
}

#'
#' Bayesian Categorical Generalized Linear Modeling
#'
bayesian_categorical_glm <- function(formula, data, method = "sample", chains = 1, parallel_chains = 1, iter = 1000, warmup = 1000, refresh = 100, format = "draws_matrix", stan_name = CATEGORICAL_GLM, stan_path = STAN_PATH) {
  stan_model <- .init_cmdstanr(stan_name, stan_path)

  f_call <- match.call(expand.dots = FALSE)
  env <- parent.frame()

  mf <- .extract_model_frame(f_call, env)

  # run stan's MCMC algorithm with extracted components
  fit <- .run_stan_algorithms(
    data = .extract_regression_objects(mf, categorical = TRUE),
    stan_model = stan_model,
    method = method,
    chains = chains,
    parallel_chains = parallel_chains,
    refresh = refresh,
    warmup = warmup, iter = iter
  )

  draws <- fit$draws(variables = "beta", format = format)

  if (method == "optimize") {
    draws <- draws[1, , drop = TRUE]
  }

  return(draws)
}

#'
#' Bayesian Gaussian Generalized Linear Modeling
#'
bayesian_gaussian_glm <- function(formula, data, method = "sample", chains = 1, parallel_chains = 1, iter = 1000, warmup = 1000, refresh = 100, format = "draws_matrix", stan_name = GAUSSIAN_GLM, stan_path = STAN_PATH) {
  stan_model <- .init_cmdstanr(stan_name, stan_path)

  f_call <- match.call(expand.dots = FALSE)
  env <- parent.frame()

  mf <- .extract_model_frame(f_call, env)

  # run stan's MCMC algorithm with extracted components
  fit <- .run_stan_algorithms(
    data = .extract_regression_objects(mf),
    stan_model = stan_model,
    method = method,
    chains = chains,
    parallel_chains = parallel_chains,
    refresh = refresh,
    warmup = warmup, iter = iter
  )

  draws <- fit$draws(variables = c("tau", "beta"), format = format)

  if (method == "optimize") {
    draws <- draws[1, , drop = TRUE]
  }

  return(draws)
}

# """
# Proportional-hazards regression models
# """

#'
#' Bayesian Exponential Proportional-Hazards Modeling
#'
bayesian_exponential_ph <- function(formula, data, method = "sample", chains = 1, parallel_chains = 1, iter = 1000, warmup = 1000, refresh = 100, max_treedepth = NULL, time_cuts = NULL_VAR, format = "draws_matrix", stan_name = EXPONENTIAL_PHREG, stan_path = STAN_PATH) {
  stan_model <- .init_cmdstanr(stan_name, stan_path)

  f_call <- match.call(expand.dots = FALSE)
  env <- parent.frame()

  mf <- .extract_model_frame(f_call, env)

  # run Stan's MCMC algorithm with extracted components
  fit <- .run_stan_algorithms(
    data = .extract_survival_objects(mf),
    stan_model = stan_model,
    method = method,
    chains = chains,
    refresh = refresh,
    warmup = warmup, iter = iter,
    max_treedepth = max_treedepth
  )

  draws <- fit$draws(variables = c("lambda", "beta"), format = "matrix")

  if (method == "optimize") {
    draws <- draws[1, , drop = TRUE]
  }

  return(draws)
}

#'
#' Bayesian Piecewise Exponential Modeling
#'
bayesian_piecewise_exponential <- function(formula, data, method = "sample", chains = 1, parallel_chains = 1, iter = 1000, warmup = 1000, refresh = 100, max_treedepth = NULL, time_cuts = NULL_VAR, format = "draws_matrix", stan_name = PIECEWISE_EXPONENTIAL_PHREG, stan_path = STAN_PATH) {
  stan_model <- .init_cmdstanr(stan_name, stan_path)

  f_call <- match.call(expand.dots = FALSE)
  env <- parent.frame()

  mf <- .extract_model_frame(f_call, env)

  # run Stan's MCMC algorithm with extracted components
  fit <- .run_stan_algorithms(
    data = .extract_split_survival_objects(mf, time_cuts),
    stan_model = stan_model,
    method = method,
    chains = chains,
    parallel_chains = parallel_chains,
    refresh = refresh,
    warmup = warmup, iter = iter,
    max_treedepth = max_treedepth
  )

  draws <- fit$draws(variables = c("lambda", "beta"), format = format)

  if (method == "optimize") {
    draws <- draws[1, , drop = TRUE]
  }

  return(draws)
}

#'
#' Bayesian Weibull Proportional-Hazards Modeling
#'
bayesian_weibull_ph <- function(formula, data, method = "sample", chains = 1, parallel_chains = 1, iter = 1000, warmup = 1000, refresh = 100, max_treedepth = NULL, time_cuts = NULL_VAR, format = "draws_matrix", stan_name = WEIBULL_PHREG, stan_path = STAN_PATH) {
  stan_model <- .init_cmdstanr(stan_name, stan_path)

  f_call <- match.call(expand.dots = FALSE)
  env <- parent.frame()

  mf <- .extract_model_frame(f_call, env)

  # run Stan's MCMC algorithm with extracted components
  fit <- .run_stan_algorithms(
    data = .extract_survival_objects(mf),
    stan_model = stan_model,
    method = method,
    chains = chains,
    parallel_chains = parallel_chains,
    refresh = refresh,
    warmup = warmup, iter = iter,
    max_treedepth = max_treedepth
  )

  draws <- fit$draws(variables = c("lambda", "nu", "beta"), format = format)

  if (method == "optimize") {
    draws <- draws[1, , drop = TRUE]
  }

  return(draws)
}

# """
# Proportional-mean regression models
# """

#'
#' Bayesian Gamma-Gompertz Proportional-Means Modeling
#'
bayesian_gamma_gompertz_pm <- function(formula, data, time_weights, method = "sample", chains = 1, parallel_chains = 1, iter = 1000, warmup = 1000, refresh = 100, max_treedepth = NULL, time_cuts = NULL_VAR, format = "draws_matrix", stan_name = GAMMA_GOMPERTZ_PMREG, stan_path = STAN_PATH) {
  stan_model <- .init_cmdstanr(stan_name, stan_path)

  f_call <- match.call(expand.dots = FALSE)
  env <- parent.frame()

  mf <- .extract_model_frame(f_call, env)

  # run Stan's MCMC algorithm with extracted components
  fit <- .run_stan_algorithms(
    data = .extract_joint_objects(mf, time_weights),
    stan_model = stan_model,
    method = method,
    chains = chains,
    parallel_chains = parallel_chains,
    refresh = refresh,
    warmup = warmup, iter = iter,
    max_treedepth = max_treedepth
  )

  draws <- fit$draws(variables = c("alpha", "eta", "xi", "beta"), format = format)

  if (method == "optimize") {
    draws <- draws[1, , drop = TRUE]
  }

  return(draws)
}

#'
#' Bayesian Gamma-Weibull Proportional-Means Modeling
#'
bayesian_gamma_weibull_pm <- function(formula, data, time_weights, method = "sample", iter = 1000, warmup = 1000, refresh = 100, time_cuts = NULL_VAR, stan_name = GAMMA_WEIBULL_PMREG, stan_path = STAN_PATH) {
  stan_model <- .init_cmdstanr(stan_name, stan_path)

  f_call <- match.call(expand.dots = FALSE)
  env <- parent.frame()

  mf <- .extract_model_frame(f_call, env)

  # run stan's MCMC algorithm with extracted components
  fit <- .run_stan_algorithms(
    data = .extract_joint_objects(mf, time_weights),
    stan_model = stan_model,
    method = method,
    chains = 1L,
    refresh = refresh,
    warmup = warmup, iter = iter
  )

  draws <- fit$draws(variables = c("alpha", "eta", "xi", "beta"), format = "matrix")

  if (method == "optimize") {
    draws <- draws[1, , drop = TRUE]
  }

  return(draws)
}

#'
#' Bayesian Gamma-Piecewise Proportional-Means Modeling
#'
bayesian_gamma_piecewise_pm <- function(formula, data, time_weights, method = "sample", chains = 1, parallel_chains = 1, iter = 1000, warmup = 1000, refresh = 100, max_treedepth = NULL, time_cuts = NULL_VAR, format = "draws_matrix", stan_name = GAMMA_PIECEWISE_PMREG, stan_path = STAN_PATH) {
  stan_model <- .init_cmdstanr(stan_name, stan_path)

  f_call <- match.call(expand.dots = FALSE)
  env <- parent.frame()

  mf <- .extract_model_frame(f_call, env)

  # run stan's MCMC algorithm with extracted components
  fit <- .run_stan_algorithms(
    data = .extract_split_joint_objects(mf, time_weights, time_cuts),
    stan_model = stan_model,
    method = method,
    chains = chains,
    parallel_chains = parallel_chains,
    refresh = refresh,
    warmup = warmup, iter = iter,
    max_treedepth = max_treedepth
  )

  draws <- fit$draws(variables = c("tau", "mu", "beta"), format = format)

  if (method == "optimize") {
    draws <- draws[1, , drop = TRUE]
  }

  return(draws)
}

# """
# Simple univariable models
# """

#'
#' Bayesian Bernoulli Modeling
#'
bayesian_bernoulli_simple <- function(formula, data, method = "sample", chains = 1, parallel_chains = 1, iter = 1000, warmup = 1000, refresh = 100, format = "draws_matrix", stan_name = BERNOULLI_SIMPLE, stan_path = STAN_PATH) {
  stan_model <- .init_cmdstanr(stan_name, stan_path)

  f_call <- match.call(expand.dots = FALSE)
  env <- parent.frame()

  mf <- .extract_model_frame(f_call, env)

  # run stan's MCMC algorithm with extracted components
  fit <- .run_stan_algorithms(
    data = .extract_univariable_objects(mf),
    stan_model = stan_model,
    method = method,
    chains = chains,
    parallel_chains = parallel_chains,
    refresh = refresh,
    warmup = warmup, iter = iter
  )

  draws <- fit$draws(variables = "p", format = format)

  if (method == "optimize") {
    draws <- draws[1, , drop = TRUE]
  }

  return(draws)
}

#'
#' Bayesian Categorical Modeling
#'
bayesian_categorical_simple <- function(formula, data, method = "sample", chains = 1, parallel_chains = 1, iter = 1000, warmup = 1000, refresh = 100, format = "draws_matrix", stan_name = CATEGORICAL_SIMPLE, stan_path = STAN_PATH) {
  stan_model <- .init_cmdstanr(stan_name, stan_path)

  f_call <- match.call(expand.dots = FALSE)
  env <- parent.frame()

  mf <- .extract_model_frame(f_call, env)

  # run stan's MCMC algorithm with extracted components
  fit <- .run_stan_algorithms(
    data = .extract_univariable_objects(mf, categorical = TRUE),
    stan_model = stan_model,
    method = method,
    chains = chains,
    parallel_chains = parallel_chains,
    refresh = refresh,
    warmup = warmup, iter = iter
  )

  draws <- fit$draws(variables = "theta", format = format)

  if (method == "optimize") {
    draws <- draws[1, , drop = TRUE]
  }

  return(draws)
}

#'
#' Bayesian Gaussian Modeling
#'
bayesian_gaussian_simple <- function(formula, data, method = "sample", chains = 1, parallel_chains = 1, iter = 1000, warmup = 1000, refresh = 100, format = "draws_matrix", stan_name = GAUSSIAN_SIMPLE, stan_path = STAN_PATH) {
  stan_model <- .init_cmdstanr(stan_name, stan_path)

  f_call <- match.call(expand.dots = FALSE)
  env <- parent.frame()

  mf <- .extract_model_frame(f_call, env)

  # run stan's MCMC algorithm with extracted components
  fit <- .run_stan_algorithms(
    data = .extract_univariable_objects(mf),
    stan_model = stan_model,
    method = method,
    chains = chains,
    parallel_chains = parallel_chains,
    refresh = refresh,
    warmup = warmup, iter = iter
  )

  draws <- fit$draws(variables = c("mu", "sigma"), format = format)

  if (method == "optimize") {
    draws <- draws[1, , drop = TRUE]
  }

  return(draws)
}
