#
# Predict methods for model fits.
#

#'
#' Predict Method for Continuation Ratio Model Fits
#'
predict_continuation_ratio <- function(pars, formula, data, index_weights) {
  f_call <- match.call(expand.dots = FALSE)
  env <- parent.frame()

  mf <- .extract_model_frame(f_call, env)

  object <- .extract_design_matrix(mf, intercept = TRUE)

  n <- object[["n"]]
  X <- object[["X"]]

  par_names <- names(pars)

  gamma <- pars[grep("gamma", par_names)]
  beta <- pars[grep("beta", par_names)]

  n_indices <- length(gamma)
  index_weights <- ifelse(index_weights > n_indices, n_indices, index_weights)

  alpha <- LaplacesDemon::invlogit(gamma[index_weights] + X %*% beta)

  return(LaplacesDemon::rbern(n, prob = alpha))
}

# """
# Generalized linear models
# """

#'
#' Predict Method for Bernoulli Generalized Linear Model Fits
#'
predict_bayesian_bernoulli_glm <- function(pars, formula, data) {
  f_call <- match.call(expand.dots = FALSE)
  env <- parent.frame()

  mf <- .extract_model_frame(f_call, env)

  object <- .extract_design_matrix(mf, intercept = TRUE)

  n <- object[["n"]]
  X <- object[["X"]]

  return(rbernglm(n, X, pars))
}

#'
#' Predict Method for Categorical Generalized Linear Model Fits
#'
predict_bayesian_categorical_glm <- function(pars, formula, data) {
  f_call <- match.call(expand.dots = FALSE)
  env <- parent.frame()

  mf <- .extract_model_frame(f_call, env)

  object <- .extract_design_matrix(mf, intercept = TRUE)

  n <- object[["n"]]
  p <- object[["p"]]
  X <- object[["X"]]

  beta <- cbind(matrix(pars, nrow = p, byrow = FALSE), 0)

  K <- ncol(beta)

  # TODO: make rcatglm() function
  # LaplacesDemon::rcat(n, p = softmax(X %*% beta))
  return(factor(apply(softmax(X %*% beta), 1, sample, x = K, size = 1, replace = FALSE), levels = 1:K))
}

#'
#' Predict Method for Gaussian Generalized Linear Model Fits
#'
predict_bayesian_gaussian_glm <- function(pars, formula, data) {
  f_call <- match.call(expand.dots = FALSE)
  env <- parent.frame()

  mf <- .extract_model_frame(f_call, env)

  object <- .extract_design_matrix(mf, intercept = TRUE)

  n <- object[["n"]]
  X <- object[["X"]]

  par_names <- names(pars)

  tau <- pars["tau"]
  beta <- pars[grep("beta", par_names)]

  return(rgaussglm(n, X, beta, tau = tau))
}

# """
# Proportional-hazards regression models
# """

#'
#' Predict Method for Bayesian Piecewise Exponential Model Fits
#'
predict_bayesian_piecewise_exponential <- function(pars, formula, data, time_cuts = NULL_VAR) {
  f_call <- match.call(expand.dots = FALSE)
  env <- parent.frame()

  mf <- .extract_model_frame(f_call, env)

  object <- .extract_design_matrix(mf)

  n <- object[["n"]]
  X <- object[["X"]]

  par_names <- names(pars)

  lambda <- pars[grep("lambda", par_names)]
  beta <- pars[grep("beta", par_names)]

  return(rpiecewiseph(n, X, beta, lambda, time_cuts))
}

#'
#' Predict Method for Bayesian Exponential Proportional-Hazards Model Fits
#'
predict_bayesian_exponential_ph <- function(pars, formula, data, time_cuts = NULL_VAR) {
  f_call <- match.call(expand.dots = FALSE)
  env <- parent.frame()

  mf <- .extract_model_frame(f_call, env)

  object <- .extract_design_matrix(mf, intercept = FALSE)

  n <- object[["n"]]
  X <- object[["X"]]

  par_names <- names(pars)

  lambda <- pars["lambda"]
  beta <- pars[grep("beta", par_names)]

  return(rexpph(n, X, beta, lambda = lambda))
}

#'
#' Predict Method for Bayesian Weibull Proportional-Hazards Model Fits
#'
predict_bayesian_weibull_ph <- function(pars, formula, data, time_cuts = NULL_VAR) {
  f_call <- match.call(expand.dots = FALSE)
  env <- parent.frame()

  mf <- .extract_model_frame(f_call, env)

  object <- .extract_design_matrix(mf, intercept = FALSE)

  n <- object[["n"]]
  X <- object[["X"]]

  par_names <- names(pars)

  lambda <- pars["lambda"]
  nu <- pars["nu"]
  beta <- pars[grep("beta", par_names)]

  return(rweibulleconph(n, X, beta, lambda, nu = nu))
}

# """
# Proportional-mean regression models
# """

#'
#' Predict Method for Bayesian Gamma-Gompertz Proportional-Means Model Fits
#'
predict_bayesian_gamma_gompertz_pm <- function(pars, formula, data, time_weights, time_cuts = NULL_VAR) {
  f_call <- match.call(expand.dots = FALSE)
  env <- parent.frame()

  mf <- .extract_model_frame(f_call, env)

  object <- .extract_joint_design_matrix(mf, time_weights, intercept = FALSE)

  n <- object[["n"]]
  t <- object[["t"]]
  X <- object[["X"]]

  par_names <- names(pars)

  alpha <- pars["alpha"]
  eta <- pars["eta"]
  xi <- pars["xi"]
  beta <- pars[grep("beta", par_names)]

  return(rggpmm(n, t, X, beta, eta, xi, alpha = alpha))
}

#'
#' Predict Method for Bayesian Gamma-Weibull Proportional-Means Model Fits
#'
predict_bayesian_gamma_weibull_pm <- function(pars, formula, data, time_weights, time_cuts = NULL_VAR) {
  f_call <- match.call(expand.dots = FALSE)
  env <- parent.frame()

  mf <- .extract_model_frame(f_call, env)

  object <- .extract_joint_design_matrix(mf, time_weights, intercept = FALSE)

  n <- object[["n"]]
  t <- object[["t"]]
  X <- object[["X"]]

  par_names <- names(pars)

  alpha <- pars["alpha"]
  eta <- pars["eta"]
  xi <- pars["xi"]
  beta <- pars[grep("beta", par_names)]

  return(rgwpmm(n, t, X, beta, eta, xi, alpha = alpha))
}

#'
#' Predict Method for Bayesian Gamma-Piecewise Proportional-Means Model Fits
#'
predict_bayesian_gamma_piecewise_pm <- function(pars, formula, data, time_weights, time_cuts = NULL_VAR) {
  f_call <- match.call(expand.dots = FALSE)
  env <- parent.frame()

  mf <- .extract_model_frame(f_call, env)

  object <- .extract_joint_design_matrix(mf, time_weights, intercept = TRUE)

  n <- object[["n"]]
  t <- object[["t"]]
  X <- object[["X"]]

  par_names <- names(pars)

  tau <- pars["tau"]
  mu <- pars[grep("mu", par_names)]
  beta <- pars[grep("beta", par_names)]

  return(rgppmm(n, t, X, beta, mu, tau, time_cuts))
}
