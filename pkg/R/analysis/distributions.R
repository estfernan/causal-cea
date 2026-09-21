#
# Custom probability distributions.
#

# """
# Bernoulli Distributions
# """

#'
#' Bernoulli Distribution: Generalized Linear Model Parameterization
#'
#' These functions provide the density and random number generation for the
#' Bernoulli distribution, given the generalized linear model parameterization
#' with a log link function and support for covariates.
#'
rbernglm <- function(n, X, beta) {
  if (is.vector(X)) {
    X <- matrix(X, ncol = 1)
  }

  p <- LaplacesDemon::invlogit(X %*% beta)
  return(LaplacesDemon::rbern(n, prob = p))
}

# """
# Exponential Distributions
# """

#'
#' Exponential Distribution: Proportional-Hazards Parameterization
#'
#' These functions provide the density and random number generation for the
#' exponential distribution, given the proportional-hazards parameterization with
#' support for covariates.
#'
rexpph <- function(n, X, beta, lambda = 1) {
  if (is.vector(X)) {
    X <- matrix(X, ncol = 1)
  }

  la_ <- lambda * exp(X %*% beta)
  return(rexp(n, rate = la_))
}

# """
# Gamma Distributions
# """

#'
#' Gamma Distribution: Proportional-Means Parameterization
#'
#' These functions provide the density and random number generation for the
#' Gamma distribution, given the proportional-means parameterization with
#' a Gompertz mean function and support for covariates.
#'
rggpmm <- function(n, t, X, beta, eta, xi, alpha = 1) {
  if (is.vector(X)) {
    X <- matrix(X, ncol = 1)
  }

  mu_0 <- function(t) eta * exp(xi * t)
  mu <- mu_0(t) * exp(X %*% beta)
  return(rgamma(n, shape = mu^2 / alpha, rate = mu / alpha))
}

#'
#' Gamma Distribution: Proportional-Means Piecewise Parameterization
#'
#' These functions provide the density and random number generation for the
#' Gamma distribution, given the proportional-means parameterization with
#' a piecewise mean function and support for covariates.
#'
rgppmm <- function(n, t, X, beta, mu, tau, time_cuts) {
  if (is.vector(X)) {
    X <- matrix(X, ncol = 1)
  }

  idx <- findInterval(t, time_cuts, all.inside = TRUE)

  lambda <- mu[idx] * exp(X %*% beta)
  return(rgamma(n, shape = lambda^2 / tau, rate = lambda / tau))
}

#'
#' Gamma Distribution: Proportional-Means Parameterization
#'
#' These functions provide the density and random number generation for the
#' Gamma distribution, given the proportional-means parameterization with
#' a Weibull mean function and support for covariates.
#'
rgwpmm <- function(n, t, X, beta, eta, xi, alpha = 1) {
  if (is.vector(X)) {
    X <- matrix(X, ncol = 1)
  }

  mu_0 <- function(t) eta * xi * t^(xi - 1)
  mu <- mu_0(t) * exp(X %*% beta)
  return(rgamma(n, shape = alpha, rate = alpha / mu))
}

# """
# Gaussian Distributions
# """

#'
#' Normal Distribution: Generalized Linear Model Parameterization
#'
#' These functions provide the density and random number generation for the
#' Normal distribution, given the generalized linear model parameterization
#' and support for covariates.
#'
rgaussglm <- function(n, X, beta, tau) {
  if (is.vector(X)) {
    X <- matrix(X, ncol = 1)
  }

  return(LaplacesDemon::rnormp(n, mean = X %*% beta, prec = tau))
}

# """
# Piecewise Constant Hazard Distributions
# """

#'
#' Piecewise Constant Hazard Distribution: Proportional-Hazards Parameterization
#'
#' These functions provide the density and random number generation for the
#' piecewise constant hazard distribution, given the proportional-hazards
#' parameterization with support for covariates.
#'
rpiecewiseph <- function(n, X, beta, lambda, time_cuts) {
  if (is.vector(X)) {
    X <- matrix(X, ncol = 1)
  }

  la_ <- exp(X %*% beta) %*% t(lambda)
  return(apply(la_, 1, mets::rpch, n = 1, breaks = time_cuts))
}

# """
# Weibull Distributions
# """

#'
#' Weibull Distribution: Alternative Econometrics Parameterization
#'
#' These functions provide the density and random number generation for the
#' Weibull distribution, given a different parameterization, often used in
#' econometrics and medical statistics.
#'
rweibullecon <- function(n, lambda, nu = 1) {
  alpha <- nu
  sigma <- lambda^(-1 / nu)
  return(rweibull(n, shape = alpha, scale = sigma))
}

#'
#' Weibull Distribution: Proportional-Hazards Parameterization
#'
#' These functions provide the density and random number generation for the
#' Weibull distribution, given the proportional-hazards parameterization with
#' support for covariates.
#'
rweibulleconph <- function(n, X, beta, lambda, nu = 1) {
  if (is.vector(X)) {
    X <- matrix(X, ncol = 1)
  }

  la_ <- lambda * exp(X %*% beta)
  return(rweibullecon(n, lambda = la_, nu = nu))
}
