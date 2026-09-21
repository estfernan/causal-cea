#
# Data generating functions to simulate observations for the numerical experiments.
#

#'
#' Simulate Longitudinal and Time-to-Event Observations
#'
simulate_observations <- function(n, theta, lambda_C = 1e-10, tau = Inf, a = NULL, target_estimand = NULL, max_iter = 250, print_diagnostics = FALSE, ...) {
  gen_treat <- is.null(a)
  calc_target <- !is.null(target_estimand)

  baseline_dinnames <- list(row = seq_len(n), col = BASELINE_COVARNAMES)
  timevarying_dinnames <- list(row = seq_len(n), col = seq_len(max_iter), slice = TIMEVARYING_COVARNAMES)

  # set data generating parameters from list
  mu_X1 <- theta[["mu_X1"]]
  tau_X1 <- theta[["tau_X1"]]
  mu_X2 <- theta[["mu_X2"]]
  tau_X2 <- theta[["tau_X2"]]
  p_X3 <- theta[["p_X3"]]
  mu_L10 <- theta[["mu_L10"]]
  tau_L10 <- theta[["tau_L10"]]
  p_L20 <- theta[["p_L20"]]
  phi_0 <- theta[["phi_0"]]
  nu_10 <- theta[["nu_10"]]
  lambda_10 <- theta[["lambda_10"]]
  psi_10 <- theta[["psi_10"]]
  nu_00 <- theta[["nu_00"]]
  lambda_00 <- theta[["lambda_00"]]
  psi_00 <- theta[["psi_00"]]
  alpha_0 <- theta[["alpha_0"]]
  eta_0 <- theta[["eta_0"]]
  xi_0 <- theta[["xi_0"]]
  beta_0 <- theta[["beta_0"]]
  beta_L1 <- theta[["beta_L1"]]
  tau_L1 <- theta[["tau_L1"]]
  beta_L2 <- theta[["beta_L2"]]
  phi <- theta[["phi"]]
  nu_1 <- theta[["nu_1"]]
  lambda_1 <- theta[["lambda_1"]]
  psi_1 <- theta[["psi_1"]]
  nu_0 <- theta[["nu_0"]]
  lambda_0 <- theta[["lambda_0"]]
  psi_0 <- theta[["psi_0"]]
  alpha <- theta[["alpha"]]
  eta <- theta[["eta"]]
  xi <- theta[["xi"]]
  beta <- theta[["beta"]]

  V <- matrix(nrow = n, ncol = max_iter)
  X <- matrix(nrow = n, ncol = P_BASELINE, dimnames = baseline_dinnames)
  L <- array(dim = c(n, max_iter, P_TIMEVARYING), dimnames = timevarying_dinnames)
  A <- matrix(nrow = n, ncol = max_iter)
  W <- matrix(nrow = n, ncol = max_iter)
  Z <- matrix(nrow = n, ncol = max_iter)
  Delta <- matrix(nrow = n, ncol = max_iter)

  if (!gen_treat) {
    A <- matrix(a, nrow = n, ncol = max_iter)
  }

  # simulate baseline confounders
  X[, 1] <- LaplacesDemon::rnormp(n, mean = mu_X1, prec = tau_X1)
  X[, 2] <- LaplacesDemon::rnormp(n, mean = mu_X2, prec = tau_X2)
  X[, 3] <- LaplacesDemon::rbern(n, prob = p_X3)

  # simulate initial time-varying confounders
  L[, 1, 1] <- LaplacesDemon::rnormp(n, mean = mu_L10, prec = tau_L10)
  L[, 1, 2] <- LaplacesDemon::rbern(n, prob = p_L20)

  if (gen_treat) {
    A[, 1] <- .simulate_binary_treatment(n, phi_0, X, L[, 1, ])
  }

  W_V <- .simulate_event_gaptime(n, psi_10, lambda_10, nu_10, X, L[, 1, ], A[, 1], int = NULL)
  W_T <- .simulate_event_gaptime(n, psi_00, lambda_00, nu_00, X, L[, 1, ], A[, 1], int = NULL)

  W_tau <- tau
  W_tilde_C <- rexp(n, rate = lambda_C)

  W_C <- pmin(W_tau, W_tilde_C)
  to_cens <- W_C < W_T

  W_obs <- pmin(W_V, W_T, W_C)
  to_claim <- W_obs == W_V

  delta <- ifelse(to_claim, REGULAR_EVENT_TYPE, to_cens * CENSOR_EVENT_TYPE)

  V[, 1] <- DEFAULT_START_TIME
  W[, 1] <- W_obs
  V[, 2] <- W[, 1]
  Delta[, 1] <- delta

  Z[, 1] <- .simulate_cost(n, W[, 1], beta_0, eta_0, xi_0, alpha_0, X, L[, 1, ], A[, 1], int = NULL)

  k <- 2L

  idx <- which(to_claim)
  n_left <- length(idx)

  while (n_left > 0) {
    L[idx, k, 1] <- .simulate_gaussian_confounder(
      n_left,
      beta_L1, tau_L1,
      X[idx, ], L[idx, k - 1, ], A[idx, k - 1]
    )

    L[idx, k, 2] <- .simulate_bernoulli_confounder(
      n_left,
      beta_L2,
      X[idx, ], L[idx, k - 1, ], A[idx, k - 1]
    )

    if (gen_treat) {
      A[idx, k] <- .simulate_binary_treatment(
        n_left, phi,
        X[idx, ], L[idx, k, ], L[idx, k - 1, ], A[idx, k - 1]
      )
    }

    W_V <- .simulate_event_gaptime(
      n_left, psi_1, lambda_1, nu_1,
      X[idx, ], L[idx, k, ], L[idx, k - 1, ], A[idx, k], A[idx, k - 1],
      int = NULL
    )

    W_T <- .simulate_event_gaptime(
      n_left, psi_0, lambda_0, nu_0,
      X[idx, ], L[idx, k, ], L[idx, k - 1, ], A[idx, k], A[idx, k - 1],
      int = NULL
    )

    W_tau <- tau - V[idx, k]
    W_tilde_C <- rexp(n_left, rate = lambda_C)

    W_C <- pmin(W_tau, W_tilde_C)
    to_cens <- W_C < W_T

    W_obs <- pmin(W_V, W_T, W_C)
    to_claim <- W_obs == W_V

    delta <- ifelse(to_claim, REGULAR_EVENT_TYPE, to_cens * CENSOR_EVENT_TYPE)

    W[idx, k] <- W_obs
    V[idx, k + 1] <- V[idx, k] + W[idx, k]
    Delta[idx, k] <- delta

    Z[idx, k] <- .simulate_cost(
      n_left, W[idx, k],
      beta, eta, xi, alpha,
      X[idx, ], L[idx, k, ], L[idx, k - 1, ], A[idx, k], A[idx, k - 1],
      int = NULL
    )

    k <- k + 1

    idx <- idx[which(to_claim)]
    n_left <- length(idx)
  }

  K_list <- rowSums(!is.na(W))

  # display the number of censored individuals at each event
  if (print_diagnostics) {
    last_idx <- max(K_list)
    full_idx <- seq_len(last_idx)

    n_dropout <- colSums(Delta == CENSOR_EVENT_TYPE, na.rm = TRUE)
    n_risk <- colSums(!is.na(Delta), na.rm = TRUE)

    n_dropout <- n_dropout[full_idx]
    n_risk <- n_risk[full_idx]

    dropout_rates <- n_dropout / n_risk

    first_dropout <- n_dropout[1] / n_risk[1]
    string <- scales::percent(first_dropout, accuracy = 0.01)
    cat("Censoring Rate at the First Event:\t", string, "\n")

    avg_dropout <- median(n_dropout[-1] / n_risk[-1], na.rm = TRUE)
    string <- scales::percent(avg_dropout, accuracy = 0.01)
    cat("Censoring Rate at the First Event:\t", string, "\n\n")
  }

  # drop entry for time of death -- not necessary for modeling
  for (i in 1:n)
  {
    V[i, K_list[i] + 1] <- NA
  }

  if (calc_target) {
    T_ <- rowSums(W, na.rm = TRUE)
    Y <- rowSums(Z, na.rm = TRUE)
    return(target_estimand(T_, Y, ...))
  }

  dtab <- .build_simulation_dataset(n, V, X, L, A, W, Z, Delta)

  return(dtab)
}

# """
# Internal Functions
# """

#'
#' Convert the Simulated Matrices to a Final Data Frame Object
#'
.build_simulation_dataset <- function(n, V, X, L, A, W, Z, Delta, lag_vars = TRUE) {
  simu_obj <- list(
    V = reshape2::melt(V, varnames = KEYS, na.rm = TRUE, value.name = "TIME"),
    W = reshape2::melt(W, varnames = KEYS, na.rm = TRUE, value.name = "GAPTIME"),
    Delta = reshape2::melt(Delta, varnames = KEYS, na.rm = TRUE, value.name = "STATUS"),
    Z = reshape2::melt(Z, varnames = KEYS, na.rm = TRUE, value.name = "COST"),
    A = reshape2::melt(A, varnames = KEYS, na.rm = TRUE, value.name = "TREATMENT"),
    L1 = reshape2::melt(L[, , 1], varnames = KEYS, na.rm = TRUE, value.name = "L1"),
    L2 = reshape2::melt(L[, , 2], varnames = KEYS, na.rm = TRUE, value.name = "L2")
  )

  merged <- Reduce(
    f = function(x, y) merge(x, y, by = KEYS, all = TRUE),
    simu_obj
  )

  dtab <- data.table::merge.data.table(
    x = data.table::as.data.table(merged),
    y = data.table::data.table(ID = 1:n, X),
    by = "ID", all = TRUE
  )

  data.table::setorder(dtab, ID, INDEX)
  data.table::setcolorder(dtab, neworder = BASELINE_COVARNAMES, after = "TREATMENT")

  if (lag_vars) {
    nms <- paste("LAG", LAG_VARS, sep = "_")
    dtab <- dtab[, (nms) := data.table::shift(.SD, n = 1), by = ID, .SDcols = LAG_VARS]
  }

  return(dtab)
}

#'
#' Generate Time-Varying Binary Confounder from a Bernoulli-GLM
#'
.simulate_bernoulli_confounder <- function(n, beta, ...) {
  H <- .design(n, ...)
  return(rbernglm(n, H, beta))
}

#'
#' Generate a Binary Treatment Assignment from a Bernoulli-GLM
#'
.simulate_binary_treatment <- function(n, beta, ...) {
  H <- .design(n, ...)
  return(rbernglm(n, H, beta))
}

#'
#' Generate a Cost from a Gamma-Gompertz Proportional-Means Model
#'
.simulate_cost <- function(n, t, beta, eta, xi, alpha, ...) {
  H <- .design(n, ...)
  return(rggpmm(n, t, H, beta, eta, xi, alpha))
}

#'
#' Generate a Waiting-Time-to-Event from a Weibull Proportional-Hazards Model
#'
.simulate_event_gaptime <- function(n, beta, lambda, nu, ...) {
  H <- .design(n, ...)
  return(rweibulleconph(n, H, beta, lambda, nu))
}

#'
#' Generate a Time-Varying Continuous Confounder from a Gaussian-GLM
#'
.simulate_gaussian_confounder <- function(n, beta, tau, ...) {
  H <- .design(n, ...)
  return(rgaussglm(n, H, beta, tau))
}
