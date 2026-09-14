# """
# File name: simulate_discrete_claim_data.R
#
# Description:
#   Function for simulating insurance claims data under the alternative
#   data-generating process with discrete-models models.
# """

simulate_discrete_claim_data <- function(
  n,
  n_intervals,
  time_cuts,
  dgp_list,
  target_estimand = NULL,
  ...,
  treatment_decision = NULL,
  time_horizon = Inf,
  censoring_rate = 1e-10
) {
  # """
  # Simulate insurance claims data under a discrete-time data-generating process.
  # """

  # setup flags
  apply_time_cutoff <- is.finite(time_horizon)
  compute_target_estimand <- !is.null(target_estimand)
  read_JSON <- is.character(dgp_list)
  simulate_treatment <- is.null(treatment_decision)

  # load parameter list from JSON file
  if (read_JSON) {
    dgp_list <- jsonlite::read_json(dgp_list, simplifyVector = TRUE)
  }

  time_mappings <- time_cuts[-(n_intervals + 1)] + diff(time_cuts) / 2

  P0 <- 3
  P <- 5

  L0 <- matrix(NA_real_, nrow = n, ncol = P0)
  V <- matrix(NA_real_, nrow = n, ncol = n_intervals)
  W <- matrix(NA_real_, nrow = n, ncol = n_intervals)
  D <- matrix(NA_integer_, nrow = n, ncol = n_intervals)
  A <- matrix(NA_integer_, nrow = n, ncol = n_intervals)
  Y <- matrix(NA_real_, nrow = n, ncol = n_intervals)
  L <- array(NA_real_, dim = c(n, n_intervals, P - P0), dimnames = NULL)

  if (!simulate_treatment) {
    A <- matrix(treatment_decision, nrow = n, ncol = n_intervals)
  }

  L0[, 1] <- rnorm(n, mean = dgp_list$mu01, sd = dgp_list$sigma01)
  L0[, 2] <- rnorm(n, mean = dgp_list$mu02, sd = dgp_list$sigma02)
  L0[, 3] <- LaplacesDemon::rbern(n, prob = dgp_list$p03)

  D[, 1] <- LaplacesDemon::rbern(n, prob = plogis(dgp_list$alpha_h[1] + L0 %*% dgp_list$psi0))
  V[, 1] <- time_mappings[1]

  X <- concat(n, 1, L0)

  L[, 1, 1] <- rnorm(n, mean = X %*% dgp_list$omegaL01, sd = dgp_list$sigmaL01)
  L[, 1, 2] <- LaplacesDemon::rbern(n, prob = plogis(X %*% dgp_list$omegaL02))

  if (simulate_treatment) {
    H <- concat(n, 1, L[, 1, ], L0)
    A[, 1] <- LaplacesDemon::rbern(n, prob = plogis(H %*% dgp_list$theta0))
  }

  H <- concat(n, L[, 1, ], L0)

  p0 <- plogis(dgp_list$theta[1] + A[, 1] * dgp_list$phiA0 + H %*% dgp_list$phiH0)
  mu0 <- exp(dgp_list$alpha_m[1] + A[, 1] * dgp_list$betaA0 + H %*% dgp_list$betaH0)

  is_pos_0 <- LaplacesDemon::rbern(n, prob = p0)

  Y[, 1] <- is_pos_0 * rgamma(n, shape = mu0^2 / dgp_list$gamma0, rate = mu0 / dgp_list$gamma0)

  idx <- D[, 1] < 1
  B <- sum(idx)

  # simulations for subsequent intervals
  for (j in 2:n_intervals) {
    if (B == 0) {
      break()
    }

    S <- concat(B, A[idx, j - 1], L[idx, j - 1, ], L0[idx, ])

    D[idx, j] <- LaplacesDemon::rbern(B, prob = plogis(dgp_list$alpha_h[j] + S %*% dgp_list$psi))
    V[idx, j] <- time_mappings[j]

    X <- concat(B, 1, A[idx, j - 1], L[idx, j - 1, ], L0[idx, ])

    L[idx, j, 1] <- rnorm(B, mean = X %*% dgp_list$omegaL1, sd = dgp_list$sigmaL1)
    L[idx, j, 2] <- LaplacesDemon::rbern(B, prob = plogis(X %*% dgp_list$omegaL2))

    # simulate or maintain the adjuvant treatment
    if (simulate_treatment) {
      H <- concat(B, 1, L[idx, j, ], A[idx, j - 1], L[idx, j - 1, ], L0[idx, ])
      a_opt <- LaplacesDemon::rbern(B, prob = plogis(H %*% dgp_list$theta))
      A[idx, j] <- ifelse(A[idx, j - 1] > 0, A[idx, j - 1], a_opt)
    }

    H <- concat(B, L[idx, j, ], A[idx, j - 1], L[idx, j - 1, ], L0[idx, ])

    p <- plogis(dgp_list$theta[j] + A[idx, j] * dgp_list$phiA + H %*% dgp_list$phiH)
    mu <- exp(dgp_list$alpha_m[j] + A[idx, j] * dgp_list$betaA + H %*% dgp_list$betaH)

    is_pos <- LaplacesDemon::rbern(B, prob = p)

    Y[idx, j] <- is_pos * rgamma(B, shape = mu^2 / dgp_list$gamma, rate = mu / dgp_list$gamma)

    idx <- ifelse(idx, D[, j] < 1, idx)
    B <- sum(idx)
  }

  if (compute_target_estimand) {
    U <- rowSums(Y, na.rm = TRUE)
    T_ <- apply(V, 1, max, na.rm = TRUE)

    if (apply_time_cutoff) {
      T_ <- pmin(T_, time_horizon)
    }

    EU <- mean(U, na.rm = TRUE)
    ET_ <- mean(T_, na.rm = TRUE)

    return(target_estimand(EU, ET_, ...))
  }

  return(build_simulated_dataset(n, P0, L0, V, W, D, L, A, Y))
}
