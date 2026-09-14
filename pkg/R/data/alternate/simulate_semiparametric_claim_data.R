# """
# File name: simulate_semiparametric_claim_data.R
#
# Description:
#   Function for simulating insurance claims data under the alternative
#   data-generating process with semiparametric models.
# """

simulate_semiparametric_claim_data <- function(
  n,
  n_claims,
  time_cuts,
  dgp_list,
  target_estimand = NULL,
  ...,
  treatment_decision = NULL,
  time_horizon = Inf,
  censoring_rate = 1e-10
) {
  # """
  # Simulate insurance claims data under a semiparametric data-generating process.
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

  P0 <- 3
  P <- 5

  L0 <- matrix(NA_real_, nrow = n, ncol = P0)
  V  <- matrix(NA_real_, nrow = n, ncol = n_claims)
  W  <- matrix(NA_real_, nrow = n, ncol = n_claims)
  D  <- matrix(NA_integer_, nrow = n, ncol = n_claims)
  A  <- matrix(NA_integer_, nrow = n, ncol = n_claims)
  Y  <- matrix(NA_real_, nrow = n, ncol = n_claims)
  L  <- array(NA_real_, dim = c(n, n_claims, P - P0), dimnames = NULL)

  if (!simulate_treatment) {
    A <- matrix(treatment_decision, nrow = n, ncol = n_claims)
  }

  L0[, 1] <- rnorm(n, mean = dgp_list$mu01, sd = dgp_list$sigma01)
  L0[, 2] <- rnorm(n, mean = dgp_list$mu02, sd = dgp_list$sigma02)
  L0[, 3] <- LaplacesDemon::rbern(n, prob = dgp_list$p03)

  w0 <- rpiecewiseph(L0, time_cuts$init, dgp_list$psi00, dgp_list$h00)
  w1 <- rpiecewiseph(L0, time_cuts$init, dgp_list$psi01, dgp_list$h01)

  C <- rexp(n, rate = censoring_rate)

  if (apply_time_cutoff) {
    C <- pmin(C, time_horizon)
  }

  W[, 1] <- pmin(w0, w1, C)
  V[, 1] <- W[, 1]
  D[, 1] <- max.col(-cbind(C, w0, w1), ties.method = "first") - 2

  X <- concat(n, 1, L0)

  L[, 1, 1] <- rnorm(n, mean = X %*% dgp_list$omegaL01, sd = dgp_list$sigmaL01)
  L[, 1, 2] <- LaplacesDemon::rbern(n, prob = plogis(X %*% dgp_list$omegaL02))

  if (simulate_treatment) {
    H <- concat(n, 1, L[, 1, ], L0)
    A[, 1] <- LaplacesDemon::rbern(n, prob = plogis(H %*% dgp_list$theta0))
  }

  H <- concat(n, L[, 1, ], L0)

  mu0 <- intensity(
    W[, 1], A[, 1],
    H,
    dgp_list$betaA0, dgp_list$betaH0,
    hazard_piecewise,
    time_cuts$init, dgp_list$m0
  )

  Y[, 1] <- rgamma(n, shape = mu0^2 / dgp_list$gamma0, rate = mu0 / dgp_list$gamma0)

  j   <- 2L
  idx <- D[, 1] > 0
  B   <- sum(idx)

  # simulations for subsequent claims
  while (B > 0) {
    if (j > n_claims) {
      stop("Insufficient storage: increase 'n_claims' to allow more events.")
    }

    S <- concat(B, A[idx, j - 1], L[idx, j - 1, ], L0[idx, ])

    w0 <- rpiecewiseph(S, time_cuts$post, dgp_list$psi0, dgp_list$h0)
    w1 <- rpiecewiseph(S, time_cuts$post, dgp_list$psi1, dgp_list$h1)

    C <- rexp(B, rate = censoring_rate)

    if (apply_time_cutoff) {
      R <- time_horizon - V[idx, j - 1]
      C <- pmin(C, R)
    }

    W[idx, j] <- pmin(w0, w1, C)
    V[idx, j] <- V[idx, j - 1] + W[idx, j]
    D[idx, j] <- max.col(-cbind(C, w0, w1), ties.method = "first") - 2

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

    mu <- intensity(
      W[idx, j], A[idx, j],
      H,
      dgp_list$betaA, dgp_list$betaH,
      hazard_piecewise,
      time_cuts$post, dgp_list$m
    )

    Y[idx, j] <- rgamma(B, shape = mu^2 / dgp_list$gamma, rate = mu / dgp_list$gamma)

    idx <- ifelse(idx, D[, j] > 0, idx)
    B   <- sum(idx)
    j   <- j + 1L
  }

  if (compute_target_estimand) {
    U   <- rowSums(Y, na.rm = TRUE)
    T_  <- rowSums(W, na.rm = TRUE)
    EU  <- mean(U, na.rm = TRUE)
    ET_ <- mean(T_, na.rm = TRUE)
    return(target_estimand(EU, ET_, ...))
  }

  return(build_simulated_dataset(n, P0, L0, V, W, D, L, A, Y))
}
