# """
# File name: fit_discrete_models.R
#
# Description:
#   Function to fit discrete-time models to simulated insurance claims data
#   using CmdStanR. Supports Bayesian posterior sampling or frequentist
#   optimization.
# """

fit_discrete_models <- function(data, time_cuts, method = c("sample", "optimize"), ...) {
  # """
  # Fit discrete-time models using CmdStanR.
  # """
  method <- match.arg(method)

  stan_file <- "pkg/stan/discrete_simulation_models.stan"
  exe_file <- tools::file_path_sans_ext(stan_file)

  stan_model <- suppressMessages(
    cmdstanr::cmdstan_model(
      stan_file,
      exe_file = exe_file,
      quiet = TRUE
    )
  )

  init_dat <- data[data$interval == 1, ]
  fu_dat <- data[data$interval > 1, ]

  S0 <- model.matrix(~ L01 + L02 + L03 - 1, data = init_dat)
  X0 <- model.matrix(~ L01 + L02 + L03, data = init_dat)
  H0 <- model.matrix(~ L1 + L2 + L01 + L02 + L03 - 1, data = init_dat)

  S <- model.matrix(~ lag_Tx + lag_L1 + lag_L2 + L01 + L02 + L03 - 1, data = fu_dat)
  X <- model.matrix(~ lag_Tx + lag_L1 + lag_L2 + L01 + L02 + L03, data = fu_dat)
  H <- model.matrix(~ L1 + L2 + lag_Tx + lag_L1 + lag_L2 + L01 + L02 + L03 - 1, data = fu_dat)

  U01 <- init_dat$L01
  U02 <- init_dat$L02
  U03 <- init_dat$L03

  d0 <- init_dat$event_type
  d <- fu_dat$event_type

  L01 <- init_dat$L1
  L02 <- init_dat$L2

  L1 <- fu_dat$L1
  L2 <- fu_dat$L2

  a0 <- init_dat$Tx
  y0 <- init_dat$cost

  a <- fu_dat$Tx
  y <- fu_dat$cost

  is_pos0 <- as.integer(y0 > 0)
  is_pos <- as.integer(y > 0)

  idxs <- fu_dat$interval

  # sample sizes and number of predictors
  n <- length(U01)
  K0 <- nrow(init_dat)
  K <- nrow(fu_dat)
  J <- length(time_cuts) - 1
  P0 <- ncol(S0)
  P <- ncol(S)
  R0 <- ncol(X0)
  R <- ncol(X)
  C0 <- ncol(H0)
  C <- ncol(H)

  # prepare data list for Stan
  stan_data <- list(
    n = n, K0 = K0, K = K,
    J = J,
    idxs = idxs,
    U01 = U01, U02 = U02, U03 = U03,
    P0 = P0, S0 = S0, d0 = d0,
    P = P, S = S, d = d,
    R0 = R0, X0 = X0, L01 = L01, L02 = L02,
    R = R, X = X, L1 = L1, L2 = L2,
    C0 = C0, H0 = H0, a0 = a0, y0 = y0, is_pos0 = is_pos0,
    C = C, H = H, a = a, y = y, is_pos = is_pos
  )

  stan_out <- stan_model[[method]](data = stan_data, ...)

  return(stan_out)
}
