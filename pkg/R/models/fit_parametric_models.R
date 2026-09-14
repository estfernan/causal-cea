# """
# File name: fit_parametric_models.R
#
# Description:
#   Function to fit correctly-specified parametric models to simulated insurance
#   claims data using CmdStanR. Supports Bayesian posterior sampling or
#   frequentist optimization.
# """

fit_parametric_models <- function(data, method = "sample", ...) {
  # """
  # Fit parametric models using CmdStanR.
  # """
  stan_file <- "pkg/stan/parametric_simulation_models.stan"
  exe_file <- tools::file_path_sans_ext(stan_file)

  stan_model <- suppressMessages(cmdstanr::cmdstan_model(stan_file, exe_file = exe_file, quiet = TRUE))

  init_dat <- data[data$claim_idx == 1, ]
  fu_dat   <- data[data$claim_idx > 1, ]

  # core settings
  n <- nrow(init_dat)
  K <- nrow(fu_dat)

  # baseline observations
  U01 <- init_dat$L01
  U02 <- init_dat$L02
  U03 <- init_dat$L03

  # baseline time-to-event observations
  S0  <- model.matrix(~ L01 + L02 + L03 - 1, data = init_dat)
  P0  <- ncol(S0)
  d00 <- as.integer(init_dat$event_type == 0)
  d01 <- as.integer(init_dat$event_type == 1)
  t0  <- init_dat$gaptime

  # follow-up time-to-event observations
  S  <- model.matrix(~ lag_Tx + lag_L1 + lag_L2 + L01 + L02 + L03 - 1, data = fu_dat)
  P  <- ncol(S)
  d0 <- as.integer(fu_dat$event_type == 0)
  d1 <- as.integer(fu_dat$event_type == 1)
  t  <- fu_dat$gaptime

  # baseline covariate observations
  X0  <- model.matrix(~ L01 + L02 + L03, data = init_dat)
  R0  <- ncol(X0)
  L01 <- init_dat$L1
  L02 <- init_dat$L2

  # follow-up covariate observations
  X  <- model.matrix(~ lag_Tx + lag_L1 + lag_L2 + L01 + L02 + L03, data = fu_dat)
  R  <- ncol(X)
  L1 <- fu_dat$L1
  L2 <- fu_dat$L2

  # baseline cost observations
  H0 <- model.matrix(~ L1 + L2 + L01 + L02 + L03 - 1, data = init_dat)
  C0 <- ncol(H0)
  a0 <- init_dat$Tx
  y0 <- init_dat$cost

  # follow-up cost observations
  H  <- model.matrix(~ L1 + L2 + lag_Tx + lag_L1 + lag_L2 + L01 + L02 + L03 - 1, data = fu_dat)
  C  <- ncol(H)
  a <- fu_dat$Tx
  y <- fu_dat$cost

  # prepare data list for Stan
  stan_data <- list(
    n = n, K = K,
    U01 = U01, U02 = U02, U03 = U03,
    P0 = P0, S0 = S0, d00 = d00, d01 = d01, t0 = t0,
    P = P, S = S, d0 = d0, d1 = d1, t = t,
    R0 = R0, X0 = X0, L01 = L01, L02 = L02,
    R = R, X = X, L1 = L1, L2 = L2,
    C0 = C0, H0 = H0, a0 = a0, y0 = y0,
    C = C, H = H, a = a, y = y
  )

  stan_out <- stan_model[[method]](data = stan_data, ...)

  return(stan_out)
}
