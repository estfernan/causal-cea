# """
# File name: fit_semiparametric_models.R
#
# Description:
#   Function to fit semiparametric models to simulated insurance claims data
#   using CmdStanR. Supports Bayesian posterior sampling or frequentist
#   optimization.
# """

fit_semiparametric_models <- function(data, time_cuts, method = c("sample", "optimize"), ...) {
  # """
  # Fit semiparametric models using CmdStanR.
  # """
  method <- match.arg(method)

  stan_file <- "pkg/stan/semiparametric_simulation_models.stan"
  exe_file <- tools::file_path_sans_ext(stan_file)

  stan_model <- cmdstanr::cmdstan_model(
    stan_file,
    exe_file = exe_file,
    quiet = TRUE
  )

  init_dat <- data[data$claim_idx == 1, ]
  fu_dat   <- data[data$claim_idx > 1, ]

  init_visit <- split_surv(init_dat, time_cuts$init, 1L, "L01 + L02 + L03 - 1")
  init_death <- split_surv(init_dat, time_cuts$init, 0L, "L01 + L02 + L03 - 1")

  fu_visit <- split_surv(fu_dat, time_cuts$post, 1L, "lag_Tx + lag_L1 + lag_L2 + L01 + L02 + L03 - 1")
  fu_death <- split_surv(fu_dat, time_cuts$post, 0L, "lag_Tx + lag_L1 + lag_L2 + L01 + L02 + L03 - 1")

  S0 <- model.matrix(~ L01 + L02 + L03 - 1, data = init_visit)
  X0 <- model.matrix(~ L01 + L02 + L03, data = init_dat)
  H0 <- model.matrix(~ L1 + L2 + L01 + L02 + L03 - 1, data = init_dat)

  S <- model.matrix(~ lag_Tx + lag_L1 + lag_L2 + L01 + L02 + L03 - 1, data = fu_visit)
  X <- model.matrix(~ lag_Tx + lag_L1 + lag_L2 + L01 + L02 + L03, data = fu_dat)
  H <- model.matrix(~ L1 + L2 + lag_Tx + lag_L1 + lag_L2 + L01 + L02 + L03 - 1, data = fu_dat)

  d00 <- init_death$event_type
  d01 <- init_visit$event_type
  D0  <- init_visit$tstop - init_visit$tstart
  I0  <- init_visit$interval - 1

  d0 <- fu_death$event_type
  d1 <- fu_visit$event_type
  D  <- fu_visit$tstop - fu_visit$tstart
  I  <- fu_visit$interval - 1

  t0 <- findInterval(init_dat$gaptime, time_cuts$init, all.inside = TRUE)
  t  <- findInterval(fu_dat$gaptime, time_cuts$post, all.inside = TRUE)

  y0  <- init_dat$cost
  a0  <- init_dat$Tx
  U01 <- init_dat$L01
  U02 <- init_dat$L02
  U03 <- init_dat$L03
  L01 <- init_dat$L1
  L02 <- init_dat$L2

  y  <- fu_dat$cost
  a  <- fu_dat$Tx
  L1 <- fu_dat$L1
  L2 <- fu_dat$L2

  # sample sizes and number of predictors
  n  <- nrow(init_dat)
  K  <- nrow(fu_dat)
  J0 <- length(time_cuts$init) - 1
  J  <- length(time_cuts$post) - 1
  P0 <- ncol(S0)
  P  <- ncol(S)
  N0 <- nrow(S0)
  N  <- nrow(S)
  C0 <- ncol(H0)
  C  <- ncol(H)
  R0 <- ncol(X0)
  R  <- ncol(X)

  # prepare data list for Stan
  stan_data <- list(
    n = n, K = K,
    J0 = J0, J = J,
    P0 = P0, P = P,
    U01 = U01, U02 = U02, U03 = U03,
    N0 = N0, N = N,
    S0 = S0, d00 = d00, d01 = d01, D0 = D0, I0 = I0,
    S = S, d0 = d0, d1 = d1, D = D, I = I,
    C0 = C0, C = C,
    H0 = H0, t0 = t0, a0 = a0, y0 = y0,
    H = H, t = t, a = a, y = y,
    R0 = R0, R = R,
    X0 = X0, L01 = L01, L02 = L02,
    X = X, L1 = L1, L2 = L2
  )

  stan_out <- stan_model[[method]](data = stan_data, ...)

  return(stan_out)
}
