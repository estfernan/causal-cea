# """
# File name: misc.R
#
# Description:
#   Functions used to facilitate operations within scripts.
# """

bayes_finite_sample_properties <- function(draws, truth, ..., alpha = 0.05, start_time = NULL) {
  # """
  # Compute the finite-sample properties of an estimator.
  # """
  lwr_prob <- alpha / 2
  upp_prob <- 1 - lwr_prob

  post_mean <- mean(draws, na.rm = TRUE)
  SE <- sd(draws, na.rm = TRUE)

  lwr_CI <- quantile(draws, p = lwr_prob, na.rm = TRUE, names = FALSE)
  upp_CI <- quantile(draws, p = upp_prob, na.rm = TRUE, names = FALSE)

  bias <- post_mean - truth
  pct_bias <- bias / abs(truth)
  MAPE <- abs(bias / truth)
  RMSE <- sqrt(bias^2 + SE^2)
  coverage <- lwr_CI < truth & truth < upp_CI
  width <- upp_CI - lwr_CI

  evalout <- data.frame(
    ...,
    truth = truth,
    mean = post_mean,
    SE = SE,
    lwr_CI = lwr_CI,
    upp_CI = upp_CI,
    bias = bias,
    pct_bias = pct_bias,
    MAPE = MAPE,
    RMSE = RMSE,
    coverage = coverage,
    CI_width = width
  )

  if (!is.null(start_time)) {
    elapsed <- difftime(time1 = Sys.time(), time2 = start_time, units = "secs")
    evalout <- cbind(evalout, runtime = elapsed)
  }

  return(evalout)
}

bootstrap_finite_sample_properties <- function(draws, truth, ..., alpha = 0.05, start_time = NULL) {
  # """
  # Compute the finite-sample properties of an estimator.
  # """
  lwr_prob <- alpha / 2
  upp_prob <- 1 - lwr_prob

  lwr_CI <- quantile(draws, p = lwr_prob, na.rm = TRUE, names = FALSE)
  upp_CI <- quantile(draws, p = upp_prob, na.rm = TRUE, names = FALSE)

  estimate <- mean(draws, na.rm = TRUE)
  SE <- sd(draws, na.rm = TRUE)

  bias <- estimate - truth
  pct_bias <- bias / abs(truth)
  MAPE <- abs(bias / truth)
  RMSE <- sqrt(bias^2 + SE^2)
  coverage <- lwr_CI < truth & truth < upp_CI
  width <- upp_CI - lwr_CI

  evalout <- data.frame(
    ...,
    truth = truth,
    mean = estimate,
    SE = SE,
    lwr_CI = lwr_CI,
    upp_CI = upp_CI,
    bias = bias,
    pct_bias = pct_bias,
    MAPE = MAPE,
    RMSE = RMSE,
    coverage = coverage,
    CI_width = width
  )

  if (!is.null(start_time)) {
    elapsed <- difftime(time1 = Sys.time(), time2 = start_time, units = "secs")
    evalout <- cbind(evalout, runtime = elapsed)
  }

  return(evalout)
}

create_discrete_time_intervals <- function(data, config) {
  # """
  # Set the time intervals for the discrete time approach.
  # """
  time_cuts <- seq(
    from = 0,
    to = max(data$time) + 0.01,
    length.out = config$n_intervals + 1
  )

  return(time_cuts)
}

create_time_cuts <- function(data, config) {
  # """
  # Set the time intervals for baseline and follow-up events.
  # """
  init_dat <- subset(data, claim_idx == 1)
  fu_dat <- subset(data, claim_idx > 1)

  init_time_cuts <- seq(
    from = 0,
    to = max(init_dat$gaptime) + 0.01,
    length.out = config$n_init_cuts + 1
  )

  time_cuts <- seq(
    from = 0,
    to = max(fu_dat$gaptime) + 0.01,
    length.out = config$n_cuts + 1
  )

  return(list(init = init_time_cuts, post = time_cuts))
}

format_simulation_output_path <- function(config, opt) {
  # """
  # Construct the file path for a single simulation replicate output.
  # """
  rate <- 100 * round(opt$rate, digits = 3)

  logdir <- file.path(
    "logs",
    opt$framework,
    opt$models,
    sprintf("n%04d", opt$n),
    sprintf("r%03d", rate)
  )

  if (!dir.exists(logdir)) {
    dir.create(logdir, recursive = TRUE)
  }

  outdir <- file.path(
    config$outdir,
    opt$framework,
    opt$models,
    sprintf("n%04d", opt$n),
    sprintf("r%03d", rate)
  )

  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }

  filename <- sprintf("rep%04d.RData", opt$iter)

  return(file.path(outdir, filename))
}

generate_person_period_dataset <- function(data, time_cuts, lag_confounders = TRUE) {
  # """
  # Create a person-period data set from the continuous-time claims data
  # """
  data$interval <- findInterval(data$time, time_cuts, all.inside = TRUE)

  collapsed <- data |>
    dplyr::summarise(
      event_type = min(event_type, na.rm = TRUE),
      Tx = max(Tx, na.rm = TRUE),
      cost = sum(cost, na.rm = TRUE),
      L1 = mean(L1, na.rm = TRUE),
      L2 = max(L2, na.rm = TRUE),
      L01 = mean(L01, na.rm = TRUE),
      L02 = mean(L02, na.rm = TRUE),
      L03 = mean(L03, na.rm = TRUE),
      .by = c(patient_id, interval)
    )

  full <- collapsed |>
    dplyr::group_by(patient_id) |>
    tidyr::complete(interval = tidyr::full_seq(1:max(interval, na.rm = TRUE), period = 1L)) |>
    tidyr::fill(event_type, Tx, L1, L2, L01, L02, L03, .direction = "down") |>
    tidyr::fill(Tx, L1, L2, L01, L02, L03, .direction = "up") |>
    dplyr::mutate(event_type = ifelse(is.na(event_type), 1L, event_type)) |>
    dplyr::ungroup(patient_id)

  full <- full |>
    dplyr::mutate(
      event_type = as.integer(event_type == 0),
      cost = ifelse(is.na(cost), 0.0, cost)
    )

  if (lag_confounders) {
    full <- full |>
      dplyr::group_by(patient_id) |>
      dplyr::mutate(dplyr::across(c(L1, L2, Tx), ~ dplyr::lag(.x, n = 1L), .names = "lag_{col}")) |>
      dplyr::ungroup()
  }

  return(full)
}

label_censoring <- function(censoring_rate) {
  # """
  # Label censoring level.
  # """
  labs <- case_when(
    censoring_rate == 1e-10 ~ "none",
    censoring_rate == 0.10 ~ "low",
    censoring_rate == 0.50 ~ "high",
  )

  return(labs)
}

label_intervals <- function(intervals) {
  # """
  # Label discrete intervals.
  # """
  labs <- case_when(
    intervals == "c010" ~ "10 intervals",
    intervals == "c025" ~ "25 intervals",
    intervals == "c050" ~ "50 intervals",
    intervals == "standard" ~ "--"
  )

  return(labs)
}

label_methods <- function(framework, models) {
  # """
  # Label simulation methods.
  # """
  labs <- case_when(
    framework == "bayes" & models == "parametric" ~ "correctly-specified Bayesian parametric models",
    framework == "bayes" & models == "semiparametric" ~ "Bayesian semiparametric models",
    framework == "bayes" & models == "misspecified" ~ "misspecified Bayesian parametric models",
    framework == "mle" & models == "discrete" ~ "maximum likelihood-based discrete-time models",
    framework == "mle" & models == "parametric" ~ "correctly-specified maximum likelihood-based models"
  )

  return(labs)
}
