#
# Posterior sampling using Stan's MCMC algorithms.
#

#'
#' Posterior Sampling for Joint Waiting-time-to-event and Longitudinal Outcomes
#'
posterior_sampler <- function(gaptime, outcome, data, covars_list, hazard = "weibull", intensity = "gamma-gompertz", init_death = TRUE, method = "sample", iter = 1000, warmup = 1000, refresh = 100, n_cuts = NULL_VAR, quantile_binning = TRUE, id_col = "ID", index_col = "INDEX", time_col = "GAPTIME", status_placeholder = "KEY", event_indicator = 1, death_indicator = 0) {
  convert_to_df <- class(data)[1] != "data.frame"
  invalid_cuts <- is.null(n_cuts)
  piecewise_setup <- (hazard == "piecewise") | (intensity == "piecewise")
  string_indicator <- (typeof(event_indicator) == "character") | (typeof(death_indicator) == "character")

  if (convert_to_df) {
    data <- as.data.frame(data)
  }

  if (string_indicator) {
    event_indicator <- paste0("'", event_indicator, "'")
    death_indicator <- paste0("'", death_indicator, "'")
  }

  f_event <- as.formula(gsub("KEY", event_indicator, deparse1(gaptime)))
  f_death <- as.formula(gsub("KEY", death_indicator, deparse1(gaptime)))

  f_init_event <- .remove_lagged_vars(f_event)
  f_init_death <- .remove_lagged_vars(f_death)
  f_init_cost <- .remove_lagged_vars(outcome)

  init_obs <- data[data[[index_col]] == DEFAULT_START_INDEX, ]
  fu_obs <- data[data[[index_col]] > DEFAULT_START_INDEX, ]

  idx_covars <- seq_along(sim_covars)

  init_time_cuts <- NULL_VAR
  fu_time_cuts <- NULL_VAR

  max_visits <- max(fu_obs[[index_col]], na.rm = TRUE)

  init_weights <- init_obs[[time_col]]
  fu_weights <- fu_obs[[time_col]]

  convert_to_numeric <- (class(init_weights) != "numeric") | (class(fu_weights) != "numeric")

  if (convert_to_numeric) {
    init_weights <- as.numeric(init_weights)
    fu_weights <- as.numeric(fu_weights)
  }

  if (piecewise_setup) {
    if (invalid_cuts) {
      stop("n_cuts_error, invalid number of cuts for the time-axis.")
    }

    init_bound <- max(init_weights, na.rm = TRUE)
    fu_bound <- max(fu_weights, na.rm = TRUE)

    init_time_cuts <- seq(from = 0, to = init_bound, length.out = n_cuts + 1)
    fu_time_cuts <- seq(from = 0, to = fu_bound, length.out = n_cuts + 1)

    # quantile-based binning to create evenly-distributed time cuts
    if (quantile_binning) {
      probs <- seq(from = 0, to = 1, length.out = n_cuts + 1)

      init_time_cuts <- quantile(init_weights, probs = probs, names = FALSE)
      fu_time_cuts <- quantile(fu_weights, probs = probs, names = FALSE)

      init_time_cuts[1] <- DEFAULT_START_TIME
      fu_time_cuts[1] <- DEFAULT_START_TIME
    }
  }

  phreg <- switch(hazard,
    "exponential" = bayesian_exponential_ph,
    "piecewise" = bayesian_piecewise_exponential,
    "weibull" = bayesian_weibull_ph
  )

  pmreg <- switch(intensity,
    "gamma-gompertz" = bayesian_gamma_gompertz_pm,
    "gamma-weibull" = bayesian_gamma_weibull_pm,
    "piecewise" = bayesian_gamma_piecewise_pm
  )

  settings <- list(
    hazard = hazard, intensity = intensity,
    init_death = init_death,
    method = method,
    iter = iter, warmup = warmup,
    max_visits = max_visits,
    init_time_cuts = init_time_cuts, fu_time_cuts = fu_time_cuts
  )

  # posterior sampling for waiting-time-to-events and longitudinal outcome
  samp_out <- .posterior_sampling_joint_components(
    init_obs, fu_obs,
    init_weights, fu_weights,
    phreg, pmreg,
    f_init_event, f_init_death, f_init_cost,
    f_event, f_death, outcome,
    init_death,
    method,
    iter, warmup, refresh,
    init_time_cuts, fu_time_cuts
  )

  # posterior sampling for confounders
  for (j in idx_covars)
  {
    f_init_covar <- covars_list[[j]][["f_init"]]
    f_covar <- covars_list[[j]][["f"]]

    spec <- covars_list[[j]][["spec"]]
    distr <- covars_list[[j]][["distr"]]
    var_col <- covars_list[[j]][["name"]]

    simplereg <- switch(distr,
      "normal" = bayesian_gaussian_simple,
      "bernoulli" = bayesian_bernoulli_simple,
      "categorical" = bayesian_categorical_simple,
      "dtph" = bayesian_bernoulli_simple
    )
    glmreg <- switch(distr,
      "normal" = bayesian_gaussian_glm,
      "bernoulli" = bayesian_bernoulli_glm,
      "categorical" = bayesian_categorical_glm,
      "dtph" = bayesian_bernoulli_glm
    )

    theta_L0 <- simplereg(
      f_init_covar,
      data = init_obs,
      method = method,
      iter = iter, warmup = warmup,
      refresh = refresh
    )

    covars_list[[j]][["init_theta"]] <- theta_L0

    # posterior sampling for time-varying confounders -- if applicable
    if (spec == "time-varying") {
      tempdat <- switch(distr,
        dtph = .construct_discrete_time_to_event_data(fu_obs, id_col, index_col, var_col),
        fu_obs
      )

      theta_L <- glmreg(
        f_covar,
        data = tempdat,
        method = method,
        iter = iter, warmup = warmup,
        refresh = refresh
      )

      covars_list[[j]][["theta"]] <- theta_L
    }
  }

  # store all MCMC results for joint components
  joint <- list(
    event_gaptime = list(
      name = "event",
      f_init = f_init_event,
      f = f_event,
      init_theta = samp_out[["theta_E0"]],
      theta = samp_out[["theta_E"]]
    ),
    death_gaptime = list(
      name = "death",
      f_init = f_init_death,
      f = f_death,
      init_theta = samp_out[["theta_D0"]],
      theta = samp_out[["theta_D"]]
    ),
    outcome = list(
      name = "outcome",
      f_init = f_init_cost,
      f = outcome,
      init_theta = samp_out[["theta_Y0"]],
      theta = samp_out[["theta_Y"]]
    )
  )

  return(list(joint = joint, covars_list = covars_list, settings = settings))
}

# """
# Internal Functions
# """

#'
#' Posterior Sampling for Waiting-time-to-event and Longitudinal Components
#'
.posterior_sampling_joint_components <- function(init_obs, fu_obs, init_weights, fu_weights, phreg, pmreg, f_init_event, f_init_death, f_init_cost, f_event, f_death, outcome, init_death, method, iter, warmup, refresh, init_time_cuts, fu_time_cuts) {
  mcmc_out <- list(
    theta_E0 = phreg(
      f_init_event,
      data = init_obs,
      method = method,
      iter = iter, warmup = warmup,
      refresh = refresh,
      time_cuts = init_time_cuts
    ),
    theta_Y0 = pmreg(
      f_init_cost,
      data = init_obs,
      time_weights = init_weights,
      method = method,
      iter = iter, warmup = warmup,
      refresh = refresh,
      time_cuts = init_time_cuts
    ),
    theta_E = phreg(
      f_event,
      data = fu_obs,
      method = method,
      iter = iter, warmup = warmup,
      refresh = refresh,
      time_cuts = fu_time_cuts
    ),
    theta_D = phreg(
      f_death,
      data = fu_obs,
      method = method,
      iter = iter, warmup = warmup,
      refresh = refresh,
      time_cuts = fu_time_cuts
    ),
    theta_Y = pmreg(
      outcome,
      data = fu_obs,
      time_weights = fu_weights,
      method = method,
      iter = iter, warmup = warmup,
      refresh = refresh,
      time_cuts = fu_time_cuts
    )
  )

  if (init_death) {
    mcmc_out[["theta_D0"]] <- phreg(
      f_init_death,
      data = init_obs,
      method = method,
      iter = iter, warmup = warmup,
      refresh = refresh,
      time_cuts = init_time_cuts
    )
  }

  return(mcmc_out)
}
