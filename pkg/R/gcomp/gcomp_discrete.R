# """
# File name: gcomp_discrete.R
#
# Description:
#   Function for running the g-computation algorithm under discrete-time
#   parametric models.
# """

gcomp_discrete <- function(data, time_cuts, target_estimand, a_opt, time_horizon = Inf, n_boot = 1000, n_samples = 1000, transform = identity, quiet = TRUE, refresh = 100, ...) {
  # """
  # Run g-computation under discrete-time parametric simulation models
  # fitted using maximum likelihood estimation.
  # """
  dgp_titles <- read_lines("config/DGP/discrete.txt")

  ids <- unique(data$patient_id)
  n_patients <- length(ids)

  n_intervals <- length(time_cuts) - 1
  n_opts <- length(a_opt)
  gcomp_out <- matrix(NA_real_, nrow = n_boot + 1L, ncol = n_opts)

  cat(sprintf("Running g-computation..."), sep = "\n")
  t0 <- Sys.time()

  stanout <- fit_discrete_models(
    data = data,
    time_cuts = time_cuts,
    method = "optimize",
    refresh = 0,
    show_messages = !quiet,
    ...
  )

  pars_list <- extract_cmdstanr_mle(stanout, dgp_titles)

  for (j in seq_len(n_opts)) {
    gcomp_out[1, j] <- simulate_discrete_claim_data(
      n = n_samples,
      n_intervals = n_intervals,
      time_cuts = time_cuts,
      dgp_list = pars_list,
      target_estimand = target_estimand,
      treatment_decision = a_opt[j],
      time_horizon = time_horizon
    )
  }

  # start bootstrap procedure for uncertainty estimation
  for (b in seq_len(n_boot)) {
    if (b == 1 || b == n_boot || b %% refresh == 0) {
      pct <- floor(b / n_boot * 100)
      cat(sprintf("Iteration: %5d / %d [%3d%%]", b, n_boot, pct), sep = "\n")
    }

    repdat <- .select_bootstrap_sample(n_patients, ids, data)

    optout <- fit_discrete_models(
      data = repdat,
      time_cuts = time_cuts,
      method = "optimize",
      refresh = 0,
      show_messages = !quiet,
      ...
    )

    boot_pars_list <- extract_cmdstanr_mle(optout, dgp_titles)

    for (j in seq_len(n_opts)) {
      gcomp_out[b + 1, j] <- simulate_discrete_claim_data(
        n = n_samples,
        n_intervals = n_intervals,
        time_cuts = time_cuts,
        dgp_list = boot_pars_list,
        target_estimand = target_estimand,
        treatment_decision = a_opt[j],
        time_horizon = time_horizon
      )
    }
  }

  elapsed <- difftime(Sys.time(), t0, units = "secs")
  cat(sprintf("Algorithm finished in %0.1f seconds.", elapsed), sep = "\n")

  colnames(gcomp_out) <- a_opt

  if (n_opts == 1) {
    gcomp_out <- as.vector(gcomp_out)
  }

  return(transform(gcomp_out))
}
