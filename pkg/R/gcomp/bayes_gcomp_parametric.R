# """
# File name: bayes_gcomp_parametric.R
#
# Description:
#   Function for running the Bayesian g-computation algorithm under
#   correctly-specified parametric models.
# """

bayes_gcomp_parametric <- function(stanout, target_estimand, a_opt, time_horizon = Inf, n_samples = 1000, n_claims = 150, transform = identity, refresh = 100) {
  # """
  # Run Bayesian g-computation using posterior draws from a fitted CmdStanR
  # object for the correctly-specified parametric simulation models.
  # """
  dgp_titles <- read_lines("config/DGP/parametric.txt")

  draws_out <- extract_cmdstanr_draws(stanout, dgp_titles)

  M <- stanout$num_chains() * stanout$metadata()$iter_sampling
  n_opts <- length(a_opt)
  gcomp_out <- matrix(NA_real_, nrow = M, ncol = n_opts)

  cat(sprintf("Running g-computation..."), sep = "\n")
  t0 <- Sys.time()

  for (m in seq_len(M)) {
    if (m == 1 || m == M || m %% refresh == 0) {
      pct <- floor(m / M * 100)
      cat(sprintf("Iteration: %5d / %d [%3d%%]", m, M, pct), sep = "\n")
    }

    pars_list <- extract_cmdstanr_draw(draws_out, m)

    for (j in seq_len(n_opts)) {
      gcomp_out[m, j] <- simulate_claim_data(
        n = n_samples,
        n_claims = n_claims,
        dgp_list = pars_list,
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
