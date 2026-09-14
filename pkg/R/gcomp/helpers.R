# """
# File name: helpers.R
#
# Description:
#   Functions used to facilitate the g-computation procedures.
# """

extract_cmdstanr_draw <- function(draws_out, m) {
  # """
  # Extract a single posterior draw as a list of parameters.
  # """
  return(sapply(draws_out, .extract_draw, iter = m))
}

extract_cmdstanr_draws <- function(fit, vars) {
  # """
  # Extract posterior draws as a named list of matrices.
  # """
  return(sapply(vars, .extract_draws, fit = fit))
}

extract_cmdstanr_mle <- function(fit, vars) {
  # """
  # Extract maximum likelihood estimates as a named list.
  # """
  return(sapply(vars, .extract_mle, fit = fit))
}

NMB_gcomp <- function(gcomp_out) {
  # """
  # Compute the net monetary benefit from the g-computation output.
  # """
  return(apply(gcomp_out, 1, diff))
}

# """
# Private functions
# """

.extract_draw <- function(x, iter) {
  # """
  # Extract a single posterior draw.
  # """
  return(x[iter, , drop = TRUE])
}

.extract_draws <- function(fit, var) {
  # """
  # Extract posterior draws.
  # """
  fit$draws(variables = var, format = "matrix")
}

.extract_mle <- function(fit, var) {
  # """
  # Extract maximum likelihood estimates as a named list.
  # """
  fit$mle(variables = var)
}

.select_patient <- function(data, id) {
  # """
  # Internal bootstrap function to select the rows for a patient.
  # """
  idxs <- data$patient_id == id
  return(data[idxs, ])
}

.select_bootstrap_sample <- function(n, ids, data) {
  # """
  # Internal function to create a bootstrap data set.
  # """
  pick <- sample(ids, size = n, replace = TRUE)
  df_list <- lapply(pick, .select_patient, data = data)
  return(do.call(rbind, df_list))
}
