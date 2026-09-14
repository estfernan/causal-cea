# """
# File name: helpers.R
#
# Description:
#   Functions used to facilitate the data-generating procedures.
# """

generate_simulated_dataset <- function(n, P0, L0, V, W, D, L, A, Y, lag_confounders = TRUE) {
  # """
  # Prepare a simulated insurance claims data set for the numerical experiments.
  # """

  # IDs and variable names
  column_ids <- c("patient_id", "claim_idx")
  lag_columns <- c("L1", "L2", "Tx")
  seq_ids <- seq_len(n)
  baseline_vars <- paste0("L0", seq_len(P0))

  base_df <- data.frame(patient_id = seq_ids, L0)
  colnames(base_df) <- c("patient_id", baseline_vars)

  simlist <- list(
    V = pivot_matrix(n, V, "time"),
    W = pivot_matrix(n, W, "gaptime"),
    D = pivot_matrix(n, D, "event_type"),
    A = pivot_matrix(n, A, "Tx"),
    Y = pivot_matrix(n, Y, "cost"),
    L1 = pivot_matrix(n, L[, , 1], "L1"),
    L2 = pivot_matrix(n, L[, , 2], "L2")
  )

  obsdat <- simlist |>
    purrr::reduce(dplyr::left_join, by = column_ids) |>
    dplyr::arrange(patient_id, claim_idx) |>
    dplyr::left_join(base_df, by = "patient_id")

  if (lag_confounders) {
    obsdat <- obsdat |>
      dplyr::group_by(patient_id) |>
      dplyr::mutate(dplyr::across(dplyr::all_of(lag_columns), ~ dplyr::lag(.x, n = 1L), .names = "lag_{col}")) |>
      dplyr::ungroup()
  }

  return(obsdat)
}

hazard_gompertz <- function(t, lambda, alpha) {
  # """
  # Baseline hazard function for the Gompertz distribution.
  # """
  return(lambda * exp(alpha * t))
}

hazard_piecewise <- function(t, time_cuts, lambda) {
  # """
  # Baseline hazard function for the piecewise constant hazard model.
  # """
  idx <- findInterval(t, time_cuts, all.inside = TRUE)
  return(lambda[idx])
}

hazard_weibull <- function(t, lambda, nu) {
  # """
  # Baseline hazard function for the Weibull distribution.
  # """
  return(lambda * nu * t^(nu - 1))
}

intensity <- function(t, a, H, betaA, betaH, base_func, ...) {
  # """
  # Compute the mean (or intensity) specification for the proportional-means model.
  # """
  m0 <- base_func(t, ...) * exp(a * betaA + H %*% betaH)
  return(m0)
}

pivot_matrix <- function(n, mat, value_name) {
  # """
  # Reshape the wide simulated matrices into a long data frame format.
  # """
  df <- as.data.frame(mat) |>
    dplyr::mutate(patient_id = seq_len(n), .before = 1) |>
    tidyr::pivot_longer(
      cols = -patient_id,
      names_to = "claim_idx",
      values_to = value_name,
      names_transform = list(claim_idx = readr::parse_number),
      values_drop_na = TRUE
    )

  return(df)
}
