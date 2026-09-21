#
# g-Computation algorithms to approximate target estimands.
#

CLUSTER_EXPORTS <- c(
  "KEYS",
  "DEFAULT_START_TIME",
  "DEFAULT_START_INDEX",
  "TREATMENT_OPTIONS",
  "CENSOR_EVENT_TYPE",
  "REGULAR_EVENT_TYPE",
  "OUTCOME_COLUMNS",
  "MST",
  "MC",
  ".bayesian_gcomputation_dtr",
  ".bayesian_gcomputation_dtr_trajectory",
  ".extract_model_frame",
  ".extract_design_matrix",
  ".extract_joint_design_matrix",
  ".simulate_baseline_confounder",
  ".simulate_timevarying_confounder",
  "predict_bayesian_piecewise_exponential",
  "predict_bayesian_gamma_gompertz_pm",
  "predict_bayesian_gamma_piecewise_pm",
  "predict_bayesian_bernoulli_glm",
  "predict_bayesian_weibull_ph",
  "predict_continuation_ratio",
  "rbernglm",
  "rpiecewiseph",
  "rggpmm",
  "rgppmm",
  "rweibullecon",
  "rweibulleconph"
)

#'
#' g-Computation Procedure
#'
gcomputation <- function(samp_out, target_estimand, a, tau = Inf, n_sim = 1000, n_within = 1, max_iter = 200, ...) {
  hazard <- samp_out[["settings"]][["hazard"]]
  intensity <- samp_out[["settings"]][["intensity"]]
  init_death <- samp_out[["settings"]][["init_death"]]
  iter <- samp_out[["settings"]][["iter"]]
  init_time_cuts <- samp_out[["settings"]][["init_time_cuts"]]
  fu_time_cuts <- samp_out[["settings"]][["fu_time_cuts"]]

  # run g-computation algorithm for each posterior draw
  gcomp_out <- pbapply::pbsapply(
    X = seq_len(iter), FUN = .bayesian_gcomputation,
    samp_out = samp_out,
    target_estimand = target_estimand,
    a = a,
    tau = tau,
    hazard = hazard, intensity = intensity,
    init_death = init_death,
    n_sim = n_sim, n_within = n_within,
    max_iter = max_iter,
    init_time_cuts = init_time_cuts, fu_time_cuts = fu_time_cuts,
    ...
  )

  return(gcomp_out)
}

#'
#' g-Computation Procedure for Dynamic Treatment Regimes
#'
gcomputation_dtr <- function(samp_out, target_estimand, g, tau = Inf, n_sim = 1000, n_within = 1, max_visits = 100, parallel = FALSE, n_cores = 4, simplify = TRUE, transpose = FALSE, ...) {
  hazard <- samp_out[["settings"]][["hazard"]]
  intensity <- samp_out[["settings"]][["intensity"]]
  init_death <- samp_out[["settings"]][["init_death"]]
  chains <- samp_out[["settings"]][["chains"]]
  iter <- samp_out[["settings"]][["iter"]]
  max_cost <- samp_out[["settings"]][["max_cost"]]
  init_time_cuts <- samp_out[["settings"]][["init_time_cuts"]]
  fu_time_cuts <- samp_out[["settings"]][["fu_time_cuts"]]

  simulate_single_dataset <- is.null(target_estimand)

  if (simulate_single_dataset) {
    chains <- 1
    iter <- 1
  }

  if (parallel) {
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)

    # run g-computation algorithm in parallel for each posterior draw
    gcomp_out <- foreach(chain_id = seq_len(chains), .export = CLUSTER_EXPORTS) %dorng% {
      curr <- sapply(
        X = seq_len(iter), FUN = .bayesian_gcomputation_dtr,
        chain_id = chain_id,
        samp_out = samp_out,
        target_estimand = target_estimand,
        g = g,
        tau = tau,
        hazard = hazard, intensity = intensity,
        init_death = init_death,
        n_sim = n_sim, n_within = n_within,
        max_visits = max_visits,
        max_cost = max_cost,
        init_time_cuts = init_time_cuts, fu_time_cuts = fu_time_cuts,
        full_data = FALSE,
        ...
      )

      if (transpose) {
        curr <- t(curr)
      }

      curr
    }

    stopCluster(cl)
  } else {
    gcomp_out <- vector(mode = "list", length = chains)

    # run g-computation algorithm for each chain and posterior draw
    for (chain_id in seq_len(chains)) {
      curr <- pbapply::pbsapply(
        X = seq_len(iter), FUN = .bayesian_gcomputation_dtr,
        chain_id = chain_id,
        samp_out = samp_out,
        target_estimand = target_estimand,
        g = g,
        tau = tau,
        hazard = hazard, intensity = intensity,
        init_death = init_death,
        n_sim = n_sim, n_within = n_within,
        max_visits = max_visits,
        max_cost = max_cost,
        init_time_cuts = init_time_cuts, fu_time_cuts = fu_time_cuts,
        full_data = FALSE,
        ...,
        simplify = simplify
      )

      if (transpose) {
        curr <- t(curr)
      }

      gcomp_out[[chain_id]] <- curr
    }
  }

  if (simulate_single_dataset) {
    gcomp_out <- gcomp_out[[1]][[1]]
  }

  return(gcomp_out)
}

#'
#' g-Computation Procedure for Dynamic Treatment Regimes and Predictive Trajectories
#'
gcomputation_dtr_trajectory <- function(samp_out, basedat, g, tau = Inf, select_chains = NULL, max_visits = 100, parallel = FALSE, n_cores = 4, combine = FALSE, ...) {
  hazard <- samp_out[["settings"]][["hazard"]]
  intensity <- samp_out[["settings"]][["intensity"]]
  init_death <- samp_out[["settings"]][["init_death"]]
  chains <- samp_out[["settings"]][["chains"]]
  iter <- samp_out[["settings"]][["iter"]]
  max_cost <- samp_out[["settings"]][["max_cost"]]
  init_time_cuts <- samp_out[["settings"]][["init_time_cuts"]]
  fu_time_cuts <- samp_out[["settings"]][["fu_time_cuts"]]

  if (is.null(select_chains)) {
    select_chains <- seq_len(chains)
  }

  if (parallel) {
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)

    # run g-computation algorithm in parallel for each posterior draw
    gcomp_out <- foreach(chain_id = select_chains, .export = CLUSTER_EXPORTS) %dorng% {
      curr <- sapply(
        X = seq_len(iter), FUN = .bayesian_gcomputation_dtr_trajectory,
        chain_id = chain_id,
        samp_out = samp_out,
        basedat = basedat,
        g = g,
        tau = tau,
        hazard = hazard, intensity = intensity,
        init_death = init_death,
        n_sim = 1,
        max_visits = max_visits,
        max_cost = max_cost,
        init_time_cuts = init_time_cuts, fu_time_cuts = fu_time_cuts,
        full_data = TRUE,
        ...
      )

      if (combine) {
        curr <- do.call(rbind, curr)
      }

      curr
    }

    stopCluster(cl)
  } else {
    gcomp_out <- vector(mode = "list", length = chains)

    # run g-computation algorithm for each chain and posterior draw
    for (chain_id in select_chains) {
      curr <- pbapply::pblapply(
        X = seq_len(iter), FUN = .bayesian_gcomputation_dtr_trajectory,
        chain_id = chain_id,
        samp_out = samp_out,
        basedat = basedat,
        g = g,
        tau = tau,
        hazard = hazard, intensity = intensity,
        init_death = init_death,
        n_sim = 1,
        max_visits = max_visits,
        max_cost = max_cost,
        init_time_cuts = init_time_cuts, fu_time_cuts = fu_time_cuts,
        full_data = TRUE,
        ...
      )

      if (combine) {
        curr <- do.call(rbind, curr)
      }

      gcomp_out[[chain_id]] <- curr
    }
  }

  if (length(select_chains) == 1) {
    gcomp_out <- gcomp_out[[select_chains]]
  }

  return(gcomp_out)
}

# """
# Internal Functions
# """

#'
#' Bayesian g-Computation for an Individual Posterior Draw
#'
.bayesian_gcomputation <- function(m, samp_out, target_estimand, a, tau, hazard, intensity, init_death, n_sim, n_within, max_iter, init_time_cuts, fu_time_cuts, ...) {
  covars_list <- samp_out[["covars_list"]]

  covar_specs <- sapply(covars_list, "[[", "spec")

  base_idx <- grep("baseline", covar_specs, fixed = TRUE)
  timevary_idx <- grep("time-varying", covar_specs, fixed = TRUE)

  predict_phreg <- switch(hazard,
    "exponential" = predict_bayesian_exponential_ph,
    "piecewise" = predict_bayesian_piecewise_exponential,
    "weibull" = predict_bayesian_weibull_ph
  )

  predict_pmreg <- switch(intensity,
    "gamma-gompertz" = predict_bayesian_gamma_gompertz_pm,
    "gamma-weibull" = predict_bayesian_gamma_weibull_pm,
    "piecewise" = predict_bayesian_gamma_piecewise_pm
  )

  ids <- seq_len(n_sim)
  idxs <- seq_len(max_iter)

  covar_names <- names(covars_list)
  col_names <- c(DEFAULT_COLUMNS, covar_names)
  base_names <- covar_names[base_idx]
  temp_names <- c(KEYS, base_names)

  simdat <- data.table::data.table(
    ID = rep(ids, each = max_iter),
    INDEX = rep(idxs, times = n_sim)
  )

  # initialize an empty 'data.table' object with the specified columns
  for (nm in col_names)
  {
    simdat[, (nm) := numeric()]
  }

  basedat <- simdat[INDEX == 1, ..temp_names]

  # simulate baseline confounders
  for (j in base_idx)
  {
    covname <- covars_list[[j]][["name"]]
    distr <- covars_list[[j]][["distr"]]

    theta <- covars_list[[j]][["init_theta"]][m, , drop = TRUE]

    basedat[, (covname) := .simulate_baseline_confounder(n_sim, distr, theta)]
    simdat[basedat, on = "ID", (covname) := get(paste0("i.", covname))]
  }

  # simulate time-varying confounders at first event
  for (j in timevary_idx)
  {
    covname <- covars_list[[j]][["name"]]
    distr <- covars_list[[j]][["distr"]]
    theta <- covars_list[[j]][["init_theta"]][m, , drop = TRUE]

    L <- .simulate_baseline_confounder(n_sim, distr, theta)

    simdat[INDEX == 1, (covname) := L]
    simdat[INDEX == 2, (paste0("LAG_", covname)) := L]
  }

  simdat[, TREATMENT := (a)]
  simdat[, LAG_TREATMENT := (a)]

  idx <- seq_len(n_sim)
  k <- 1L

  b_curr_index <- simdat[["INDEX"]] == k
  b_next_index <- simdat[["INDEX"]] == k + 1

  b_left <- simdat[["ID"]] %in% idx

  b_curr <- b_curr_index & b_left
  b_next <- b_next_index & b_left

  curr <- simdat[b_curr, ]

  f_event <- samp_out[["joint"]][["event_gaptime"]][["f_init"]][-2]
  theta <- samp_out[["joint"]][["event_gaptime"]][["init_theta"]][m, , drop = TRUE]
  W_V <- predict_phreg(theta, f_event, data = curr, time_cuts = init_time_cuts)

  f_death <- samp_out[["joint"]][["death_gaptime"]][["f_init"]][-2]
  theta <- samp_out[["joint"]][["death_gaptime"]][["init_theta"]][m, , drop = TRUE]
  W_T <- predict_phreg(theta, f_death, data = curr, time_cuts = init_time_cuts)

  if (!init_death) {
    W_T <- rep(Inf, times = n_sim)
  }

  W_tau <- tau
  to_cutoff <- W_tau < W_T

  W_obs <- pmin(W_V, W_T, W_tau)
  to_claim <- W_obs == W_V

  delta <- ifelse(to_claim, REGULAR_EVENT_TYPE, to_cutoff * CENSOR_EVENT_TYPE)

  simdat[b_curr, TIME := DEFAULT_START_TIME]
  simdat[b_curr, GAPTIME := W_obs]
  simdat[b_curr, STATUS := delta]
  simdat[b_next, TIME := W_obs]

  curr <- simdat[b_curr, ]

  f_init <- samp_out[["joint"]][["outcome"]][["f_init"]][-2]
  theta <- samp_out[["joint"]][["outcome"]][["init_theta"]][m, , drop = TRUE]

  Z <- predict_pmreg(theta, f_init, data = curr, time_weights = curr[["GAPTIME"]], time_cuts = init_time_cuts)

  simdat[b_curr, NEXT_INPATIENT_COST := Z]
  simdat[b_curr, BALANCE := Z]

  k <- k + 1

  idx <- which(to_claim)
  n_left <- length(idx)

  while (n_left > 0) {
    b_prev_index <- simdat[["INDEX"]] == k - 1
    b_curr_index <- simdat[["INDEX"]] == k
    b_next_index <- simdat[["INDEX"]] == k + 1

    b_left <- simdat[["ID"]] %in% idx

    b_prev <- b_prev_index & b_left
    b_curr <- b_curr_index & b_left
    b_next <- b_next_index & b_left

    curr <- simdat[b_curr, ]

    # simulate time-varying confounders at follow-up events
    for (j in timevary_idx)
    {
      covname <- covars_list[[j]][["name"]]
      distr <- covars_list[[j]][["distr"]]

      f_covar <- covars_list[[j]][["f"]][-2]
      theta <- covars_list[[j]][["theta"]][m, , drop = TRUE]

      L <- .simulate_timevarying_confounder(curr, distr, theta, f_covar)

      simdat[b_curr, (covname) := L]
      simdat[b_next, (paste0("LAG_", covname)) := L]

      curr <- simdat[b_curr, ]
    }

    f_event <- samp_out[["joint"]][["event_gaptime"]][["f"]][-2]
    theta <- samp_out[["joint"]][["event_gaptime"]][["theta"]][m, , drop = TRUE]
    W_V <- predict_phreg(theta, f_event, data = curr, time_cuts = fu_time_cuts)

    f_death <- samp_out[["joint"]][["death_gaptime"]][["f"]][-2]
    theta <- samp_out[["joint"]][["death_gaptime"]][["theta"]][m, , drop = TRUE]
    W_T <- predict_phreg(theta, f_death, data = curr, time_cuts = fu_time_cuts)

    W_tau <- tau - simdat[b_curr, TIME]

    to_cutoff <- W_tau < W_T

    W_obs <- pmin(W_V, W_T, W_tau)
    to_claim <- W_obs == W_V

    delta <- ifelse(to_claim, REGULAR_EVENT_TYPE, to_cutoff * CENSOR_EVENT_TYPE)

    simdat[b_curr, GAPTIME := W_obs]
    simdat[b_curr, STATUS := delta]
    simdat[b_next, TIME := simdat[b_curr, GAPTIME] + simdat[b_curr, TIME]]

    curr <- simdat[b_curr, ]

    f_outcome <- samp_out[["joint"]][["outcome"]][["f"]][-2]
    theta <- samp_out[["joint"]][["outcome"]][["theta"]][m, , drop = TRUE]

    Z <- predict_pmreg(theta, f_outcome, data = curr, time_weights = curr[["GAPTIME"]], time_cuts = fu_time_cuts)

    simdat[b_curr, NEXT_INPATIENT_COST := Z]
    simdat[b_curr, BALANCE := simdat[b_curr, NEXT_INPATIENT_COST] + simdat[b_prev, BALANCE]]

    k <- k + 1

    idx <- idx[which(to_claim)]
    n_left <- length(idx)
  }

  # check if events exceed storage
  if (k > max_iter) {
    stop("InsufficientStorage: increase 'max_store' to store all observations.")
  }

  simdat <- simdat[!is.na(TIME), ]

  eol <- simdat[, lapply(.SD, max, na.rm = TRUE), by = ID, .SDcols = OUTCOME_COLUMNS]

  T_ <- eol[["TIME"]]
  Y <- eol[["BALANCE"]]

  if (n_within > 1) {
    T_ <- zoo::rollapply(T_, n_within, mean, by = n_within)
    Y <- zoo::rollapply(Y, n_within, mean, by = n_within)
  }

  return(target_estimand(T_, Y, ...))
}

#'
#' Bayesian g-Computation with Dynamic Treatment Regimes for an Individual Posterior Draw
#'
.bayesian_gcomputation_dtr_v0 <- function(m, samp_out, target_estimand, g, tau, hazard, intensity, init_death, n_sim, n_within, max_visits, visit_shift, init_time_cuts, fu_time_cuts, ...) {
  covars_list <- samp_out[["covars_list"]]

  covar_names <- sapply(covars_list, "[[", "name")
  covar_distr <- sapply(covars_list, "[[", "distr")
  covar_specs <- sapply(covars_list, "[[", "spec")

  base_idx <- grep("baseline", covar_specs, fixed = TRUE)
  timevary_idx <- grep("time-varying", covar_specs, fixed = TRUE)

  predict_phreg <- switch(hazard,
    "piecewise" = predict_bayesian_piecewise_exponential,
    "weibull" = predict_bayesian_weibull_ph
  )

  predict_pmreg <- switch(intensity,
    "gamma-gompertz" = predict_bayesian_gamma_gompertz_pm,
    "piecewise" = predict_bayesian_gamma_piecewise_pm
  )

  ids <- seq_len(n_sim)
  idxs <- seq_len(max_visits + visit_shift)

  DEFAULT_COLUMNS <- c("TIME", "GAPTIME", "STATUS", "NEXT_INPATIENT_COST", "BALANCE", "ADJUVANT_TREATMENT")

  covar_names <- names(covars_list)
  col_names <- c(DEFAULT_COLUMNS, covar_names)
  base_names <- covar_names[base_idx]
  temp_names <- c(KEYS, base_names)

  idxs <- factor(idxs)

  simdat <- data.table::data.table(
    ID = rep(ids, each = max_visits + visit_shift),
    INDEX = rep(idxs, times = n_sim)
  )

  # initialize an empty 'data.table' object with the specified columns
  for (nm in col_names)
  {
    if (nm %in% covar_names) {
      if (covar_distr[[nm]] == "categorical") {
        simdat[, (nm) := factor()]
      } else {
        simdat[, (nm) := numeric()]
      }
    } else {
      if (nm == "TREATMENT" | nm == "ADJUVANT_TREATMENT") {
        simdat[, (nm) := factor()]
      } else {
        simdat[, (nm) := numeric()]
      }
    }
  }

  basedat <- simdat[INDEX == 1, ..temp_names]

  # simulate baseline confounders
  for (j in base_idx)
  {
    covname <- covars_list[[j]][["name"]]
    distr <- covars_list[[j]][["distr"]]
    theta <- covars_list[[j]][["init_theta"]][m, , drop = TRUE]

    basedat[, (covname) := .simulate_baseline_confounder(n_sim, distr, theta)]
    simdat[basedat, on = "ID", (covname) := get(paste0("i.", covname))]
  }

  # simulate time-varying confounders at first event
  for (j in timevary_idx)
  {
    covname <- covars_list[[j]][["name"]]
    distr <- covars_list[[j]][["distr"]]
    theta <- covars_list[[j]][["init_theta"]][m, , drop = TRUE]

    L <- .simulate_baseline_confounder(n_sim, distr, theta)

    simdat[INDEX == 1, (covname) := L]
    simdat[INDEX == 2, (paste0("LAG_", covname)) := L]

    # TODO: check if this is still necessary
    # if (distr == "DTPH") {
    #   simdat[ID %in% which(L == 1), (covname) := 1]
    #   simdat[ID %in% which(L == 1), (paste0("LAG_", covname)) := 1]
    # }
  }

  idx <- seq_len(n_sim)
  k <- 1L

  b_curr_index <- simdat[["INDEX"]] == k
  b_next_index <- simdat[["INDEX"]] == k + 1

  b_left <- simdat[["ID"]] %in% idx

  b_curr <- b_curr_index & b_left
  b_next <- b_next_index & b_left

  simdat[b_curr, TIME := DEFAULT_START_TIME]

  a <- g(simdat[b_curr, ])

  # dynamically assign a treatment
  simdat[b_curr, ADJUVANT_TREATMENT := (a)]
  simdat[b_next, LAG_ADJUVANT_TREATMENT := (a)]

  curr <- simdat[b_curr, ]

  f_event <- samp_out[["joint"]][["event_gaptime"]][["f_init"]][-2]
  theta <- samp_out[["joint"]][["event_gaptime"]][["init_theta"]][m, , drop = TRUE]
  W_V <- predict_phreg(theta, f_event, data = curr, time_cuts = init_time_cuts)

  f_death <- samp_out[["joint"]][["death_gaptime"]][["f_init"]][-2]
  theta <- samp_out[["joint"]][["death_gaptime"]][["init_theta"]][m, , drop = TRUE]
  W_T <- predict_phreg(theta, f_death, data = curr, time_cuts = init_time_cuts)

  if (!init_death) {
    W_T <- rep(Inf, times = n_sim)
  }

  W_tau <- tau
  to_cutoff <- W_tau < W_T

  W_obs <- pmin(W_V, W_T, W_tau)
  to_claim <- W_obs == W_V

  delta <- ifelse(to_claim, REGULAR_EVENT_TYPE, to_cutoff * CENSOR_EVENT_TYPE)

  simdat[b_curr, GAPTIME := W_obs]
  simdat[b_curr, STATUS := delta]
  simdat[b_next, TIME := W_obs]
  simdat[b_next, LAG_GAPTIME := W_obs]

  curr <- simdat[b_curr, ]

  f_init <- samp_out[["joint"]][["outcome"]][["f_init"]][-2]
  theta <- samp_out[["joint"]][["outcome"]][["init_theta"]][m, , drop = TRUE]

  Z <- predict_pmreg(theta, f_init, data = curr, time_weights = curr[["GAPTIME"]], time_cuts = init_time_cuts)
  Z <- ifelse(curr[["UPCOMING_TREATMENT"]] == 1, 0, Z)

  simdat[b_curr, NEXT_INPATIENT_COST := Z]
  simdat[b_curr, BALANCE := Z]
  simdat[b_next, LAG_NEXT_INPATIENT_COST := Z]

  k <- k + 1

  idx <- which(to_claim)
  n_left <- length(idx)

  while (n_left > 0) {
    b_prev_index <- simdat[["INDEX"]] == k - 1
    b_curr_index <- simdat[["INDEX"]] == k
    b_next_index <- simdat[["INDEX"]] == k + 1

    b_left <- simdat[["ID"]] %in% idx

    b_prev <- b_prev_index & b_left
    b_curr <- b_curr_index & b_left
    b_next <- b_next_index & b_left

    curr <- simdat[b_curr, ]

    # simulate time-varying confounders at follow-up events
    for (j in timevary_idx)
    {
      covname <- covars_list[[j]][["name"]]
      distr <- covars_list[[j]][["distr"]]

      if (distr == "DTPH") {
        next()
      }

      f_covar <- covars_list[[j]][["f"]][-2]
      theta <- covars_list[[j]][["theta"]][m, , drop = TRUE]

      L <- .simulate_timevarying_confounder(curr, distr, theta, f_covar)

      simdat[b_curr, (covname) := L]
      simdat[b_next, (paste0("LAG_", covname)) := L]

      curr <- simdat[b_curr, ]
    }

    a <- g(curr)

    # dynamically assign a treatment
    simdat[b_curr, ADJUVANT_TREATMENT := (a)]
    simdat[b_next, LAG_ADJUVANT_TREATMENT := (a)]

    curr <- simdat[b_curr, ]

    # simulate scheduling assignment
    for (j in timevary_idx)
    {
      covname <- covars_list[[j]][["name"]]
      distr <- covars_list[[j]][["distr"]]

      if (distr != "DTPH") {
        next()
      }

      f_covar <- covars_list[[j]][["f"]][-2]
      theta <- covars_list[[j]][["theta"]][m, , drop = TRUE]

      f_covar <- as.formula(paste(as.character(f_covar), collapse = " "))

      L <- .simulate_timevarying_confounder(curr, distr, theta, f_covar)

      tmp <- simdat[b_left, sum(.SD, na.rm = TRUE) > 0, by = ID, .SD = covname]
      L_new <- ifelse(tmp[["V1"]], 0, L)

      simdat[b_curr, (covname) := L_new]
      simdat[b_next, (paste0("LAG_", covname)) := L_new]

      # TODO: check if this is still necessary
      # if (distr == "DTPH") {
      #   idx_new <- simdat[b_curr, get("ID")][which(L == 1)]
      #
      #   simdat[ID %in% idx_new & as.numeric(INDEX) > k, (covname) := 1]
      #   simdat[ID %in% idx_new & as.numeric(INDEX) > k, (paste0("LAG_", covname)) := 1]
      # }

      curr <- simdat[b_curr, ]
    }

    f_event <- samp_out[["joint"]][["event_gaptime"]][["f"]][-2]
    theta <- samp_out[["joint"]][["event_gaptime"]][["theta"]][m, , drop = TRUE]
    W_V <- predict_phreg(theta, f_event, data = curr, time_cuts = fu_time_cuts)

    f_death <- samp_out[["joint"]][["death_gaptime"]][["f"]][-2]
    theta <- samp_out[["joint"]][["death_gaptime"]][["theta"]][m, , drop = TRUE]
    W_T <- predict_phreg(theta, f_death, data = curr, time_cuts = fu_time_cuts)

    W_tau <- tau - simdat[b_curr, TIME]

    to_cutoff <- W_tau < W_T

    W_obs <- pmin(W_V, W_T, W_tau)
    to_claim <- W_obs == W_V

    delta <- ifelse(to_claim, REGULAR_EVENT_TYPE, to_cutoff * CENSOR_EVENT_TYPE)

    simdat[b_curr, GAPTIME := W_obs]
    simdat[b_curr, STATUS := delta]
    simdat[b_next, LAG_GAPTIME := W_obs]
    simdat[b_next, TIME := simdat[b_curr, GAPTIME] + simdat[b_curr, TIME]]

    curr <- simdat[b_curr, ]

    f_outcome <- samp_out[["joint"]][["outcome"]][["f"]][-2]
    theta <- samp_out[["joint"]][["outcome"]][["theta"]][m, , drop = TRUE]

    Z <- predict_pmreg(theta, f_outcome, data = curr, time_weights = curr[["GAPTIME"]], time_cuts = fu_time_cuts)
    Z <- ifelse(curr[["UPCOMING_TREATMENT"]] == 1, 0, Z)

    simdat[b_curr, NEXT_INPATIENT_COST := Z]
    simdat[b_curr, BALANCE := simdat[b_curr, NEXT_INPATIENT_COST] + simdat[b_prev, BALANCE]]
    simdat[b_next, LAG_NEXT_INPATIENT_COST := Z]

    k <- k + 1

    drop <- is.infinite(Z)

    if (sum(drop) > 0) {
      simdat <- simdat[!(ID %in% idx[which(drop)]), ]
    }

    idx <- idx[which(to_claim & !drop)]
    n_left <- length(idx)
  }

  # check if events exceed storage
  if (k > max_visits + visit_shift) {
    print("InsufficientStorage: increase 'max_store' to store all observations.")
  }

  simdat <- simdat[!is.na(TIME), ]

  eol <- simdat[, lapply(.SD, max, na.rm = TRUE), by = ID, .SDcols = OUTCOME_COLUMNS]

  T_ <- eol[["TIME"]]
  Y <- eol[["BALANCE"]]

  if (n_within > 1) {
    T_ <- zoo::rollapply(T_, n_within, mean, by = n_within)
    Y <- zoo::rollapply(Y, n_within, mean, by = n_within)
  }

  return(target_estimand(T_, Y, ...))
}

#'
#' Bayesian g-Computation with Dynamic Treatment Regimes for an Individual Posterior Draw
#'
.bayesian_gcomputation_dtr <- function(m, chain_id, samp_out, target_estimand, g, tau, hazard, intensity, init_death, n_sim, n_within, max_visits, max_cost, init_time_cuts, fu_time_cuts, full_data, ...) {
  simulate_single_dataset <- is.null(target_estimand)

  covars_list <- samp_out[["covars_list"]]

  covar_names <- sapply(covars_list, "[[", "name")
  covar_distr <- sapply(covars_list, "[[", "distr")
  covar_specs <- sapply(covars_list, "[[", "spec")

  base_idx <- grep("baseline", covar_specs, fixed = TRUE)
  timevary_idx <- grep("time-varying", covar_specs, fixed = TRUE)

  predict_phreg <- switch(hazard,
    "piecewise" = predict_bayesian_piecewise_exponential,
    "weibull" = predict_bayesian_weibull_ph
  )

  predict_pmreg <- switch(intensity,
    "gamma-gompertz" = predict_bayesian_gamma_gompertz_pm,
    "piecewise" = predict_bayesian_gamma_piecewise_pm
  )

  ids <- seq_len(n_sim)
  idxs <- seq_len(max_visits)

  DEFAULT_COLUMNS <- c("TIME", "GAPTIME", "STATUS", "DEATH", "UPCOMING_COST", "BALANCE", "TREATMENT_READY", "ADJUVANT_TREATMENT")

  covar_names <- names(covars_list)
  col_names <- c(DEFAULT_COLUMNS, covar_names)
  base_names <- covar_names[base_idx]
  temp_names <- c(KEYS, base_names)

  simdat <- data.table::data.table(
    ID = rep(ids, each = max_visits),
    INDEX = rep(idxs, times = n_sim)
  )

  # initialize an empty 'data.table' object with the specified columns
  for (nm in col_names)
  {
    if (nm %in% covar_names) {
      if (covar_distr[[nm]] == "categorical") {
        simdat[, (nm) := factor()]
      } else {
        simdat[, (nm) := numeric()]
      }
    } else {
      if (nm == "ADJUVANT_TREATMENT") {
        simdat[, (nm) := factor()]
      } else {
        simdat[, (nm) := numeric()]
      }
    }
  }

  basedat <- simdat[INDEX == 1, ..temp_names]

  # simulate baseline confounders
  for (j in base_idx)
  {
    covname <- covars_list[[j]][["name"]]
    distr <- covars_list[[j]][["distr"]]
    theta <- covars_list[[j]][["init_theta"]][m, chain_id, , drop = TRUE]

    if (simulate_single_dataset) {
      theta <- apply(covars_list[[j]][["init_theta"]], 3, mean)
    }

    basedat[, (covname) := .simulate_baseline_confounder(n_sim, distr, theta)]
    simdat[basedat, on = "ID", (covname) := get(paste0("i.", covname))]
  }

  # simulate time-varying confounders at first event
  for (j in timevary_idx)
  {
    covname <- covars_list[[j]][["name"]]
    distr <- covars_list[[j]][["distr"]]
    f_covar <- covars_list[[j]][["f_init"]][-2]
    theta <- covars_list[[j]][["init_theta"]][m, chain_id, , drop = TRUE]

    if (simulate_single_dataset) {
      theta <- apply(covars_list[[j]][["init_theta"]], 3, mean)
    }

    L <- .simulate_timevarying_confounder(basedat, distr, theta, f_covar)

    simdat[INDEX == 1, (covname) := L]
    simdat[INDEX == 2, (paste0("LAG_", covname)) := L]
  }

  simdat[INDEX == 1, TREATMENT_READY := 0]
  simdat[INDEX == 2, LAG_TREATMENT_READY := 0]

  idx <- seq_len(n_sim)
  k <- 1L

  b_curr_index <- simdat[["INDEX"]] == k
  b_next_index <- simdat[["INDEX"]] == k + 1

  b_left <- simdat[["ID"]] %in% idx

  b_curr <- b_curr_index & b_left
  b_next <- b_next_index & b_left

  simdat[b_curr, TIME := DEFAULT_START_TIME]

  a <- g(simdat[b_curr, ])

  # dynamically assign a treatment
  simdat[b_curr, ADJUVANT_TREATMENT := (a)]
  simdat[b_next, LAG_ADJUVANT_TREATMENT := (a)]

  curr <- simdat[b_curr, ]

  f_event <- samp_out[["joint"]][["event_gaptime"]][["f_init"]][-2]
  theta <- samp_out[["joint"]][["event_gaptime"]][["init_theta"]][m, chain_id, , drop = TRUE]

  if (simulate_single_dataset) {
    theta <- apply(samp_out[["joint"]][["event_gaptime"]][["init_theta"]], 3, mean)
  }

  W_V <- predict_phreg(theta, f_event, data = curr, time_cuts = init_time_cuts)

  f_death <- samp_out[["joint"]][["death_gaptime"]][["f_init"]][-2]
  theta <- samp_out[["joint"]][["death_gaptime"]][["init_theta"]][m, chain_id, , drop = TRUE]

  if (simulate_single_dataset) {
    theta <- apply(samp_out[["joint"]][["death_gaptime"]][["init_theta"]], 3, mean)
  }

  W_T <- predict_phreg(theta, f_death, data = curr, time_cuts = init_time_cuts)

  if (!init_death) {
    W_T <- rep(Inf, times = n_sim)
  }

  W_tau <- tau
  to_cutoff <- W_tau < W_T

  W_obs <- pmin(W_V, W_T, W_tau)
  to_claim <- W_obs == W_V

  delta <- ifelse(to_claim, REGULAR_EVENT_TYPE, to_cutoff * CENSOR_EVENT_TYPE)

  simdat[b_curr, GAPTIME := W_obs]
  simdat[b_curr, STATUS := delta]
  simdat[b_curr, DEATH := (delta == 0)]
  simdat[b_next, TIME := W_obs]
  simdat[b_next, LAG_GAPTIME := W_obs]

  curr <- simdat[b_curr, ]

  f_init <- samp_out[["joint"]][["outcome"]][["f_init"]][-2]
  theta <- samp_out[["joint"]][["outcome"]][["init_theta"]][m, chain_id, , drop = TRUE]

  if (simulate_single_dataset) {
    theta <- apply(samp_out[["joint"]][["outcome"]][["init_theta"]], 3, mean)
  }

  Z <- predict_pmreg(theta, f_init, data = curr, time_weights = curr[["GAPTIME"]], time_cuts = init_time_cuts)

  simdat[b_curr, UPCOMING_COST := Z]
  simdat[b_curr, BALANCE := Z]
  simdat[b_next, LAG_UPCOMING_COST := Z]

  k <- k + 1

  idx <- which(to_claim)
  n_left <- length(idx)

  while (n_left > 0) {
    b_prev_index <- simdat[["INDEX"]] == k - 1
    b_curr_index <- simdat[["INDEX"]] == k
    b_next_index <- simdat[["INDEX"]] == k + 1

    b_left <- simdat[["ID"]] %in% idx

    b_prev <- b_prev_index & b_left
    b_curr <- b_curr_index & b_left
    b_next <- b_next_index & b_left

    curr <- simdat[b_curr, ]

    # simulate time-varying confounders at follow-up events
    for (j in timevary_idx)
    {
      covname <- covars_list[[j]][["name"]]
      distr <- covars_list[[j]][["distr"]]

      f_covar <- covars_list[[j]][["f"]][-2]
      theta <- covars_list[[j]][["theta"]][m, chain_id, , drop = TRUE]

      if (simulate_single_dataset) {
        theta <- apply(covars_list[[j]][["theta"]], 3, mean)
      }

      L <- .simulate_timevarying_confounder(curr, distr, theta, f_covar)

      simdat[b_curr, (covname) := L]
      simdat[b_next, (paste0("LAG_", covname)) := L]

      curr <- simdat[b_curr, ]
    }

    f_initiate <- samp_out[["joint"]][["treatment_initiation"]][["f"]][-2]
    theta <- samp_out[["joint"]][["treatment_initiation"]][["theta"]][m, chain_id, , drop = TRUE]

    if (simulate_single_dataset) {
      theta <- apply(samp_out[["joint"]][["treatment_initiation"]][["theta"]], 3, mean)
    }

    S <- predict_continuation_ratio(theta, f_initiate, data = curr, index_weights = curr[["INDEX"]])

    # check if treatment has been initiated already
    S <- ifelse(curr[["LAG_TREATMENT_READY"]] == 0, S, 1)

    simdat[b_curr, TREATMENT_READY := S]
    simdat[b_next, LAG_TREATMENT_READY := S]

    curr <- simdat[b_curr, ]

    a <- g(curr)

    # dynamically assign a treatment
    simdat[b_curr, ADJUVANT_TREATMENT := (a)]
    simdat[b_next, LAG_ADJUVANT_TREATMENT := (a)]

    curr <- simdat[b_curr, ]

    f_event <- samp_out[["joint"]][["event_gaptime"]][["f"]][-2]
    theta <- samp_out[["joint"]][["event_gaptime"]][["theta"]][m, chain_id, , drop = TRUE]

    if (simulate_single_dataset) {
      theta <- apply(samp_out[["joint"]][["event_gaptime"]][["theta"]], 3, mean)
    }

    W_V <- predict_phreg(theta, f_event, data = curr, time_cuts = fu_time_cuts)

    f_death <- samp_out[["joint"]][["death_gaptime"]][["f"]][-2]
    theta <- samp_out[["joint"]][["death_gaptime"]][["theta"]][m, chain_id, , drop = TRUE]

    if (simulate_single_dataset) {
      theta <- apply(samp_out[["joint"]][["death_gaptime"]][["theta"]], 3, mean)
    }

    W_T <- predict_phreg(theta, f_death, data = curr, time_cuts = fu_time_cuts)

    W_tau <- tau - simdat[b_curr, TIME]

    to_cutoff <- W_tau < W_T

    W_obs <- pmin(W_V, W_T, W_tau)
    to_claim <- W_obs == W_V

    delta <- ifelse(to_claim, REGULAR_EVENT_TYPE, to_cutoff * CENSOR_EVENT_TYPE)

    simdat[b_curr, GAPTIME := W_obs]
    simdat[b_curr, STATUS := delta]
    simdat[b_curr, DEATH := (delta == 0)]
    simdat[b_next, LAG_GAPTIME := W_obs]
    simdat[b_next, TIME := simdat[b_curr, GAPTIME] + simdat[b_curr, TIME]]

    curr <- simdat[b_curr, ]

    f_outcome <- samp_out[["joint"]][["outcome"]][["f"]][-2]
    theta <- samp_out[["joint"]][["outcome"]][["theta"]][m, chain_id, , drop = TRUE]

    if (simulate_single_dataset) {
      theta <- apply(samp_out[["joint"]][["outcome"]][["theta"]], 3, mean)
    }

    Z <- predict_pmreg(theta, f_outcome, data = curr, time_weights = curr[["GAPTIME"]], time_cuts = fu_time_cuts)

    simdat[b_curr, UPCOMING_COST := Z]
    simdat[b_curr, BALANCE := simdat[b_prev, BALANCE] + Z]
    simdat[b_next, LAG_UPCOMING_COST := Z]

    k <- k + 1

    drop <- is.infinite(Z) | simdat[b_curr, BALANCE] > max_cost

    if (sum(drop) > 0) {
      simdat <- simdat[!(ID %in% idx[which(drop)]), ]
    }

    idx <- idx[which(to_claim & !drop)]
    n_left <- length(idx)
  }

  if (full_data) {
    # TODO: need to implement for posterior trajectories
    browser()
  }

  simdat <- simdat[!is.na(TIME), ]

  if (simulate_single_dataset) {
    return(simdat)
  }

  eol <- simdat[, lapply(.SD, max, na.rm = TRUE), by = ID, .SDcols = OUTCOME_COLUMNS]

  T_ <- eol[["TIME"]]
  Y <- eol[["BALANCE"]]

  if (n_within > 1) {
    T_ <- zoo::rollapply(T_, n_within, mean, by = n_within)
    Y <- zoo::rollapply(Y, n_within, mean, by = n_within)
  }

  return(target_estimand(T_, Y, ...))
}

#'
#' Bayesian g-Computation with Dynamic Treatment Regimes for an Individual Predictive Trajectory
#'
.bayesian_gcomputation_dtr_trajectory <- function(m, chain_id, samp_out, basedat, g, tau, hazard, intensity, init_death, n_sim, max_visits, max_cost, init_time_cuts, fu_time_cuts, full_data, ...) {
  covars_list <- samp_out[["covars_list"]]

  covar_names <- sapply(covars_list, "[[", "name")
  covar_distr <- sapply(covars_list, "[[", "distr")
  covar_specs <- sapply(covars_list, "[[", "spec")

  base_idx <- grep("baseline", covar_specs, fixed = TRUE)
  timevary_idx <- grep("time-varying", covar_specs, fixed = TRUE)

  predict_phreg <- switch(hazard,
    "piecewise" = predict_bayesian_piecewise_exponential,
    "weibull" = predict_bayesian_weibull_ph
  )

  predict_pmreg <- switch(intensity,
    "gamma-gompertz" = predict_bayesian_gamma_gompertz_pm,
    "piecewise" = predict_bayesian_gamma_piecewise_pm
  )

  ids <- seq_len(n_sim)
  idxs <- seq_len(max_visits)

  DEFAULT_COLUMNS <- c("TIME", "GAPTIME", "STATUS", "DEATH", "UPCOMING_COST", "BALANCE", "TREATMENT_READY", "ADJUVANT_TREATMENT")

  covar_names <- names(covars_list)
  col_names <- c(DEFAULT_COLUMNS, covar_names)
  base_names <- covar_names[base_idx]
  temp_names <- c(KEYS, base_names)

  simdat <- data.table::data.table(
    ID = rep(ids, each = max_visits),
    INDEX = rep(idxs, times = n_sim)
  )

  # initialize an empty 'data.table' object with the specified columns
  for (nm in col_names)
  {
    if (nm %in% covar_names) {
      if (covar_distr[[nm]] == "categorical") {
        simdat[, (nm) := factor()]
      } else {
        simdat[, (nm) := numeric()]
      }
    } else {
      if (nm == "ADJUVANT_TREATMENT") {
        simdat[, (nm) := factor()]
      } else {
        simdat[, (nm) := numeric()]
      }
    }
  }

  # simulate baseline confounders
  for (j in base_idx)
  {
    covname <- covars_list[[j]][["name"]]
    simdat[basedat, on = "ID", (covname) := get(paste0("i.", covname))]
  }

  # simulate time-varying confounders at first event
  for (j in timevary_idx)
  {
    covname <- covars_list[[j]][["name"]]
    distr <- covars_list[[j]][["distr"]]
    f_covar <- covars_list[[j]][["f_init"]][-2]
    theta <- covars_list[[j]][["init_theta"]][m, chain_id, , drop = TRUE]

    L <- .simulate_timevarying_confounder(basedat, distr, theta, f_covar)

    simdat[INDEX == 1, (covname) := L]
    simdat[INDEX == 2, (paste0("LAG_", covname)) := L]
  }

  simdat[INDEX == 1, TREATMENT_READY := 0]
  simdat[INDEX == 2, LAG_TREATMENT_READY := 0]

  idx <- seq_len(n_sim)
  k <- 1L

  b_curr_index <- simdat[["INDEX"]] == k
  b_next_index <- simdat[["INDEX"]] == k + 1

  b_left <- simdat[["ID"]] %in% idx

  b_curr <- b_curr_index & b_left
  b_next <- b_next_index & b_left

  simdat[b_curr, TIME := DEFAULT_START_TIME]

  a <- g(simdat[b_curr, ])

  # dynamically assign a treatment
  simdat[b_curr, ADJUVANT_TREATMENT := (a)]
  simdat[b_next, LAG_ADJUVANT_TREATMENT := (a)]

  curr <- simdat[b_curr, ]

  f_event <- samp_out[["joint"]][["event_gaptime"]][["f_init"]][-2]
  theta <- samp_out[["joint"]][["event_gaptime"]][["init_theta"]][m, chain_id, , drop = TRUE]
  W_V <- predict_phreg(theta, f_event, data = curr, time_cuts = init_time_cuts)

  f_death <- samp_out[["joint"]][["death_gaptime"]][["f_init"]][-2]
  theta <- samp_out[["joint"]][["death_gaptime"]][["init_theta"]][m, chain_id, , drop = TRUE]
  W_T <- predict_phreg(theta, f_death, data = curr, time_cuts = init_time_cuts)

  if (!init_death) {
    W_T <- rep(Inf, times = n_sim)
  }

  W_tau <- tau
  to_cutoff <- W_tau < W_T

  W_obs <- pmin(W_V, W_T, W_tau)
  to_claim <- W_obs == W_V

  delta <- ifelse(to_claim, REGULAR_EVENT_TYPE, to_cutoff * CENSOR_EVENT_TYPE)

  simdat[b_curr, GAPTIME := W_obs]
  simdat[b_curr, STATUS := delta]
  simdat[b_curr, DEATH := (delta == 0)]
  simdat[b_next, TIME := W_obs]
  simdat[b_next, LAG_GAPTIME := W_obs]

  curr <- simdat[b_curr, ]

  f_init <- samp_out[["joint"]][["outcome"]][["f_init"]][-2]
  theta <- samp_out[["joint"]][["outcome"]][["init_theta"]][m, chain_id, , drop = TRUE]
  Z <- predict_pmreg(theta, f_init, data = curr, time_weights = curr[["GAPTIME"]], time_cuts = init_time_cuts)

  simdat[b_curr, UPCOMING_COST := Z]
  simdat[b_curr, BALANCE := Z]
  simdat[b_next, LAG_UPCOMING_COST := Z]

  k <- k + 1

  idx <- which(to_claim)
  n_left <- length(idx)

  while (n_left > 0) {
    b_prev_index <- simdat[["INDEX"]] == k - 1
    b_curr_index <- simdat[["INDEX"]] == k
    b_next_index <- simdat[["INDEX"]] == k + 1

    b_left <- simdat[["ID"]] %in% idx

    b_prev <- b_prev_index & b_left
    b_curr <- b_curr_index & b_left
    b_next <- b_next_index & b_left

    curr <- simdat[b_curr, ]

    # simulate time-varying confounders at follow-up events
    for (j in timevary_idx)
    {
      covname <- covars_list[[j]][["name"]]
      distr <- covars_list[[j]][["distr"]]

      f_covar <- covars_list[[j]][["f"]][-2]
      theta <- covars_list[[j]][["theta"]][m, chain_id, , drop = TRUE]

      L <- .simulate_timevarying_confounder(curr, distr, theta, f_covar)

      simdat[b_curr, (covname) := L]
      simdat[b_next, (paste0("LAG_", covname)) := L]

      curr <- simdat[b_curr, ]
    }

    f_initiate <- samp_out[["joint"]][["treatment_initiation"]][["f"]][-2]
    theta <- samp_out[["joint"]][["treatment_initiation"]][["theta"]][m, chain_id, , drop = TRUE]
    S <- predict_continuation_ratio(theta, f_initiate, data = curr, index_weights = curr[["INDEX"]])

    # check if treatment has been initiated already
    S <- ifelse(curr[["LAG_TREATMENT_READY"]] == 0, S, 1)

    simdat[b_curr, TREATMENT_READY := S]
    simdat[b_next, LAG_TREATMENT_READY := S]

    curr <- simdat[b_curr, ]

    a <- g(curr)

    # dynamically assign a treatment
    simdat[b_curr, ADJUVANT_TREATMENT := (a)]
    simdat[b_next, LAG_ADJUVANT_TREATMENT := (a)]

    curr <- simdat[b_curr, ]

    f_event <- samp_out[["joint"]][["event_gaptime"]][["f"]][-2]
    theta <- samp_out[["joint"]][["event_gaptime"]][["theta"]][m, chain_id, , drop = TRUE]
    W_V <- predict_phreg(theta, f_event, data = curr, time_cuts = fu_time_cuts)

    f_death <- samp_out[["joint"]][["death_gaptime"]][["f"]][-2]
    theta <- samp_out[["joint"]][["death_gaptime"]][["theta"]][m, chain_id, , drop = TRUE]
    W_T <- predict_phreg(theta, f_death, data = curr, time_cuts = fu_time_cuts)

    W_tau <- tau - simdat[b_curr, TIME]

    to_cutoff <- W_tau < W_T

    W_obs <- pmin(W_V, W_T, W_tau)
    to_claim <- W_obs == W_V

    delta <- ifelse(to_claim, REGULAR_EVENT_TYPE, to_cutoff * CENSOR_EVENT_TYPE)

    simdat[b_curr, GAPTIME := W_obs]
    simdat[b_curr, STATUS := delta]
    simdat[b_curr, DEATH := (delta == 0)]
    simdat[b_next, LAG_GAPTIME := W_obs]
    simdat[b_next, TIME := simdat[b_curr, GAPTIME] + simdat[b_curr, TIME]]

    curr <- simdat[b_curr, ]

    f_outcome <- samp_out[["joint"]][["outcome"]][["f"]][-2]
    theta <- samp_out[["joint"]][["outcome"]][["theta"]][m, chain_id, , drop = TRUE]
    Z <- predict_pmreg(theta, f_outcome, data = curr, time_weights = curr[["GAPTIME"]], time_cuts = fu_time_cuts)

    simdat[b_curr, UPCOMING_COST := Z]
    simdat[b_curr, BALANCE := simdat[b_prev, BALANCE] + Z]
    simdat[b_next, LAG_UPCOMING_COST := Z]

    k <- k + 1

    drop <- is.infinite(Z) | simdat[b_curr, BALANCE] > max_cost

    if (sum(drop) > 0) {
      simdat <- simdat[!(ID %in% idx[which(drop)]), ]
    }

    idx <- idx[which(to_claim & !drop)]
    n_left <- length(idx)
  }

  simdat[, BALANCE := data.table::shift(BALANCE, fill = 0, type = "lag"), by = ID]
  simdat[, EVENT := data.table::shift(STATUS, fill = 1, type = "lag"), by = ID]

  simdat <- simdat[!is.na(TIME), ]

  W <- data.table::shift(simdat[["GAPTIME"]], fill = 0, type = "lag")
  Y <- data.table::shift(simdat[["UPCOMING_COST"]], fill = 0, type = "lag")

  T_ <- simdat[["TIME"]]
  U <- simdat[["BALANCE"]]
  dt <- simdat[["EVENT"]]

  J <- length(T_)

  return(cbind(m, j = 1:J, W = W, Y = Y, T_, U, dt))
}

#'
#' Generate a Baseline Confounder
#'
.simulate_baseline_confounder <- function(n, distr, theta) {
  if (distr == "bernoulli") {
    x <- LaplacesDemon::rbern(n, prob = theta)
  } else if (distr == "normal") {
    x <- rnorm(n, mean = theta[1], sd = theta[2])
  } else if (distr == "categorical") {
    x <- factor(LaplacesDemon::rcat(n, p = theta), levels = 1:length(theta))
  }

  return(x)
}

#'
#' Generate a Time-varying Confounder
#'
.simulate_timevarying_confounder <- function(data, distr, theta, f_covar) {
  if (distr == "bernoulli") {
    x <- predict_bayesian_bernoulli_glm(theta, f_covar, data = data)
  } else if (distr == "normal") {
    x <- predict_bayesian_gaussian_glm(theta, f_covar, data = data)
  } else if (distr == "categorical") {
    x <- predict_bayesian_categorical_glm(theta, f_covar, data = data)
  }

  return(x)
}
