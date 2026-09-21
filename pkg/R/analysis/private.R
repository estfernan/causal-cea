#
# Internal procedures used within project functions.
#

# """
# Design Matrices
# """

#'
#' Construct Subject-level for Discrete-Time Hazard Modeling
#'
.construct_discrete_time_to_event_data_individual <- function(data, index_column, var_column) {
  W <- data[[var_column]]

  pos <- which(W > 0)
  N <- length(W)

  keep <- seq(from = 1, to = ifelse(sum(W) > 0, pos, N), by = 1)

  return(data[keep, ])
}

#'
#' Construct Observations for Discrete-Time Hazard Modeling
#'
.construct_discrete_time_to_event_data <- function(data, id_column, index_column, var_column) {
  ids <- data[[id_column]]

  datlist <- lapply(
    split(data, ids),
    .construct_discrete_time_to_event_data_individual,
    index_column = index_column,
    var_column = var_column
  )

  tte <- do.call(rbind, datlist)

  return(tte[tte[[index_column]] > 1, ])
}

#'
#' Construct Observations for Discrete-Time Hazard Modeling
#'
setup_discrete_time_to_event_data <- function(data, id_column, index_column, var_column) {
  ids <- data[[id_column]]

  datlist <- lapply(
    split(data, ids),
    .construct_discrete_time_to_event_data_individual,
    index_column = index_column,
    var_column = var_column
  )

  return(do.call(rbind, datlist))
}

#'
#' Construct Design Matrices or Vectors
#'
.design <- function(n, ..., int = 1) {
  func <- ifelse(n > 1, cbind, c)
  return(func(int, ...))
}

# """
# Model Formulae
# """

#'
#' Remove Lagged Variables from Model Formula
#'
.remove_lagged_vars <- function(formula, prefix = "^LAG_") {
  covars <- attr(terms(formula), "term.labels")
  keep <- covars[!grepl(prefix, covars)]

  lhs <- as.character(rlang::f_lhs(formula))

  if ("Surv" %in% lhs) {
    lhs <- paste0(lhs[1], "(", lhs[2], ", ", lhs[3], ")")
  }

  terms <- paste(keep, collapse = " + ")
  newfm <- paste(lhs, terms, sep = " ~ ")

  return(as.formula(newfm))
}

# """
# Model Frames
# """

#'
#' Construct a Design Matrix from Model Objects
#'
.extract_design_matrix <- function(x, ..., intercept = TRUE) {
  X <- model.matrix(x, data = x)

  if (!intercept) {
    X <- X[, -1, drop = FALSE]
  }

  n <- nrow(X)
  p <- ncol(X)

  return(list(n = n, p = p, X = X, ...))
}

#'
#' Construct a Design Matrix and Time Outcomes from Model Objects
#'
.extract_joint_design_matrix <- function(x, t, ..., intercept = FALSE) {
  X <- model.matrix(x, data = x)

  if (!intercept) {
    X <- X[, -1, drop = FALSE]
  }

  n <- nrow(X)
  p <- ncol(X)

  return(list(n = n, p = p, t = t, X = X, ...))
}

#'
#' Extract Components for Joint Modeling
#'
.extract_joint_objects <- function(x, t, ..., intercept = FALSE, non_zero_check = TRUE, epsilon = 1e-5) {
  idx_missing <- attr(x, "na.action")

  y <- model.response(x, type = "double")
  X <- model.matrix(x, data = x)

  if (non_zero_check) {
    y[y <= 0] <- epsilon
  }

  if (!intercept) {
    X <- X[, -1, drop = FALSE]
  }

  drop_times <- length(idx_missing) > 0

  if (drop_times) {
    t <- t[-idx_missing]
  }

  n <- nrow(X)
  p <- ncol(X)

  return(list(n = n, p = p, y = y, t = t, X = X, ...))
}

#'
#' Extract the Model Frame from a Function Call
#'
.extract_model_frame <- function(f_call, env) {
  args <- c("formula", "data", "subset", "weights", "na.action", "offset")

  idx <- match(args, names(f_call), nomatch = 0L)

  f_call <- f_call[c(1L, idx)]
  f_call[[1]] <- quote(stats::model.frame)

  mf <- eval(f_call, envir = env)

  return(mf)
}

#'
#' Construct Model Objects for Regression Analysis
#'
.extract_regression_objects <- function(x, ..., categorical = FALSE, intercept = TRUE) {
  y <- model.response(x, type = "any")
  X <- model.matrix(x, data = x)

  if (!intercept) {
    X <- X[, -1, drop = FALSE]
  }

  n <- nrow(X)
  p <- ncol(X)

  datalist <- list(n = n, p = p, X = X, y = y, ...)

  if (categorical) {
    K <- nlevels(y)
    datalist <- c(datalist, K = K)
  }

  return(datalist)
}

#'
#' Construct Model for Discrete-time Hazard Analysis
#'
.extract_discrete_tte_objects <- function(x, index_weights, ..., intercept = FALSE) {
  y <- model.response(x, type = "double")
  X <- model.matrix(x, data = x)

  pos <- attr(x, "na.action")

  if (length(pos) > 0) {
    index_weights <- index_weights[-pos]
  }

  n <- nrow(X)
  p <- ncol(X)

  K <- max(index_weights)

  return(list(n = n, p = p, K = K, X = X, y = y, idx = index_weights, ...))
}

#'
#' Construct Model Sub-Objects for Piecewise Joint Outcome Analysis
#'
.extract_split_joint_objects <- function(x, t, time_cuts, ..., intercept = FALSE, non_zero_check = TRUE, epsilon = 1e-5) {
  y <- model.response(x, type = "double")
  X <- model.matrix(x, data = x)

  if (non_zero_check) {
    y[y <= 0] <- epsilon
  }

  pos <- attr(x, "na.action")

  if (length(pos) > 0) {
    t <- t[-pos]
  }

  # TODO: Remove this feature
  # if (!intercept) {
  #   X <- X[, -1, drop = FALSE]
  # }

  idx <- findInterval(t, time_cuts, all.inside = TRUE)

  n <- nrow(X)
  p <- ncol(X)

  K <- length(time_cuts) - 1

  return(list(n = n, p = p, K = K, y = y, t = t, X = X, idx = idx, ...))
}

#'
#' Construct Model Sub-Objects for Piecewise Survival Analysis
#'
.extract_split_survival_objects <- function(x, time_cuts, ..., intercept = FALSE) {
  y <- model.response(x, type = "double")
  X <- model.matrix(x, data = x)

  # TODO: Remove this feature
  # if (!intercept) {
  #   X <- X[, -1, drop = FALSE]
  # }

  t <- y[, 1, drop = TRUE]
  dt <- y[, 2, drop = TRUE]

  surv_out <- survival::Surv(t, dt)

  df <- survival::survSplit(
    surv_out ~ .,
    data = data.frame(surv_out, X, check.names = FALSE),
    subset = NULL,
    cut = time_cuts,
    id = "id", episode = "interval",
    event = "delta"
  )

  H <- df[["surv_out"]][, "stop"] - df[["surv_out"]][, "start"]
  dt <- df[["surv_out"]][, "status"]

  idx <- df[["interval"]] - 1

  X <- df[, colnames(X)]

  n <- nrow(X)
  p <- ncol(X)

  K <- length(time_cuts) - 1

  return(list(n = n, p = p, K = K, H = H, X = X, delta = dt, idx = idx, ...))
}

#'
#' Extract Components for Survival Modeling
#'
.extract_survival_objects <- function(x, ..., intercept = FALSE) {
  y <- model.response(x, type = "double")
  X <- model.matrix(x, data = x)

  if (!intercept) {
    X <- X[, -1, drop = FALSE]
  }

  t <- y[, 1, drop = TRUE]
  dt <- y[, 2, drop = TRUE]

  n <- nrow(X)
  p <- ncol(X)

  return(list(n = n, p = p, t = t, delta = dt, X = X, ...))
}

#'
#' Construct Model Objects for Univariable Analysis
#'
.extract_univariable_objects <- function(x, ..., categorical = FALSE) {
  y <- model.response(x, type = "any")
  n <- length(y)

  datalist <- list(n = n, y = y, ...)

  if (categorical) {
    K <- nlevels(y)
    datalist <- c(datalist, K = K)
  }

  return(datalist)
}

# """
# Stan
# """

#'
#' Create a New Stan Model Object using CmdStanR
#'
.init_cmdstanr <- function(filename, stan_path, ..., parent_dir = getwd()) {
  path <- file.path(parent_dir, stan_path)
  exe_file <- file.path(path, filename)
  stan_file <- paste0(exe_file, STAN_EXTENSION)

  # Windows cannot handle compilation paths with spaces due to processx limitations
  is_windows <- .Platform$OS.type == "windows"
  has_spaces <- grepl(" ", exe_file, fixed = TRUE)

  if (is_windows && has_spaces) {
    mod <- cmdstanr::cmdstan_model(
      stan_file = stan_file,
      dir = tempdir(),
      ...
    )
  } else {
    mod <- cmdstanr::cmdstan_model(
      stan_file = stan_file,
      exe_file = exe_file,
      ...
    )
  }

  return(mod)
}

#'
#' Create a New Stan Model Object using CmdStanR
#'
# .init_cmdstanr <- function(filename, stan_path, ..., parent_dir = getwd()) {
#   path <- file.path(parent_dir, stan_path)
#
#   exe_file <- file.path(path, filename)
#   stan_file <- paste0(exe_file, STAN_EXTENSION)
#
#   mod <- cmdstanr::cmdstan_model(stan_file = stan_file, exe_file = exe_file, ...)
#
#   return(mod)
# }

#'
#' Run Stan's MCMC or Optimization Algorithms
#'
.run_stan_algorithms <- function(data, stan_model, method, chains, parallel_chains, refresh, warmup, iter, max_treedepth = NULL) {
  stanargs <- list(data = data, refresh = refresh)

  if (method == "laplace") {
    stanargs <- c(stanargs, init = 0, sig_figs = 8, draws = iter)
  }

  if (method == "optimize") {
    stanargs <- c(stanargs, algorithm = "lbfgs", init = 0, jacobian = TRUE)
  }

  if (method == "sample") {
    stanargs <- c(stanargs, chains = chains, parallel_chains = parallel_chains, iter_warmup = warmup, iter_sampling = iter, max_treedepth = max_treedepth)
  }

  if (method == "variational") {
    stanargs <- c(stanargs, tol_rel_obj = 0.01, iter = 5000, output_samples = 1000)
  }

  return(do.call(stan_model[[method]], args = stanargs))
}
