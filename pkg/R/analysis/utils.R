#
# Utility functions used in project scripts.
#

#'
#' Save and crop a ggplot
#'
saveplot <- function(filename, plot = ggplot2::last_plot(),
                     device = NULL,
                     path = NULL,
                     scale = 1, width = NA, height = NA,
                     units = c("in", "cm", "mm", "px"),
                     dpi = 300,
                     create.dir = FALSE, quiet = TRUE, ...) {
  figpath <- file.path(path, filename)

  ggplot2::ggsave(
    figpath,
    plot = plot,
    device = device,
    path = NULL,
    scale = scale, width = width, height = height,
    units = units,
    dpi = dpi,
    create.dir = create.dir,
    ...
  )

  return(knitr::plot_crop(figpath, quiet = quiet))
}

# """
# Base
# """

#'
#' Number of Unique Elements
#'
length_unique <- function(x, ...) {
  return(length(unique(x, ...)))
}

# """
# Environment Variables
# """

#'
#' Get Environment Variables
#'
Sys_getenv <- function(x = NULL, unset = "", names = NA, to = "string", ...) {
  x <- Sys.getenv(x, unset = unset, names = names)

  if (to == "integer") {
    x <- strtoi(x, ...)
  }

  return(x)
}

# """
# Files
# """

#'
#' Use C-style String Formatting for Simulation Files
#'
simfile <- function(hazard, intensity, method, rep, lambda_C, n_obs, outdir = "") {
  lambda_C <- 1000 * round(lambda_C, digits = 3)
  fmt <- "%s-%s-%s-%03d-%04d-%04d.RData"
  outfile <- sprintf(fmt, hazard, intensity, method, rep, lambda_C, n_obs)

  if (outdir != EMPTY_STRING) {
    outfile <- file.path(outdir, outfile)
  }

  return(outfile)
}

# """
# Operations
# """

#'
#' Finite-sample properties of an estimator
#'
estimator_properties <- function(draws, truth, alpha = 0.05, ...) {
  upp_prob <- 1 - alpha / 2
  lwr_prob <- 1 - upp_prob

  other_args <- list(...)
  n_args <- length(other_args)

  post_mean <- mean(draws, na.rm = TRUE)
  post_se <- sd(draws, na.rm = TRUE)
  lwr_CI <- quantile(draws, p = lwr_prob, na.rm = TRUE, names = FALSE)
  upp_CI <- quantile(draws, p = upp_prob, na.rm = TRUE, names = FALSE)
  bias <- post_mean - truth
  pct_bias <- bias / abs(truth)
  coverage <- lwr_CI < truth & truth < upp_CI
  width <- upp_CI - lwr_CI

  props <- c(
    truth = truth,
    mean = post_mean, SE = post_se,
    lwr_CI = lwr_CI, upp_CI = upp_CI,
    bias = bias,
    pct_bias = pct_bias,
    RMSE = sqrt(bias^2 + post_se^2),
    coverage = coverage, CI_width = width
  )

  if (n_args > 0) {
    tr_props <- t(props)
    props <- data.frame(other_args, tr_props)
  }

  return(props)
}

#'
#' Softmax function
#'
softmax <- function(x, na_rm = FALSE) {
  is_matrix <- !is.null(dim(x))

  if (is_matrix) {
    bounds <- apply(x, 1, max)
    df <- exp(sweep(x, 1, bounds))
    tots <- rowSums(df, na.rm = na_rm)
    return(df / tots)
  }

  bound <- max(x, na.rm = na_rm)
  df <- exp(x - bound)
  tot <- sum(df, na.rm = na_rm)
  return(df / tot)
}

#'
#' Quantile-based binning to create time cuts
#'
quantile_binning <- function(times, n_cuts)
{
  probs <- seq(from = 0, to = 1, length.out = n_cuts + 1)
  time_cuts <- quantile(times, probs = probs, names = FALSE)
  time_cuts[1] <- DEFAULT_START_TIME
  return(time_cuts)
}

# """
# Stan
# """

#'
#' Compile Stan Programs
#'
compile_stanmodels <- function(stan_path = STAN_PATH, force_recompile = TRUE, ...) {
  files <- list.files(
    path = stan_path,
    pattern = STAN_EXTENSION,
    full.names = TRUE, recursive = TRUE
  )

  stanmodels <- sapply(
    X = files,
    FUN = cmdstanr::cmdstan_model,
    force_recompile = force_recompile,
    ...
  )

  return(invisible(stanmodels))
}
