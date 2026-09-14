# """
# File name: distributions.R
#
# Description:
#   Functions for implementing probability distributions.
# """

rpiecewiseph <- function(X, time_cuts, beta, lambda) {
  # """
  # Sample survival times from a piecewise constant hazard model.
  # """
  reformat <- is.vector(X)

  # ensure matrix form
  if (reformat) {
    X <- matrix(X, nrow = 1L)
  } else {
    X <- as.matrix(X)
  }

  n <- nrow(X)
  e_lp <- exp(drop(X %*% beta))

  h <- outer(e_lp, lambda, FUN = "*")

  times <- vapply(
    X = seq_len(n),
    FUN = .rpch1,
    FUN.VALUE = numeric(1L),
    h = h,
    time_cuts = time_cuts
  )

  return(times)
}

# """
# Private functions
# """

.rpch1 <- function(i, h, time_cuts) {
  # """
  # Internal piecewise constant hazard distribution for a single matrix row.
  # """
  return(mets::rpch(1L, lambda = h[i, ], breaks = time_cuts))
}
