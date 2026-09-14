# """
# File name: math.R
#
# Description:
#   Functions to compute mathematical and statistical quantities.
# """

MV <- function(U, T_, kappa = 1) {
  # """
  # Compute the monetary value for a given cost and survival time.
  # """
  return(kappa * T_ - U)
}
