# """
# File name: internal.R
#
# Description:
#   Functions used to facilitate operations within others.
# """

concat <- function(n, ...) {
  # """
  # Concatenate input vectors or matrices into a single structure.
  # """
  func <- ifelse(n > 1, cbind, c)
  return(func(...))
}
