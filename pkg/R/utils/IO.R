# """
# File name: IO.R
#
# Description:
#   Functions for processing inputs and outputs, such as loading, saving, etc.
# """

read_lines <- function(file, ..., transform = identity) {
  # """
  # Read text lines from a file
  # """
  out <- readLines(file, ...)
  return(transform(out))
}
