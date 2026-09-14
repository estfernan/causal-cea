#!/usr/bin/env Rscript

# """
# File name: init.R
#
# Description:
#   Loads libraries, constant variables and all functions in the source directory.
# """

library(cmdstanr)
library(dplyr)
library(ggplot2)
library(optparse)
library(readr)
library(scales)

source_all <- function(dir, recursive = TRUE) {
  # """
  # Read all R code in a folder.
  # """
  files <- list.files(dir, pattern = "R", full.names = TRUE, recursive = recursive)
  sapply(files, source)
  return(invisible(files))
}

# setup flags and variables
target_file <- "outputs/target_estimands.RData"

# load all functions in source directory
source_all("pkg/R", recursive = TRUE)

# load important simulation objects
if (file.exists(target_file)) {
  load(target_file)
}

# store default command line options for the simulation scripts
option_list <- list(
  make_option("--n", type = "integer", default = 1000),
  make_option("--rate", type = "double", default = 0.1),
  make_option("--iter", type = "integer", default = 1)
)
