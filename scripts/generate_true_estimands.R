#!/usr/bin/env Rscript

# """
# File name: generate_true_estimands.R
#
# Description:
#   Compute the "true" target estimand (net monetary benefit) for the simulation
#   study under two treatment strategies: always treated and never treated.
#   The truth is estimated by simulating a very large synthetic dataset (n_truth)
#   under each strategy.
# """

source("init.R")

set.seed(106, kind = "L'Ecuyer-CMRG")

opt <- yaml::read_yaml("config/simulations.yaml")
seeds <- read_lines("config/seeds/random_seeds.txt", transform = as.integer)

# set individual configuration settings
config <- opt$settings
specs <- opt$specs
dgp_list <- opt$DGP

outpath <- file.path(config$outdir, "target_estimands.RData")

MV1 <- simulate_claim_data(
  n = config$n_truth,
  n_claims = config$n_claims,
  dgp_list = dgp_list,
  target_estimand = MV,
  treatment_decision = 1L,
  time_horizon = config$time_horizon
)

MV0 <- simulate_claim_data(
  n = config$n_truth,
  n_claims = config$n_claims,
  dgp_list = dgp_list,
  target_estimand = MV,
  treatment_decision = 0L,
  time_horizon = config$time_horizon
)

truth <- MV1 - MV0

save(seeds, config, specs, dgp_list, MV1, MV0, truth, file = outpath)
