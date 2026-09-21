#!/usr/bin/env Rscript

# """
# File name: run_mle_parametric.R
#
# Description:
#   Runs one independent replication of the simulation study using
#   correctly-specified parametric models fitted using maximum likelihood
#   estimation.
# """

source("init.R")

opt <- parse_args(OptionParser(option_list = option_list))

set.seed(seeds[opt$iter], kind = "L'Ecuyer-CMRG")

# script metadata
opt$framework <- "mle"
opt$models <- "parametric"

obsdat <- simulate_claim_data(
  n = opt$n,
  n_claims = config$n_claims,
  dgp_list = dgp_list,
  censoring_rate = opt$rate
)

toc <- Sys.time()

draws <- gcomp_parametric(
  data = obsdat,
  target_estimand = MV,
  a_opt = c(0, 1),
  time_horizon = config$time_horizon,
  n_boot = specs$n_boot,
  n_samples = specs$n_samples,
  n_claims = config$n_claims,
  transform = NMB_gcomp,
  refresh = 100,
  init = 0,
  algorithm = "lbfgs"
)

metrics <- bootstrap_finite_sample_properties(
  draws, truth,
  framework = opt$framework,
  models = opt$models,
  n_patients = opt$n,
  censoring_rate = opt$rate,
  rep = opt$iter,
  wald_ci = FALSE,
  start_time = toc
)

outpath <- format_simulation_output_path(config, opt)
save(draws, metrics, file = outpath)
