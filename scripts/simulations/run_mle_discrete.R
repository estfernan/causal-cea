#!/usr/bin/env Rscript

# """
# File name: run_mle_discrete.R
#
# Description:
#   Runs one independent replication of the simulation study using
#   discrete-time models fitted using maximum likelihood estimation.
# """

source("init.R")

opt <- parse_args(OptionParser(option_list = option_list))

set.seed(seeds[opt$iter], kind = "L'Ecuyer-CMRG")

# script metadata
opt$framework <- "mle"
opt$models <- "discrete"

obsdat <- simulate_claim_data(
  n = opt$n,
  n_claims = config$n_claims,
  dgp_list = dgp_list,
  censoring_rate = opt$rate,
  lag_confounders = FALSE
)

time_cuts <- create_discrete_time_intervals(obsdat, config)
ppdat <- generate_person_period_dataset(obsdat, time_cuts)

toc <- Sys.time()

draws <- gcomp_discrete(
  data = ppdat,
  time_cuts = time_cuts,
  target_estimand = MV,
  a_opt = c(0, 1),
  time_horizon = config$time_horizon,
  n_boot = specs$n_boot,
  n_samples = specs$n_samples,
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
  start_time = toc
)

outpath <- format_simulation_output_path(config, opt)
save(draws, metrics, file = outpath)
