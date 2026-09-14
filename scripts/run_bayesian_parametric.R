#!/usr/bin/env Rscript

# """
# File name: run_bayesian_parametric.R
#
# Description:
#   Runs one independent replication of the simulation study using
#   correctly-specified Bayesian parametric models.
# """

source("init.R")

opt <- parse_args(OptionParser(option_list = option_list))

set.seed(seeds[opt$iter], kind = "L'Ecuyer-CMRG")

# script metadata
opt$framework <- "bayes"
opt$models <- "parametric"

obsdat <- simulate_claim_data(
  n = opt$n,
  n_claims = config$n_claims,
  dgp_list = dgp_list,
  censoring_rate = opt$rate
)

toc <- Sys.time()

stanout <- fit_parametric_models(
  data = obsdat,
  method = "sample",
  seed = seeds[opt$iter],
  refresh = 100,
  chains = 1,
  iter_warmup = specs$n_warmup,
  iter_sampling = specs$n_iter
)

draws <- bayes_gcomp_parametric(
  stanout = stanout,
  target_estimand = MV,
  a_opt = c(0, 1),
  time_horizon = config$time_horizon,
  n_samples = specs$n_samples,
  n_claims = config$n_claims,
  transform = NMB_gcomp,
  refresh = 100
)

metrics <- bayes_finite_sample_properties(
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
