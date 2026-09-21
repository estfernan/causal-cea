#!/usr/bin/env Rscript

# """
# File name: run_bayesian_semiparametric.R
#
# Description:
#   Runs one independent replication of the simulation study using
#   Bayesian semiparametric models.
# """

library(survival)

source("init.R")

opt <- parse_args(OptionParser(option_list = option_list))

set.seed(seeds[opt$iter], kind = "L'Ecuyer-CMRG")

# script metadata
opt$framework <- "bayes"
opt$models <- "semiparametric"

obsdat <- simulate_claim_data(
  n = opt$n,
  n_claims = config$n_claims,
  dgp_list = dgp_list,
  censoring_rate = opt$rate
)

time_cuts <- create_time_cuts(obsdat, config)

toc <- Sys.time()

stanout <- fit_semiparametric_models(
  data = obsdat,
  time_cuts = time_cuts,
  method = "sample",
  seed = seeds[opt$iter],
  refresh = 100,
  chains = 1,
  iter_warmup = specs$n_warmup,
  iter_sampling = specs$n_iter
)

draws <- bayes_gcomp_semiparametric(
  stanout = stanout,
  time_cuts = time_cuts,
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
