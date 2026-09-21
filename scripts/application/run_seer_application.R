#!/usr/bin/env Rscript

# """
# File name: run_seer_application.R
#
# Description:
#   Execution script for the SEER-Medicare application analysis.
#   Fits confounder, survival, treatment, and cost models under parametric or
#   semiparametric specifications and executes g-computation.
# """

library(doParallel)
library(doRNG)
library(modeest)
library(survival)

# ------------------------------------------------------------------------------
# 1. Configuration and Paths
# ------------------------------------------------------------------------------

type <- "parametric" # Options: "semiparametric", "parametric"
n_cuts <- 10
chains <- 4
warmup <- 1000
iter <- 1000
refresh <- 100

set.seed(106, kind = "L'Ecuyer-CMRG")

output_dir <- file.path("outputs", "analysis", type, "MCMC")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 2. Source Dependencies
# ------------------------------------------------------------------------------

source("pkg/R/analysis/constants.R")
source("pkg/R/analysis/distributions.R")
source("pkg/R/analysis/DTRs.R")
source("pkg/R/analysis/g-computation.R")
source("pkg/R/analysis/models.R")
source("pkg/R/analysis/predict.R")
source("pkg/R/analysis/private.R")
source("pkg/R/analysis/utils.R")

# ------------------------------------------------------------------------------
# 3. Data Loading
# ------------------------------------------------------------------------------

load("data/obsdat.RData")

init_obs <- subset(obsdat, INDEX == 1)
fu_obs <- subset(obsdat, INDEX > 1)

# ------------------------------------------------------------------------------
# 4. Baseline and Time-Varying Confounder Models
# ------------------------------------------------------------------------------

cat("Fitting confounder models...\n")

f_comorbidity <- HIGH_COMORBIDITY ~
  AGE + RACE_WBO + MARRIED + MEDIAN_INCOME + METROPOLITAN + NODES_EXAMINED +
  FIGO_STAGE_IB + LAG_HIGH_COMORBIDITY + LAG_ANY_OTHER_TX +
  LAG_TREATMENT_READY + LAG_ADJUVANT_TREATMENT + LAG_GAPTIME +
  LAG_UPCOMING_COST

f_other <- ANY_OTHER_TX ~
  AGE + RACE_WBO + MARRIED + MEDIAN_INCOME + METROPOLITAN + NODES_EXAMINED +
  FIGO_STAGE_IB + LAG_HIGH_COMORBIDITY + LAG_ANY_OTHER_TX +
  LAG_TREATMENT_READY + LAG_ADJUVANT_TREATMENT + LAG_GAPTIME +
  LAG_UPCOMING_COST

covars_dict <- list(
  AGE = list(name = "AGE", spec = "baseline", distr = "normal"),
  RACE_WBO = list(name = "RACE_WBO", spec = "baseline", distr = "categorical"),
  MARRIED = list(name = "MARRIED", spec = "baseline", distr = "bernoulli"),
  MEDIAN_INCOME = list(name = "MEDIAN_INCOME", spec = "baseline", distr = "categorical"),
  METROPOLITAN = list(name = "METROPOLITAN", spec = "baseline", distr = "bernoulli"),
  NODES_EXAMINED = list(name = "NODES_EXAMINED", spec = "baseline", distr = "bernoulli"),
  FIGO_STAGE_IB = list(name = "FIGO_STAGE_IB", spec = "baseline", distr = "bernoulli"),
  HIGH_COMORBIDITY = list(name = "HIGH_COMORBIDITY", spec = "time-varying", distr = "bernoulli", f = f_comorbidity),
  ANY_OTHER_TX = list(name = "ANY_OTHER_TX", spec = "time-varying", distr = "bernoulli", f = f_other)
)

for (j in seq_along(covars_dict)) {
  covar <- covars_dict[[j]][["name"]]
  distr <- covars_dict[[j]][["distr"]]
  f <- covars_dict[[j]][["f"]]
  spec <- covars_dict[[j]][["spec"]]

  if (spec == "baseline") {
    f_init <- paste(covar, "~", 1)

    simplereg <- switch(distr,
      "normal" = bayesian_gaussian_simple,
      "bernoulli" = bayesian_bernoulli_simple,
      "categorical" = bayesian_categorical_simple
    )

    theta_L0 <- simplereg(
      f_init,
      data = init_obs,
      chains = chains,
      parallel_chains = chains,
      iter = iter,
      warmup = warmup,
      refresh = refresh,
      format = "draws_array"
    )
  } else if (spec == "time-varying") {
    f_init <- .remove_lagged_vars(f)

    glmreg <- switch(distr,
      "normal" = bayesian_gaussian_glm,
      "bernoulli" = bayesian_bernoulli_glm,
      "categorical" = bayesian_categorical_glm
    )

    theta_L0 <- glmreg(
      f_init,
      data = init_obs,
      chains = chains,
      parallel_chains = chains,
      iter = iter,
      warmup = warmup,
      refresh = refresh,
      format = "draws_array"
    )
  }

  covars_dict[[j]][["f_init"]] <- f_init
  covars_dict[[j]][["init_theta"]] <- theta_L0

  if (spec == "time-varying") {
    theta_L <- glmreg(
      f,
      data = fu_obs,
      chains = chains,
      parallel_chains = chains,
      iter = iter,
      warmup = warmup,
      refresh = refresh,
      format = "draws_array"
    )

    covars_dict[[j]][["theta"]] <- theta_L
  }
}

save(covars_dict, file = file.path(output_dir, "covars_dict.RData"))

# ------------------------------------------------------------------------------
# 5. Outcome & Gaptime Models
# ------------------------------------------------------------------------------

cat("Fitting survival, treatment initiation, and cost outcome models...\n")

f_init_death <- Surv(GAPTIME, DEATH) ~
  AGE + RACE_WBO + MARRIED + MEDIAN_INCOME + METROPOLITAN + NODES_EXAMINED +
  FIGO_STAGE_IB + HIGH_COMORBIDITY + ANY_OTHER_TX + TREATMENT_READY +
  ADJUVANT_TREATMENT - 1

f_death <- Surv(GAPTIME, DEATH) ~
  AGE + RACE_WBO + MARRIED + MEDIAN_INCOME + METROPOLITAN + NODES_EXAMINED +
  FIGO_STAGE_IB + HIGH_COMORBIDITY + ANY_OTHER_TX + TREATMENT_READY +
  ADJUVANT_TREATMENT + LAG_HIGH_COMORBIDITY + LAG_ANY_OTHER_TX +
  LAG_TREATMENT_READY + LAG_ADJUVANT_TREATMENT + LAG_GAPTIME +
  LAG_UPCOMING_COST - 1

f_init_event <- Surv(GAPTIME, STATUS) ~
  AGE + RACE_WBO + MARRIED + MEDIAN_INCOME + METROPOLITAN + NODES_EXAMINED +
  FIGO_STAGE_IB + HIGH_COMORBIDITY + ANY_OTHER_TX + TREATMENT_READY +
  ADJUVANT_TREATMENT - 1

f_event <- Surv(GAPTIME, STATUS) ~
  AGE + RACE_WBO + MARRIED + MEDIAN_INCOME + METROPOLITAN + NODES_EXAMINED +
  FIGO_STAGE_IB + HIGH_COMORBIDITY + ANY_OTHER_TX + TREATMENT_READY +
  ADJUVANT_TREATMENT + LAG_HIGH_COMORBIDITY + LAG_ANY_OTHER_TX +
  LAG_TREATMENT_READY + LAG_ADJUVANT_TREATMENT + LAG_GAPTIME +
  LAG_UPCOMING_COST - 1

f_init_cost <- UPCOMING_COST ~
  AGE + RACE_WBO + MARRIED + MEDIAN_INCOME + METROPOLITAN + NODES_EXAMINED +
  FIGO_STAGE_IB + HIGH_COMORBIDITY + ANY_OTHER_TX + TREATMENT_READY +
  ADJUVANT_TREATMENT + DEATH - 1

f_cost <- UPCOMING_COST ~
  AGE + RACE_WBO + MARRIED + MEDIAN_INCOME + METROPOLITAN + NODES_EXAMINED +
  FIGO_STAGE_IB + HIGH_COMORBIDITY + ANY_OTHER_TX + TREATMENT_READY +
  ADJUVANT_TREATMENT + DEATH + LAG_HIGH_COMORBIDITY + LAG_ANY_OTHER_TX +
  LAG_TREATMENT_READY + LAG_ADJUVANT_TREATMENT + LAG_GAPTIME +
  LAG_UPCOMING_COST - 1

init_weights <- init_obs[["GAPTIME"]]
fu_weights <- fu_obs[["GAPTIME"]]

if (type == "semiparametric") {
  init_time_cuts <- quantile_binning(init_weights, n_cuts)
  fu_time_cuts <- quantile_binning(fu_weights, n_cuts)

  theta_D0 <- bayesian_piecewise_exponential(
    formula = f_init_death, data = init_obs, chains = chains,
    parallel_chains = chains, iter = iter, warmup = warmup,
    refresh = refresh, max_treedepth = 12, time_cuts = init_time_cuts,
    format = "draws_array"
  )

  theta_D <- bayesian_piecewise_exponential(
    formula = f_death, data = fu_obs, chains = chains,
    parallel_chains = chains, iter = iter, warmup = warmup,
    refresh = refresh, time_cuts = fu_time_cuts, format = "draws_array"
  )

  theta_E0 <- bayesian_piecewise_exponential(
    formula = f_init_event, data = init_obs, chains = chains,
    parallel_chains = chains, iter = iter, warmup = warmup,
    refresh = refresh, time_cuts = init_time_cuts, format = "draws_array"
  )

  theta_E <- bayesian_piecewise_exponential(
    formula = f_event, data = fu_obs, chains = chains,
    parallel_chains = chains, iter = iter, warmup = warmup,
    refresh = refresh, time_cuts = fu_time_cuts, format = "draws_array"
  )

  theta_Y0 <- bayesian_gamma_piecewise_pm(
    formula = f_init_cost, data = init_obs, time_weights = init_weights,
    chains = chains, parallel_chains = chains, iter = iter, warmup = warmup,
    refresh = refresh, time_cuts = init_time_cuts, format = "draws_array"
  )

  theta_Y <- bayesian_gamma_piecewise_pm(
    formula = f_cost, data = fu_obs, time_weights = fu_weights,
    chains = chains, parallel_chains = chains, iter = iter, warmup = warmup,
    refresh = refresh, time_cuts = fu_time_cuts, format = "draws_array"
  )
} else if (type == "parametric") {
  init_time_cuts <- NULL
  fu_time_cuts <- NULL

  theta_D0 <- bayesian_weibull_ph(
    formula = f_init_death, data = init_obs, chains = chains,
    parallel_chains = chains, iter = iter, warmup = warmup,
    refresh = refresh, format = "draws_array"
  )

  theta_D <- bayesian_weibull_ph(
    formula = f_death, data = fu_obs, chains = chains,
    parallel_chains = chains, iter = iter, warmup = warmup,
    refresh = refresh, format = "draws_array"
  )

  theta_E0 <- bayesian_weibull_ph(
    formula = f_init_event, data = init_obs, chains = chains,
    parallel_chains = chains, iter = iter, warmup = warmup,
    refresh = refresh, format = "draws_array"
  )

  theta_E <- bayesian_weibull_ph(
    formula = f_event, data = fu_obs, chains = chains,
    parallel_chains = chains, iter = iter, warmup = warmup,
    refresh = refresh, format = "draws_array"
  )

  theta_Y0 <- bayesian_gamma_gompertz_pm(
    formula = f_init_cost, data = init_obs, time_weights = init_weights,
    chains = chains, parallel_chains = chains,
    iter = iter, warmup = warmup,
    refresh = refresh, format = "draws_array"
  )

  theta_Y <- bayesian_gamma_gompertz_pm(
    formula = f_cost, data = fu_obs, time_weights = fu_weights,
    chains = chains, parallel_chains = chains,
    iter = iter, warmup = warmup,
    refresh = refresh, format = "draws_array"
  )
}

save(f_init_death, theta_D0, init_time_cuts, file = file.path(output_dir, "theta_D0.RData"))
save(f_death, theta_D, fu_time_cuts, file = file.path(output_dir, "theta_D.RData"))
save(f_init_event, theta_E0, file = file.path(output_dir, "theta_E0.RData"))
save(f_event, theta_E, file = file.path(output_dir, "theta_E.RData"))
save(f_init_cost, theta_Y0, file = file.path(output_dir, "theta_Y0.RData"))
save(f_cost, theta_Y, file = file.path(output_dir, "theta_Y.RData"))

# Treatment Continuation Model
f_treatment <- TREATMENT_READY ~
  AGE + RACE_WBO + MARRIED + MEDIAN_INCOME + METROPOLITAN + NODES_EXAMINED +
  FIGO_STAGE_IB + HIGH_COMORBIDITY + ANY_OTHER_TX + LAG_HIGH_COMORBIDITY +
  LAG_ANY_OTHER_TX + LAG_TREATMENT_READY + LAG_ADJUVANT_TREATMENT +
  LAG_GAPTIME + LAG_UPCOMING_COST - 1

fu_obs_trt <- setup_discrete_time_to_event_data(
  data = fu_obs,
  id_column = "ID",
  index_column = "INDEX",
  var_column = "TREATMENT_READY"
)

theta_S <- bayesian_continuation_ratio(
  formula = f_treatment,
  data = fu_obs_trt,
  index_weights = fu_obs_trt[["INDEX"]],
  chains = chains,
  parallel_chains = chains,
  iter = iter,
  warmup = warmup,
  refresh = refresh,
  format = "draws_array"
)

save(f_treatment, theta_S, file = file.path(output_dir, "theta_S.RData"))

# ------------------------------------------------------------------------------
# 6. G-Computation Execution and Trajectory Simulation
# ------------------------------------------------------------------------------

cat("Running g-computation and dynamic treatment trajectory simulations...\n")

joint <- list(
  treatment_initiation = list(name = "initiate", f = f_treatment, theta = theta_S),
  event_gaptime = list(name = "event", f_init = f_init_event, f = f_event, init_theta = theta_E0, theta = theta_E),
  death_gaptime = list(name = "death", f_init = f_init_death, f = f_death, init_theta = theta_D0, theta = theta_D),
  outcome = list(name = "outcome", f_init = f_init_cost, f = f_cost, init_theta = theta_Y0, theta = theta_Y)
)

samp_out <- list(
  joint = joint,
  covars_list = covars_dict,
  settings = list(
    hazard = ifelse(type == "semiparametric", "piecewise", "weibull"),
    intensity = ifelse(type == "semiparametric", "piecewise", "gamma-gompertz"),
    init_death = TRUE,
    method = "sample",
    chains = chains,
    iter = iter,
    init_time_cuts = init_time_cuts,
    fu_time_cuts = fu_time_cuts
  )
)

save(samp_out, file = file.path(output_dir, sprintf("samp_out_%s.RData", type)))

time_cutoff <- 3

gcomp_EBRT <- gcomputation_dtr(
  samp_out,
  target_estimand = SM, g = g_EBRT, tau = time_cutoff,
  n_sim = 1000, parallel = TRUE, n_cores = 4, transpose = TRUE
)

gcomp_VBT <- gcomputation_dtr(
  samp_out,
  target_estimand = SM, g = g_VBT, tau = time_cutoff,
  n_sim = 1000, parallel = TRUE, n_cores = 4, transpose = TRUE
)

gcomp_EBRT_all <- do.call("rbind", gcomp_EBRT)
gcomp_VBT_all <- do.call("rbind", gcomp_VBT)

save(time_cutoff, gcomp_EBRT_all, gcomp_VBT_all, file = file.path(output_dir, "gcomp_out.RData"))

basedat <- data.frame(
  ID = 1,
  INDEX = 1,
  AGE = mean(obsdat[["AGE"]], na.rm = TRUE),
  RACE_WBO = mlv(obsdat[["RACE_WBO"]], na.rm = TRUE),
  MARRIED = mlv(obsdat[["MARRIED"]], na.rm = TRUE),
  MEDIAN_INCOME = mlv(obsdat[["MEDIAN_INCOME"]], na.rm = TRUE),
  METROPOLITAN = mlv(obsdat[["METROPOLITAN"]], na.rm = TRUE),
  NODES_EXAMINED = mlv(obsdat[["NODES_EXAMINED"]], na.rm = TRUE),
  FIGO_STAGE_IB = mlv(obsdat[["FIGO_STAGE_IB"]], na.rm = TRUE)
)

pred_EBRT <- gcomputation_dtr_trajectory(
  samp_out, basedat,
  g = g_EBRT, tau = time_cutoff, select_chains = 1, combine = TRUE
)

pred_VBT <- gcomputation_dtr_trajectory(
  samp_out, basedat,
  g = g_VBT, tau = time_cutoff, select_chains = 1, combine = TRUE
)

save(pred_EBRT, pred_VBT, file = file.path(output_dir, "gcomp_trajectories.RData"))

cat("SEER-Medicare application and trajectory pipeline completed successfully!\n")
