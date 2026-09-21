#!/usr/bin/env Rscript

# File name: generate_synthetic_obsdat.R
#
# Description:
#   Generates synthetic independent test data to test run_seer_application.R
#   without requiring complex longitudinal DGPs.

set.seed(106)

n_patients <- 100
n_intervals <- 3
n_obs <- n_patients * n_intervals

# Exact factor levels from the original SEER-Medicare dataset
race_levels <- c("1", "2", "3")
income_levels <- c("1", "2", "3", "4")
trt_levels <- c("1= None", "2= EBRT", "3= VBT", "4= Both")

ids <- as.integer(rep(seq_len(n_patients), each = n_intervals))
indices <- as.integer(rep(seq_len(n_intervals), times = n_patients))

obsdat <- data.frame(
  ID = ids,
  PATIENT_ID = ids,
  INDEX = indices,
  TIME = as.numeric(rep(0, n_obs)),
  GAPTIME = round(rexp(n_obs, rate = 1 / 12) + 0.1, 2),
  STATUS = as.numeric(rbinom(n_obs, 1, 0.25)),
  DEATH = as.numeric(rbinom(n_obs, 1, 0.15)),
  UPCOMING_COST = round(rgamma(n_obs, shape = 2, scale = 1.5), 2),
  BALANCE = as.numeric(rep(0, n_obs)),
  TREATMENT_READY = as.numeric(rbinom(n_obs, 1, 0.5)),
  ADJUVANT_TREATMENT = factor(sample(trt_levels, n_obs, replace = TRUE), levels = trt_levels),
  AGE = round(rnorm(n_obs, mean = 68, sd = 6)),
  RACE_WBO = factor(sample(race_levels, n_obs, replace = TRUE), levels = race_levels),
  MARRIED = as.numeric(rbinom(n_obs, 1, 0.6)),
  MEDIAN_INCOME = factor(sample(income_levels, n_obs, replace = TRUE), levels = income_levels),
  METROPOLITAN = as.numeric(rbinom(n_obs, 1, 0.75)),
  NODES_EXAMINED = as.numeric(rbinom(n_obs, 1, 0.5)),
  FIGO_STAGE_IB = as.numeric(rbinom(n_obs, 1, 0.4)),
  HIGH_COMORBIDITY = as.numeric(rbinom(n_obs, 1, 0.3)),
  ANY_OTHER_TX = as.numeric(rbinom(n_obs, 1, 0.2)),
  LAG_HIGH_COMORBIDITY = as.integer(rbinom(n_obs, 1, 0.3)),
  LAG_ANY_OTHER_TX = as.integer(rbinom(n_obs, 1, 0.2)),
  LAG_TREATMENT_READY = as.numeric(rbinom(n_obs, 1, 0.5)),
  LAG_ADJUVANT_TREATMENT = factor(sample(trt_levels, n_obs, replace = TRUE), levels = trt_levels),
  LAG_GAPTIME = round(rexp(n_obs, rate = 1 / 12) + 0.1, 2),
  LAG_UPCOMING_COST = round(rgamma(n_obs, shape = 2, scale = 1.5), 2)
)

# Logical consistency check
obsdat$TREATMENT_INITIATED <- ifelse(obsdat$ADJUVANT_TREATMENT == "1= None", FALSE, TRUE)

# Verify representation of all 4 levels in baseline and follow-up
stopifnot(
  all(trt_levels %in% unique(obsdat$ADJUVANT_TREATMENT[obsdat$INDEX == 1])),
  all(trt_levels %in% unique(obsdat$ADJUVANT_TREATMENT[obsdat$INDEX > 1]))
)

dir.create("data", recursive = TRUE, showWarnings = FALSE)
save(obsdat, file = "data/obsdat.RData")
