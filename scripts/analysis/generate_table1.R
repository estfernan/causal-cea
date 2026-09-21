#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(scales)
library(readr)

source("init.R")

# 1. Load intermediate simulation output files
files <- list.files(
  path = "outputs",
  pattern = "^rep\\d{4}\\.RData$",
  full.names = TRUE,
  recursive = TRUE
)

if (length(files) == 0) {
  stop("No outputs found. Run simulations first.")
}

# tabout <- data.frame()
# for (file in files) {
#   load(file)
#   tabout <- rbind(tabout, metrics)
# }

tabout <- data.frame()
for (file in files) {
  load(file)

  # Extract interval tag from path if it exists, otherwise set as NA or "standard"
  if (grepl("/c\\d{3}/", file)) {
    interval_tag <- regmatches(file, regexpr("c\\d{3}", file))
    metrics$interval_setting <- interval_tag
  } else {
    metrics$interval_setting <- "standard"
  }

  tabout <- rbind(tabout, metrics)
}

# Compute performance metrics for Table 1
infotab <- tabout %>%
  group_by(censoring_rate, framework, models, interval_setting) %>%
  summarise(
    n_reps   = n(),
    avg_est  = mean(mean, na.rm = TRUE),
    pct_bias = mean(pct_bias, na.rm = TRUE),
    SE       = mean(SE, na.rm = TRUE),
    RMSE     = mean(RMSE, na.rm = TRUE),
    CI_width = mean(CI_width, na.rm = TRUE),
    coverage = mean(coverage, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  mutate(
    censoring_level = label_censoring(censoring_rate),
    intervals = label_intervals(interval_setting),
    method_label = label_methods(framework, models)
  )

table1_formatted <- infotab %>%
  transmute(
    `Censoring Rate` = percent(censoring_rate),
    `Method` = method_label,
    `Interval Setting` = intervals,
    `Estimate` = number(avg_est, accuracy = 0.001),
    `% Bias` = percent(pct_bias, accuracy = 0.01),
    `SE` = number(SE, accuracy = 0.001),
    `RMSE` = number(RMSE, accuracy = 0.001),
    `Width` = number(CI_width, accuracy = 0.001),
    `CR (%)` = percent(coverage, accuracy = 0.1)
  ) %>%
  arrange(`Censoring Rate`, Method)

# Export Table 1
write_csv(table1_formatted, "results/table1.csv")
