#!/usr/bin/env Rscript

# File name: generate_fig2.R
#
# Description:
#   Generates publication figure for posterior predictive cost
#   trajectories for a sample patient profile across DTRs (EBRT vs. VBT).

library(dplyr)
library(ggplot2)
library(latex2exp)
library(scales)

source("pkg/R/analysis/utils.R")

type <- "semiparametric"
input_path <- file.path("outputs", "analysis", type, "MCMC", "gcomp_trajectories.RData")

load(input_path)

cb_palette <- c("1= EBRT" = "#0072B2", "2= VBT" = "#E69F00")
time_cutoff <- 3

gcompdat <- bind_rows(
  data.frame(g = "1= EBRT", pred_EBRT),
  data.frame(g = "2= VBT", pred_VBT)
)

time_cuts <- quantile(
  gcompdat[["T_"]],
  probs = seq(from = 0, to = 1, by = 0.05),
  names = FALSE
)

# Convert treatment options to LaTeX code for labels
gcompdat <- gcompdat %>%
  mutate(
    g = factor(
      x = g,
      levels = c("1= EBRT", "2= VBT"),
      labels = c(TeX("EBRT, $d_1$"), TeX("VBT, $d_2$"))
    )
  )

# Compute posterior quantities using time cuts
postdat <- gcompdat %>%
  filter(j > 1) %>%
  mutate(
    idx = findInterval(T_, time_cuts, all.inside = TRUE),
    .after = j
  ) %>%
  summarise(
    lwr_Y = quantile(Y, p = 0.025, na.rm = TRUE),
    upp_Y = quantile(Y, p = 0.975, na.rm = TRUE),
    Y = mean(Y, na.rm = TRUE),
    .by = c(g, idx)
  ) %>%
  arrange(g, idx) %>%
  mutate(
    m = 0,
    lwr_U = cumsum(lwr_Y),
    upp_U = cumsum(upp_Y),
    U = cumsum(Y),
    T_ = time_cuts[idx + 1],
    .by = g,
    .keep = "none"
  )

p1 <- ggplot(gcompdat, mapping = aes(x = T_, y = U, group = m)) +
  geom_step(color = "#999999", alpha = 0.35, linewidth = 0.6) +
  geom_line(
    data = postdat,
    mapping = aes(x = T_, y = U, color = g),
    linewidth = 1.1
  ) +
  geom_ribbon(
    data = postdat,
    mapping = aes(ymin = lwr_U, ymax = upp_U, fill = g),
    alpha = 0.1,
    color = NA
  ) +
  facet_wrap(~g, labeller = "label_parsed") +
  scale_color_manual(
    name = NULL,
    values = c("#0072B2", "#E69F00"),
    labels = c(TeX("EBRT, $d_1$"), TeX("VBT, $d_2$"))
  ) +
  scale_fill_manual(
    name = NULL,
    values = c("#0072B2", "#E69F00"),
    labels = c(TeX("EBRT, $d_1$"), TeX("VBT, $d_2$"))
  ) +
  scale_y_continuous(labels = dollar) +
  labs(
    x = "Time Post-Hysterectomy in Years",
    y = "Accrued Cost in Hundred Thousand USD"
  ) +
  theme_bw(base_size = 12) +
  theme(
    aspect.ratio = 1,
    axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
    legend.position = "none"
  )

results_dir <- "results"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

saveplot(
  filename = "fig2.pdf",
  plot = p1,
  path = results_dir,
  width = 8,
  height = 6,
  dpi = 600
)
