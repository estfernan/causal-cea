#!/usr/bin/env Rscript

# """
# File name: generate_fig1.R
#
# Description:
#   Generates publication figures for Net Monetary Benefit (NMB),
#   Restricted Mean Survival Time (RMST), and Mean Cost (MC) posterior
#   distributions from g-computation output.
# ""

library(dplyr)
library(ggplot2)
library(ggpubr)
library(latex2exp)
library(purrr)
library(scales)
library(tidyr)

source("pkg/R/analysis/constants.R")
source("pkg/R/analysis/utils.R")

type <- "semiparametric"
input_path <- file.path("outputs", "analysis", type, "MCMC", "gcomp_out.RData")

load(input_path)

cb_palette <- c("1= EBRT" = "#0072B2", "2= VBT" = "#E69F00")
kappa_values <- seq(from = 0, to = 5, by = 0.5)

# Combine posterior draws for both EBRT and VBT
gcompdat <- bind_rows(
  data.frame(
    TREATMENT = "1= EBRT",
    DRAW = seq_len(nrow(gcomp_EBRT_all)),
    gcomp_EBRT_all
  ),
  data.frame(
    TREATMENT = "2= VBT",
    DRAW = seq_len(nrow(gcomp_VBT_all)),
    gcomp_VBT_all
  )
)

NMB <- gcompdat %>%
  crossing(KAPPA = kappa_values) %>%
  mutate(MV = KAPPA * MST - MC) %>%
  pivot_wider(
    id_cols = c(KAPPA, DRAW),
    names_from = TREATMENT,
    values_from = MV
  ) %>%
  mutate(NMB = `1= EBRT` - `2= VBT`) %>%
  summarise(
    post_mean = mean(NMB, na.rm = TRUE),
    post_lwr = quantile(NMB, p = 0.025, na.rm = TRUE),
    post_upp = quantile(NMB, p = 0.975, na.rm = TRUE),
    .by = KAPPA
  )

p1 <- ggplot(
  NMB,
  mapping = aes(x = KAPPA, y = post_mean, ymin = post_lwr, ymax = post_upp)
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.7,
    color = "gray50"
  ) +
  geom_point(color = "#0072B2", size = 2.5) +
  geom_errorbar(color = "#0072B2", linewidth = 0.9, width = 0.1) +
  scale_x_continuous(labels = dollar) +
  scale_y_continuous(labels = dollar) +
  labs(
    x = TeX("WTP in Hundred Thousand USD per Year, $\\kappa$"),
    y = TeX("NMB in Hundred Thousand USD per Year, $\\Psi_{d_1, d_2, \\kappa}$")
  ) +
  theme_bw(base_size = 12) +
  theme(
    aspect.ratio = 1,
    axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0))
  )

postdat <- gcompdat %>%
  summarise(MST = mean(MST), MC = mean(MC), .by = TREATMENT)

p2 <- ggplot(
  gcompdat,
  mapping = aes(x = MST, color = TREATMENT, fill = TREATMENT)
) +
  geom_density(alpha = 0.15, linewidth = 1) +
  geom_vline(
    data = postdat,
    mapping = aes(xintercept = MST, color = TREATMENT),
    linetype = "dashed",
    linewidth = 0.9
  ) +
  scale_color_manual(
    name = NULL,
    values = cb_palette,
    labels = c(TeX("EBRT, $d_1$"), TeX("VBT, $d_2$"))
  ) +
  scale_fill_manual(
    name = NULL,
    values = cb_palette,
    labels = c(TeX("EBRT, $d_1$"), TeX("VBT, $d_2$"))
  ) +
  labs(x = "RMST in Years", y = "Posterior Density") +
  theme_bw(base_size = 12) +
  theme(
    aspect.ratio = 1,
    axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
    legend.position = "top"
  )

p3 <- ggplot(
  gcompdat,
  mapping = aes(x = MC, color = TREATMENT, fill = TREATMENT)
) +
  geom_density(alpha = 0.15, linewidth = 1) +
  geom_vline(
    data = postdat,
    mapping = aes(xintercept = MC, color = TREATMENT),
    linetype = "dashed",
    linewidth = 0.9
  ) +
  scale_color_manual(
    name = NULL,
    values = cb_palette,
    labels = c(TeX("EBRT, $d_1$"), TeX("VBT, $d_2$"))
  ) +
  scale_fill_manual(
    name = NULL,
    values = cb_palette,
    labels = c(TeX("EBRT, $d_1$"), TeX("VBT, $d_2$"))
  ) +
  labs(x = "MC in Hundred Thousand USD", y = "Posterior Density") +
  theme_bw(base_size = 12) +
  theme(
    aspect.ratio = 1,
    axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
    legend.position = "top"
  )

results_dir <- "results"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

saveplot(
  filename = "fig1.pdf",
  plot = ggarrange(p1, p2, p3, nrow = 1, align = "hv"),
  path = results_dir,
  width = 14,
  height = 5,
  dpi = 600
)
