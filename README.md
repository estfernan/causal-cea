# Bayesian Inference for the Effects of Dynamic Treatments on Joint Longitudinal-Survival Outcomes: Applications in Cost-Effectiveness Analyses

This repository contains the application and simulation codebase, data-generating processes, model estimation scripts, and pipelines for the manuscript.

## Directory Structure

- `/bash`: Reference SLURM batch scripts for HPC execution.
- `/config`: YAML settings, random seed files, and DGP parameter configurations.
- `/logs`: SLURM log destinations.
- `/outputs`: Intermediate MCMC sampling outputs and replication results.
- `/pkg`: Core package containing R functions (`/pkg/R`) and Stan source models (`/pkg/stan`).
- `/results`: Final exported summary tables and publication figures.
- `/scripts`: Execution pipelines categorized into modular subdirectories:
  - `/scripts/application`: SEER-Medicare empirical application pipeline.
  - `/scripts/reporting`: Scripts for generating paper figures and summary tables.
  - `/scripts/simulations`: Data-generating processes and Monte Carlo simulation runners.

## Step-by-Step Execution Guide for Simulations

Due to the computational requirements of fitting the models used in the manuscript across 1,000 Monte Carlo replications, the pipeline is designed to be modular rather than executed via a single wrapper script. This allows users to inspect, test, or re-run individual components (e.g., debugging of a single replication vs. scaling via HPC cluster arrays).

### 1. Generate True Target Estimands

Compute the true population target estimands prior to evaluating performance metrics:

``` r
Rscript scripts/generate_true_estimands.R
```

### 2. Running Simulations

-   **Option A (HPC Cluster):** Use the provided templates in the `/bash` folder to submit array jobs via SLURM (e.g., `sbatch bash/run_bayes_semiparametric.sh`).
-   **Option B (Single Replication):** You can execute individual replication scripts locally by specifying the iteration index, sample size, and censoring rate:

``` bash
Rscript scripts/run_bayes_semiparametric.R --n 1000 --rate 0.10 --iter 1
```

### 3. Summarize Results & Generate Outputs

Once replication `.RData` files populate the `/outputs` directory, compile Table 1 directly into the `/results` folder by executing:

``` r
Rscript scripts/summarize_simulation_results.R
```

## Step-by-Step Execution Guide for SEER-Medicare Application

### 1. Model Fitting & G-Computation

Fit the confounder, survival, treatment initiation, and cost outcome models, followed by g-computation simulations:

``` bash
Rscript scripts/application/run_seer_application.R
```

Model specification (`parametric` or `semiparametric`) can be toggled inside `scripts/application/run_seer_application.R`. MCMC draws and trajectory arrays will be written to `outputs/analysis/[type]/MCMC/`.

### 2. Generate Application Figures and Tables

Once model fitting completes, execute the reporting scripts to build the final publication figures and tables in `/results`:

``` bash
# Generate posterior density plots (Figure 1)
Rscript scripts/reporting/generate_fig1.R

# Generate posterior predictive cost trajectories (Figure 2)
Rscript scripts/reporting/generate_fig2.R

# Generate qpplication summary table (Table 2)
Rscript scripts/reporting/generate_table2.R
```

## Software Dependencies & Session Info

The codebase was developed and tested under **R version 4.5.2**. Core dependencies include `CmdStanR`, `optparse`, and `tidyverse`. Below is the execution environment session information (`sessionInfo()`):

``` r
# R version 4.5.2 (2025-10-31)
# Platform: aarch64-apple-darwin20
# Running under: macOS Tahoe 26.5.1
# 
# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# time zone: America/New_York
# tzcode source: internal
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] scales_1.4.0   readr_2.1.6    optparse_1.7.5 ggplot2_4.0.0  dplyr_1.1.4   
# [6] cmdstanr_0.9.0
# 
# loaded via a namespace (and not attached):
#  [1] gtable_0.3.6         jsonlite_2.0.0       compiler_4.5.2      
#  [4] tidyselect_1.2.1     R6_2.6.1             generics_0.1.4      
#  [7] distributional_0.5.0 knitr_1.50           backports_1.5.0     
# [10] checkmate_2.3.3      tibble_3.3.0         tzdb_0.5.0          
# [13] pillar_1.11.1        RColorBrewer_1.1-3   posterior_1.6.1     
# [16] rlang_1.1.6          getopt_1.20.4        xfun_0.54           
# [19] S7_0.2.0             cli_3.6.5            withr_3.0.2         
# [22] magrittr_2.0.4       ps_1.9.1             grid_4.5.2          
# [25] processx_3.8.6       rstudioapi_0.17.1    hms_1.1.4           
# [28] lifecycle_1.0.4      vctrs_0.6.5          evaluate_1.0.5      
# [31] glue_1.8.0           tensorA_0.36.2.1     farver_2.1.2        
# [34] sessioninfo_1.2.3    abind_1.4-8          tools_4.5.2         
# [37] pkgconfig_2.0.3
```
