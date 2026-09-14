#!/bin/bash

#SBATCH --job-name=run_mle_parametric_n1000_r010
#SBATCH --ntasks=1
#SBATCH --mem=3G
#SBATCH --time=05:00:00
#SBATCH --array=1-1000%100
#SBATCH --output=logs/mle/parametric/n1000/r010/%A_%a.out
#SBATCH --error=logs/mle/parametric/n1000/r010/%A_%a.err

module load r/4.4.0-yycctsj

# script arguments
N=1000
RATE="0.10"

Rscript --vanilla scripts/run_mle_parametric.R \
  --n $N \
  --rate $RATE \
  --iter $SLURM_ARRAY_TASK_ID
