#!/bin/bash
#SBATCH -J manhattan 
#SBATCH --time=0:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --array=0-599
#SBATCH -o logs/manhattan/%A_%a.out

set -euo pipefail

source ./utils.sh
read EMBEDDING AGE < <(inverse_map ${SLURM_ARRAY_TASK_ID})

Rscript plots.R $EMBEDDING $AGE 
