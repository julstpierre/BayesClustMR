#!/bin/bash
#SBATCH --account=def-mlegault
#SBATCH --job-name=mr_clust_single
#SBATCH --time=05:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=../results/logs/%x_%j.out

cd ~/projects/def-mlegault/gf591137/BayesClustMR/simulations/scripts/
module load r/4.5.0

echo "Running $1, prior=$4, chunks $2 to $3"

Rscript sim_estim.R "$1" "$2" "$3" "$4"

# DO NOT UNCOMMENT
# sbatch run_single_job.sh D_S1 351 400 right_diffuse