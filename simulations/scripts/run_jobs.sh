#!/bin/bash
#SBATCH --account=def-mlegault
#SBATCH --job-name=mr_clust_jobs
#SBATCH --time=05:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=../results/logs/%x_%A_%a.out
#SBATCH --array=0-599

cd ~/links/projects/def-mlegault/gf591137/BayesClustMR/simulations/scripts/
module load r/4.5.0

SCENARIOS=("D_S1" "D_S2" "D_S3")
PRIORS=("right_confident" "right_diffuse" "uniform" "wrong_diffuse" "wrong_confident")
N_CHUNKS=40
CHUNK_SIZE=25

N_SCENARIOS=${#SCENARIOS[@]}
N_PRIORS=${#PRIORS[@]}

# Decompose SLURM_ARRAY_TASK_ID into scenario, prior, and chunk indices
SCENARIO_IDX=$(( SLURM_ARRAY_TASK_ID / (N_PRIORS * N_CHUNKS) ))
REMAINDER=$(( SLURM_ARRAY_TASK_ID % (N_PRIORS * N_CHUNKS) ))
PRIOR_IDX=$(( REMAINDER / N_CHUNKS ))
CHUNK_IDX=$(( REMAINDER % N_CHUNKS ))

SCENARIO=${SCENARIOS[$SCENARIO_IDX]}
PRIOR=${PRIORS[$PRIOR_IDX]}
START_CHUNK=$(( CHUNK_IDX * CHUNK_SIZE + 1 ))
END_CHUNK=$(( (CHUNK_IDX + 1) * CHUNK_SIZE ))

echo "Running $SCENARIO, prior=$PRIOR, chunks $START_CHUNK to $END_CHUNK (array task $SLURM_ARRAY_TASK_ID)"

Rscript sim_estim.R "$SCENARIO" "$START_CHUNK" "$END_CHUNK" "$PRIOR"
