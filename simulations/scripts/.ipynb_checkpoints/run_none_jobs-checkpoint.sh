#!/bin/bash
#SBATCH --account=def-mlegault
#SBATCH --job-name=mr_clust_none_jobs
#SBATCH --time=20:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=../results/logs/%x_%A_%a.out
#SBATCH --array=0-299

cd ~/links/projects/def-mlegault/gf591137/BayesClustMR/simulations/scripts/
module load r/4.5.0

SCENARIOS=("D_S1" "D_S2" "D_S3")
PRIOR="none"
N_CHUNKS=100
CHUNK_SIZE=10
N_SCENARIOS=${#SCENARIOS[@]}

# Decompose SLURM_ARRAY_TASK_ID into scenario and chunk indices only
SCENARIO_IDX=$(( SLURM_ARRAY_TASK_ID / N_CHUNKS ))
CHUNK_IDX=$(( SLURM_ARRAY_TASK_ID % N_CHUNKS ))

SCENARIO=${SCENARIOS[$SCENARIO_IDX]}
START_CHUNK=$(( CHUNK_IDX * CHUNK_SIZE + 1 ))
END_CHUNK=$(( (CHUNK_IDX + 1) * CHUNK_SIZE ))

echo "Running $SCENARIO, prior=$PRIOR, chunks $START_CHUNK to $END_CHUNK (array task $SLURM_ARRAY_TASK_ID)"

Rscript sim_estim.R "$SCENARIO" "$START_CHUNK" "$END_CHUNK" "$PRIOR"