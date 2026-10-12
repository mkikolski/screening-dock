#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --job-name=dock_prod
#SBATCH --output=logs/dock_prod_%A_%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --partition=phd
#SBATCH --array=0-7
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=XXXX

set -euo pipefail

PROJECT="${DOCKPIPE_PROJECT:-/storage/users/XXXX/Gen_network/dockpipe}"
CONFIG="${DOCKPIPE_CONFIG:-$PROJECT/configs/panel500.yaml}"
ENV_NAME="${DOCKPIPE_ENV:-dockpipe}"
MAMBA_ROOT="${MAMBA_ROOT_PREFIX:-/storage/users/XXXX/micromamba}"

module purge

export MAMBA_ROOT_PREFIX="$MAMBA_ROOT"
eval "$("$MAMBA_ROOT/bin/micromamba" shell hook --shell bash)"
micromamba activate "$ENV_NAME"

export PYTHONNOUSERSITE=1
export PYTHONPATH="$PROJECT${PYTHONPATH:+:$PYTHONPATH}"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

cd "$PROJECT"
mkdir -p logs

SHARD_INDEX="${SLURM_ARRAY_TASK_ID:-0}"
SHARD_COUNT="${SLURM_ARRAY_TASK_COUNT:-1}"
WORKERS="${SLURM_CPUS_PER_TASK:-12}"

echo "Shard $SHARD_INDEX of $SHARD_COUNT, $WORKERS workers"

python run_pipeline.py dock \
  --config "$CONFIG" \
  --workers "$WORKERS" \
  --shard-index "$SHARD_INDEX" \
  --shard-count "$SHARD_COUNT"

echo "Shard $SHARD_INDEX complete."
