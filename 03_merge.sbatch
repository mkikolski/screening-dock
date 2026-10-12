#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --job-name=dock_merge
#SBATCH --output=logs/dock_merge_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --partition=phd

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

cd "$PROJECT"
python merge_results.py --config "$CONFIG"
echo "Merge complete."
