#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --job-name=dock_prep
#SBATCH --output=logs/dock_prep_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --partition=phd
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

echo "=== environment ==="
python -c "import sys, rdkit, meeko, vina; \
print('python', sys.version.split()[0]); \
print('rdkit ', rdkit.__version__); \
print('meeko ', meeko.__version__); \
print('vina  ', vina.__version__)"
echo "==================="

WORKERS="${SLURM_CPUS_PER_TASK:-12}"

python run_pipeline.py prep-receptors --config "$CONFIG"

python run_pipeline.py redock --config "$CONFIG"

python run_pipeline.py prep-ligands --config "$CONFIG" --workers "$WORKERS"

echo "Prep complete. Submit 02_dock.sbatch next."
