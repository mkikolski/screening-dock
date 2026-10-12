#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --job-name=dockpipe
#SBATCH --output=logs/dockpipe_%j.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=24G
#SBATCH --time=16:00:00
#SBATCH --partition=phd
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mikolaj.kikolski@pwr.edu.pl

set -euo pipefail

# Must be the directory CONTAINING dockpipe/, not dockpipe/ itself: the stages
# below run run_pipeline.py, which sits at the repo root.
PROJECT="${DOCKPIPE_PROJECT:-/storage/users/kikolski/mateusz-dock}"
CONFIG="${DOCKPIPE_CONFIG:-$PROJECT/configs/panel500.yaml}"
ENV_NAME="${DOCKPIPE_ENV:-dockpipe}"
MAMBA_ROOT="${MAMBA_ROOT_PREFIX:-/storage/users/kikolski/micromamba}"

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

WORKERS="${SLURM_CPUS_PER_TASK:-12}"

# Code provenance: the commit every stage below runs from. The same values are
# written into each stage manifest by dockpipe.common.environment_provenance().
COMMIT="$(git -C "$PROJECT" rev-parse HEAD 2>/dev/null || echo unknown)"
DIRTY="$(git -C "$PROJECT" status --porcelain 2>/dev/null | head -c1)"
echo "=== dockpipe job ${SLURM_JOB_ID:-local} ==="
echo "commit : $COMMIT${DIRTY:+ (DIRTY working tree)}"
python -c "import sys, rdkit, meeko, vina; \
print('python :', sys.version.split()[0]); \
print('rdkit  :', rdkit.__version__); \
print('meeko  :', meeko.__version__); \
print('vina   :', vina.__version__)"
echo "workers: $WORKERS"
echo "config : $CONFIG"
echo "=========================================="

# One job, five job steps. Each srun is a separate step so `sacct -j <id>` shows
# per-stage timing and exit codes, and set -e aborts the rest on failure.
srun --job-name=prep-receptors \
  python run_pipeline.py prep-receptors --config "$CONFIG"

srun --job-name=redock \
  python run_pipeline.py redock --config "$CONFIG"

srun --job-name=prep-ligands \
  python run_pipeline.py prep-ligands --config "$CONFIG" --workers "$WORKERS"

srun --job-name=dock \
  python run_pipeline.py dock --config "$CONFIG" --workers "$WORKERS"

srun --job-name=merge \
  python merge_results.py --config "$CONFIG"

echo "=== pipeline complete (commit $COMMIT) ==="
