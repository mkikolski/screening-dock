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

# Code provenance: the commit every stage below runs from. The same values, plus
# package versions, resolved tool paths and $LOADEDMODULES, are written into
# every stage manifest by dockpipe.common.environment_provenance().
COMMIT="$(git -C "$PROJECT" rev-parse HEAD 2>/dev/null || echo unknown)"
DIRTY="$(git -C "$PROJECT" status --porcelain 2>/dev/null | head -c1)"

echo "=== dockpipe job ${SLURM_JOB_ID:-local} on $(hostname) ==="
echo "commit : $COMMIT${DIRTY:+ (DIRTY working tree -- result not reproducible)}"
echo "config : $CONFIG"
echo "workers: $WORKERS"
echo "modules: ${LOADEDMODULES:-<none>}"

# Fail before the queue time is spent, not at stage 1. The Python `vina`
# package is what docking.py imports -- a `module load vina` CLI does not
# satisfy it. Receptor prep additionally shells out to a separate binary.
python - <<'PY'
import importlib, shutil, sys
bad = []
for mod in ("rdkit", "meeko", "molscrub", "vina", "pandas", "pdbfixer", "openmm"):
    try:
        m = importlib.import_module(mod)
        print(f"{mod:<9}: {getattr(m, '__version__', 'unknown')}")
    except Exception as e:
        print(f"{mod:<9}: MISSING ({e.__class__.__name__})")
        bad.append(mod)
try:
    import spyrmsd
    print(f"{'spyrmsd':<9}: {spyrmsd.__version__}")
except Exception:
    print(f"{'spyrmsd':<9}: absent -- redocking RMSD will NOT be "
          "symmetry-corrected (inflated for symmetric ligands)")
prep = [t for t in ("mk_prepare_receptor.py", "prepare_receptor")
        if shutil.which(t)]
print(f"{'rec-prep':<9}: {prep[0] + ' -> ' + shutil.which(prep[0]) if prep else 'MISSING'}")
if not prep:
    bad.append("mk_prepare_receptor.py/prepare_receptor")
if bad:
    sys.exit("missing dependencies: " + ", ".join(bad))
PY
echo "=========================================="

# A resolvable record of the whole env, next to the results.
micromamba list --export > "logs/env_${SLURM_JOB_ID:-local}.txt" 2>/dev/null || true

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
