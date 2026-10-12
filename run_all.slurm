#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --job-name=dockpipe
#SBATCH --output=logs/dockpipe_%j.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
# This cluster: partition `standard` (the only one), 2 nodes x 32 CPUs x ~31 GB,
# MaxTime=UNLIMITED. 32 cpus takes a whole node; drop to 16 to share one.
#SBATCH --cpus-per-task=32
#SBATCH --mem=24G
# No cluster time limit, so this is only a safety net. Raise it rather than risk
# a timeout: the dock stage has no checkpointing and restarts from zero.
#SBATCH --time=48:00:00
#SBATCH --partition=standard
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mikolaj.kikolski@pwr.edu.pl

set -euo pipefail

# Must be the directory CONTAINING dockpipe/, not dockpipe/ itself: the stages
# below run run_pipeline.py, which sits at the repo root.
PROJECT="${DOCKPIPE_PROJECT:-/storage/users/kikolski/screening-dock}"
CONFIG="${DOCKPIPE_CONFIG:-$PROJECT/configs/panel500.yaml}"
ENV_NAME="${DOCKPIPE_ENV:-dockpipe}"
# ENV_NAME may be a name in $CONDA_ROOT/envs, or an absolute prefix path if the
# env lives in your own storage.
CONDA_ROOT="${DOCKPIPE_CONDA_ROOT:-/opt/modules/miniforge3}"

# --export=NONE means Lmod's shell function is not inherited either, and a
# non-login bash does not source /etc/profile.d. Without this, `module` is
# simply not a command inside the job.
if ! command -v module >/dev/null 2>&1; then
  for init in /etc/profile.d/lmod.sh /etc/profile.d/modules.sh \
              /usr/share/lmod/lmod/init/bash; do
    # shellcheck disable=SC1090
    [[ -r "$init" ]] && source "$init" && break
  done
fi

# Purge only if there is a module system; without this guard `set -e` would
# abort the job on a host that has none, which the conda path does not need.
if command -v module >/dev/null 2>&1; then
  module purge
elif [[ -n "${DOCKPIPE_MODULE:-}" ]]; then
  echo "DOCKPIPE_MODULE=$DOCKPIPE_MODULE but no module command found" >&2
  exit 1
fi

# Two ways to get the environment. Set DOCKPIPE_MODULE to use the Lmod
# modulefile (modulefiles/dockpipe/phase1.lua); otherwise activate the conda env
# directly from the miniforge install. Loading the module also supplies
# DOCKPIPE_PROJECT, but $PROJECT was resolved above, so an explicit
# DOCKPIPE_PROJECT still wins.
if [[ -n "${DOCKPIPE_MODULE:-}" ]]; then
  module use "${DOCKPIPE_MODULEPATH:-$PROJECT/modulefiles}"
  module load "$DOCKPIPE_MODULE"
else
  # conda.sh defines the shell function; `conda activate` does not work in a
  # non-interactive shell without it.
  source "$CONDA_ROOT/etc/profile.d/conda.sh"
  conda activate "$ENV_NAME"
fi

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
for mod in ("rdkit", "gemmi", "meeko", "molscrub", "vina", "pandas",
            "pdbfixer", "openmm"):
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

# A reinstallable record of the whole env, next to the logs. Which of these
# works depends on how the env was provided, so try each.
{
  conda list --export 2>/dev/null \
    || python -m pip freeze 2>/dev/null \
    || true
} > "logs/env_${SLURM_JOB_ID:-local}.txt"

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
