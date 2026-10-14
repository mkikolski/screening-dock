#!/usr/bin/env bash
# ALL (SLURM's default) so the job inherits a usable environment -- with
# NONE even LOGNAME is unset and Lmod's own init script dies on it.
# Reproducibility does not depend on this: the job purges modules,
# activates one explicit env prefix, dumps the resolved environment to
# logs/env_<jobid>.txt, and records the commit and package versions in
# every manifest. What DOES leak in is sanitized below.
#SBATCH --export=ALL
#SBATCH --job-name=dockpipe
# Relative to the SUBMIT directory, so submit from the project root:
#   cd /storage/users/kikolski/screening-dock && sbatch run_all.slurm
# sbatch opens this file before the script runs, so logs/ must already
# exist -- it ships with the repo via logs/.gitkeep.
#SBATCH --output=logs/dockpipe_%j.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
# This cluster: partition `standard` (the only one), 2 nodes x 32 CPUs x ~31 GB,
# MaxTime=UNLIMITED. 32 cpus takes a whole node; drop to 16 to share one.
#SBATCH --cpus-per-task=32
#SBATCH --mem=24G
# Measured: ~30 min per docking at exhaustiveness 32 on this panel, so 500
# ligand states x 5 receptors is 39-78 h at 32 workers. MaxTime=UNLIMITED
# here, so this is only a safety net and is deliberately well clear of the
# estimate. Both long stages checkpoint, so rerunning the same command
# resumes rather than starting over.
#SBATCH --time=168:00:00
#SBATCH --partition=standard
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mikolaj.kikolski@pwr.edu.pl

set -euo pipefail

# Must be the directory CONTAINING dockpipe/, not dockpipe/ itself: the stages
# below run run_pipeline.py, which sits at the repo root.
PROJECT="${DOCKPIPE_PROJECT:-/storage/users/kikolski/screening-dock}"
CONFIG="${DOCKPIPE_CONFIG:-$PROJECT/configs/panel500.yaml}"
CONDA_ROOT="${DOCKPIPE_CONDA_ROOT:-/opt/modules/miniforge3}"

# Variables that would let the submitting shell shadow the env's own
# packages: a stale PYTHONPATH or an active virtualenv silently wins over
# conda activate, and LD_LIBRARY_PATH breaks compiled extensions in ways that
# surface hours later as an import error. Dropped before anything is set up.
unset PYTHONPATH PYTHONHOME PYTHONSTARTUP LD_PRELOAD LD_LIBRARY_PATH
unset VIRTUAL_ENV CONDA_PREFIX CONDA_DEFAULT_ENV

# Third-party init scripts run with `set -u` off and `set -e` relaxed. Neither
# Lmod's init nor conda.sh is written to survive those options, and an unbound
# variable inside someone else's script must not kill a 40-hour job.
set +eu

# A login shell supplies these; a batch job is not guaranteed one, and Lmod's
# init reads LOGNAME unguarded. Derived here, after set +eu, because getent can
# legitimately fail and must not abort anything.
export LOGNAME="${LOGNAME:-${SLURM_JOB_USER:-$(id -un)}}"
export USER="${USER:-$LOGNAME}"
export HOME="${HOME:-$(getent passwd "$LOGNAME" | cut -d: -f6)}"

# Found rather than hardcoded: the env has lived in two different places
# during setup, and a wrong guess here fails the job after the queue wait.
# Always a full prefix path -- `conda activate dockpipe` only resolves names
# inside conda's own envs_dirs, which none of these are. A real env has both
# bin/python and conda-meta/; a half-created directory has neither, which is
# exactly what "Not a conda environment" reports.
if [[ -n "${DOCKPIPE_ENV:-}" ]]; then
  ENV_CANDIDATES=("$DOCKPIPE_ENV")
else
  ENV_CANDIDATES=(
    "$HOME/.conda/envs/dockpipe"
    "/storage/users/$LOGNAME/envs/dockpipe"
    "$CONDA_ROOT/envs/dockpipe"
  )
fi

ENV_NAME=""
for cand in "${ENV_CANDIDATES[@]}"; do
  if [[ -x "$cand/bin/python" && -d "$cand/conda-meta" ]]; then
    ENV_NAME="$cand"
    break
  fi
done

if [[ -z "$ENV_NAME" ]]; then
  echo "No usable dockpipe conda env found. Tried:" >&2
  printf '  %s\n' "${ENV_CANDIDATES[@]}" >&2
  echo "Set DOCKPIPE_ENV=/path/to/env, or create one:" >&2
  echo "  conda env create -p <path> -f $PROJECT/environment.yml" >&2
  exit 1
fi
echo "env    : $ENV_NAME"

if ! command -v module >/dev/null 2>&1; then
  for init in /etc/profile.d/lmod.sh /etc/profile.d/modules.sh \
              /usr/share/lmod/lmod/init/bash; do
    # shellcheck disable=SC1090
    [[ -r "$init" ]] && source "$init" && break
  done
fi

if command -v module >/dev/null 2>&1; then
  module purge
fi

# Two ways to get the environment. Set DOCKPIPE_MODULE to use the Lmod
# modulefile (modulefiles/dockpipe/phase1.lua); otherwise activate the conda env
# directly from the miniforge install. Loading the module also supplies
# DOCKPIPE_PROJECT, but $PROJECT was resolved above, so an explicit
# DOCKPIPE_PROJECT still wins.
if [[ -n "${DOCKPIPE_MODULE:-}" ]]; then
  if ! command -v module >/dev/null 2>&1; then
    echo "DOCKPIPE_MODULE=$DOCKPIPE_MODULE but no module command found" >&2
    exit 1
  fi
  module use "${DOCKPIPE_MODULEPATH:-$PROJECT/modulefiles}"
  module load "$DOCKPIPE_MODULE"
else
  # conda.sh defines the shell function; `conda activate` does not work in a
  # non-interactive shell without it.
  source "$CONDA_ROOT/etc/profile.d/conda.sh"
  conda activate "$ENV_NAME"
fi

set -eu

# Verify by hand what `set -e` was not watching for above. Everything after
# this point is ours again and aborts on error as normal.
if ! command -v python >/dev/null 2>&1; then
  echo "No python on PATH after environment setup." >&2
  echo "  CONDA_ROOT=$CONDA_ROOT" >&2
  echo "  ENV_NAME=$ENV_NAME" >&2
  echo "  CONDA_PREFIX=${CONDA_PREFIX:-<unset>}" >&2
  exit 1
fi
echo "python : $(command -v python)"

export PYTHONNOUSERSITE=1
# Set exactly, never appended: an inherited PYTHONPATH is the thing most
# likely to make this job import a different copy of the code than the commit
# recorded below.
export PYTHONPATH="$PROJECT"
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
# Required by three of the five receptors: Meeko has no residue template
# for NAP/NDP, so the cofactor is typed by obabel and appended.
obabel = shutil.which("obabel")
print(f"{'obabel':<9}: {obabel or 'MISSING'}")
if not obabel:
    bad.append("obabel")
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

# --allow-fail is required, not optional: gyrb_ecoli_4zvi and
# fabi_saureus_3gr6 fail the 2 A top-pose gate for documented reasons (Vina
# mis-ranks a pose it did find -- see README). Without the flag the gate
# sys.exit()s and set -e kills the job before any docking starts.
srun --job-name=redock \
  python run_pipeline.py redock --config "$CONFIG" --allow-fail

srun --job-name=prep-ligands \
  python run_pipeline.py prep-ligands --config "$CONFIG" --workers "$WORKERS"

srun --job-name=dock \
  python run_pipeline.py dock --config "$CONFIG" --workers "$WORKERS"

srun --job-name=merge \
  python merge_results.py --config "$CONFIG"

echo "=== pipeline complete (commit $COMMIT) ==="
