# dockpipe

AutoDock Vina screening pipeline: SMILES panel → protomers/3D (molscrub) → PDBQT
(Meeko) → docking against a configurable receptor panel → merged scores with
ligand efficiency and a selectivity axis.

## Layout

The Python modules form a package. `run_pipeline.py` imports `dockpipe.common`,
and `receptors.py` / `docking.py` use relative imports (`from .common import`),
so they must sit inside a `dockpipe/` directory with the entrypoints beside it:

```
<PROJECT>/
  run_pipeline.py          # CLI: prep-receptors | redock | prep-ligands | dock | all
  merge_results.py         # shard CSVs → merged/per-ligand/wide tables + summary.json
  configs/panel500.yaml
  dockpipe/
    __init__.py
    common.py              # config dataclasses, logging, provenance manifests
    receptors.py           # PDB fetch, split, protonate, PDBQT, box derivation
    ligands.py             # SMILES → molscrub states → Meeko PDBQT
    docking.py             # Vina driver, redocking control, scores CSV
  run_all.slurm            # whole pipeline, one job, one srun step per stage
```

`run_all.slurm` sets `PYTHONPATH="$PROJECT"`; that resolves `dockpipe.*` only
when `$PROJECT` is the directory *containing* `dockpipe/`.

## Stages

1. **`prep-receptors`** — fetch each PDB from RCSB (or PDB-REDO with
   `--pdb-redo`), filter to `keep_chains`, route HETATMs to either the receptor
   (`keep_hetatms`, i.e. cofactors) or the reference ligand
   (`redock_ligand_resname`), protonate with PDBFixer at the entry's `ph`, write
   receptor PDBQT via Meeko's `mk_prepare_receptor.py` (falling back to the ADFR
   suite's `prepare_receptor`), and derive the docking box from the reference
   ligand centroid. Writes `receptors/receptor_index.json`.
2. **`redock`** — re-docks each crystal ligand into its own site and compares to
   the crystal pose. RMSD is symmetry-corrected via spyrmsd when importable,
   otherwise naive and order-dependent. Gates the run at 2.0 Å; pass
   `--allow-fail` to continue past a failure.
3. **`prep-ligands`** — molscrub produces protonation states and 3D coords,
   Meeko writes one PDBQT per state (`<name>_<n>.pdbqt`). One `Scrub` is built
   per worker process. Writes `ligands/ligand_index.json`.
4. **`dock`** — every ligand state × every receptor through the Vina Python API.
   Vina maps are cached per worker and keyed on receptor, so `chunksize=1` keeps
   workers reusing their maps instead of recomputing them.
5. **`merge_results.py`** — concatenates shard CSVs, dedups on
   (ligand, receptor), collapses states to the best score per molecule, adds
   `ligand_efficiency` (−score / heavy atoms) and a `dhfr_selectivity` column
   (human minus bacterial; positive favours the bacterial target).

`all` runs stages 1–4 in order.

## Running

Locally, one stage at a time:

```bash
python run_pipeline.py prep-receptors --config configs/panel500.yaml
python run_pipeline.py redock         --config configs/panel500.yaml
python run_pipeline.py prep-ligands   --config configs/panel500.yaml --workers 32
python run_pipeline.py dock           --config configs/panel500.yaml --workers 32
python merge_results.py               --config configs/panel500.yaml
```

Validate the config first — seconds, and safe on a login node:

```bash
python preflight.py --config configs/panel500.yaml
```

It fetches each PDB, splits it, derives the box, and checks that `ligands_smi`
exists and that `output_root` has a writable ancestor. Non-zero exit means do
not submit yet. It does **not** check SLURM or the environment; the job script
does that itself and aborts before stage 1 if anything is missing.

Then a single submission:

```bash
sbatch run_all.slurm
squeue -u "$USER"            # the one job
sacct -j <jobid> --format=JobID,JobName,Elapsed,State,ExitCode
```

All five stages run as `srun` steps inside that one job, in order. `set -e`
means a failing stage — including the redocking gate — aborts the rest instead
of docking against a broken receptor; `sacct` shows which step failed, since
each step carries its stage name. The node's CPUs are used by the worker `Pool`
within a stage, so there is no job array and no sharding to coordinate.

### Cluster resources

This cluster has **one** partition, `standard` (also the default), with 2 nodes
× 32 CPUs × ~31 GB and `MaxTime=UNLIMITED`. The header asks for a whole node
(`--cpus-per-task=32`) for 48 h. `WORKERS` comes from `SLURM_CPUS_PER_TASK`, so
changing that one directive rescales the worker pool.

Re-check before assuming any of that still holds:

```bash
sinfo -s
scontrol show partition | grep -E 'PartitionName|MaxTime|MaxMemPerNode|State'
scontrol show config | grep -i selecttype   # linear => whole nodes regardless
```

`Requested partition configuration not available right now` means the ask can
never be satisfied as configured — a nonexistent partition name, or a `--time`,
`--mem` or `--cpus-per-task` above the partition's limits — not that the queue
is busy. `sinfo -p <part> -R` shows drained nodes. Command-line flags override
the `#SBATCH` directives, so alternatives can be tried without editing:

```bash
sbatch --partition=<other> --cpus-per-task=16 --time=04:00:00 run_all.slurm
```

Two cautions at 32 workers: each worker caches its own Vina maps, so if the job
is OOM-killed lower `--cpus-per-task` rather than raise `--mem` (the node has
only ~31 GB); and a 32-CPU request takes a whole node out of the two available,
so `--cpus-per-task=16` is the neighbourly option when someone else is active.

### How long it will take

Do not guess — measure. `scores*.csv` carries a `seconds` column per docking:

```bash
python -c "import pandas as pd; print(pd.read_csv('results/panel500_v1/docking/scores.csv').seconds.describe())"
```

Estimated hours ≈ `p90_seconds × n_states × n_receptors / workers / 3600`. Use
p90, not the median: Vina's cost tracks the ligand's torsion count and the tail
dominates. The current panel is ~500 molecules of median 41 heavy atoms, so
roughly 3750–5000 dockings at `exhaustiveness: 32`.

If it will not fit, in order of preference: raise `--time` (no cluster limit);
split by receptor into independent jobs (`--receptor`, and `merge_results.py`
already globs and dedups the shards); re-shard with `--shard-index` /
`--shard-count`; or halve `exhaustiveness` and re-check the redocking gate.
**The dock stage has no checkpointing** — a timeout loses every docking, so
prefer a long `--time` or smaller independent jobs over a tight limit.

Override paths with the `DOCKPIPE_PROJECT`, `DOCKPIPE_CONFIG`, `DOCKPIPE_ENV`,
`DOCKPIPE_CONDA_ROOT`, `DOCKPIPE_MODULE` and `DOCKPIPE_MODULEPATH` environment
variables.

Useful CLI flags: `--receptor <name>` restricts a run to one target,
`--ligands` overrides `ligands_smi`, `--workers` overrides `n_workers`,
`--shard-index` / `--shard-count` partition the ligand-state list (kept for
splitting a rerun by hand; `run_all.slurm` does not use them).

Outputs land under `<output_root>/<run_name>/`: `receptors/`, `ligands/`,
`docking/` (per-receptor poses, `scores*.csv`, `redock_control.json`,
`summary.json`), `manifests/` (package versions, git commit, SLURM ids and
SHA-256 checksums per stage), `logs/`.

## Selecting target PDBs

Edit the `receptors:` list in the config. Each entry:

```yaml
- name: dhfr_ecoli_7mym          # output directory and CSV key
  pdb_id: 7MYM                   # auto-downloaded; or local_pdb: /path/to.pdb
  keep_chains: [A]
  keep_hetatms: [NDP]            # cofactors kept as part of the receptor
  redock_ligand_resname: TOP     # extracted as the redocking reference
  ph: 7.4
  notes: >
    Free text, copied into the receptor record.
  box:
    reference_ligand_resname: TOP  # usually the same resname
    size: [20.0, 20.0, 20.0]
```

That list is the whole mechanism: `prepare_receptors` iterates it and `dock`
iterates the resulting index. Adding a target means adding an entry; a receptor
that raises during prep is recorded as `status: failed` and skipped downstream
rather than killing the run.

**The box needs a center from somewhere.** Either set
`box.reference_ligand_resname` (center = centroid of those HETATM atoms) or set
`box.center: [x, y, z]` explicitly, which bypasses ligand-derived boxes. With
both null, `derive_box` raises `No reference-ligand atoms found` and the
receptor fails. The two GyrB entries in `panel500.yaml` are in exactly that
state — placeholders to fill in after inspecting the entries. Omitting
`box.size` falls back to the ligand extent plus 8 Å and `2 × padding`.

Cofactor choices are load-bearing, not cosmetic: the `notes` fields record why
NADPH is kept for DHFR (TMP binding is cooperative; removing it changes the M20
loop pocket) and NADP for FabI (triclosan stacks on the nicotinamide ring).
Check for structural waters bridging the pocket before dropping all HETATMs —
add `HOH` to `keep_hetatms` if one sits in the site.

## Exhaustiveness and Vina settings

Set under `vina:` in the config; there is no CLI override.

```yaml
vina:
  scoring: vina        # sf_name passed to Vina
  exhaustiveness: 32   # search effort
  n_poses: 9
  seed: 2137           # fixed, so runs reproduce
  cpu: 1               # one thread per Vina; parallelism comes from the pool
  min_rmsd: 1.0        # pose-clustering cutoff
```

The value flows `PipelineConfig.from_yaml` → `VinaSpec` → `DockTask` →
`v.dock(exhaustiveness=...)` in `docking.py`. If the `vina:` block is absent the
`VinaSpec` defaults apply (exhaustiveness 32, seed 42).

**Choosing a value.** Vina's own default is 8; 32 is the usual choice for
ranking a panel rather than refining a single pose. Cost is roughly linear, so
64 doubles wall time for a modest reproducibility gain. Two things should drive
the number up:

- **Box volume.** `derive_box` warns above ~27000 Å³ because a larger search
  space needs proportionally more sampling to cover. The 20–22 Å cubes here are
  8000–10648 Å³, comfortably inside that.
- **Ligand flexibility.** Panels with many rotatable bonds need more effort to
  converge; a panel of rigid fragments needs less.

The empirical check is the redocking gate plus a seed sweep: if re-running a
handful of ligands with different `seed` values moves the top score by more than
a few tenths of a kcal/mol, exhaustiveness is too low for that panel.
`cpu: 1` is deliberate — Vina's internal threading and the worker `Pool` would
otherwise oversubscribe the node, which is also why run_all.slurm exports
`OMP_NUM_THREADS=1`.

## Creating the environment

The miniforge install at `/opt/modules/miniforge3` is not on `PATH` and has no
module of its own, so `conda` has to be bootstrapped by sourcing its init
script. That is all `run_all.slurm` needs at run time, and it is how you create
the env once:

`envs/` there is group-writable (`drwxrwxr-x`), so the env can be created in
place under its own name, which is where the modulefile already looks. Confirm
you are in the owning group first — the mode only grants write to that group,
and `other` has `r-x`:

```bash
stat -c '%A %U:%G' /opt/modules/miniforge3/envs   # note the group
id -nG                                            # must contain it
```

If it does:

```bash
source /opt/modules/miniforge3/etc/profile.d/conda.sh   # defines `conda`
export CONDA_PKGS_DIRS=/storage/users/kikolski/conda_pkgs

conda env create -f environment.yml          # -> envs/dockpipe, per `name:`
conda activate dockpipe
```

In practice conda fell back to `~/.conda/envs/dockpipe` rather than
`miniforge3/envs`, which means the group write bit did not apply to this
account. That works, but verify two things before a long job — a conda env is
several GB, and `$HOME` is not guaranteed to be mounted on the compute nodes:

```bash
du -sh ~/.conda/envs/dockpipe        # typically 3-6 GB
quota -s 2>/dev/null || df -h ~      # is there room
srun -n1 -t 00:01:00 ls ~/.conda/envs/dockpipe/bin/python   # visible on a node?
```

If that `srun` fails, the env is unusable from the job and has to be rebuilt on
shared storage:

```bash
conda env create -f environment.yml -p /storage/users/kikolski/envs/dockpipe
```

then set `DOCKPIPE_ENV` to that path and update `env_root` in
`modulefiles/dockpipe/phase1.lua` to match. `conda activate dockpipe` by name
keeps working either way, because conda searches all of `envs_dirs`.

`CONDA_PKGS_DIRS` is still worth setting: only `envs/` was opened up, and the
package cache lives next to it at `miniforge3/pkgs`. Without a writable cache
conda falls back to `~/.conda/pkgs`, and these packages run to several GB
against a home quota that usually cannot take them.

If you are *not* in that group, build the env in your own storage instead and
pass the prefix — `DOCKPIPE_ENV` accepts either a name or a path:

```bash
conda env create -f environment.yml -p /storage/users/kikolski/envs/dockpipe
DOCKPIPE_ENV=/storage/users/kikolski/envs/dockpipe sbatch run_all.slurm
```

In that case also update `env_root` in `modulefiles/dockpipe/phase1.lua` (there
is a commented-out line for it) before using the module path.

## Environment via Lmod

`modulefiles/dockpipe/phase1.lua` wraps the conda env as a module, following the
group's convention (PATH prepend plus `CONDA_*`, no `conda activate`). To use it
instead of activating the env directly, set `DOCKPIPE_MODULE`:

```bash
DOCKPIPE_MODULE=dockpipe/phase1 sbatch run_all.slurm
# or interactively:
module use /storage/users/kikolski/screening-dock/modulefiles
module load dockpipe/phase1
```

With `DOCKPIPE_MODULE` unset the script instead sources
`$DOCKPIPE_CONDA_ROOT/etc/profile.d/conda.sh` (default
`/opt/modules/miniforge3`) and runs `conda activate "$DOCKPIPE_ENV"`, so both
work and nothing breaks before the module is installed. `DOCKPIPE_ENV` may be a
name under `$CONDA_ROOT/envs` or an absolute prefix path if the env lives in
your own storage. `DOCKPIPE_MODULEPATH` overrides where `module use` looks, for
when an admin installs the file on the system `MODULEPATH` instead.

The modulefile sets `DOCKPIPE_PROJECT`, so a `module load` alone points the
pipeline at the right checkout. It deliberately does **not** set thread limits:
Vina runs single-threaded with parallelism from the worker pool, which is a
property of the job, not the environment.

Two things to know:

- **A `vina` module is not a substitute.** `docking.py` uses the Vina Python
  API (`from vina import Vina`); a CLI-only `module load vina` does not satisfy
  that import, and the job would die at the docking stage after prep had already
  run. The bindings come from `pip install vina` inside the env.
- **`module` may not exist inside the job.** `--export=NONE` stops Lmod's shell
  function from being inherited and a non-login bash does not source
  `/etc/profile.d`, so the script sources the Lmod init itself before calling
  `module purge`. Any `module load` must come after that purge, since modules
  loaded on the login node do not survive into the job.

`family("conda_env")` in the modulefile prevents it loading alongside another
conda-env module, which would fight over `PATH` and `CONDA_PREFIX`. The guard is
only mutual if the other modules declare the same family.

## Environment

`environment.yml` builds the `dockpipe` conda env: Python 3.11, RDKit 2024.03.5,
OpenMM + PDBFixer from conda-forge, then meeko ≥ 0.6.0, molscrub, vina ≥ 1.2.5
and spyrmsd from pip. spyrmsd is optional at runtime but without it redocking
RMSD is not symmetry-corrected, which inflates it for symmetric ligands.

## Redocking validation

Run `redock` before any screen. It re-docks each co-crystal ligand into its own
box and reports RMSD against the crystal pose. Bond orders come from the
chemical component's ideal SDF, fetched from RCSB -- a PDB entry stores no
connectivity for ligands, so without that template RDKit guesses bonds by
proximity and the control measures a saturated analogue rather than the real
ligand.

Result for `configs/panel500.yaml` (exhaustiveness 32, 9 poses):

| receptor | top pose | best pose | rank | verdict |
|---|---|---|---|---|
| gyrb_saureus_3u2d | 0.52 | 0.52 | 1 | PASS |
| gyrb_ecoli_4zvi   | 9.89 | 1.86 | 4 | FAIL (ranking) |
| dhfr_ecoli_7mym   | 0.89 | 0.85 | 2 | PASS |
| dhfr_human_2w3a   | 0.41 | 0.41 | 1 | PASS |
| fabi_saureus_3gr6 | 5.90 | 1.70 | 2 | FAIL (ranking) |

All five reproduce the native pose within 2 A somewhere in the top 9, so the
boxes, receptors and cofactors are correct. The two failures are Vina ranking
a non-native pose first, not a setup error, and neither has a configuration
fix (see the per-receptor notes in the config). Since a compound's reported
score is the best of its 9 poses, this matters little for ranking compounds
*within* a receptor; treat absolute cross-receptor comparisons involving
gyrb_ecoli_4zvi and fabi_saureus_3gr6 as weaker evidence.

Box sizes are deliberately larger than the crystal ligands need: the control
has to use the same box as the screen, and the panel is median 41 heavy atoms
against triclosan's 17.

## Provenance: which code produced a result

The repo is under git, and every stage records the commit it ran from. Two
places to look:

- The job log header (`logs/dockpipe_<jobid>.log`) prints `commit : <sha>` once,
  plus `(DIRTY working tree)` if the submitted tree had uncommitted changes.
- Each stage manifest under `<output_root>/<run_name>/manifests/*.json` carries
  `provenance.git_commit` and `provenance.git_dirty` alongside package versions,
  hostname, SLURM ids and SHA-256 checksums of that stage's outputs.

The hash comes from the *source tree*, not the working directory —
`dockpipe.common._safe` runs git with `cwd=CODE_ROOT` (the parent of
`dockpipe/`), so it stays correct no matter where the job is launched from.

A `git_dirty: true` manifest means the result is not reproducible from any
commit. Commit before submitting a run whose numbers you intend to keep:

```bash
git commit -am "panel500: <what changed>" && sbatch run_all.slurm
```

To go from a result back to its code: `git checkout $(jq -r
.provenance.git_commit manifests/dock.json)`.

## Other files

- `exp_c5_selected_500.smi` — the 500-ligand panel the config points at, two
  fields per line (SMILES, name). Tracked in the repo deliberately: it travels
  with the code, so the commit hash in a manifest identifies the panel that
  produced the scores as well as the code that ran them. Override per run with
  `--ligands`.
`.env` is gitignored on purpose: nothing in the pipeline loads one, and a
`env > .env` dump describes whichever machine made it — sourcing a laptop's
copy on a compute node would point `PATH` and `CONDA_PREFIX` at paths that do
not exist there. Configuration lives in `configs/*.yaml` plus the `DOCKPIPE_*`
overrides.

Removed: `analyze_500.py` (prototype that `ligands.py` replaced),
`draft.slurm` (Meeko-CLI sketch that never reached docking),
`01_prep.sbatch` / `02_dock.sbatch` / `03_merge.sbatch` (superseded by
`run_all.slurm`), and a stray `.env` shell dump. All but the `.env`, which was
never tracked, are in git history if ever needed.
