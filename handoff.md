# Handoff — 2026-10-26

**State: run complete, analysed.** Job 1994 finished 2026-10-18 22:44 on
`serwer1`, commit `9770a4e`. Nothing is pending.

## Result

The panel does not work as an antibacterial set, and the evidence is clean.

| comparison | favours bacterial | ≥0.5 kcal/mol | ≥1.0 |
|---|---|---|---|
| *E. coli* vs human DHFR (like-for-like) | **0 / 492** | 0 | 0 |
| best of any bacterial target vs human DHFR | 36 / 492 (7%) | 9 | 2 |

Those rows are not equally trustworthy. **The first is the one to quote:**
same fold, same cofactor, both receptors passed redocking (0.89 and 0.40 Å),
identical ligands — the only cross-target comparison this data genuinely
supports. Not one molecule out of 492 prefers the bacterial enzyme. The median
is −1.39 kcal/mol and even the 75th percentile is −1.07, so the entire
distribution sits on the wrong side.

The second row compares FabI and GyrB against human DHFR across different
folds, which Vina is far weaker at. It is also concentrated in the
**least-validated receptor**: nearly every near-miss is `fabi_saureus_3gr6`,
whose redock top pose was 5.91 Å (native at rank 2). So the only selectivity
signal in the panel sits on the target whose pose ranking we trust least.

Ligand efficiency is low throughout. Best is 0.352 — on the *human* target.
That is `...-00005047`, the smallest molecule in the panel at 21 heavy atoms,
which tops LE on all five receptors and is itself slightly human-preferring.
Excluding it, the next best is ~0.27 against a conventional hit floor of ~0.3.

Best single candidate in 492: `...-00004571` — FabI, −8.68 kcal/mol, 34 heavy
atoms, LE 0.26, beats human DHFR by 1.00. Still under the LE floor, and on the
weakest-validated target.

**It is systematic, not a bad sampling run.** Per generator seed the
selectivity margins are indistinguishable — medians −0.82 / −0.86 / −0.86,
permutation test on the widest gap p=0.521 (10k shuffles). The tail counts
differ (17/164, 6/134, 13/194) but a proportion test on the widest pair gives
p≈0.06, consistent with noise.

### What follows from it

The ligand set is named `exp_c5_admet_mpo_guided` — generated against an
ADMET/MPO objective with no selectivity term. The result is exactly what that
predicts: drug-like molecules that bind human DHFR better than any bacterial
target in the panel. If the goal is antibacterial, **human DHFR has to be an
explicit anti-target inside the generator's objective**, not a filter applied
afterwards. Filtering 492 molecules produced zero survivors; it will keep
producing zero.

### Correction worth keeping

I asserted twice during this work that Vina's score-vs-size bias would inflate
the rankings here. **It does not apply to this panel.** Measured: Pearson
+0.04 to +0.14 across the five receptors, i.e. slightly *positive* — bigger
molecules score marginally worse. The size range is narrow (median 41, min 21,
max 57) and the set is chemically homogeneous, which washes the correlation
out. Raw score is usable for ranking within this panel; LE still reorders it
because of the one small outlier. Don't reintroduce the caveat without
re-measuring.

## Reproducing / extending

```bash
cd /storage/users/kikolski/screening-dock
source /opt/modules/miniforge3/etc/profile.d/conda.sh
conda activate /storage/users/kikolski/envs/dockpipe

./check_run.sh                  # run status digest, read-only, pure bash
python analyze.py --top 15      # selectivity + LE + size bias + per-seed
python analyze.py --le-min 0.25 --margin-min 0.25   # looser thresholds
```

`analyze.py` writes `docking/analysis/{per_molecule,shortlist,top_by_le}.csv`.
It needs only pandas and pyyaml — no rdkit — because the heavy-atom counts are
already in `scores_by_ligand.csv`.

To rerun the screen, `sbatch run_all.slurm`. Both long stages checkpoint, so a
resubmit resumes from the last finished item rather than restarting:

- `results/panel500_v1/ligands/.ckpt_ligands.jsonl`
- `results/panel500_v1/docking/.ckpt_<receptor>.jsonl`

`python test_checkpoint.py` guards that logic; run it if you touch either
stage.

## What the run actually cost

| | |
|---|---|
| per docking | **104.8 min** (median 67 min, p90 130 min, max 5.3 h) |
| total | 134 h of the 168 h limit — 80% used |
| dockings | 2460 = 492 molecules × 5 receptors, **0 failures** |
| molecules | 492 of 500; 8 lost to `scrub_failed` |
| states | exactly 1 per molecule |

My pre-run estimate was 30 min per docking — **wrong by 3.4×**, because it came
from `head -10` of a sorted SMILES file. 48 h would have died at ~36% done and
120 h at ~90%. If you rerun with more ligands, shard across both nodes
(`--shard-count 2`) rather than trusting a time limit again; `merge_results.py`
already globs `scores*.csv` and dedups on `(ligand, receptor)`.

## Redocking validation: 3/5 PASS, as accepted

| receptor | top pose | best pose | rank | verdict |
|---|---|---|---|---|
| gyrb_saureus_3u2d | 0.53 | 0.53 | 1 | PASS |
| gyrb_ecoli_4zvi | 9.89 | 1.47 | 5 | FAIL (ranking) |
| dhfr_ecoli_7mym | 0.89 | 0.89 | 1 | PASS |
| dhfr_human_2w3a | 0.40 | 0.40 | 1 | PASS |
| fabi_saureus_3gr6 | 5.91 | 1.71 | 2 | FAIL (ranking) |

All five reproduce the native pose within 2 Å somewhere in the top 9, so the
boxes, receptors and cofactors are correct. The two failures are Vina ranking a
non-native pose first, and neither has a configuration fix — both were checked
rather than assumed:

- **3GR6** — triclosan stacks face-to-face on the NAD+ nicotinamide and Vina
  under-rewards pi-stacking. Not a missing water: the closest HOH is 7.07 Å and
  nothing lies within 4.5 Å.
- **4ZVI** — only a marginal 3.33 Å water contact exists, and 3U2D passes at
  0.53 Å with its four proper 2.72 Å waters excluded, so a missing bridging
  water does not explain it either.

**Don't tighten the boxes** to "fix" these. The control must use the same box as
the screen, and the panel is median 41 heavy atoms against triclosan's 17.
22³ = 10,648 Å³, well under Vina's ~27,000 Å³ degradation point.

The `GLU A:230` / `THR A:392` RuntimeErrors in the job log are **not failures**
— they are the first receptor-prep attempt with rebuilt heavy atoms, which then
succeeded without them. Both GyrB receptors produced full 492-row sets.

## Cluster facts, hard-won

The env must be **writable by you** *and* **mounted on the compute nodes**.
Only `/storage` is both. This cost four failed submits.

| location | writable | on compute nodes |
|---|---|---|
| `/opt/modules/miniforge3/envs` | no — `general:general`, we are not in it | yes |
| `~/.conda/envs` | yes | **no** |
| `/storage/users/$USER/envs` | yes | yes |

`$HOME` is the trap: an env there works interactively on the login node and
then does not exist inside a batch job, which reads as a broken env rather than
a missing mount. The check that catches it:

```bash
srun -p standard -n1 ls /storage/users/$USER/envs/dockpipe/bin/python
```

Others, all now fixed in the repo:

- `logs/` must exist at submit time — sbatch opens `--output` *before* running
  the script, so `mkdir -p logs` inside it is too late. Tracked via
  `logs/.gitkeep`.
- `--export=NONE` leaves the environment so bare that Lmod's own init dies
  reading `LOGNAME`. Now `--export=ALL`, with `PYTHONPATH`, `PYTHONHOME`,
  `LD_LIBRARY_PATH`, `VIRTUAL_ENV`, `CONDA_PREFIX` unset explicitly so the
  submitting shell cannot shadow the env.
- `conda env create` needs `-p <full path>`, never `-n`: `-n` writes to
  `envs_dirs[0]` = `~/.conda/envs`, back onto the unshared filesystem.
- `conda.sh` and Lmod's init are not `set -eu` clean; the setup block relaxes
  both, restores them, then verifies `python` is on `PATH` by hand.
- Partition `standard` is the only one: 2 nodes × 32 CPUs, `MaxTime=UNLIMITED`.

## Three silent-wrong-answer bugs from this work

1. **The redock control was measuring the wrong molecule.** A PDB entry stores
   no connectivity for ligands, and `ligand_lines_to_pdbqt` was always called
   with `template_smiles=None`, so RDKit perceived bonds by proximity and made
   every one single; `SanitizeMol` + `AddHs` then saturated the aromatics. Each
   reference ligand went into Vina with 6–8 bonds wrong. Bond orders now come
   from the chemical component's ideal SDF, fetched and cached per ligand. 3U2D
   went 2.46 → 0.53 Å on that alone; the two DHFR "passes" beforehand were luck.
2. **Only pose 1 was scored** — `_heavy_atom_coords_pdbqt` broke at the first
   `ENDMDL`, so a native pose at rank 2–5 read as a failure.
3. **`PYTHONPATH` was appended, not set** — under `--export=ALL` a stale
   `PYTHONPATH` or an active virtualenv could have made the job import
   different code than the commit it records in every manifest.

Also worth remembering: meeko ≥0.6 treats the receptor as a polymer needing a
template per residue and ships none for NAP/NDP, so cofactors are prepared
separately and typed rigid by obabel (Vina reads only the atom-type column, so
bond orders and partial charges do not need recovering from coordinates).
Receptor prep drops exactly one residue across the whole panel: TYR A:14 in
3U2D, 25–27 Å from the box centre, backbone+CB only in the crystal — verified
harmless. Per-receptor reports are in
`results/.../receptors/<name>/<name>.prep.log`.

## Open items

- `--cpus-per-task=32` takes a whole node and there are only 2. If OOM-killed,
  lower cpus rather than raise `--mem`.
- GyrB `keep_hetatms: []` — the 2 MG in 3U2D were never confirmed to sit
  outside the pyrrolamide site. 3U2D redocks at 0.53 Å without them, so this is
  not urgent.
- The 8 `scrub_failed` molecules are listed in
  `results/panel500_v1/manifests/ligands.json` under `failures`, if the
  denominator ever needs defending.
