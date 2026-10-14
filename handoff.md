# Handoff — 2026-10-14

**State: submitted and running.** All work is pushed (`origin/master` ==
`HEAD`). The env is at `/storage/users/kikolski/envs/dockpipe`.

Nothing needs doing unless something below says so.

## While it runs

The job log is `logs/dockpipe_<jobid>.log` in the project root. The first ~2
min should show, in order:

```
env    : /storage/users/kikolski/envs/dockpipe
python : /storage/users/kikolski/envs/dockpipe/bin/python
commit : <40-char sha>
rdkit    : 2024.03.5          <- the dependency gate; all must be present
...
obabel   : /storage/.../bin/obabel
```

Then prep-receptors (~20 s), redock (~10 min, 3/5 PASS expected),
prep-ligands, then docking for hours.

Two things worth checking once, in the first hour:

```bash
# 1. how many ligand states molscrub produced -> which runtime row applies
python -c "import json; print(len(json.load(open(
  'results/panel500_v1/ligands/ligand_index.json'))))"

# 2. the live rate, printed every 25 dockings
grep "per docking" logs/dockpipe_*.log | tail -3
```

| ligand states | dockings | wall at 32 workers |
|---|---|---|
| 500 | 2,500 | 39 h |
| 750 | 3,750 | 59 h |
| 1,000 | 5,000 | 78 h |

`--time` is 168 h, so all three fit. The 30 min/docking figure came from
`head -10` of the SMILES file, so if that file is sorted the sample is biased —
trust check #2 over the table.

If the state count is near 1,000 and node 2 is free, `--shard-count 2` across a
2-task array halves the wall time; `merge_results.py` already globs
`scores*.csv` and dedups on `(ligand, receptor)`.

## If it dies

**Resubmit the same script.** Both long stages checkpoint and resume from the
last finished item:

- `results/panel500_v1/ligands/.ckpt_ligands.jsonl`
- `results/panel500_v1/docking/.ckpt_<receptor>.jsonl`

`scores.csv` is rewritten after every receptor and `ligand_index.json` after
the ligand stage, so even a killed run leaves usable output for `merge`.
Recorded failures are *not* retried — the resume line reports how many there
were, since every failure mode here is deterministic.

`python test_checkpoint.py` guards this logic; run it if you touch either
stage.

## When it finishes

```bash
ls results/panel500_v1/docking/scores_merged.csv   # merge ran
python -c "
import csv, statistics as st
r = list(csv.DictReader(open('results/panel500_v1/docking/scores_merged.csv')))
s = [float(x['seconds']) for x in r if x['seconds']]
print(f'rows={len(r)} failed={sum(x[\"status\"]!=\"ok\" for x in r)}')
print(f'median={st.median(s):.0f}s p90={sorted(s)[int(.9*len(s))-1]:.0f}s')"
```

Outputs: `scores_merged.csv` (every docking), `scores_by_ligand.csv`,
`scores_wide.csv` (ligand × receptor matrix). `seconds` is a real column, so
the per-docking distribution is available for the write-up.

## Redocking: 3/5 PASS, and that is the accepted state

| receptor | top pose | best pose | rank | verdict |
|---|---|---|---|---|
| gyrb_saureus_3u2d | 0.52 | 0.52 | 1 | PASS |
| gyrb_ecoli_4zvi | 9.89 | 1.86 | 4 | FAIL (ranking) |
| dhfr_ecoli_7mym | 0.89 | 0.85 | 2 | PASS |
| dhfr_human_2w3a | 0.41 | 0.41 | 1 | PASS |
| fabi_saureus_3gr6 | 5.90 | 1.70 | 2 | FAIL (ranking) |

All five reproduce the native pose within 2 Å somewhere in the top 9, so the
boxes, receptors and cofactors are correct. The two failures are Vina ranking
a non-native pose first. Neither has a configuration fix, and both were
checked rather than assumed:

- **3GR6** — triclosan stacks face-to-face on the NAD+ nicotinamide and Vina
  under-rewards pi-stacking. Not a missing water: the closest HOH is 7.07 Å,
  nothing lies within 4.5 Å.
- **4ZVI** — only a marginal 3.33 Å water contact exists, and 3U2D passes at
  0.52 Å with its four proper 2.72 Å waters excluded, so a missing bridging
  water does not explain it either.

Since a compound's score is the best of its 9 poses, this matters little for
ranking compounds *within* a receptor. Treat absolute cross-receptor
comparisons involving 4ZVI and 3GR6 as weaker evidence. Per-receptor notes are
in `configs/panel500.yaml`; the table is also in `README.md`.

**Don't tighten the boxes** to "fix" these. The control must use the same box
as the screen, and the panel is median 41 heavy atoms against triclosan's 17.
22³ = 10,648 Å³, well under Vina's ~27,000 Å³ degradation point.

## Cluster facts, hard-won

The env must be **writable by you** *and* **mounted on the compute nodes**.
Only `/storage` is both. This cost four failed submits.

| location | writable | on compute nodes |
|---|---|---|
| `/opt/modules/miniforge3/envs` | no — `general:general`, we are not in it | yes |
| `~/.conda/envs` | yes | **no** |
| `/storage/users/$USER/envs` | yes | yes |

`$HOME` is the trap: an env there works interactively on the login node and
then does not exist inside a batch job, which reads as a broken env rather
than a missing mount. The pre-flight check that catches it:

```bash
srun -p standard -n1 ls /storage/users/$USER/envs/dockpipe/bin/python
```

Other things that bit, all now fixed in the repo:

- `logs/` must exist at submit time — sbatch opens `--output` *before* running
  the script, so `mkdir -p logs` inside it is too late. Tracked via
  `logs/.gitkeep`.
- `--export=NONE` leaves the environment so bare that Lmod's own init script
  dies reading `LOGNAME`. Now `--export=ALL`, with `PYTHONPATH`, `PYTHONHOME`,
  `LD_LIBRARY_PATH`, `VIRTUAL_ENV`, `CONDA_PREFIX` etc. unset explicitly so the
  submitting shell cannot shadow the env.
- `conda env create` needs `-p <full path>`, never `-n`: `-n` writes to
  `envs_dirs[0]` = `~/.conda/envs`, i.e. back onto the unshared filesystem.
- `conda.sh` and Lmod's init are not `set -eu` clean; the setup block relaxes
  both, restores them, then verifies `python` is on `PATH` by hand.
- Partition is `standard` (the only one), 2 nodes × 32 CPUs,
  `MaxTime=UNLIMITED`.

## What was fixed this session

| commit | |
|---|---|
| `b808955` | cofactors prepared separately — meeko has no NAP/NDP template |
| `857492d` | leftover `het` in the preflight print |
| `7d06be9` | retry receptor prep without rebuilt heavy atoms |
| `2966b5f` | redock control: real bond orders, all 9 poses |
| `ad6cfe8` | recorded the redock verdicts |
| `13f143e` | dock stage checkpointing |
| `f436e9b` | four defects that would have cost the run |
| `125392a` | ship `logs/` so sbatch can open its output file |
| `80db7bc` | survive the environment `--export=NONE` gives the job |
| `f1efe0c` | `--export=ALL` + sanitize what leaks in |
| `9770a4e` | env on `/storage`, the only location that works |

Three silent-wrong-answer bugs worth remembering:

1. **The redock control was measuring the wrong molecule.** A PDB entry stores
   no connectivity for ligands, and `ligand_lines_to_pdbqt` was always called
   with `template_smiles=None`, so RDKit perceived bonds by proximity and made
   every one single; `SanitizeMol` + `AddHs` then saturated the aromatics. Each
   reference ligand went into Vina with 6–8 bonds wrong. Bond orders now come
   from the chemical component's ideal SDF, fetched and cached per ligand.
   3U2D went 2.46 → 0.52 Å on this alone; the two DHFR "passes" before the fix
   were luck.
2. **Only pose 1 was being scored** — `_heavy_atom_coords_pdbqt` broke at the
   first `ENDMDL`, so a native pose at rank 2–4 read as a failure.
3. **`PYTHONPATH` was appended, not set** — under `--export=ALL` a stale
   `PYTHONPATH` or an active virtualenv could have made the job import a
   different copy of the code than the commit it records in every manifest.

Receptor prep drops exactly one residue across the whole panel: TYR A:14 in
3U2D, 25–27 Å from the box center, backbone+CB only in the crystal. Verified
harmless. The dropped-residue report is in each
`results/.../receptors/<name>/<name>.prep.log`.

## Open items

- `--cpus-per-task=32` takes a whole node and there are only 2. If OOM-killed,
  lower cpus rather than raise `--mem`.
- GyrB `keep_hetatms: []` — the 2 MG in 3U2D were never confirmed to sit
  outside the pyrrolamide site. 3U2D redocks at 0.52 Å without them, so this
  is not urgent.
- `prep-ligands` failures (unparseable SMILES, scrub failures) are listed in
  `results/panel500_v1/manifests/ligands.json` under `failures` — worth a look
  for how much of the 500 actually made it through.
