# Handoff — 2026-11-03

**State: run complete and analysed; a LaTeX report exists but still needs real
data.** Job 1994 finished 2026-10-18 22:44 on `serwer1`, commit `9770a4e`.

**Pick up at step 1 below.** The PDF in `report/` currently builds from a
synthetic fixture, so the numbers in it are not yours.

## Next steps

### 1. Build the real PDF (~10 min, blocks nothing else)

```bash
# cluster
cd /storage/users/kikolski/screening-dock/results/panel500_v1
tar czf /tmp/report_data.tgz -C docking analysis/per_molecule.csv \
    redock_control.json summary.json -C ../manifests ligands.json

# locally
cd ~/mateusz-dock && git pull
mkdir -p report/data && tar xzf /path/to/report_data.tgz -C report/data
cd report && make venv && make      # make venv once: system python3 lacks pandas
```

Confirm it is real data, not the fixture: the near-miss table must show **9**
rows and the selectivity table **0 / 492**. The fixture gives 7 and 8.

### 2. FabL, if the isoform argument matters (~27 h, unattended)

`target_choice.md` §4.4 lists FabL as priority 4, and §5's FabI+FabL
cross-isoform comparison is the stated justification for including FabI at
all — but FabL was never docked, so that argument has no data behind it.

Residue names verified against the entry (192 NDP atoms / 4 chains = 48 each,
68 TCL / 4 = 17 each, 1.80 Å, tetramer A–D). Add to `configs/panel500.yaml`:

```yaml
  - name: fabl_bsubtilis_3oid
    pdb_id: 3OID
    keep_chains: [A]
    keep_hetatms: [NDP] # verified NDP in this entry: 48 atoms/chain
    redock_ligand_resname: TCL # verified: triclosan, 17 atoms
    ph: 7.4
    notes: >
      B. subtilis FabL, the triclosan-resistance isoform. Reversibly
      inhibited, unlike FabI's slow-binding ternary complex, so ordinary
      docking occupancy is the right model. 1.80 A, R-free 0.181 --
      highest-quality structure in the panel. Tetramer A-D; chain A only.
    box:
      reference_ligand_resname: TCL
      size: [22.0, 22.0, 22.0]
```

```bash
python preflight.py --config configs/panel500.yaml   # catches a wrong resname
sbatch run_all.slurm
```

No special-casing needed. `prep-ligands` finds all 492 recorded and returns
instantly, the five existing receptors report "already complete, nothing to
dock", and only FabL's 492 dockings run. Then re-pull and `make` — the report's
title, receptor count and target list are all generated, so nothing goes stale.

### 3. LLE, optional (~15 min of work)

`target_choice.md` §7 asks for triage by LE **and** LLE; only LE exists.
RDKit's Crippen estimator gives logP from the SMILES already in
`scores_by_ligand.csv`, so no new data is needed — but LLE is properly defined
against *measured* potency, and substituting a docking score makes it a
rescaled LE rather than a true LLE. Label it as such if it gets added.

## The report

`report/` builds a 4–5 page PDF: `make venv` once, then `make`. Sources are
tracked (`summary.tex`, `make_report.py`, `Makefile`); the data, figures and
PDF are not.

Nothing in the prose is typed by hand. `make_report.py` reads `analyze.py`'s
output and emits `numbers.tex` (one macro per quantity), four complete
tabulars and three matplotlib figures; `summary.tex` is prose plus `\input`.
**If a number looks wrong, fix the generator, not the `.tex`.** A macro used
but never generated fails the build rather than printing silently.

Per the brief the report states what was measured and makes no recommendation.

## Reconciliation with `target_choice.md` (2026-11-03)

Reading the target rationale against the results changed what the report says:

- **Framing.** The panel was generated target-agnostically under an ADMET
  surrogate reward, with receptors chosen retrospectively. Without that stated,
  "0 of 492" reads as *a designed antibacterial set failed* instead of *a
  target-agnostic set shows no incidental DHFR selectivity*. That paragraph is
  now in the summary and is the single most important thing in the document.
- **Gate redefined.** §8 said targets failing redocking are dropped, not
  carried with a caveat — which the report was doing for 4ZVI and 3GR6. The
  gate is now defined as recovery within the nine retained poses (all five
  pass), with the pose's rank reported as a confidence weight instead.
  `target_choice.md` §3 records this.
- **Corrections pushed back into the document:** 7MYM's cofactor is **NAP, not
  NDP** as written (silent failure — the wrong name drops the cofactor and
  docks into an empty pocket); the panel is `exp_c5_selected_500`, not EXP-C2;
  §7's claimed Vina size bias is **absent here** (+0.04 to +0.14, measured);
  §4.1's GyrB bridging-water pitfall was tested (3U2D redocks at 0.53 Å with
  its four waters excluded).
- **Still outstanding by the document's own criteria:** no actives-vs-decoys
  enrichment test (ROC AUC, EF1%, BEDROC), so all targets remain provisional;
  FabL not docked.

## Summary for the supervisor

Talking points, in the order they hold up to questioning. Detail and numbers
are in the sections below.

- **The screen is done and the setup is validated.** 492 of 500 molecules
  against 5 receptors — 4 bacterial (*E. coli* DHFR, *E. coli* and *S. aureus*
  GyrB, *S. aureus* FabI) plus human DHFR as the anti-target. 2460 dockings,
  zero failures. Redocking recovers the native pose within 2 Å in the top 9 for
  all five receptors; 3 of 5 also rank it first.
- **The headline is a clean negative: 0 of 492 molecules prefer *E. coli* DHFR
  over human DHFR.** Median −1.39 kcal/mol, and the 75th percentile is still
  −1.07 — the whole distribution is on the wrong side. This is the one
  cross-target comparison the data supports (same fold, same cofactor, both
  receptors validated, identical ligands), which is what makes it quotable.
- **The only selectivity signal is FabI, and it is weak.** 9 molecules beat
  human DHFR by ≥0.5 kcal/mol, almost all at FabI — the receptor whose redock
  put the native pose at rank 2 behind a 5.91 Å top pose. Best single candidate
  is `...-00004571` (FabI, −8.68 kcal/mol, LE 0.26), still under a 0.3 ligand
  efficiency floor.
- **It is systematic, not a bad sampling run.** The three generator seeds are
  statistically indistinguishable (medians −0.82 / −0.86 / −0.86, permutation
  p=0.521).
- **The recommendation:** the set was generated `admet_mpo_guided`, with no
  selectivity term, and the result is what that predicts. Human DHFR needs to
  be an explicit anti-target *inside* the generator objective. Post-hoc
  filtering of 492 molecules produced zero survivors and will keep producing
  zero.
- **Cost, for planning:** 134 h on one 32-CPU node at exhaustiveness 32,
  ~16 min per molecule across all five receptors.

Worth saying out loud if asked what could be wrong: the FabI and *E. coli* GyrB
receptors fail the strict top-pose redocking gate, so their absolute scores are
weaker evidence than the three that pass. That does not touch the DHFR result,
which is the one the conclusion rests on.

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

### When does anything need re-running?

The commit is recorded in the job log and in every stage manifest, so this is
a one-command answer rather than a guess:

```bash
git diff --name-only 9770a4e..HEAD -- dockpipe/ run_pipeline.py \
    merge_results.py configs/ environment.yml
```

Non-empty means the pipeline changed and the results no longer match the code.
Empty means they still do. Checked on 2026-10-26: only `analyze.py`,
`check_run.sh` and `handoff.md` have changed since the run, so the 2460 results
stand. Re-run `analyze.py` only when `analyze.py` itself changes.

### Throughput, if you rescreen

Point the config at a new set (`ligands_smi`, `run_name`) and
`sbatch run_all.slurm`. Measured cost on one 32-CPU node at exhaustiveness 32:

| molecules | wall time |
|---|---|
| 500 | 134 h (measured) |
| 1,000 | ~11 days |
| 2,000 | ~23 days, or ~11 sharded across both nodes |

Exhaustiveness 32 does not scale past roughly 1,000 compounds here. Vina's
default is 8; for *ranking* a library that is a defensible 4x saving, with 32
kept for a small confirmation set. Don't raise it — pose ranking is already the
limiting factor on two receptors and more sampling will not fix a scoring-
function preference.

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

## Redocking validation: 5/5 pass the gate, 3/5 rank the pose first

The gate is recovery of the crystal pose within 2 Å among the nine retained
poses — see the reconciliation section above for why it is defined that way.
All five pass. The `verdict` column below is `analyze.py`'s *top-pose* test,
kept because the rank is a useful per-receptor confidence weight, but it is no
longer the pass/fail criterion.

| receptor | top pose | best pose | rank | verdict |
|---|---|---|---|---|
| gyrb_saureus_3u2d | 0.53 | 0.53 | 1 | rank 1 |
| gyrb_ecoli_4zvi | 9.89 | 1.47 | 5 | pose at rank 5 |
| dhfr_ecoli_7mym | 0.89 | 0.89 | 1 | rank 1 |
| dhfr_human_2w3a | 0.40 | 0.40 | 1 | rank 1 |
| fabi_saureus_3gr6 | 5.91 | 1.71 | 2 | pose at rank 2 |

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
- **If FabI becomes the focus, validate that receptor first.** It is where the
  only selectivity signal is, and it is also the receptor whose pose ranking is
  least reliable (native pose at rank 2, 5.91 Å top pose). Building a campaign
  on `...-00004571` without resolving that is building on the weakest part of
  the panel.
- FabL and the enrichment test are in **Next steps**, not here: both are
  preconditions the project set for itself, not optional extras.
- The 8 `scrub_failed` molecules are listed in
  `results/panel500_v1/manifests/ligands.json` under `failures`, if the
  denominator ever needs defending.
