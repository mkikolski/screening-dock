# Handoff — 2026-11-17

**State: three follow-up runs are coded and committed, nothing submitted yet.**
The original screen (job 1994, commit `9770a4e`, finished 2026-10-18 22:44 on
`serwer1`) is complete and written up in `report/summary.pdf` from real data.
The runs below exist because review found that report's headline confounded.

**Pick up at step 0: `git push`.** 15 commits are local-only and the cluster
pulls from the remote, so nothing run there has the new code until you do.

## Why these runs exist

The report said 0 of 492 molecules beat human DHFR. The number is real but its
reading was never established, because the DHFR pair was not the controlled
comparison the report claimed:

- **7MYM** is 3.04 Å and holds **NAP** (NADP⁺, oxidised).
- **2W3A** is 1.50 Å and holds **NDP** (NADPH, reduced).

The pair differed in resolution *and* in cofactor oxidation state at once.
Three things say the gap is a receptor artefact rather than selectivity:

1. Four receptors — two folds, three proteins, both box sizes, cofactor and
   none — cluster between −6.27 and −6.54 kcal/mol. **2W3A alone is −7.69.**
2. The gap is **flat across ligand size**: +1.32 / +1.39 / +1.41 across
   heavy-atom terciles, r(gap, heavy atoms) = 0.05. A pocket-fit effect would
   grow with size. A constant offset is a baseline.
3. **Trimethoprim runs the other way**: −7.638 at 7MYM against −7.423 at 2W3A.
   Of 493 ligands docked into both, the single exception is the only one with
   crystal evidence in both structures.

Resolution also rank-correlates with mean score across the five receptors
(ρ 0.70, or 0.40 with 2W3A dropped), pointing the same way.

6XG5 separates the offset from the chemistry. 7MYM stays in the panel so the
offset is measured rather than assumed away.

| | 7MYM (keep) | **6XG5 (new)** | 2W3A (human) |
|---|---|---|---|
| resolution | 3.04 Å | **1.90 Å** | 1.50 Å |
| cofactor | NAP (NADP⁺) | **NDP (NADPH)** | NDP (NADPH) |
| ligand | TOP | TOP | TOP |
| chains | A, B, C | **A only** | A, B |

6XG5 is *E. coli* K-12 *folA*, 165 residues (160–165 an expression tag), also
holding GOL and CL, both dropped. Manna et al., *Nat Commun* **12**:2949
(2021), doi:10.1038/s41467-021-23191-z.

## Runbook

### 0. Push, then pull on the cluster

```bash
cd ~/mateusz-dock && git push -u origin master     # master has no upstream yet

# on the cluster
cd /storage/users/kikolski/screening-dock && git pull
conda activate /storage/users/kikolski/envs/dockpipe
```

### 1. Check the molscrub signature (1 min, blocks step 2)

molscrub 0.8.0 documents `num_etkdg_attempts` and `use_random_coords` as CLI
flags but does not publish its Python signature. `_init_worker` raises rather
than dropping them silently, so confirm the names first:

```bash
python -c "import inspect, molscrub; print(inspect.signature(molscrub.Scrub.__init__))"
```

Both present → go. Different → edit the two keys in `dockpipe/ligands.py:52-56`;
the `RuntimeError` lists what is actually accepted.

### 2. Panel run: 6XG5 plus the 8 recovered molecules (~27 h)

```bash
R=/storage/users/kikolski/screening-dock/results/panel500_v1
cp $R/ligands/.ckpt_ligands.jsonl $R/ligands/.ckpt_ligands.jsonl.bak
cp $R/manifests/receptors.json   $R/manifests/receptors.json.bak

cd /storage/users/kikolski/screening-dock
DOCKPIPE_RETRY_FAILED=1 sbatch run_all.slurm
```

One pass does both: the 8 are re-prepped first, then 6XG5 docks the full set.
The five existing receptors report complete and are skipped.

**Watch for this line early in the log.** It is the difference between
recovering the 8 and silently not:

```
retry_failed: 8 non-ok records in .../.ckpt_ligands.jsonl will be redone
```

Missing → the env var did not reach the script and the losses stay lost.

### 3. GyrB reference set (~6 h)

```bash
python fetch_gyrb_actives.py --selftest      # do not skip
python fetch_gyrb_actives.py --out-dir /storage/users/kikolski/screening-dock
# then READ gyrb_actives.csv before docking
DOCKPIPE_CONFIG=$PWD/configs/gyrb_knowns.yaml sbatch run_all.slurm
```

**ChEMBL is unreliable and the script is built around that.** Measured
2026-11-17: it 500s in bursts lasting a minute or two, fast-failing in ~5s,
then recovers and answers in ~0.3s. Its own `/status.json` 500s too, so there
is nothing to health-check. Consequences baked in:

- `limit=1000` (the documented maximum) times out; 100 and 50 both 500 even
  when healthy. Page size is **25**.
- 8 attempts with exponential backoff, ~4.5 min of patience per page.
- Pages cache to `<out-dir>/.chembl_cache`. **A rerun resumes** — a full fetch
  is 26 pages and replays from cache in 0.00s. Delete the directory to force
  fresh data.
- A dead (target, potency type) pair is logged and skipped, not fatal. The
  first version raised `SystemExit` and discarded the IC50 pages it already
  had when Ki started 500ing.

Verified against the live API: **540 potency records, 456 distinct compounds,
0.05 nM to 1.08 mM**, fetched in 10 s once the service was responding.

**Targets are pinned ChEMBL IDs, not UniProt lookups.** `CHEMBL1826`
(*E. coli*) and `CHEMBL1921666` (*S. aureus*). P66936 — the accession
BindingDB shows for *S. aureus* GyrB — maps to **zero** ChEMBL targets, so the
original accession lookup returned an empty *S. aureus* half while appearing
to work. Re-derive with
`target.json?pref_name__icontains=gyrase%20subunit%20B`.

**How inclusion is decided.** The assay is the evidence, not the scaffold: 540
of 540 records carry descriptions like "Inhibition of ATPase activity of
*S. aureus* GyrB", which measures the ATP site directly — a fluoroquinolone
cannot inhibit GyrB ATPase. So:

1. An off-site scaffold (fluoroquinolone, NBTI) excludes outright, whatever
   the assay says.
2. A known ATP-site scaffold includes.
3. Otherwise the assay decides.

`--strict-chemotype` requires an allowlist match instead, which drops most
modern GyrB series (pyrimidinoindoles, thiazolopyridines) and is usually the
wrong choice. Every row's `mechanism_evidence` column records which rule
applied, and `gyrb_excluded.csv` holds the rejects with reasons — a good
compound in there means the filter is too strict.

**`--max-compounds` defaults to 60**, the most potent after filtering. Uncapped
is ~450 compounds, which is ~900 dockings and **~50 h on two receptors** —
more than the panel run it exists to calibrate.

The selftest asserts the parts that fail silently: that a fluoroquinolone stays
excluded even under an ATPase assay and even when it also carries a urea, that
assay evidence admits an unlisted scaffold, that a whole-cell MIC assay does
not, and that `--strict-chemotype` overrides all of it.

Separate `run_name` (`gyrb_knowns_v1`), so it can run alongside step 2 if the
partition allows two jobs.

**State the free result first when reporting this.** The co-crystal ligands
redock at −8.004 (08B, 3U2D) and −7.970 (4S4, 4ZVI). The panel's best GyrB
scores are −8.75 and −8.33. The generated molecules already out-score the real
inhibitors, which is a calibration warning rather than a hit.

## When the results land

### Pull them

`results/` and `report/data/` are both gitignored, so result data never travels
through the repo. The login host is recorded nowhere here — fill it in below
once known.

```bash
cd ~/mateusz-dock && mkdir -p report/data
R=<user>@<login-host>:/storage/users/kikolski/screening-dock/results/panel500_v1
rsync -av "$R"/docking/{analysis/per_molecule.csv,scores_wide.csv,redock_control.json,summary.json} \
          "$R"/manifests/ligands.json  report/data/
cd report && make venv && make      # make venv once: system python3 lacks pandas
```

Files land flat in `report/data/`, which is what `make_report.py` expects
(`data/per_molecule.csv`, not `data/analysis/...`).

### Checks, in order of what they would invalidate

```bash
# 1. The 5 existing receptors must be unchanged. Receptor prep is not
#    checkpointed, so all six are re-prepared; prep is deterministic, so
#    verify rather than assume. A change here invalidates every comparison
#    with the first run.
python - <<'PY'
import json
old = {r['name']: r.get('checksums') for r in json.load(open('manifests/receptors.json.bak'))}
new = {r['name']: r.get('checksums') for r in json.load(open('manifests/receptors.json'))}
for k in old: print(k, 'UNCHANGED' if old[k] == new[k] else '*** CHANGED ***')
PY

# 2. 6XG5 kept the right cofactor. This trap already bit once: the config
#    named NDP where the entry held NAP, and the splitter dropped it silently.
grep NDP receptors/dhfr_ecoli_6xg5/dhfr_ecoli_6xg5_cofactor.pdbqt | wc -l   # expect 48
grep -c 'GOL\|CL ' receptors/dhfr_ecoli_6xg5/dhfr_ecoli_6xg5_clean.pdb      # expect 0

# 3. 6XG5 redocks its own ligand, and the other five still do
python -c "import json;[print(r['receptor'],round(r['rmsd_best'],2),r['best_pose_rank']) for r in json.load(open('docking/redock_control.json'))]"

# 4. The denominator moved: 500 of 500 if all 8 recovered, 492+n otherwise
python -c "import json;d=json.load(open('manifests/ligands.json'));print(d['n_ok'],'of',d['n_input'])"

# 5. dhfr_selectivity means what it says -- both columns must be present
head -1 docking/scores_wide.csv
```

### The number this was all for

Median `selectivity_dhfr_ecoli_6xg5` against median
`selectivity_dhfr_ecoli_7mym` in `scores_wide.csv`. Positive favours the
bacterial enzyme; the original run gave −1.39 for 7MYM.

| 6XG5 median | 7MYM median | reading |
|---|---|---|
| ≈ 0 | ≈ −1.39 | **The headline was receptor preparation.** The 0/492 was largely 2W3A's baseline, not selectivity. The report's conclusion needs rewriting, not updating. |
| ≈ −1.39 | ≈ −1.39 | **The negative result survives a matched structure**, and is much stronger than before: same cofactor, 0.40 Å resolution gap, offset still there. |
| in between | ≈ −1.39 | Part artefact, part real. Report both and quantify the split rather than picking one. |
| anything at all, but 6XG5 redocks badly | — | Discount it. A receptor that cannot recover its own crystal pose cannot settle this. Run check 3 before reading the medians. |

Recompute these too, because a fifth bacterial receptor changes them for
structural rather than chemical reasons: `best_bacterial`, `margin_vs_human`,
and the "36 of 492 beat human at some target" figure. Report the old and new
panels side by side instead of silently replacing.

### Report

Agreed deliverable: regenerate `summary.pdf` with 6XG5 as the primary axis and
7MYM retained as the offset control, **plus a one-page note listing which
numbers moved and why**. The cofactor corrections to Table 1 and Methods are
already committed — put them in that note rather than letting them pass
silently, since they propped up the "same fold with the same cofactor" claim.

## Still optional: LLE (~15 min of work, no docking)

`target_choice.md` §7 asks for triage by LE **and** LLE; only LE exists.
RDKit's Crippen estimator gives logP from the SMILES already in
`scores_by_ligand.csv`, so no new data is needed — but LLE is properly defined
against *measured* potency, and substituting a docking score makes it a
rescaled LE rather than a true LLE. Label it as such if it gets added.

## FabL: dropped, deliberately (2026-11-03)

FabL (3OID) was planned as priority 4 and is **not being run**. It is written
up in `report/summary.tex` §"A possible extension: FabL" with the reasons it
would be worth adding, and that is where it stays.

Why it is still worth describing rather than deleting: FabL is the triclosan
*resistance* determinant in *B. subtilis* (`fabL` knockouts ~250-fold more
sensitive, `fabI` knockouts not), it is inhibited reversibly rather than
through FabI's slow-binding ternary complex, so plain docking occupancy is
actually the right model for it, and FabI-against-FabL would be the same
cross-isoform question that the DHFR pair asks about host selectivity. 3OID at
1.80 A is also the best structure in the candidate set.

If it is ever run, the verified config block is in git history at `521b97d`
(`handoff.md` step 2): chain A, `keep_hetatms: [NDP]` (48 atoms/chain), redock
on TCL (17 atoms/chain), 22 A box. Checkpointing means it is a plain
`sbatch run_all.slurm` — the five existing receptors report "already complete"
and only FabL's 492 dockings run, about 27 h.

## The report

`report/` builds a 7-page PDF: `make venv` once, then `make`. Sources are
tracked (`summary.tex`, `make_report.py`, `Makefile`); the data, figures and
PDF are not. `report/notes.md` holds the per-section speaker notes — the
detail deliberately kept out of the PDF.

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
  §7's claimed Vina size bias is **absent here** (+0.04 to +0.13, measured);
  §4.1's GyrB bridging-water pitfall was tested (3U2D redocks at 0.53 Å with
  its four waters excluded).
- **Still outstanding by the document's own criteria:** no actives-vs-decoys
  enrichment test (ROC AUC, EF1%, BEDROC), so all targets remain provisional.
  FabL is no longer outstanding; it is dropped (see above).
- **7MYM is 3.04 A**, failing §3's ≤ 2.5 A criterion — §8's resolution column
  was "—" pending the check, and the check was never done before the run. It
  redocks at 0.89 A ranked first, so the site is fine, but swap it for a
  higher-resolution E. coli DHFR:TMP:NADPH entry on any re-run.

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
+0.04 to +0.13 across the five receptors, i.e. slightly *positive* — bigger
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
- The enrichment test (actives vs decoys) is a precondition the project set
  for itself and has not been run, so every target is formally provisional.
  The GyrB reference run above is **not** that test: it is a positive-control
  comparison with no decoys, so it gives no ROC AUC, EF1% or BEDROC. The full
  test needs property-matched decoys and roughly 1–10 days of compute
  depending on set size, and was deliberately scoped separately.
  FabL is dropped, not pending — if asked, the answer is that the isoform
  argument is described in the report as an extension and FabI stands on
  triclosan and afabicin instead.
- The 8 `scrub_failed` molecules are being recovered by the step 2 run above.
  Root cause: `n_scrub_attempts: 5` retried an identical deterministic `Scrub`
  call five times while molscrub's own `num_etkdg_attempts` and
  `use_random_coords` were never set — its error message named both knobs.
  All eight have a quaternary carbon bonded directly to a trisubstituted
  nitrogen, which over-determines the distance bounds around that bond; two
  also carry strained fused bicyclics. They are 38–48 heavy atoms against a
  panel median of 41, so this was never a size problem. If some still fail
  after the retry, record that: it bounds what the generator can produce that
  is embeddable at all. The originals remain in
  `results/panel500_v1/manifests/ligands.json` under `failures`.
- **If a rerun ever looks like a no-op, check this first.** `read_checkpoint`
  treats a recorded failure as done unless `retry_failed=True`, on the stated
  reasoning that failures are deterministic. That reasoning expired when the
  `Scrub` call changed. `DOCKPIPE_RETRY_FAILED=1` is the opt-in, and the
  `retry_failed:` log line is the proof it took effect.
- **`primary_bacterial_dhfr` is now mandatory with two bacterial DHFRs.**
  `merge_results.py` used to pick that column positionally, so
  `dhfr_ecoli_6xg5` sorting before `dhfr_ecoli_7mym` would have repointed
  `dhfr_selectivity` with no error and no warning. It now raises instead.
  Adding a third bacterial DHFR means setting this key, not trusting the sort.
