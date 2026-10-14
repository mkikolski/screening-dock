# Handoff — 2026-10-14

State: validated and ready to submit. Nothing is pushed yet.

## Do this first

```bash
# 1. local
git push

# 2. cluster -- ONE-TIME: create the env in the shared miniforge install.
#    $HOME does not travel to the compute node, so an env under ~/.conda is
#    unusable from a batch job even though it works interactively.
cd /storage/users/kikolski/screening-dock && git pull
conda env create -p /opt/modules/miniforge3/envs/dockpipe -f environment.yml

# 3. submit
sbatch run_all.slurm
squeue -u $USER
```

`run_all.slurm` runs all five stages itself: prep-receptors, redock
(`--allow-fail` is already in the file), prep-ligands, dock, merge. No manual
edits are needed any more.

Both long stages checkpoint, so **if anything dies, resubmit the same script**
and it resumes from the last finished item:

- `results/panel500_v1/ligands/.ckpt_ligands.jsonl`
- `results/panel500_v1/docking/.ckpt_<receptor>.jsonl`

`scores.csv` is rewritten after every receptor and `ligand_index.json` after
the ligand stage, so a killed run still leaves usable output for `merge`.

Runtime at 32 workers, from the measured ~30 min per docking:

| ligand states | dockings | wall |
|---|---|---|
| 500 | 2,500 | 39 h |
| 750 | 3,750 | 59 h |
| 1,000 | 5,000 | 78 h |

`--time` is 168 h. The progress line prints live min/docking and hours
remaining every 25 dockings -- check it in the first hour, since the 30 min
figure came from `head -10` of the SMILES file and that sample may be biased.

If the state count is near 1,000 and node 2 is free, `--shard-count 2` across
a 2-task array halves it; `merge_results.py` already globs `scores*.csv` and
dedups on `(ligand, receptor)`.

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
  0.52 Å with its four proper 2.72 Å waters excluded. So a missing bridging
  water does not explain it either.

Since a compound's score is the best of its 9 poses, this matters little for
ranking compounds *within* a receptor. Treat absolute cross-receptor
comparisons involving 4ZVI and 3GR6 as weaker evidence. Per-receptor notes are
in `configs/panel500.yaml`; the table is also in `README.md`.

**Don't tighten the boxes** to "fix" these. The control must use the same box
as the screen, and the panel is median 41 heavy atoms against triclosan's 17.
22³ = 10,648 Å³, well under Vina's ~27,000 Å³ degradation point.

## What was fixed this session

| commit | |
|---|---|
| `b808955` | cofactors prepared separately — meeko has no NAP/NDP template |
| `857492d` | leftover `het` in the preflight print |
| `7d06be9` | retry receptor prep without rebuilt heavy atoms |
| `2966b5f` | redock control: real bond orders, all 9 poses |
| `ad6cfe8` | recorded the redock verdicts |
| `13f143e` | dock stage checkpointing |

Two of these were silent-wrong-answer bugs worth remembering:

1. **The redock control was measuring the wrong molecule.** A PDB entry stores
   no connectivity for ligands, and `ligand_lines_to_pdbqt` was always called
   with `template_smiles=None`, so RDKit perceived bonds by proximity and made
   every one single; `SanitizeMol` + `AddHs` then saturated the aromatics. Each
   reference ligand went into Vina with 6–8 bonds wrong. Bond orders now come
   from the chemical component's ideal SDF (fetched and cached per ligand).
   3U2D went 2.46 → 0.52 Å on this alone. The two DHFR "passes" before the fix
   were luck.
2. **Only pose 1 was being scored** — `_heavy_atom_coords_pdbqt` broke at the
   first `ENDMDL`, so a native pose at rank 2–4 read as a failure.

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
- After the run: `seconds` is a column in `scores.csv`, so the real
  per-docking distribution is available for writing up.
