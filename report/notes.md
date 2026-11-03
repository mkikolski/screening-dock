# Meeting notes — what is *not* in summary.pdf

Keyed to the report's sections. The PDF states the facts and shows the plots.
These are the things to say out loud, including the ones that only matter if
someone asks.

For the ordered argument rather than the per-section detail, see
`handoff.md` § "Summary for the supervisor".

---

## Summary

- 500 generated, **492 docked**. The 8 losses are molscrub failures, listed in
  `results/panel500_v1/manifests/ligands.json` under `failures` if the
  denominator ever needs defending.
- 2460 dockings, **zero failures**. Job 1994, 134 h wall on one 32-core node,
  mean 104.8 min per docking.
- The headline is a negative result, and a clean one: the distribution does not
  straddle zero, it sits entirely on one side of it.
- **Left out deliberately:** any recommendation. The obvious implication —
  selectivity has to be in the generator's objective rather than discovered
  afterwards — is for the meeting, not the document.

## Targets

- **GyrB and FabI were your suggestions. The DHFR pair is the one I added**,
  and it is the only comparison in the panel that controls for anything.
- Why it controls: same ligands, same fold, same cofactor, two receptors
  differing in the pocket. A FabI-vs-human-DHFR score difference confounds the
  comparison with Vina's calibration across unrelated sites. A DHFR-vs-DHFR
  difference does not.
- **The NADPH cooperativity, if it comes up.** Trimethoprim binds E·NADPH more
  tightly than free E — about −2.9 kcal/mol of coupling in *L. casei*. Two
  causes: the diaminopyrimidine contacts the nicotinamide ring directly, and
  NADPH orders the Met20 loop onto the inhibitor. The NMR work found it
  entropy-driven with enthalpy opposing, and found that interligand contact
  *cannot* explain the human/bacterial difference — so the selectivity comes
  from loop dynamics. Krucinska 2022 (our 7MYM citation) is the reverse
  demonstration: plasmid DfrA resistance enzymes lose the cooperativity via
  proline substitutions that restrict Met20. **A rigid receptor cannot model
  any of this.** We capture the geometry the cofactor creates, not the
  coupling. Retaining NADPH is still right, but because it is a pocket wall,
  not because we reproduce the effect.
- **7MYM is 3.04 Å**, which fails the ≤ 2.5 Å criterion in `target_choice.md`
  §3. That column was left as "—" pending a check that never happened before
  the run. It redocks at 0.89 Å ranked first, so it does not appear to have
  cost anything, but it is the first thing to swap.
- **2W3A's RCSB citation record is misleading** — it carries an unrelated
  *M. avium* DHFR title and "to be published". The entry itself is verified
  human DHFR, 187 residues, NADPH + trimethoprim, 1.50 Å. Cited as a
  deposition for that reason.
- **The NAP/NDP bug, if the methods are questioned.** The config named NDP as
  the cofactor for 7MYM and 3GR6, which actually hold NAP. The splitter dropped
  it **silently** — no error — so the panel would have docked into a
  cofactor-less pocket and the whole selectivity axis would have been
  meaningless. Caught before production.
- The bridging water was tested, not assumed: 3U2D has four ordered waters
  within 3.5 Å of the ligand, closest at 2.72 Å, and redocks at 0.53 Å with all
  of them removed.
- Any row of the excluded-targets table can be expanded on. The three method
  exclusions (covalent, metal, charged nucleotide) are the defensible ones.
- **FabL:** ~27 h of wall time and nothing else. The config block is verified
  (chain A, NDP 48 atoms/chain, redock on TCL 17 atoms/chain, 22 Å box) and
  checkpointing means it is a plain `sbatch` — the five existing receptors
  report complete and only FabL's 492 dockings run.

## Result — selectivity

- **0 of 492.** Median −1.39 kcal/mol in favour of the human enzyme, upper
  quartile −1.07. Nothing in the upper tail crosses zero.
- Widening to any bacterial target: 36/492 (7%) beat human DHFR, 9 by
  ≥ 0.5 kcal/mol, 2 by ≥ 1.0.
- **7 of those 9 are FabI, 2 are GyrB *S. aureus*.** The only positive signal
  sits on the receptor with the weakest validation and the contested biology.
  Worth saying before anyone picks a molecule off that table.

## Result — ligand efficiency

- Best LE in the panel 0.352, at **human** DHFR. Best bacterial 0.332. Both
  from the smallest molecule in the set (21 heavy atoms).
- Nothing that prefers a bacterial target clears LE 0.30.
- **The size bias was measured and is absent**: Pearson +0.04 to +0.13 across
  the five receptors, i.e. marginally *positive*. `target_choice.md` §7
  asserted the bias as the reason to triage by LE. That reasoning does not
  apply to this panel. Likely because the heavy-atom range is narrow — median
  41, 21 to 57.
- **LLE was not computed.** Substituting a docking score for measured potency
  makes it a rescaled LE, so it should not carry the LLE label. Its one real
  use would be testing whether the FabI margins are a lipophilicity artefact —
  FabI's pocket is greasy and triclosan is a chlorinated diphenyl ether, so
  "the 36 are just the greasiest molecules in the panel" is a live alternative
  explanation. About 15 min of work, no new docking.

## Result — generator seeds

- Three seeds, n = 164 / 134 / 194. Medians −0.82 / −0.86 / −0.86.
- Permutation test on the widest pair: p = 0.522 over 10,000 shuffles. The
  seeds are indistinguishable on this measure, so the result is systematic
  rather than one bad sampling run.
- **If pressed on the p-value:** it is not reproducible to three decimals
  without care. |obs| is 0.041 and about 40 of 10,000 shuffled gaps land on it
  exactly, because Vina scores are discretised to 3 dp. Ties now count toward
  the null on both code paths. This moved the figure from 0.519/0.521 to 0.522
  and changes nothing at p ≈ 0.5.

## Methods and validation

- **The gate was redefined after the first run** to best-of-nine rather than
  top-scored pose. Under the top-pose reading, 4ZVI and 3GR6 would have been
  dropped despite recovering the crystal pose at 1.47 Å and 1.71 Å. Rank is now
  recorded as a per-receptor confidence weight. Be ready to defend this as a
  scoring-function limitation rather than a moved goalpost.
- **The redock control initially measured the wrong molecule.** With no
  template SMILES, bond orders were perceived from coordinates and every
  reference ligand had 6–8 bonds wrong. Fixed by pulling each component's
  idealised SDF from RCSB. 3U2D went 2.46 Å → 0.53 Å. The three "failures"
  before that fix were an artefact of the control, not of the docking.
- Provenance: every stage manifest carries the pipeline commit, package
  versions, resolved tool paths and SHA-256 checksums of the inputs.
- Both long stages checkpoint to append-only JSONL, so a killed job resumes
  instead of restarting. That is what makes FabL cheap.

## Scope

- **The enrichment test is the real gap.** `target_choice.md` §8 says every
  target is provisional until it passes actives-vs-decoys (ROC AUC, EF1%,
  BEDROC), and none has. Expect this question — the criterion is the project's
  own. It needs ChEMBL actives per target plus property-matched decoys and
  another docking run, so it is real work, not an afternoon.
- FAS-II in Gram-positives: Brinster 2009 against Parsons & Rock 2011.
  Afabicin reaching phase 2 cuts against the strongest form of the objection,
  but the debate is open and it lands on FabI.

## If asked

- **Why only 500 molecules?** 104.8 min per docking × 492 × 5 = 134 h. A
  10,000-molecule panel is about 46 days on one node, so it needs sharding
  across a job array rather than a longer wall clock.
- **Why no `module load vina`?** The pipeline uses the Vina *Python* API, which
  a module providing only the CLI cannot satisfy.
- **Cluster constraint worth knowing:** compute nodes do not mount `$HOME`, so
  the conda environment has to live under `/storage`.
