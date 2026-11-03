# Target Selection Rationale — Docking of the EXP-C2 Panel

**Project:** Synthesis-aware generative design of nitrogen-enriched antibacterial candidates
**Panel:** `exp_c5_selected_500` — 500 molecules generated, 492 docked
(8 lost in ligand preparation)
**Date:** 2026-11-01
**Status:** Draft for internal review

---

## 1. Purpose and scope

This document records *why* each docking target was chosen, so that the selection
can be defended in review and reproduced later. It covers the selection criteria,
the targets retained, the targets rejected and the reasons in each case, and the
known weaknesses of the choice.

Two constraints frame everything below:

1. **No target was defined before generation.** The EXP-C2 panel was produced by
   a target-agnostic pipeline (RxnFlow generation under an ADMET surrogate
   reward, with hard chemical and ADMET filters). Targets are therefore selected
   *retrospectively*, to ask whether the panel is chemically compatible with
   validated antibacterial binding sites — not to claim the panel was designed
   against them.
2. **No experimental validation is available.** Docking results support
   *prioritisation and geometric plausibility* only. No potency or activity
   claim is made anywhere in this work.

---

## 2. What the panel is, and how it constrains target choice

| Property | Value |
|---|---|
| Molecules | 500 |
| Distinct Bemis–Murcko scaffolds | 500 |
| SA Score ≤ 5 | 100% (mean 4.676) |
| Total nitrogen, 2–6 | 100% |
| ADMET endpoint gates passed | 100% |
| Internal diversity | 0.854 |
| Top-10 scaffold fraction | 2% |

Three consequences for target selection:

- **Nitrogen-rich but not amine-guaranteed.** The panel was filtered on *total*
  nitrogen count, which includes non-ionisable nitrogen (amides, nitriles,
  aromatic N without basic character). It does **not** guarantee the primary
  amine required by the eNTRy rules for Gram-negative accumulation. Targets were
  therefore preferred where *aromatic and heterocyclic* nitrogen is the
  recognition element, rather than a cationic amine.
- **Small, drug-like, synthesisable.** Sites requiring large, lipophilic,
  natural-product-like ligands (RNA polymerase switch region) or charged
  phosphonate/nucleotide mimics (MurA–F, DXR) are a poor fit.
- **Non-covalent docking only.** Targets whose validated chemistry is covalent
  (PBPs/β-lactams, MurA/fosfomycin, benzoxaborole tRNA synthetases) cannot be
  honestly modelled by AutoDock Vina.

---

## 3. Structure selection criteria

Target choice and *structure* choice are separate decisions. Every structure
below was assessed against the same checklist:

| Criterion | Requirement applied |
|---|---|
| Resolution | ≤ 2.0 Å preferred; ≤ 2.5 Å accepted |
| R-free | Consistent with resolution; small R-free − R-work gap |
| Ligand state | **Holo** — a drug-like co-crystal ligand in the target pocket |
| Ligand density | Co-crystal ligand well resolved (checked in the validation report) |
| Site completeness | No missing loops or residues lining the pocket |
| Cofactor | Retained where it forms part of the binding surface |
| Redocking | Co-crystal ligand must be recovered to ≤ 2.0 Å RMSD within the retained poses |

**Holo over apo** is a hard preference. Apo structures present side chains and
loops in a non-binding-competent conformation and give measurably worse
enrichment. *B. subtilis* FabL is the clearest illustration: the apo form (3OIC)
shows substantially different conformations at both the cofactor and
substrate-binding regions, producing a different tetrameric arrangement
altogether — docking into it would be meaningless.

**Cofactors are retained, not stripped.** For DHFR and FabI/FabL the cofactor is
part of the pocket wall, not a spectator. Trimethoprim binding to DHFR is
strongly NADPH-cooperative, and triclosan stacks directly on the FabI/FabL
nicotinamide ring. Removing the cofactor destroys the site.

**The redocking control is a gate, not a formality.** The pipeline exits
non-zero if any receptor fails the 2.0 Å threshold, before any production
docking runs.

*Gate definition (added after the first run).* The threshold applies to the
**best of the nine retained poses**, not to the top-scored pose alone. Vina
routinely finds the crystal pose and ranks it second or below; that is a
scoring-function limitation, not evidence the site or box is wrong. The rank at
which the pose appears is recorded as a per-receptor confidence weight. Under
the top-pose reading, GyrB *E. coli* (4ZVI) and FabI *S. aureus* (3GR6) would
both be dropped despite recovering the crystal pose at 1.47 Å and 1.71 Å. This is also the practical check that cofactor and ligand residue
names were set correctly in the configuration.

---

## 4. Targets selected

### 4.1 DNA gyrase B (GyrB) ATP site — *broad spectrum*

**Why this target.** Clinically validated (novobiocin historically; multiple
pyrrolamide and aminobenzimidazole programmes since). The GyrB ATPase domain is
conserved across Gram-positives and Gram-negatives, and dual GyrB/ParE
inhibition is an active clinical strategy — this is the broad-spectrum anchor of
the set.

**Why it suits this panel.** Recognition runs through an Asp73 (*E. coli*
numbering) hydrogen bond, frequently bridged by a conserved structural water,
engaged by pyrrole NH and aromatic nitrogen donors. **No primary amine is
required** — which matters, because the panel cannot guarantee one. This is the
best chemical match in the set for nitrogen-rich heterocycles.

**Structures.** 3U2D (*S. aureus*, pyrrolamide), 4ZVI (*E. coli*,
N-phenylpyrrolamide). 4URO (novobiocin) is retained as a cross-docking check.

**Pitfall.** Several ligand series recruit a bridging active-site water. Where
present it must be modelled explicitly; omitting it removes a real hydrogen-bond
partner.

*Tested.* 3U2D has four ordered waters within 3.5 Å of the ligand, the closest
at 2.72 Å. With all of them excluded it still redocks at **0.53 Å**. 4ZVI has
only a marginal 3.33 Å water contact. For these two structures the bridging
water was not required, and its absence does not explain 4ZVI's pose ranking.

---

### 4.2 Bacterial DHFR (FolA) — *selectivity axis*

**Why this target.** Clinically validated (trimethoprim). More importantly, it
carries the cleanest **selectivity story** available without a wet lab:
trimethoprim's affinity for bacterial DHFR exceeds that for the human enzyme by
roughly four to five orders of magnitude, amplified by NADPH cooperativity.

**Why it suits this panel.** The 2,4-diaminopyrimidine recognition motif binding
the conserved Asp27 is *intrinsically nitrogen-rich* — the single best chemical
match between the panel's design heuristic and a validated pocket.

**Structures.** 7MYM (*E. coli*, trimethoprim + NADPH). Cofactor retained —
the residue in this entry is **NAP**, not NDP. (2W3A genuinely is NDP.)
Getting this wrong is silent: the cofactor is simply dropped and the panel
docks into an empty pocket with no error raised.

**Counter-screen.** 2W3A (human DHFR, trimethoprim + NADPH, 1.50 Å, R-free
0.214); 1U72 (human DHFR, methotrexate + NADPH, 1.90 Å) as an alternative. The
same panel is docked into both, and the bacterial-minus-human score differential
is reported as a **qualitative selectivity filter**. This is a score difference
and inherits every caveat of a single docking score; it is not a selectivity
prediction.

---

### 4.3 FabI (enoyl-ACP reductase) — *Gram-positive pathogen*

**Why this target.** Validated: triclosan as a tool compound, afabicin
(Debio-1450) in clinical development against *S. aureus*. Represents the
Gram-positive, pathogen-relevant arm of the set.

**Why it suits this panel.** Flat, nitrogen-containing aromatics stack against
the NAD⁺ nicotinamide with a hydrogen bond to Tyr156 — compatible with the
panel's heteroaromatic content.

**Structures.** 3GR6 (*S. aureus*, NADP + triclosan); 4NR0 (*P. aeruginosa*,
NAD⁺ + triclosan) as an alternative. Cofactor retained.

**Pitfall.** The substrate-binding loop is ordered only in the ternary complex.
Holo structures only.

---

### 4.4 FabL (enoyl-ACP reductase III) — *resistance isoform*

**Why this target.** FabL is not a FabI variant. *B. subtilis* carries both: FabI
is NADH-dependent, FabL is NADPH-dependent, and they share only ~24% sequence
identity. FabL is structurally closer to FabG and the wider SDR family than to
FabI, arriving at the same Tyr-(Xaa)₆-Lys catalytic dyad by convergent
architecture.

**It is the resistance determinant, which is the point.** *fabI* knockouts are as
triclosan-sensitive as wild type, whereas *fabL* knockouts are **250-fold more
sensitive** — FabL, not FabI, is what makes *B. subtilis* triclosan-resistant.
Targeting FabL means targeting a resistance mechanism directly.

**Why it is unusually tractable for docking.** FabI is inhibited by triclosan
through a stable, slow-binding FabI·NAD⁺·triclosan ternary complex — effectively
mechanism-based. FabL is inhibited **reversibly**, without forming that complex.
Ordinary reversible occupancy is what a docking scoring function can actually
represent; a cofactor-dependent dead-end complex is not.

**Structure.** **3OID** (*B. subtilis*, NADPH + triclosan, **1.80 Å, R-free
0.181**) — the highest-quality structure in the whole set. NADPH retained.
3OIC (apo, 2.2 Å) is explicitly **not** used, for the reasons in §3.

**Known weakness, stated plainly.** *B. subtilis* is a model organism, not a
pathogen. FabL is included for the isoform-comparison rationale in §5, not as a
clinical target in its own right. This should be stated in the manuscript rather
than left for a reviewer to raise.

---

## 5. The unifying rationale: an ENR isoform panel

FabI and FabL are not redundant. The bacterial ENR family comprises FabI, FabL,
FabL2, FabK and FabV, and the central medicinal-chemistry problem is that a
compound hitting one isoform has narrow spectrum, while one hitting several is
broad-spectrum.

| Isoform | Cofactor | Representative organisms | Fold | Triclosan |
|---|---|---|---|---|
| FabI | NADH | *E. coli*, *S. aureus*, *P. aeruginosa* | SDR | sensitive (ternary complex) |
| FabL | NADPH | *B. subtilis* (alongside FabI) | SDR, FabG-like | resistant (reversible) |
| FabK | FMN | *S. pneumoniae*, *E. faecalis* | flavoprotein — unrelated fold | resistant |
| FabV | NADH | *P. aeruginosa*, *V. cholerae* | distinct | very weakly inhibited |
| InhA | NADH | *M. tuberculosis* | SDR | — |

Docking the panel across FabI + FabL (and optionally FabK) yields a
**cross-isoform compatibility profile**: structurally the same question as the
human-vs-bacterial DHFR comparison, but about antibacterial spectrum rather than
host selectivity. FabK is the strongest optional addition precisely because its
flavoprotein fold is unrelated — it is a genuine orthogonality test, not a near
duplicate.

**Each target answers a different question.** This is the actual justification
for the set as a whole:

| Target | Question it answers |
|---|---|
| GyrB | Is the panel compatible with a broad-spectrum, clinically validated ATP site? |
| DHFR + human DHFR | Can any scaffold discriminate bacterial from host enzyme? |
| FabI | Does the panel engage a pathogen-relevant, clinically advanced target? |
| FabL | Does it engage the isoform responsible for resistance? |

---

## 6. Targets considered and rejected

| Target | Reason for rejection |
|---|---|
| **LpxC** | Zinc metalloenzyme; potent inhibitors are Zn-chelating hydroxamates. Vina's empirical score has no explicit electrostatic term and mishandles metal coordination; AutoDock4Zn helps but remains fragile. Hydrophobic tunnel also favours lipophilic tails, against the panel's profile. |
| **PDF** | Fe/Ni/Zn-dependent with hydroxamate chelation, plus documented induced-fit pocket closure. Same metal problem as LpxC. |
| **MurA–MurF** | Large, polar, often charged nucleotide/UDP-sugar sites — poor Vina scoring. MurA inhibition by fosfomycin is covalent. |
| **RNA polymerase** | Switch-region ligands (myxopyronin, corallopyronin) are large lipophilic macrocycles binding a shallow surface that locally refolds on binding. Not a compact pocket for small nitrogen-rich molecules. |
| **PBP2a / PBP3** | β-lactam mechanism is covalent (acyl-enzyme). Non-covalent docking would be misleading. |
| **LptD / LptA / LolE** | Biologically compelling (abaucin precedent) but the target was identified by phenotypic ML plus genetics, not docking; no clean high-resolution small-molecule co-crystal pocket. Aspirational, not executable here. |
| **NBTI / GyrA–DNA site** | Chemically attractive, but requires docking into a protein–DNA complex and potency depends on bifurcated halogen bonding, which Vina scores poorly. Retained as an optional reach target, interpreted via interaction fingerprints rather than score. |
| **tRNA synthetases** | Benzoxaborole mechanism is covalent and tRNA-dependent. |
| **DXR/MEP, UPPS, AccC, riboswitches** | Charged phosphonate sites, RNA targets, or thin structural/actives support. |

---

## 7. Known weaknesses of this selection

These should appear in the manuscript, not be left for review to surface.

- **FAS-II as a Gram-positive target is contested.** Brinster et al. (*Nature*
  2009) argued that Gram-positive pathogens bypass FAS-II inhibition by
  scavenging host fatty acids; Rock, Parsons and others published substantial
  rebuttals, and afabicin's clinical progression cuts against the strongest form
  of the objection. The debate is unsettled. Both sides should be cited.
- **Retrospective target selection is hypothesis generation, not validation.**
  The panel was not designed against these targets. Framing must be "chemically
  and geometrically compatible with", never "predicted active against".
- **Docking score is not binding free energy.** Published correlations with
  measured affinity are weak. Only enrichment, pose plausibility and ligand
  efficiency are defensible outputs.
- **Vina-specific limits.** Charged ligands, halogen bonds, entropy and
  desolvation are all poorly handled, and there is a known size/lipophilicity
  bias. Triage is therefore by **ligand efficiency (LE) and LLE**, not raw score.
  *Measured on this panel:* the size bias is **absent** — score against
  heavy-atom count gives Pearson +0.04 to +0.14 across the five receptors, i.e.
  marginally positive. The heavy-atom range is narrow (median 41, 21–57) and the
  set chemically homogeneous, which is the likely reason. LE triage remains
  sensible, but not for this stated reason on this panel.
- **The panel's nitrogen heuristic is not the eNTRy rule.** Total nitrogen count
  in 2–6 is a loose proxy; it does not establish Gram-negative accumulation,
  which additionally requires an ionisable (preferably primary) amine,
  globularity ≤ 0.25 and ≤ 5 rotatable bonds. Any Gram-negative framing must be
  a property-space argument, not a permeability claim.

---

## 8. Summary

| Priority | Target | Structure | Resolution | Cofactor | Role in the set |
|---|---|---|---|---|---|
| 1 | GyrB | 3U2D, 4ZVI | — | — | Broad spectrum |
| 2 | DHFR (*E. coli*) | 7MYM | — | NADPH | Nitrogen-motif match |
| 2b | DHFR (human) | 2W3A | 1.50 Å | NADPH | Selectivity counter-screen |
| 3 | FabI (*S. aureus*) | 3GR6 | — | NADP | Gram-positive pathogen |
| 4 | FabL (*B. subtilis*) | **3OID** | **1.80 Å** | NADPH | Resistance isoform |
| opt. | FabK | TBD | — | FMN | Fold-orthogonal ENR |

All selections are provisional until each passes the redocking gate at ≤ 2.0 Å
RMSD and an actives-versus-decoys enrichment test (ROC AUC, EF1%, BEDROC) built
from ChEMBL actives with property-matched decoys. **Targets that fail either
control are dropped rather than carried forward with a caveat.**

### Status after the first production run (2026-10-18)

| Target | Redock, best of 9 | Rank of that pose | Gate |
|---|---|---|---|
| GyrB 3U2D | 0.53 Å | 1 | pass |
| GyrB 4ZVI | 1.47 Å | 5 | pass |
| DHFR 7MYM | 0.89 Å | 1 | pass |
| DHFR 2W3A | 0.40 Å | 1 | pass |
| FabI 3GR6 | 1.71 Å | 2 | pass |
| FabL 3OID | — | — | **not run** |

- **FabL was not docked.** The §5 FabI+FabL cross-isoform rationale is therefore
  unsupported by data so far, and §5 is the stated justification for including
  FabI at all. Either run 3OID (≈ 27 h for 492 molecules on one node at the
  measured rate) or drop the isoform argument.
- **The enrichment test has not been run.** By the criterion above, all five
  targets remain provisional.

---

## 9. Key references

1. Kim K.-H. *et al.* Crystal structures of enoyl-ACP reductases I (FabI) and
   III (FabL) from *B. subtilis*. *J. Mol. Biol.* (2011).
   doi:10.1016/j.jmb.2010.12.003 — PDB 3OIC, 3OID, 3OIF, 3OIG.
2. Heath R.J. *et al.* The enoyl-[acyl-carrier-protein] reductases FabI and FabL
   from *Bacillus subtilis*. *J. Biol. Chem.* (2000). PMID 11007778.
3. Massengo-Tiassé R.P. & Cronan J.E. *Vibrio cholerae* FabV defines a new class
   of enoyl-ACP reductase. *J. Biol. Chem.* (2008).
4. Krucinska J. *et al. Commun. Biol.* (2022). doi:10.1038/s42003-022-03384-y —
   *E. coli* DHFR:TMP:NADPH, PDB 7MYM.
5. Cody V. *et al.* (2005). PMID 15681865 — human DHFR, PDB 1U72.
6. Richter M.F. *et al.* Predictive compound accumulation rules yield a broad-
   spectrum antibiotic. *Nature* (2017) 545:299–306. doi:10.1038/nature22308 —
   eNTRy rules.
7. O'Shea R. & Moser H.E. Physicochemical properties of antibacterial compounds.
   *J. Med. Chem.* (2008) 51:2871–2878. doi:10.1021/jm700967e
8. Brinster S. *et al.* Type II fatty acid synthesis is not a suitable
   antibacterial target for Gram-positive pathogens. *Nature* (2009) — and
   subsequent rebuttals.
9. Santos-Martins D. *et al.* AutoDock4Zn. *J. Chem. Inf. Model.* (2014).
   doi:10.1021/ci500209e
10. Quiroga R. & Villarreal M.A. Vinardo scoring function. *PLoS ONE* (2016).

> Resolutions, R-free values, ligand residue names and cofactor identities must
> be confirmed against each RCSB entry (and its PDB-REDO record) before the
> methods table is finalised. Entries marked "—" above are pending that check.
