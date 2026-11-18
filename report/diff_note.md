# What changed in `summary.pdf`, and why

Companion to the regenerated report. Previous version: commit `9770a4e`,
SLURM job 1994, 492 molecules × 5 receptors.

**Basis for the numbers below.** Interim read of the running job at 467 of
498 molecules scored at 6XG5. The medians are stable at that sample size;
every figure here is regenerated from the completed run before circulation.

---

## 1. A sixth receptor, and why the old pair was not a controlled comparison

The previous report compared 7MYM (3.04 Å, NADP⁺) against human 2W3A
(1.50 Å, NADPH) and described it as "the same fold with the same cofactor".
Those two structures differ in resolution **and** in cofactor oxidation state
at the same time, so the comparison was not controlled and that sentence was
not true.

6XG5 was added: *E. coli* K-12 DHFR (*folA*) with trimethoprim and NADPH at
1.90 Å, single chain A (Manna et al., *Nat Commun* 12:2949, 2021). Against
2W3A it matches the cofactor exactly and cuts the resolution gap from 1.54 Å
to 0.40 Å. 7MYM is deliberately retained as the control on the offset itself.

## 2. The offset test — the substantive new result

Measured over the molecules scored at both bacterial structures:

| | 7MYM (mismatched) | 6XG5 (matched) |
|---|---|---|
| median selectivity vs human | −1.388 | −1.274 |
| molecules favouring bacterial | 0 / 467 | 4 / 467 |

**Shift: +0.114 kcal/mol, or 8% of the gap.**

Matching resolution and cofactor oxidation state moves the result by 8%. The
remaining 92% survives a structure matched to the human enzyme on both
counts. The negative result is therefore **stronger** than the previous
report claimed, not weaker: the most obvious way for it to have been a
receptor artefact has been tested and ruled out.

What is still unexplained is why a 1.90 Å structure of the same enzyme, with
the same inhibitor and the same reduced cofactor, scores every molecule in the
panel more weakly than the human orthologue. Ligand size is already excluded
(the score/heavy-atom correlation is marginally positive at every receptor).
The remaining candidate is the pocket, which would make the offset a property
of the generated chemotype — large and lipophilic — rather than of the method.
An actives-versus-decoys enrichment run at both structures would settle it.
This screen does not.

## 3. The headline selectivity axis moved

`dhfr_selectivity` is now score(2W3A) − score(6XG5); it was score(2W3A) −
score(7MYM). Both are reported per receptor as `selectivity_<receptor>`, so
neither is hidden.

## 4. Table 1 cofactor rows corrected

Table 1 previously listed both 7MYM and 2W3A as "trimethoprim + NADPH"; 7MYM
holds NAP (NADP⁺, oxidised). Methods said the opposite — "NADP⁺ in the two
DHFR structures" — which made 2W3A wrong there instead. Between them the two
errors propped up the claim of "the same fold with the same cofactor" in the
Summary and in §2.2.

`panel500.yaml` recorded the oxidation states correctly in every receptor
block throughout. The error was prose-only and never reached the docking.

## 5. Denominator: 492 → 498 of 500

All eight losses were ETKDG conformer-embedding failures, not chemistry
filters. Root cause: `n_scrub_attempts: 5` retried an *identical*
deterministic call five times, while molscrub's own `num_etkdg_attempts` and
`use_random_coords` were left at defaults — the retry budget bought nothing.
Both are now set.

Six of the eight were recovered. Two still fail to embed, which bounds what
the generator produces that is embeddable at all and is reported as such.
None of the six ranks among the panel's strongest at any target, so the loss
was not concealing a hit population.

## 6. Three silent failures removed from the pipeline

- `read_checkpoint` counted failed records as done, so re-running ligand prep
  would have skipped the eight losses with no message. Now gated behind
  `--retry-failed`.
- `analyze.py` and `merge_results.py` each chose the bacterial DHFR by column
  order. Adding 6XG5 sorts it ahead of 7MYM, so `dhfr_selectivity` would have
  silently changed meaning and every downstream number would have inherited
  it. Both now take it from `primary_bacterial_dhfr` and fail loudly when it
  is unset and ambiguous. `analyze.py` is the one that feeds the report.
- `make_report.py` labelled the headline column by column order too, which
  would have mislabelled it without failing. It now identifies the primary
  from the data.

## 7. Figure 1 overlays both bacterial DHFRs

The selectivity histogram now shows 6XG5 and 7MYM on the same axis and the
same bins, so the two distributions can be read against each other directly.

---

## What did not change

The conclusion. No molecule in the set shows useful selectivity for a
bacterial DHFR over the human enzyme, and the panel's best GyrB scores still
exceed both co-crystal ligands (08B −8.004, 4S4 −7.970), which remains a
calibration warning rather than a hit.

---

## 8. Reference-set potencies: three are bounds, not measurements

Checking the strongest reference compound (CHEMBL4165664, flagged because
its chemotype was `unclassified`) found its ChEMBL record is **Ki < 0.05 nM**,
not Ki = 0.05 nM — a reference compound quoted below the assay's limit of
quantitation in Zidar et al., *Eur J Med Chem* 154:117 (2018).

`fetch_gyrb_actives.py` never requested `standard_relation`, so it imported
the bound as an equality. Because the selection sorts ascending by potency, a
censored record goes straight to the top of both the potency ranking and the
chemotype round-robin. Three of the ten are affected (`<0.05`, `<1.0`,
`<8.0` nM).

The compounds themselves are sound. CHEMBL4165664 is measured in a direct
single-protein *E. coli* GyrB assay (ChEMBL confidence 9) and is dual-active
against topoisomerase IV, which is the expected profile for an ATP-site
binder, since ParE is the GyrB homologue. No re-docking is needed.

What changed is the correlation. Excluding the three bounds flips the rank
correlation between log potency and score from **+0.13 / +0.43** to
**−0.29 / −0.14** — the right direction rather than the wrong one, though
still weak at seven compounds. The censored values had been inverting it.

The fetcher now requests `standard_relation`, records it per compound, and
skips censored records when ranking. A re-run would therefore return a
slightly different ten.
