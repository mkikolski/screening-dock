#!/usr/bin/env python3
"""Build a DHFR reference set that measures the selectivity axis, not just potency.

Why this is not a copy of fetch_gyrb_actives.py
-----------------------------------------------
The GyrB reference set answered "does Vina rank real inhibitors above our
generated panel?" -- a one-sided question, so compounds were chosen to span
chemotypes at good potency.

The DHFR claim in the report is different, and stronger: it asserts that no
generated molecule scores better against bacterial DHFR than against human
DHFR. That is a claim about a *difference* of two scores. A reference set of
potent DHFR inhibitors cannot test it. What tests it is a set of compounds
whose bacterial-vs-human selectivity is already *measured*, spanning both
directions, docked through the same pipeline:

    selectivity_vina = score(human) - score(E. coli)
    selectivity_expt = log10( Ki_human / Ki_E.coli )

If those two correlate, the report's selectivity axis is measuring selectivity
and the negative result is about chemistry. If the Vina differential is flat
across compounds that really do span 10^5-fold measured selectivity, then the
axis cannot resolve selectivity at all and the headline result is uninformative
rather than negative. Either outcome is worth more than the number of knowns
that beat the panel.

This is why the set is selected by stratifying on measured log10 selectivity
rather than by potency or chemotype.

Mechanism filtering is assay-level, not structural
--------------------------------------------------
The GyrB script needed SMARTS because "DNA gyrase" annotations mix ATP-site
binders with GyrA/DNA cleavage-complex poisons on the same protein complex, so
the annotation alone cannot tell them apart. DHFR has one active site and both
targets here are SINGLE PROTEIN, so every record is already site-specific. The
real contamination is a different thing: ChEMBL's human DHFR target carries 192
ADMET (assay_type A) and 42 toxicity (T) records alongside 1917 binding (B)
ones. Those are not DHFR inhibition and would corrupt the human potency.
Measured against the live service 2026-11-22. Hence ASSAY_TYPE below.

Every compound that ships therefore has, at BOTH enzymes, an uncensored nM
binding constant from a single-protein assay. That is stronger mechanism
evidence than any structural guess, and it is why nothing is dropped on
structure here -- the 2,4-diamino motif is recorded as an annotation only.

Targets are pinned as ChEMBL IDs, not resolved from accessions, and the two
near misses are excluded explicitly because both would silently poison the set:

  CHEMBL2627    E. coli DHFR "type 1" = P00382, the Tn7 dfrA1 plasmid enzyme
                that confers trimethoprim resistance. Its SAR runs opposite to
                the chromosomal enzyme by construction. 7MYM and 6XG5 are both
                chromosomal folA (P0ABQ4), so this target is a mechanism
                mismatch of the same kind as docking a fluoroquinolone into the
                GyrB ATP box.
  CHEMBL1681620 S. aureus DHFR. No S. aureus DHFR structure is in the panel.

Run (rdkit optional -- without it the annotation columns are blank and the
.smi is still correct, so the selection can be reviewed off-cluster):

    python fetch_dhfr_actives.py --selftest
    python fetch_dhfr_actives.py --out-dir /storage/users/kikolski/screening-dock

Writes dhfr_actives.smi (pipeline input), dhfr_actives.csv (provenance, both
potencies, measured selectivity) and dhfr_rejected.csv (everything considered
and why it was not taken). configs/dhfr_knowns.yaml points at the .smi.
"""
from __future__ import annotations

import argparse
import csv
import json
import logging
import math
import statistics
import sys
from pathlib import Path

# The HTTP layer is reused rather than rewritten: it already carries the page
# size, retry budget and on-disk cache that were tuned against this same
# service for the GyrB set (ChEMBL 500s in bursts; limit=1000 and limit=50 both
# fail, 25 works). Importing it also means one place to fix if that changes.
from fetch_gyrb_actives import (  # noqa: E402
    FetchError,
    PAGE,
    _get,
    usable_records,
)

log = logging.getLogger("fetch_dhfr")

# Chromosomal folA only. The accessions are recorded for provenance; nothing is
# looked up by them (an accession lookup is what returned a silently empty
# S. aureus half for GyrB).
TARGETS = {
    "ecoli": {"chembl": "CHEMBL1809", "uniprot": "P0ABQ4",
              "organism": "Escherichia coli",
              "receptors": ("dhfr_ecoli_7mym", "dhfr_ecoli_6xg5")},
    "human": {"chembl": "CHEMBL202", "uniprot": "P00374",
              "organism": "Homo sapiens",
              "receptors": ("dhfr_human_2w3a",)},
}

EXCLUDED_TARGETS = {
    "CHEMBL2627": "E. coli DHFR type 1 (P00382, Tn7 dfrA1) -- trimethoprim-"
                  "resistant plasmid enzyme, not the chromosomal folA docked here",
    "CHEMBL1681620": "S. aureus DHFR -- no S. aureus DHFR structure in the panel",
}

POTENCY_TYPES = ("IC50", "Ki", "Kd")
# Binding assays only. See the module docstring: A and T records on the human
# target are ADMET and toxicity measurements, not DHFR inhibition.
ASSAY_TYPE = "B"

# Sanity bounds on a median potency, in nM.
#   MAX_ANY   a "measurement" of 500 mM is an extrapolation from an inactive
#             curve, not a constant. Two such records sit at the human-selective
#             extreme of the raw intersection (10^8-10^9 nM) and would otherwise
#             define the ladder's endpoint.
#   MIN_BEST  a compound inactive at both enzymes carries no selectivity
#             information, however reliably it was measured.
MAX_ANY = 1e6
MIN_BEST = 1e4

# Recorded, never filtered on. The classical antifolate pharmacophore is the
# 2,4-diaminopyrimidine ring, which is also present in the diaminopteridines
# (methotrexate, aminopterin), the diaminoquinazolines (trimetrexate,
# piritrexim) and the diaminopyrimidines proper (trimethoprim, iclaprim). A
# compound that lacks it is still shipped if the assay evidence holds --
# pemetrexed is 2-amino-4-oxo and is a genuine, if weak, human DHFR inhibitor.
DIAMINOPYRIMIDINE = "[NX3;H1,H2]c1nc([NX3;H1,H2])ncc1"

CSV_COLS = [
    "name", "chembl_id", "pref_name", "max_phase", "smiles",
    "inchikey", "heavy_atoms", "motif",
    "ecoli_nM", "ecoli_n", "ecoli_types",
    "human_nM", "human_n", "human_types",
    "log10_selectivity", "fold_selectivity", "direction",
    "ladder_bin", "selected_as", "assay_confidence",
]


# --- potency ---------------------------------------------------------------

def binding_records(target_chembl_id: str) -> list[dict]:
    """Every uncensored nM binding record for one target.

    Unlike the GyrB fetcher this does not stop at a potency cutoff: the
    selection needs the whole measured range at both enzymes, because a
    compound that is *weak* at human DHFR is exactly what makes it an
    informative bacterial-selective reference.
    """
    rows: list[dict] = []
    for stype in POTENCY_TYPES:
        offset = 0
        while True:
            try:
                data = _get("activity", target_chembl_id=target_chembl_id,
                            standard_type=stype, standard_units="nM",
                            assay_type=ASSAY_TYPE,
                            only="molecule_chembl_id,canonical_smiles,"
                                 "standard_type,standard_relation,"
                                 "standard_value,standard_units,"
                                 "assay_chembl_id,assay_description",
                            limit=PAGE, offset=offset)
            except FetchError as exc:
                # One dead potency type must not discard the pages already in
                # hand -- that is how the GyrB run lost its Ki half once.
                log.error("skipping %s/%s at offset %d: %s",
                          target_chembl_id, stype, offset, exc)
                break
            batch = data.get("activities", [])
            usable, censored = usable_records(batch)
            if censored:
                log.info("  %s/%s: %d censored record(s) dropped at offset %d",
                         target_chembl_id, stype, censored, offset)
            rows.extend(usable)
            if offset == 0:
                log.info("  %s/%s: %s records available", target_chembl_id,
                         stype, data.get("page_meta", {}).get("total_count"))
            if len(batch) < PAGE:
                break
            offset += PAGE
    return rows


def summarise(rows: list[dict]) -> dict[str, dict]:
    """{molecule: {nM, n, types, smiles, assays}} using the MEDIAN.

    Median, not minimum. Trimethoprim has 34 E. coli and 47 human records
    spanning roughly a log unit each; taking the most potent of each would
    build the selectivity ratio out of the two most favourable outliers in
    opposite directions and inflate it.
    """
    by_mol: dict[str, dict] = {}
    for r in rows:
        m = by_mol.setdefault(r["molecule_chembl_id"],
                              {"vals": [], "types": set(), "assays": set(),
                               "smiles": r.get("canonical_smiles")})
        m["vals"].append(float(r["standard_value"]))
        m["types"].add(r["standard_type"])
        if r.get("assay_chembl_id"):
            m["assays"].add(r["assay_chembl_id"])
    return {
        k: {"nM": statistics.median(v["vals"]), "n": len(v["vals"]),
            "types": "/".join(sorted(v["types"])), "smiles": v["smiles"],
            "assays": sorted(v["assays"])}
        for k, v in by_mol.items()
    }


def pair(ecoli: dict, human: dict) -> tuple[list[dict], list[dict]]:
    """(compounds measured at both enzymes, records of what was rejected)."""
    kept, rejected = [], []
    for mol in sorted(set(ecoli) & set(human)):
        e, h = ecoli[mol], human[mol]
        row = {"chembl_id": mol, "smiles": e["smiles"] or h["smiles"],
               "ecoli_nM": e["nM"], "ecoli_n": e["n"], "ecoli_types": e["types"],
               "human_nM": h["nM"], "human_n": h["n"], "human_types": h["types"],
               "assays": e["assays"] + h["assays"]}
        if max(e["nM"], h["nM"]) > MAX_ANY:
            row["why"] = f"median > {MAX_ANY:.0g} nM at one enzyme -- " \
                         "extrapolated from an inactive curve, not a constant"
            rejected.append(row)
            continue
        if min(e["nM"], h["nM"]) > MIN_BEST:
            row["why"] = f"inactive at both enzymes (best {min(e['nM'], h['nM']):.0f}" \
                         f" nM > {MIN_BEST:.0g}) -- carries no selectivity signal"
            rejected.append(row)
            continue
        if not row["smiles"]:
            row["why"] = "no structure in ChEMBL"
            rejected.append(row)
            continue
        row["log10_selectivity"] = math.log10(h["nM"] / e["nM"])
        kept.append(row)
    for mol in sorted(set(ecoli) ^ set(human)):
        side = "E. coli" if mol in ecoli else "human"
        rejected.append({"chembl_id": mol, "why": f"measured at {side} only -- "
                                                  "no measured selectivity"})
    return kept, rejected


# --- selection -------------------------------------------------------------

# Bin edges in log10(Ki_human / Ki_ecoli). Positive = bacterial-selective,
# which is the direction the report's dhfr_selectivity column is built to find.
# Deliberately coarse and asymmetric: the pool is dense near zero and thin past
# a factor of 100, so equal-width bins would spend the whole budget on
# non-selective compounds and leave the informative tails with one point each.
LADDER = [
    (2.0, 99.0, "bacterial >100x"),
    (1.0, 2.0, "bacterial 10-100x"),
    (0.5, 1.0, "bacterial 3-10x"),
    (-0.5, 0.5, "non-selective"),
    (-1.0, -0.5, "human 3-10x"),
    (-2.0, -1.0, "human 10-100x"),
    (-99.0, -2.0, "human >100x"),
]


def ladder_bin(sel: float) -> str:
    for lo, hi, label in LADDER:
        if lo <= sel < hi:
            return label
    return "out of range"


def select(pool: list[dict], per_bin: int, max_compounds: int) -> list[dict]:
    """Stratify on measured selectivity, preferring the best-measured compound.

    Two rules, in order:

    1. A named drug beats an unnamed ChEMBL ID at equal measurement weight.
       The set has to be *readable* in a report: "trimethoprim sits here and
       methotrexate sits there" is an argument a reader can check against what
       they already know, where "CHEMBL3827532 sits here" is not. The named
       antifolates also happen to carry the most measurements by a wide margin
       (methotrexate n=127, trimethoprim n=81).
    2. Within that, more measurements first, since the x-axis of the whole
       experiment is the measured ratio and a one-shot IC50 makes a noisy
       abscissa.

    Note what is NOT a rule: potency. A set chosen for potency collapses onto
    the non-selective bin, because the most potent DHFR inhibitors ever made
    are the tight-binding antifolates that hit both enzymes equally.
    """
    def weight(r):
        named = bool(r.get("pref_name"))
        phase = r.get("max_phase") or 0
        return (named, float(phase), r["ecoli_n"] + r["human_n"])

    buckets: list[list[dict]] = []
    for lo, hi, label in LADDER:
        bucket = sorted((r for r in pool if ladder_bin(r["log10_selectivity"]) == label),
                        key=weight, reverse=True)
        for r in bucket[:per_bin]:
            r["ladder_bin"] = label
            r["selected_as"] = "named drug" if r.get("pref_name") else "ladder fill"
        buckets.append(bucket[:per_bin])
        if not bucket:
            log.warning("ladder bin %-18s is EMPTY -- the selectivity axis will "
                        "have a gap there", label)

    # Round-robin, NOT sort-then-slice. per_bin * len(LADDER) can exceed
    # max_compounds, and slicing a selectivity-sorted list then amputates whole
    # bins from one end: at per_bin=3 it dropped the entire "human >100x" bin
    # and said nothing. A gap at one end of the axis is the one defect this set
    # cannot have, since a correlation fitted to a truncated range is exactly
    # the artefact the reference set exists to rule out.
    chosen: list[dict] = []
    for rank in range(per_bin):
        for bucket in buckets:
            if rank < len(bucket) and len(chosen) < max_compounds:
                chosen.append(bucket[rank])
    dropped = sum(len(b) for b in buckets) - len(chosen)
    if dropped:
        log.warning("max_compounds=%d dropped %d candidate(s); every bin keeps "
                    "at least its best-measured compound", max_compounds, dropped)
    chosen.sort(key=lambda r: -r["log10_selectivity"])
    return chosen


def annotate(rows: list[dict]) -> None:
    """Heavy atoms, InChIKey and the 2,4-diamino motif flag, if rdkit is here.

    Optional on purpose. The .smi needs only a name and a SMILES, so the whole
    selection can be built and reviewed on a machine without rdkit and the
    columns filled in later on the cluster.
    """
    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem import Descriptors
    except ImportError:
        log.warning("rdkit absent -- heavy_atoms/inchikey/motif left blank "
                    "(the .smi and the potency columns are unaffected)")
        for r in rows:
            r.setdefault("heavy_atoms", "")
            r.setdefault("inchikey", "")
            r.setdefault("motif", "")
        return
    RDLogger.DisableLog("rdApp.*")
    patt = Chem.MolFromSmarts(DIAMINOPYRIMIDINE)
    assert patt is not None, f"bad SMARTS: {DIAMINOPYRIMIDINE}"
    for r in rows:
        mol = Chem.MolFromSmiles(r["smiles"])
        if mol is None:
            r["heavy_atoms"] = r["inchikey"] = ""
            r["motif"] = "UNPARSEABLE"
            log.error("rdkit cannot parse %s: %s", r["chembl_id"], r["smiles"])
            continue
        r["heavy_atoms"] = Descriptors.HeavyAtomCount(mol)
        r["inchikey"] = Chem.MolToInchiKey(mol)
        r["motif"] = ("2,4-diaminopyrimidine" if mol.HasSubstructMatch(patt)
                      else "other folate-site")


# --- build -----------------------------------------------------------------

def build(out_dir: Path, per_bin: int, max_compounds: int) -> None:
    pool_by_target = {}
    for key, spec in TARGETS.items():
        log.info("fetching %s (%s, %s)", key, spec["chembl"], spec["uniprot"])
        pool_by_target[key] = summarise(binding_records(spec["chembl"]))
        log.info("  %d distinct molecules with an exact nM binding value",
                 len(pool_by_target[key]))

    pool, rejected = pair(pool_by_target["ecoli"], pool_by_target["human"])
    log.info("%d compounds measured at both enzymes (%d rejected)",
             len(pool), len(rejected))
    if not pool:
        raise SystemExit("No compound has an exact binding value at both "
                         "enzymes. Check ASSAY_TYPE and the target IDs.")

    names = fetch_names([r["chembl_id"] for r in pool])
    for r in pool:
        r["pref_name"], r["max_phase"] = names.get(r["chembl_id"], (None, None))

    chosen = select(pool, per_bin, max_compounds)
    annotate(chosen)

    # An unnamed compound needs a stable, filesystem-safe ligand name, and the
    # pipeline keys every checkpoint and pose file on it. Rank by selectivity so
    # the name itself orders the ladder: dhfr_known_001 is the most
    # bacterial-selective compound in the set.
    for i, r in enumerate(chosen, 1):
        r["name"] = f"dhfr_known_{i:03d}"
        r["fold_selectivity"] = f"{10 ** abs(r['log10_selectivity']):.1f}"
        r["direction"] = ("bacterial" if r["log10_selectivity"] > 0.5
                          else "human" if r["log10_selectivity"] < -0.5
                          else "neither")
        r["log10_selectivity"] = f"{r['log10_selectivity']:+.2f}"
        r["ecoli_nM"] = f"{r['ecoli_nM']:.4g}"
        r["human_nM"] = f"{r['human_nM']:.4g}"
        r["assay_confidence"] = ""

    out_dir.mkdir(parents=True, exist_ok=True)
    smi = out_dir / "dhfr_actives.smi"
    smi.write_text("".join(f"{r['smiles']}\t{r['name']}\n" for r in chosen),
                   encoding="utf-8")

    with open(out_dir / "dhfr_actives.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=CSV_COLS, extrasaction="ignore")
        w.writeheader()
        w.writerows(chosen)

    with open(out_dir / "dhfr_rejected.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["chembl_id", "why", "ecoli_nM",
                                           "human_nM", "ecoli_n", "human_n"],
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rejected)

    (out_dir / "dhfr_actives_provenance.json").write_text(json.dumps({
        "targets": TARGETS, "excluded_targets": EXCLUDED_TARGETS,
        "assay_type": ASSAY_TYPE, "potency_types": POTENCY_TYPES,
        "max_any_nM": MAX_ANY, "min_best_nM": MIN_BEST,
        "ladder": [list(x) for x in LADDER],
        "n_paired": len(pool), "n_selected": len(chosen),
    }, indent=2), encoding="utf-8")

    print(f"\n{len(chosen)} compounds -> {smi}")
    print(f"{'name':<16}{'log10':>7} {'E.coli nM':>11} {'human nM':>11}  "
          f"{'bin':<18} compound")
    for r in chosen:
        print(f"{r['name']:<16}{r['log10_selectivity']:>7} {r['ecoli_nM']:>11} "
              f"{r['human_nM']:>11}  {r['ladder_bin']:<18} "
              f"{r['pref_name'] or r['chembl_id']}")
    print(f"\n{len(chosen) * 3} dockings (3 DHFR receptors) "
          f"= {len(chosen) * 3 / 18.4:.1f} h at 18.4 dockings/h")


def fetch_names(ids: list[str]) -> dict[str, tuple]:
    """{chembl_id: (pref_name, max_phase)}, in pages the service will answer."""
    out: dict[str, tuple] = {}
    for i in range(0, len(ids), 20):
        chunk = ids[i:i + 20]
        try:
            data = _get("molecule", molecule_chembl_id__in=",".join(chunk),
                        only="molecule_chembl_id,pref_name,max_phase", limit=PAGE)
        except FetchError as exc:
            # Names are presentation, not evidence. Losing them must not lose
            # the potencies, so this degrades to unnamed rather than aborting.
            log.error("no names for %d compounds: %s", len(chunk), exc)
            continue
        for m in data.get("molecules", []):
            out[m["molecule_chembl_id"]] = (m.get("pref_name"), m.get("max_phase"))
    return out


# --- selftest --------------------------------------------------------------

def selftest() -> None:
    """Checks the logic that decides what ships, offline."""
    assert ladder_bin(3.44) == "bacterial >100x"
    assert ladder_bin(1.33) == "bacterial 10-100x"
    assert ladder_bin(0.11) == "non-selective"
    assert ladder_bin(-0.92) == "human 3-10x"
    # Bins are half-open [lo, hi), so they close at the more-selective edge:
    # pemetrexed's exactly 100-fold counts as 10-100x, not >100x.
    assert ladder_bin(-2.00) == "human 10-100x"
    assert ladder_bin(-2.01) == "human >100x"
    assert ladder_bin(2.00) == "bacterial >100x"
    # Every bin edge lands in exactly one bin.
    for lo, hi, _ in LADDER:
        for probe in (lo, (lo + hi) / 2):
            hits = [l for a, b, l in LADDER if a <= probe < b]
            assert len(hits) == 1, f"{probe} matched {hits}"

    # summarise takes the median, not the minimum: a single optimistic assay
    # must not set the potency.
    rows = [{"molecule_chembl_id": "M", "canonical_smiles": "CC",
             "standard_value": v, "standard_type": "Ki",
             "assay_chembl_id": f"A{i}"}
            for i, v in enumerate(("1", "10", "1000"))]
    s = summarise(rows)["M"]
    assert s["nM"] == 10.0 and s["n"] == 3, s

    # The two sanity bounds, and the pairing itself.
    e = {"GOOD": {"nM": 7.0, "n": 34, "types": "Ki", "smiles": "CC", "assays": []},
         "WILD": {"nM": 6.3e8, "n": 1, "types": "IC50", "smiles": "CC", "assays": []},
         "DEAD": {"nM": 5e4, "n": 1, "types": "IC50", "smiles": "CC", "assays": []},
         "LONELY": {"nM": 5.0, "n": 2, "types": "Ki", "smiles": "CC", "assays": []}}
    h = {"GOOD": {"nM": 19070.0, "n": 47, "types": "Ki", "smiles": "CC", "assays": []},
         "WILD": {"nM": 12.9, "n": 2, "types": "Ki", "smiles": "CC", "assays": []},
         "DEAD": {"nM": 6e4, "n": 1, "types": "IC50", "smiles": "CC", "assays": []}}
    kept, rej = pair(e, h)
    assert [r["chembl_id"] for r in kept] == ["GOOD"], kept
    assert abs(kept[0]["log10_selectivity"] - 3.4357) < 1e-3, kept[0]
    why = {r["chembl_id"]: r["why"] for r in rej}
    assert "extrapolated" in why["WILD"] and "inactive at both" in why["DEAD"]
    assert "E. coli only" in why["LONELY"], why["LONELY"]

    # A named phase-4 drug outranks a better-measured unnamed compound, and an
    # empty bin warns instead of silently shrinking the set.
    pool = [{"chembl_id": "A", "log10_selectivity": 3.4, "pref_name": "TRIMETHOPRIM",
             "max_phase": 4, "ecoli_n": 34, "human_n": 47, "smiles": "CC"},
            {"chembl_id": "B", "log10_selectivity": 3.0, "pref_name": None,
             "max_phase": None, "ecoli_n": 90, "human_n": 90, "smiles": "CC"}]
    got = select(pool, per_bin=1, max_compounds=9)
    assert [r["chembl_id"] for r in got] == ["A"], got
    assert got[0]["selected_as"] == "named drug"

    # max_compounds truncates, and the surviving set stays sorted.
    pool2 = [{"chembl_id": f"C{i}", "log10_selectivity": s, "pref_name": None,
              "max_phase": None, "ecoli_n": 1, "human_n": 1, "smiles": "CC"}
             for i, s in enumerate((3.0, 1.5, 0.0, -1.5, -3.0))]
    got2 = select(pool2, per_bin=1, max_compounds=3)
    assert len(got2) == 3 and [r["log10_selectivity"] for r in got2] == [3.0, 1.5, 0.0]

    # The cap must not amputate one end of the ladder. With 7 populated bins,
    # per_bin=3 and a cap of 18, every bin must still be represented -- the
    # earlier sort-then-slice silently dropped all of "human >100x".
    pool3 = [{"chembl_id": f"D{i}-{j}", "log10_selectivity": s, "pref_name": None,
              "max_phase": None, "ecoli_n": j, "human_n": 1, "smiles": "CC"}
             for i, s in enumerate((3.0, 1.5, 0.7, 0.0, -0.7, -1.5, -3.0))
             for j in range(3)]
    got3 = select(pool3, per_bin=3, max_compounds=18)
    assert len(got3) == 18, len(got3)
    bins3 = {r["ladder_bin"] for r in got3}
    assert len(bins3) == 7, sorted(bins3)
    assert max(r["log10_selectivity"] for r in got3) == 3.0
    assert min(r["log10_selectivity"] for r in got3) == -3.0

    # rdkit is optional: annotate() must leave usable blanks, not raise.
    rs = [{"smiles": "CC1=CN=C(N)N=C1N", "chembl_id": "X"}]
    annotate(rs)
    assert "heavy_atoms" in rs[0] and "motif" in rs[0]

    assert "CHEMBL2627" in EXCLUDED_TARGETS and "CHEMBL1809" not in EXCLUDED_TARGETS
    print("selftest OK")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out-dir", type=Path, default=Path("."))
    p.add_argument("--per-bin", type=int, default=3,
                   help="compounds per selectivity bin (7 bins; default 3)")
    p.add_argument("--max-compounds", type=int, default=18,
                   help="hard cap; 3 dockings each (default 18 = ~2.9 h)")
    p.add_argument("--selftest", action="store_true")
    p.add_argument("-v", "--verbose", action="store_true")
    a = p.parse_args()
    logging.basicConfig(level=logging.DEBUG if a.verbose else logging.INFO,
                        format="%(levelname)s %(message)s")
    if a.selftest:
        selftest()
        return
    build(a.out_dir, a.per_bin, a.max_compounds)


if __name__ == "__main__":
    main()
