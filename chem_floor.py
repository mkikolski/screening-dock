#!/usr/bin/env python3
"""How much of the generated panel is in the right chemical space at all.

The ADMET comparison in admet_table.py found that the panel's top-scoring
molecules are physicochemically indistinguishable from the rest of the panel --
so docking score is not degrading properties -- while the whole panel sits
outside the space any known inhibitor of these targets occupies. Two numbers
separate those facts and neither follows from the medians:

  1. What fraction of the 500 generated molecules clears a floor that every
     reference compound clears trivially?
  2. Do any of the molecules that beat human DHFR at some target -- the
     report's shortlist -- clear it?

The second is the one that matters. If the shortlist is empty after the floor,
then the molecules the report puts forward cannot be progressed, and that
belongs next to the shortlist rather than in a footnote.

The floor
---------
At least one aromatic ring, and at most one carboxylic acid. Deliberately a
floor and not an optimisation target: it is set where the reference compounds
are unanimous, not where they are typical. All ten GyrB inhibitors have
aromatic rings (mean 3.8) and none carries a single carboxylate, so a molecule
failing this is not a borderline case -- it is outside the envelope of every
compound known to inhibit the target. The antifolates sit the same side of it,
since 2,4-diaminopyrimidines are aromatic.

Both thresholds are arguable, so --min-aromatic-rings and --max-acids exist and
the sensitivity is printed across neighbouring values. What is not arguable is
the direction.

The substructure counting and the SMILES reading are imported from
admet_table.py rather than rewritten: those SMARTS are covered by its selftest
(a polyacid counts as three acids, an amide nitrogen is not a base), and a
second copy would drift.

    python chem_floor.py --selftest
    python chem_floor.py
"""
from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path

from admet_table import descriptors, read_smi

log = logging.getLogger("chem_floor")


def passes(row: dict, min_ar: int, max_acid: int) -> tuple[bool, str]:
    """(clears the floor, why not).

    A blank descriptor means rdkit could not parse the SMILES. Counted as a
    failure but reported separately: folding it into "fails the floor" would
    overstate the result.
    """
    if row.get("aromatic_rings") == "" or row.get("n_acid") == "":
        return False, "unparseable"
    ar, acid = int(row["aromatic_rings"]), int(row["n_acid"])
    why = []
    if ar < min_ar:
        why.append(f"aromatic_rings {ar} < {min_ar}")
    if acid > max_acid:
        why.append(f"acids {acid} > {max_acid}")
    return (not why), "; ".join(why)


def survival(rows: list[dict], min_ar: int, max_acid: int) -> dict:
    """Counts, plus which condition did the rejecting. Stamps passes/fail_reason."""
    out = {"n": len(rows), "pass": 0, "fail_aromatic": 0, "fail_acid": 0,
           "fail_both": 0, "unparseable": 0}
    for r in rows:
        ok, why = passes(r, min_ar, max_acid)
        r["passes"] = int(ok)
        r["fail_reason"] = why
        if ok:
            out["pass"] += 1
        elif why == "unparseable":
            out["unparseable"] += 1
        elif "aromatic_rings" in why and "acids" in why:
            out["fail_both"] += 1
        elif "aromatic_rings" in why:
            out["fail_aromatic"] += 1
        else:
            out["fail_acid"] += 1
    return out


def pct(n: int, d: int) -> str:
    return "n/a" if not d else f"{100 * n / d:.1f}%"


def show(label: str, s: dict) -> None:
    print(f"  {label:<42} {s['pass']:>4} / {s['n']:<4} = {pct(s['pass'], s['n']):>6}")
    if s["pass"] < s["n"]:
        bits = [f"no aromatic ring {s['fail_aromatic']}",
                f"too many acids {s['fail_acid']}", f"both {s['fail_both']}"]
        if s["unparseable"]:
            bits.append(f"unparseable {s['unparseable']}")
        print(f"  {'':<42} failed on: " + ", ".join(bits))


def build(args) -> None:
    panel_smi = read_smi(args.panel_smi)
    scored = {r["ligand_base"]: r
              for r in csv.DictReader(args.panel_csv.open(encoding="utf-8"))}
    log.info("%d generated molecules in %s, %d scored in %s",
             len(panel_smi), args.panel_smi.name, len(scored), args.panel_csv.name)

    rows = [{"id": name, "set": "panel", "smiles": smi,
             "best_bacterial": scored.get(name, {}).get("best_bacterial", ""),
             "margin_vs_human": scored.get(name, {}).get("margin_vs_human", ""),
             "scored": name in scored}
            for name, smi in panel_smi.items()]

    knowns: list[dict] = []
    for d in args.knowns_dir:
        for a in (sorted(d.glob("*_actives.csv")) if d.is_dir() else []):
            for r in csv.DictReader(a.open(encoding="utf-8")):
                knowns.append({"id": r["name"],
                               "set": a.name.replace("_actives.csv", ""),
                               "smiles": r.get("smiles", ""),
                               "best_bacterial": "", "margin_vs_human": "",
                               "scored": False})
    if not knowns:
        log.warning("no reference compounds in %s -- the floor's calibration "
                    "cannot be checked",
                    ", ".join(str(d) for d in args.knowns_dir))

    if not descriptors(rows + knowns):
        raise SystemExit("rdkit is required here: the whole analysis is a "
                         "substructure count. Use .venv/bin/python.")

    # Analysis 2's population. margin_vs_human = score(human) - score(best
    # bacterial), so POSITIVE means the bacterial target is favoured. Taking
    # this the other way would select the most human-selective molecules and
    # invert the question.
    shortlist = [r for r in rows if r["margin_vs_human"] not in ("", None)
                 and float(r["margin_vs_human"]) > 0]
    scored_rows = [r for r in rows if r["scored"]]

    print(f"\nFloor: at least {args.min_aromatic_rings} aromatic ring(s) and at "
          f"most {args.max_acids} carboxylic acid(s)\n")
    for label, subset in (("all generated molecules", rows),
                          ("  of those, the scored ones", scored_rows),
                          ("  of those, beating human DHFR (margin>0)", shortlist)):
        show(label, survival(subset, args.min_aromatic_rings, args.max_acids))
    for s in sorted({k["set"] for k in knowns}):
        show(f"reference compounds ({s})",
             survival([k for k in knowns if k["set"] == s],
                      args.min_aromatic_rings, args.max_acids))

    if shortlist:
        survival(shortlist, args.min_aromatic_rings, args.max_acids)
        surv = [r for r in shortlist if r["passes"]]
        print(f"\nShortlist detail ({len(shortlist)} molecules beat human DHFR "
              f"at some target):")
        if surv:
            print(f"  {len(surv)} clear the floor:")
            for r in sorted(surv, key=lambda x: float(x["margin_vs_human"]),
                            reverse=True):
                print(f"    {r['id'][-22:]}  margin {float(r['margin_vs_human']):+.3f}"
                      f"  score {r['best_bacterial']}"
                      f"  rings {r['aromatic_rings']} acids {r['n_acid']}")
        else:
            print("  NONE clears the floor.")
            print("  best five by margin, and why each fails:")
            for r in sorted(shortlist, key=lambda x: -float(x["margin_vs_human"]))[:5]:
                print(f"    {r['id'][-22:]}  margin {float(r['margin_vs_human']):+.3f}"
                      f"  rings {r['aromatic_rings']} acids {r['n_acid']}"
                      f"  -- {r['fail_reason']}")

    print("\nSensitivity to where the floor is set:")
    print(f"  {'min aromatic rings':>20} {'max acids':>10} {'scored panel':>14}"
          f" {'shortlist':>11}")
    for ar in (1, 2):
        for ac in (0, 1, 2):
            a = survival(scored_rows, ar, ac)
            b = survival(shortlist, ar, ac)
            print(f"  {ar:>20} {ac:>10} {pct(a['pass'], a['n']):>14}"
                  f" {pct(b['pass'], b['n']):>11}")
    # The sensitivity sweep left the flags on its last cell, so restamp at the
    # real floor before the rows are written out.
    survival(rows + knowns, args.min_aromatic_rings, args.max_acids)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    cols = ["id", "set", "smiles", "aromatic_rings", "n_acid", "passes",
            "fail_reason", "best_bacterial", "margin_vs_human"]
    with args.out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows + knowns)
    print(f"\nPer-molecule verdicts -> {args.out}")


def selftest() -> None:
    mk = lambda ar, ac: {"aromatic_rings": ar, "n_acid": ac}
    assert passes(mk(1, 0), 1, 1)[0]
    assert passes(mk(4, 0), 1, 1)[0]
    assert passes(mk(1, 1), 1, 1)[0]          # boundary, both inclusive
    ok, why = passes(mk(0, 0), 1, 1)
    assert not ok and why == "aromatic_rings 0 < 1", why
    ok, why = passes(mk(2, 3), 1, 1)
    assert not ok and why == "acids 3 > 1", why
    ok, why = passes(mk(0, 4), 1, 1)
    assert not ok and "aromatic_rings" in why and "acids" in why, why
    assert passes({"aromatic_rings": "", "n_acid": ""}, 1, 1) == (False, "unparseable")

    s = survival([mk(1, 0), mk(0, 0), mk(2, 5), mk(0, 9),
                  {"aromatic_rings": "", "n_acid": ""}], 1, 1)
    assert s == {"n": 5, "pass": 1, "fail_aromatic": 1, "fail_acid": 1,
                 "fail_both": 1, "unparseable": 1}, s

    # survival() must stamp the flags it reports: the CSV is written from them
    # and a stale flag would disagree with the printed counts.
    rows = [mk(1, 0), mk(0, 3)]
    survival(rows, 1, 1)
    assert rows[0]["passes"] == 1 and rows[1]["passes"] == 0
    assert "aromatic_rings" in rows[1]["fail_reason"]
    assert pct(1, 4) == "25.0%" and pct(0, 0) == "n/a"

    # The real SMARTS end to end, on molecules from each side of the floor.
    probes = [
        ("COc1cc(Cc2cnc(N)nc2N)cc(OC)c1OC", True),   # trimethoprim
        ("CN(Cc1cnc2nc(N)nc(N)c2n1)c1ccc(C(=O)NC(CCC(=O)O)C(=O)O)cc1", False),
        ("OC(=O)C1COCC(C(=O)O)(N(C)CC(=O)O)C1", False),   # panel-like triacid
    ]
    rs = [{"smiles": s, "id": f"p{i}"} for i, (s, _) in enumerate(probes)]
    if descriptors(rs):
        for r, (smi, want) in zip(rs, probes):
            got, why = passes(r, 1, 1)
            assert got is want, (smi, r["aromatic_rings"], r["n_acid"], why)
        print("selftest OK")
    else:
        print("selftest OK (rdkit absent: SMARTS probes skipped)")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--panel-smi", type=Path, default=Path("exp_c5_selected_500.smi"))
    p.add_argument("--panel-csv", type=Path,
                   default=Path("report/data/per_molecule.csv"))
    p.add_argument("--knowns-dir", type=Path, nargs="*",
                   default=[Path("report/data_knowns"), Path("report/data_dhfr")])
    p.add_argument("--min-aromatic-rings", type=int, default=1)
    p.add_argument("--max-acids", type=int, default=1)
    p.add_argument("--out", type=Path, default=Path("report/admet/chem_floor.csv"))
    p.add_argument("--selftest", action="store_true")
    a = p.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    if a.selftest:
        selftest()
        return
    build(a)


if __name__ == "__main__":
    main()
