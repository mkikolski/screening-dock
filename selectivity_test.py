#!/usr/bin/env python3
"""Does the bacterial-minus-human score difference measure selectivity?

The report's headline is a claim about a difference of two Vina scores: that no
generated molecule scores better against bacterial DHFR than against human
DHFR. Nothing in the panel itself can show that difference carries information.
So 18 DHFR inhibitors whose selectivity is already MEASURED were docked through
the same pipeline, spanning +3.44 to -2.18 in log10(Ki_human / Ki_E.coli) --
about 400,000-fold -- and this script regresses one on the other.

Three tests, because the first is underpowered on its own:

  1. Correlation across all 18. Direct but weak at this n.
  2. Bacterial-selective group vs human-selective group, by permutation. More
     power than a correlation over a noisy range, and it is the comparison the
     report actually needs.
  3. Calibration of the panel against where the real bacterial-selective drugs
     land. This replaces counting molecules that cross zero, which turns out to
     be the wrong threshold: most genuinely selective drugs do not cross it
     either, so "0 of 492" measured the threshold rather than the molecules.

Sign convention, taken from analyze.py and verified against the data:
selectivity = score(human) - score(bacterial), so POSITIVE favours bacterial.
Ranking it the other way returns the most human-selective compounds.

    python selectivity_test.py --selftest
    python selectivity_test.py
"""
from __future__ import annotations

import argparse
import csv
import math
import random
import statistics as st
from pathlib import Path

RECEPTORS = ("dhfr_ecoli_6xg5", "dhfr_ecoli_7mym")
# log10 selectivity beyond which a compound counts as taking a side. 0.5 is
# about three-fold, comfortably outside assay scatter, and it leaves the
# genuinely non-selective antifolates (methotrexate +0.32, aminopterin +0.03,
# trimetrexate -0.09) out of both groups where they belong.
GROUP_CUT = 0.5
# kcal/mol per log10 unit of affinity ratio at 298 K: 2.303 RT.
KCAL_PER_LOG = 1.364


def rank(v: list[float]) -> list[float]:
    """Average ranks, so ties do not bias the Spearman coefficient."""
    order = sorted(range(len(v)), key=lambda i: v[i])
    out = [0.0] * len(v)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
            j += 1
        avg = (i + j) / 2 + 1
        for k in range(i, j + 1):
            out[order[k]] = avg
        i = j + 1
    return out


def pearson(a: list[float], b: list[float]) -> float:
    ma, mb = st.mean(a), st.mean(b)
    num = sum((x - ma) * (y - mb) for x, y in zip(a, b))
    den = math.sqrt(sum((x - ma) ** 2 for x in a) * sum((y - mb) ** 2 for y in b))
    return num / den if den else float("nan")


def spearman(a: list[float], b: list[float]) -> float:
    return pearson(rank(a), rank(b))


def permutation_p(a: list[float], b: list[float], n: int = 200_000,
                  seed: int = 2137) -> tuple[float, float]:
    """(observed difference of means, two-sided p).

    A permutation test rather than a t-test: n = 9 against 6 with no reason to
    expect normal residuals. Seeded so the p-value in the report is stable.
    """
    obs = st.mean(a) - st.mean(b)
    pool = list(a) + list(b)
    na = len(a)
    rnd = random.Random(seed)
    hits = 0
    for _ in range(n):
        rnd.shuffle(pool)
        if abs(st.mean(pool[:na]) - st.mean(pool[na:])) >= abs(obs) - 1e-12:
            hits += 1
    return obs, (hits + 1) / (n + 1)


def load(dhfr_dir: Path, panel_csv: Path):
    act = {r["name"]: r for r in
           csv.DictReader((dhfr_dir / "dhfr_actives.csv").open(encoding="utf-8"))}
    kn = {r["ligand_base"]: r for r in
          csv.DictReader((dhfr_dir / "per_molecule.csv").open(encoding="utf-8"))}
    rows = []
    for name, a in act.items():
        s = kn.get(name)
        if not s:
            continue
        row = {"compound": a.get("pref_name") or a["chembl_id"],
               "measured_log10": float(a["log10_selectivity"]),
               "heavy_atoms": int(a["heavy_atoms"]) if a.get("heavy_atoms") else ""}
        for rec in RECEPTORS:
            col = f"selectivity_{rec}"
            row[f"vina_{rec}"] = float(s[col]) if col in s else None
        row["group"] = ("bacterial" if row["measured_log10"] > GROUP_CUT else
                        "human" if row["measured_log10"] < -GROUP_CUT else "neither")
        rows.append(row)
    missing = set(act) - set(kn)

    # Which bacterial DHFR is the panel's dhfr_selectivity actually built on?
    # Read from the data, never assumed from the config. The panel and the
    # reference set are separate runs that can carry different receptor sets --
    # the panel's per_molecule.csv had no 6XG5 column at all while the
    # reference run did -- and comparing the panel at one receptor against the
    # knowns at another is the same cross-receptor error already removed from
    # merge_results.py and analyze.py. Found here only because this script
    # made it too, and quoted a panel-at-7MYM figure under a 6XG5 heading.
    prows = list(csv.DictReader(panel_csv.open(encoding="utf-8")))
    panel = [float(r["dhfr_selectivity"]) for r in prows]
    basis = None
    if prows:
        human = next((c for c in prows[0] if "dhfr" in c and "human" in c), None)
        for rec in (c for c in prows[0] if "dhfr" in c and "human" not in c
                    and not c.startswith("selectivity")):
            if human and all(
                    abs((float(r[human]) - float(r[rec])) - float(r["dhfr_selectivity"]))
                    < 1e-6 for r in prows):
                basis = rec
                break
    return sorted(rows, key=lambda r: -r["measured_log10"]), panel, missing, basis


def report(rows, panel, missing, basis, out: Path) -> None:
    if basis is None:
        print("WARNING could not tell which bacterial DHFR the panel's "
              "dhfr_selectivity is built on; panel comparisons skipped")
    else:
        print(f"panel dhfr_selectivity is built on {basis}")
    if missing:
        print(f"WARNING {len(missing)} reference compound(s) have no score: "
              f"{', '.join(sorted(missing))}")
    print(f"n = {len(rows)} reference compounds, {len(panel)} generated molecules\n")
    print(f"{'compound':<16}{'measured':>10}" +
          "".join(f"{r.split('_')[-1]:>11}" for r in RECEPTORS))
    for r in rows:
        print(f"{r['compound'][:15]:<16}{r['measured_log10']:>+10.2f}" +
              "".join(f"{r[f'vina_{rec}']:>+11.3f}" for rec in RECEPTORS))

    exp = [r["measured_log10"] for r in rows]
    for rec in RECEPTORS:
        v = [r[f"vina_{rec}"] for r in rows]
        bact = [r[f"vina_{rec}"] for r in rows if r["group"] == "bacterial"]
        hum = [r[f"vina_{rec}"] for r in rows if r["group"] == "human"]
        obs, p = permutation_p(bact, hum)
        ref = st.mean(bact)
        print(f"\n=== {rec} ===")
        print(f"  correlation with measured selectivity   "
              f"Pearson {pearson(exp, v):+.3f}   Spearman {spearman(exp, v):+.3f}")
        print(f"  bacterial-selective knowns  n={len(bact)}  mean {st.mean(bact):+.3f}")
        print(f"  human-selective knowns      n={len(hum)}  mean {st.mean(hum):+.3f}")
        print(f"  difference {obs:+.3f} kcal/mol, correct direction, "
              f"permutation p = {p:.3f}")
        print(f"  knowns with positive selectivity: "
              f"{sum(1 for x in v if x > 0)}/{len(v)}  <- zero is the wrong threshold")
        # Like-for-like only. Quoting the panel at one bacterial DHFR under
        # another receptor's heading is the whole error this guard exists for.
        if basis == rec:
            reach = sum(1 for x in panel if x > ref)
            print(f"  generated panel: median {st.median(panel):+.3f}, "
                  f"{ref - st.mean(panel):+.3f} kcal/mol below the "
                  f"bacterial-selective mean")
            print(f"  panel molecules reaching that mean: {reach}/{len(panel)} "
                  f"({100 * reach / len(panel):.1f}%)")
        else:
            print(f"  panel comparison SKIPPED: the panel's dhfr_selectivity is "
                  f"built on {basis or 'an unknown receptor'}, not {rec}. "
                  f"Re-run merge_results.py then analyze.py so the panel carries "
                  f"a selectivity_{rec} column.")
        span_real = (max(exp) - min(exp)) * KCAL_PER_LOG
        span_vina = max(v) - min(v)
        print(f"  dynamic range: {span_real:.1f} kcal/mol of real ddG mapped onto "
              f"{span_vina:.2f} kcal/mol ({span_real / span_vina:.1f}x compression)")

    out.parent.mkdir(parents=True, exist_ok=True)
    cols = ["compound", "measured_log10", *[f"vina_{r}" for r in RECEPTORS],
            "heavy_atoms", "group"]
    with out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"\nPer-compound values -> {out}")


def selftest() -> None:
    assert rank([10.0, 20.0, 30.0]) == [1.0, 2.0, 3.0]
    assert rank([5.0, 5.0, 9.0]) == [1.5, 1.5, 3.0]          # ties averaged
    assert abs(pearson([1, 2, 3], [2, 4, 6]) - 1.0) < 1e-12
    assert abs(pearson([1, 2, 3], [6, 4, 2]) + 1.0) < 1e-12
    assert math.isnan(pearson([1, 1, 1], [1, 2, 3]))         # zero variance
    # Spearman is rank-based, so a monotone but curved relation is still 1.0
    assert abs(spearman([1, 2, 3, 4], [1, 4, 9, 16]) - 1.0) < 1e-12

    # Separated groups give a small p; identical groups give p near 1.
    obs, p = permutation_p([5.0, 5.1, 4.9], [1.0, 1.1, 0.9], n=2000)
    assert abs(obs - 4.0) < 1e-9 and p < 0.15, (obs, p)
    obs, p = permutation_p([1.0, 2.0, 3.0], [1.0, 2.0, 3.0], n=2000)
    assert abs(obs) < 1e-12 and p > 0.9, (obs, p)
    # Seeded, so the reported p does not drift between runs.
    assert permutation_p([3.0, 4.0], [1.0, 2.0], n=500) == \
           permutation_p([3.0, 4.0], [1.0, 2.0], n=500)

    # The group cut must leave the non-selective antifolates out of both sides.
    for v, want in ((3.44, "bacterial"), (0.98, "bacterial"), (0.32, "neither"),
                    (-0.09, "neither"), (-0.85, "human"), (-2.18, "human")):
        got = ("bacterial" if v > GROUP_CUT else
               "human" if v < -GROUP_CUT else "neither")
        assert got == want, (v, got, want)
    # The basis detection, which is what stops a cross-receptor comparison.
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        t = Path(td)
        (t / "dhfr_actives.csv").write_text(
            "name,chembl_id,pref_name,log10_selectivity,heavy_atoms\n"
            "k1,C1,ALPHA,+3.00,21\nk2,C2,BETA,-3.00,22\n", encoding="utf-8")
        (t / "per_molecule.csv").write_text(
            "ligand_base,selectivity_dhfr_ecoli_6xg5,selectivity_dhfr_ecoli_7mym\n"
            "k1,0.40,0.03\nk2,-2.00,-1.90\n", encoding="utf-8")
        # Panel whose dhfr_selectivity is human - 7MYM: basis must read 7MYM,
        # never 6XG5, even though a 6XG5 column is present.
        p = t / "panel.csv"
        p.write_text(
            "ligand_base,dhfr_ecoli_7mym,dhfr_ecoli_6xg5,dhfr_human_2w3a,"
            "dhfr_selectivity\n"
            "m1,-6.0,-6.5,-7.4,-1.4\nm2,-6.2,-6.9,-7.0,-0.8\n", encoding="utf-8")
        rows, panel, missing, basis = load(t, p)
        assert basis == "dhfr_ecoli_7mym", basis
        assert len(rows) == 2 and panel == [-1.4, -0.8] and not missing
        assert [r["group"] for r in rows] == ["bacterial", "human"]

        # A panel whose column matches nothing must report no basis rather than
        # guessing one, so the comparison is skipped instead of being wrong.
        p.write_text(
            "ligand_base,dhfr_ecoli_7mym,dhfr_human_2w3a,dhfr_selectivity\n"
            "m1,-6.0,-7.4,+9.9\n", encoding="utf-8")
        assert load(t, p)[3] is None

        # A reference compound with no score is reported, not silently dropped.
        (t / "per_molecule.csv").write_text(
            "ligand_base,selectivity_dhfr_ecoli_6xg5,selectivity_dhfr_ecoli_7mym\n"
            "k1,0.40,0.03\n", encoding="utf-8")
        assert load(t, p)[2] == {"k2"}
    print("selftest OK")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dhfr-dir", type=Path, default=Path("report/data_dhfr"))
    p.add_argument("--panel-csv", type=Path,
                   default=Path("report/data/per_molecule.csv"))
    p.add_argument("--out", type=Path,
                   default=Path("report/admet/selectivity_test.csv"))
    p.add_argument("--selftest", action="store_true")
    a = p.parse_args()
    if a.selftest:
        selftest()
        return
    rows, panel, missing, basis = load(a.dhfr_dir, a.panel_csv)
    if not rows:
        raise SystemExit(f"No reference compound in {a.dhfr_dir} has a score.")
    report(rows, panel, missing, basis, a.out)


if __name__ == "__main__":
    main()
