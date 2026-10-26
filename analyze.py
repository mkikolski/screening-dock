#!/usr/bin/env python3
"""Post-run analysis: selectivity, ligand efficiency, and the size bias.

Reads docking/scores_by_ligand.csv (written by merge_results.py) and answers
the three questions the raw score table does not:

1. Does the panel favour the bacterial targets or the human one? The run's
   headline result is that it does not, and this quantifies per molecule how
   badly.
2. Which compounds are worth looking at? Not the ones with the best Vina
   score -- Vina's score grows with molecular size, so on a panel of median
   41 heavy atoms the top of that ranking is largely just the biggest
   molecules. Ligand efficiency (-score / heavy atoms) is the honest ranking.
3. How strong is that size bias here? Reported as a correlation per receptor,
   so the caveat is a number rather than an assertion.

    python analyze.py [--config configs/panel500.yaml] [--le-min 0.30]
                      [--margin-min 0.5]

Writes docking/analysis/{per_molecule,shortlist,top_by_le}.csv and prints a
report. Read-only with respect to everything the pipeline produced.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from dockpipe.common import PipelineConfig

HUMAN_KEY = "human"


def load(dock_dir: Path) -> pd.DataFrame:
    path = dock_dir / "scores_by_ligand.csv"
    if not path.exists():
        raise SystemExit(f"Missing {path}. Run merge_results.py first.")
    df = pd.read_csv(path)
    missing = {"ligand_base", "receptor", "best_score", "ligand_efficiency"} - set(df)
    if missing:
        raise SystemExit(f"{path} lacks columns: {sorted(missing)}")
    return df


def pivot(df: pd.DataFrame, values: str) -> pd.DataFrame:
    out = df.pivot_table(index="ligand_base", columns="receptor", values=values)
    # Otherwise "receptor" prints as a stray header line above the real columns.
    out.columns.name = None
    return out


def fmt(v: float) -> str:
    return f"{v:.2f}"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default="configs/panel500.yaml")
    ap.add_argument("--le-min", type=float, default=0.30,
                    help="ligand-efficiency floor for the shortlist")
    ap.add_argument("--margin-min", type=float, default=0.5,
                    help="kcal/mol a molecule must beat the human target by")
    ap.add_argument("--top", type=int, default=15)
    args = ap.parse_args()

    cfg = PipelineConfig.from_yaml(args.config)
    dock_dir = Path(cfg.output_root) / cfg.run_name / "docking"
    out_dir = dock_dir / "analysis"
    out_dir.mkdir(parents=True, exist_ok=True)

    df = load(dock_dir)
    scores = pivot(df, "best_score")
    le = pivot(df, "ligand_efficiency")

    human = [c for c in scores.columns if HUMAN_KEY in c]
    bact = [c for c in scores.columns if HUMAN_KEY not in c]
    if not human:
        raise SystemExit("No human receptor in the panel; nothing to compare against.")
    human_col = human[0]

    heavy = (df.groupby("ligand_base")["heavy_atoms"].first()
             if "heavy_atoms" in df else None)

    print(f"\n{len(scores)} molecules x {len(scores.columns)} receptors")
    print(f"human target: {human_col}")
    print(f"bacterial   : {', '.join(bact)}")

    # ---- 1. size bias, so the caveat is a number -------------------------
    if heavy is not None and heavy.notna().any():
        print("\n--- Vina score vs molecular size ---")
        print("  negative r = bigger molecules score 'better'; this is why the")
        print("  ranking below is by ligand efficiency, not by raw score")
        print(f"  heavy atoms: median={heavy.median():.0f} "
              f"min={heavy.min():.0f} max={heavy.max():.0f}")
        for rec in scores.columns:
            pair = pd.concat([scores[rec], heavy], axis=1).dropna()
            if len(pair) > 2:
                print(f"  {rec:<22} pearson={pair.corr().iloc[0, 1]:+.2f}  "
                      f"spearman={pair.corr(method='spearman').iloc[0, 1]:+.2f}")

    # ---- 2. selectivity --------------------------------------------------
    # Positive = the bacterial target scores better (more negative) than human.
    per_mol = scores.copy()
    per_mol["best_bacterial"] = scores[bact].min(axis=1)
    per_mol["best_bacterial_target"] = scores[bact].idxmin(axis=1)
    per_mol["margin_vs_human"] = scores[human_col] - per_mol["best_bacterial"]

    dhfr_bact = [c for c in bact if "dhfr" in c]
    if dhfr_bact:
        # Like-for-like: same fold, same cofactor, both passed redocking, so
        # this is the one cross-target comparison the data really supports.
        per_mol["dhfr_selectivity"] = scores[human_col] - scores[dhfr_bact[0]]

    if heavy is not None:
        per_mol["heavy_atoms"] = heavy
    per_mol["le_at_best_bacterial"] = [
        le.at[m, t] if pd.notna(t) else float("nan")
        for m, t in zip(per_mol.index, per_mol["best_bacterial_target"])
    ]

    print(f"\n--- Selectivity vs {human_col} ---")
    for label, col in [("any bacterial target", "margin_vs_human"),
                       (f"{dhfr_bact[0] if dhfr_bact else 'n/a'} (like-for-like)",
                        "dhfr_selectivity")]:
        if col not in per_mol:
            continue
        v = per_mol[col].dropna()
        if v.empty:
            continue
        n = len(v)
        print(f"  {label}: median={v.median():+.2f} kcal/mol")
        print(f"    favours bacterial at all : {(v > 0).sum():4d}/{n} "
              f"({100 * (v > 0).mean():.0f}%)")
        print(f"    by >= {args.margin_min} kcal/mol      : "
              f"{(v >= args.margin_min).sum():4d}/{n}")
        print(f"    by >= 1.0 kcal/mol       : {(v >= 1.0).sum():4d}/{n}")
        print(f"    quartiles: {v.quantile(.25):+.2f} / {v.median():+.2f} "
              f"/ {v.quantile(.75):+.2f}")

    # ---- 3. the shortlist ------------------------------------------------
    sel = per_mol[(per_mol["margin_vs_human"] >= args.margin_min)
                  & (per_mol["le_at_best_bacterial"] >= args.le_min)]
    sel = sel.sort_values("le_at_best_bacterial", ascending=False)

    print(f"\n--- Shortlist: beats {human_col} by >= {args.margin_min} "
          f"kcal/mol AND LE >= {args.le_min} ---")
    if sel.empty:
        print("  EMPTY -- no molecule here is both selective for a bacterial")
        print("  target and efficient. The closest candidates, to show how far")
        print("  off they are:")
        near = per_mol.sort_values("margin_vs_human", ascending=False).head(args.top)
        print(near[["best_bacterial_target", "margin_vs_human",
                    "le_at_best_bacterial"]].to_string(float_format=fmt))
    else:
        print(f"  {len(sel)} molecule(s):")
        cols = ["best_bacterial_target", "best_bacterial", "margin_vs_human",
                "le_at_best_bacterial"]
        if "heavy_atoms" in sel:
            cols.append("heavy_atoms")
        print(sel[cols].head(args.top).to_string(float_format=fmt))

    # ---- 4. best per receptor, by LE not by score ------------------------
    print(f"\n--- Top {args.top} by ligand efficiency, per receptor ---")
    keep = [c for c in ("ligand_base", "best_score", "ligand_efficiency",
                        "heavy_atoms") if c in df]
    tops = []
    for rec in sorted(scores.columns):
        grp = df[df["receptor"] == rec].nlargest(args.top, "ligand_efficiency")[keep]
        out = grp.copy()
        out.insert(0, "receptor", rec)
        tops.append(out)
        print(f"\n  {rec}")
        print(grp.to_string(index=False, float_format=lambda v: f"{v:.3f}"))

    # ---- 5. per generator seed ------------------------------------------
    seeds = per_mol.index.to_series().str.extract(r"seed-(\d+)", expand=False)
    if seeds.notna().any():
        per_mol["seed"] = seeds.values
        print("\n--- Per generator seed ---")
        print(f"  {'seed':<6}{'n':>5}{'margin med':>12}{'margin IQR':>18}"
              f"{'>0':>5}{'LE med':>9}{'LE max':>9}")
        for sd, g in per_mol.groupby("seed"):
            m = g["margin_vs_human"].dropna()
            l = g["le_at_best_bacterial"].dropna()
            iqr = f"{m.quantile(.25):+.2f}..{m.quantile(.75):+.2f}"
            print(f"  {sd:<6}{len(g):>5}{m.median():>+12.2f}{iqr:>18}"
                  f"{(m > 0).sum():>5}{l.median():>9.2f}{l.max():>9.2f}")

        # Is the best-vs-worst seed gap real or sampling noise? Permutation
        # test on the median difference. numpy only: scipy is not a declared
        # dependency of this project, so it must not be relied on here.
        meds = per_mol.groupby("seed")["margin_vs_human"].median().sort_values()
        if len(meds) > 1:
            hi, lo = meds.index[-1], meds.index[0]
            a = per_mol.loc[per_mol["seed"] == hi, "margin_vs_human"].dropna().to_numpy()
            b = per_mol.loc[per_mol["seed"] == lo, "margin_vs_human"].dropna().to_numpy()
            obs = float(np.median(a) - np.median(b))
            pool = np.concatenate([a, b])
            rng = np.random.default_rng(0)
            n_a, n_shuf = len(a), 10000
            gaps = np.empty(n_shuf)
            for i in range(n_shuf):
                q = rng.permutation(pool)
                gaps[i] = np.median(q[:n_a]) - np.median(q[n_a:])
            pval = float((np.abs(gaps) >= abs(obs)).mean())
            verdict = "distinguishable" if pval < 0.05 else "NOT distinguishable"
            print(f"  seed {hi} vs seed {lo}: median gap {obs:+.2f} kcal/mol, "
                  f"permutation p={pval:.3f} ({n_shuf} shuffles) -> {verdict}")

    # ---- write-out -------------------------------------------------------
    per_mol.sort_values("margin_vs_human", ascending=False).to_csv(
        out_dir / "per_molecule.csv")
    sel.to_csv(out_dir / "shortlist.csv")
    pd.concat(tops, ignore_index=True).to_csv(out_dir / "top_by_le.csv",
                                              index=False)
    print(f"\nWrote {out_dir}/per_molecule.csv, shortlist.csv, top_by_le.csv")


if __name__ == "__main__":
    main()
