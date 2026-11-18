#!/usr/bin/env python3
"""Interim read of a dock run, straight from the per-receptor checkpoints.

merge_results.py reads docking/scores*.csv, which cmd_dock rewrites only after
each receptor finishes. Mid-run that file is missing every receptor the loop
has not reached yet -- including dhfr_human_2w3a, without which merge skips its
whole selectivity block and reports no selectivity at all, with no warning.
The .ckpt_<receptor>.jsonl files are appended and flushed per docking, so they
are the only complete view while a job is still running.

Read-only, stdlib only: opens nothing but the checkpoints.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from statistics import median

STATE_RE = re.compile(r"^(?P<base>.+)_(?P<state>\d+)$")


def base_name(ligand: str) -> str:
    """Strip molscrub's _<state> suffix, as merge_results.base_name does."""
    m = STATE_RE.match(ligand)
    return m.group("base") if m else ligand


def load(dock_dir: Path) -> dict[str, dict[str, float]]:
    """{receptor: {molecule: best score over its states}}."""
    out: dict[str, dict[str, float]] = {}
    for ck in sorted(dock_dir.glob(".ckpt_*.jsonl")):
        scores: dict[str, float] = {}
        receptor = None
        for line in ck.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if not line:
                continue
            rec = json.loads(line)
            if not rec.get("receptor"):      # not a docking checkpoint
                break
            receptor = rec["receptor"]
            if rec.get("status") != "ok" or rec.get("best_score") is None:
                continue
            b, s = base_name(rec["ligand"]), float(rec["best_score"])
            if b not in scores or s < scores[b]:
                scores[b] = s
        if receptor:
            out[receptor] = scores
    return out


def report(data: dict[str, dict[str, float]], primary: str, control: str,
           human: str) -> None:
    print("=== per receptor, status ok ===")
    print(f"{'receptor':<26}{'n':>6}{'median':>9}{'best':>9}")
    for name in sorted(data):
        s = list(data[name].values())
        if not s:
            print(f"{name:<26}{0:>6}")
            continue
        print(f"{name:<26}{len(s):>6}{median(s):>9.2f}{min(s):>9.2f}")

    absent = [r for r in (primary, control, human) if not data.get(r)]
    if absent:
        print(f"\nNo DHFR comparison yet: nothing recorded for {', '.join(absent)}.")
    else:
        common = set(data[primary]) & set(data[control]) & set(data[human])
        print(f"\n=== DHFR selectivity, {len(common)} molecules scored at all three ===")
        print("selectivity = human - bacterial; positive favours the bacterial target")
        meds = {}
        for b in (control, primary):
            sel = [data[human][m] - data[b][m] for m in sorted(common)]
            meds[b] = median(sel)
            print(f"  vs {b:<24} median {meds[b]:+7.3f}   "
                  f"favouring bacterial {sum(v > 0 for v in sel):>4}/{len(sel)}")
        print(f"  paired shift, {primary} minus {control}: "
              f"{meds[primary] - meds[control]:+7.3f}")

    # The recovered ligands: present at some receptor the loop has already
    # passed, absent from the human receptor it has not reached yet.
    scored = set().union(*(set(v) for v in data.values())) if data else set()
    extra = sorted(scored - set(data.get(human, {})))
    if extra and data.get(human):
        order = sorted(data)
        print(f"\n=== {len(extra)} molecules not yet scored at {human} ===")
        print("molecule".ljust(30) + "".join(f"{r.rsplit('_', 1)[-1]:>9}"
                                             for r in order))
        for m in extra:
            cells = "".join(f"{data[r][m]:>9.2f}" if m in data[r] else f"{'--':>9}"
                            for r in order)
            print(f"{m:<30}{cells}")


def selftest() -> None:
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        p = Path(d)

        def ck(name: str, rows: list[dict]) -> None:
            (p / f".ckpt_{name}.jsonl").write_text(
                "".join(json.dumps(r) + "\n" for r in rows), encoding="utf-8")

        ck("dhfr_human_2w3a", [
            {"ligand": "m1_0", "receptor": "dhfr_human_2w3a",
             "status": "ok", "best_score": -8.0}])
        ck("dhfr_ecoli_6xg5", [
            {"ligand": "m1_0", "receptor": "dhfr_ecoli_6xg5",
             "status": "ok", "best_score": -7.9},
            {"ligand": "m1_1", "receptor": "dhfr_ecoli_6xg5",
             "status": "ok", "best_score": -6.0},     # worse state, must lose
            {"ligand": "m2_0", "receptor": "dhfr_ecoli_6xg5",
             "status": "ok", "best_score": -7.0},     # recovered ligand
            {"ligand": "m3_0", "receptor": "dhfr_ecoli_6xg5",
             "status": "dock_failed", "best_score": None}])
        ck("ligands", [{"ligand": "m1", "status": "ok"}])   # not a dock ckpt
        data = load(p)

    assert set(data) == {"dhfr_human_2w3a", "dhfr_ecoli_6xg5"}, data.keys()
    assert data["dhfr_ecoli_6xg5"]["m1"] == -7.9, "best state per molecule wins"
    assert "m3" not in data["dhfr_ecoli_6xg5"], "non-ok records dropped"
    assert set(data["dhfr_ecoli_6xg5"]) - set(data["dhfr_human_2w3a"]) == {"m2"}
    print("selftest ok")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dock-dir", default=".",
                    help="directory holding .ckpt_<receptor>.jsonl")
    ap.add_argument("--primary", default="dhfr_ecoli_6xg5")
    ap.add_argument("--control", default="dhfr_ecoli_7mym")
    ap.add_argument("--human", default="dhfr_human_2w3a")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()
    if args.selftest:
        selftest()
        return
    data = load(Path(args.dock_dir))
    if not data:
        raise SystemExit(f"No docking checkpoints in {args.dock_dir}")
    report(data, args.primary, args.control, args.human)


if __name__ == "__main__":
    main()
