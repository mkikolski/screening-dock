#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import logging
import re
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger

from dockpipe.common import PipelineConfig, setup_logging, write_manifest

RDLogger.DisableLog("rdApp.*")
log = logging.getLogger("merge")

STATE_RE = re.compile(r"^(?P<base>.+)_(?P<state>\d+)$")

# Every scores*.csv this script writes, so the shard glob can skip them. Add to
# this list whenever a new one is written into dock_dir -- the glob picks up
# anything matching, and an output read back as an input is silent: the rows
# have no ligand or receptor column, so they merge as NaN rather than failing.
OUTPUT_NAMES = frozenset({"scores_merged.csv", "scores_by_ligand.csv",
                          "scores_wide.csv"})

# What a real shard must contain. Checked per file as a backstop to
# OUTPUT_NAMES, so a future output named scores_<something>.csv that nobody
# remembered to list is skipped rather than silently merged.
SHARD_COLUMNS = frozenset({"ligand", "receptor", "best_score", "status"})


def base_name(ligand: str) -> str:
    m = STATE_RE.match(ligand)
    return m.group("base") if m else ligand


def heavy_atoms(smiles: str) -> int | None:
    mol = Chem.MolFromSmiles(smiles)
    return mol.GetNumHeavyAtoms() if mol else None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()

    cfg = PipelineConfig.from_yaml(args.config)
    root = Path(cfg.output_root) / cfg.run_name
    dock_dir = root / "docking"
    setup_logging("INFO", root / "logs" / "merge.log")

    # Every CSV this script writes into dock_dir also matches scores*.csv, so
    # each one has to be excluded by name or the next run reads its own output
    # back in as an input. scores_wide.csv was added without being listed here
    # and did exactly that: its 492 rows carry no ligand or receptor column, so
    # they all arrived as NaN and drop_duplicates collapsed them into a single
    # row -- which is why the merge reported 2989 rows against the 2988 real
    # dockings, and why base_name() then died on a float.
    shards = [s for s in sorted(dock_dir.glob("scores*.csv"))
              if s.name not in OUTPUT_NAMES]
    if not shards:
        raise SystemExit(f"No shard CSVs found in {dock_dir}")
    log.info("Merging %d shard files", len(shards))

    frames = []
    for s in shards:
        d = pd.read_csv(s)
        absent = SHARD_COLUMNS - set(d.columns)
        if absent:
            log.warning("%s is not a shard CSV (no %s) -- skipping; it looks "
                        "like a previous output of this script", s.name,
                        ", ".join(sorted(absent)))
            continue
        d["_shard"] = s.name
        log.info("  %s: %d rows", s.name, len(d))
        frames.append(d)
    if not frames:
        raise SystemExit(
            f"None of the {len(shards)} scores*.csv in {dock_dir} has the "
            f"columns of a docking shard ({', '.join(sorted(SHARD_COLUMNS))}). "
            "Run the dock stage before merging.")
    df = pd.concat(frames, ignore_index=True)

    # A row with no ligand or receptor cannot be keyed on, and base_name()
    # raises TypeError on the NaN pandas leaves there -- which is how this was
    # found: panel500_v1 merged to 2989 rows where six receptors x 498 ligand
    # states is 2988, so exactly one row was junk. cmd_dock rewrites scores.csv
    # after every receptor, so a job ending mid-write leaves a truncated final
    # line.
    #
    # Dropped, but never quietly: the count and the file are logged, and more
    # than one bad row in a single shard aborts. One partial line is all an
    # interrupted write can leave behind; several mean the file is structurally
    # wrong, and merging it would under-report the panel with nothing else to
    # show for it.
    bad = df[df["ligand"].isna() | df["receptor"].isna()]
    if not bad.empty:
        per_shard = bad["_shard"].value_counts()
        for name, n in per_shard.items():
            log.warning("%s: %d row(s) with no ligand or receptor", name, n)
        if per_shard.max() > 1:
            raise SystemExit(
                f"{per_shard.idxmax()} has {per_shard.max()} unkeyable rows, "
                "more than an interrupted write can leave. Inspect it rather "
                "than merging: the result would be missing dockings with no "
                "other sign.")
        log.warning("dropping %d truncated row(s); %d remain",
                    len(bad), len(df) - len(bad))
        df = df.drop(bad.index)

    df = df.drop(columns=["_shard"])
    before = len(df)
    df = df.drop_duplicates(subset=["ligand", "receptor"], keep="first")
    log.info("%d rows after dedup (%d duplicate(s) removed)",
             len(df), before - len(df))

    lig_index = root / "ligands" / "ligand_index.json"
    smiles_by_name: dict[str, str] = {}
    if lig_index.exists():
        for rec in json.loads(lig_index.read_text()):
            smiles_by_name[rec["name"]] = rec["smiles"]

    df["ligand_base"] = df["ligand"].map(base_name)
    df["smiles"] = df["ligand_base"].map(smiles_by_name)

    merged = dock_dir / "scores_merged.csv"
    df.to_csv(merged, index=False)

    ok = df[df["status"] == "ok"].copy()
    if ok.empty:
        raise SystemExit("No successful docking results to summarize.")

    best = (ok.sort_values("best_score")
              .groupby(["ligand_base", "receptor"], as_index=False)
              .first())

    best["heavy_atoms"] = best["smiles"].map(
        lambda s: heavy_atoms(s) if isinstance(s, str) else None
    )
    best["ligand_efficiency"] = -best["best_score"] / best["heavy_atoms"]

    by_lig = dock_dir / "scores_by_ligand.csv"
    best.to_csv(by_lig, index=False)

    wide = best.pivot_table(index="ligand_base", columns="receptor",
                            values="best_score")
    bact = sorted(c for c in wide.columns if "dhfr" in c and "human" not in c)
    human = [c for c in wide.columns if "dhfr_human" in c]
    if bact and human:
        # One column per bacterial DHFR, so every pair stays visible and
        # comparable instead of collapsing into a single number.
        for b in bact:
            wide[f"selectivity_{b}"] = wide[human[0]] - wide[b]

        # Then one headline column. This used to be wide[bact[0]], picked by
        # sort order: adding a second bacterial DHFR silently repointed
        # dhfr_selectivity at the new receptor, with no error and no warning,
        # and every downstream number inherited it. Name it in the config.
        primary = cfg.primary_bacterial_dhfr
        if primary is None:
            if len(bact) > 1:
                raise SystemExit(
                    f"{len(bact)} bacterial DHFR receptors in the panel "
                    f"({', '.join(bact)}) but primary_bacterial_dhfr is unset. "
                    "Set it in the config: dhfr_selectivity must not be chosen "
                    "by column order."
                )
            primary = bact[0]
        elif primary not in bact:
            raise SystemExit(
                f"primary_bacterial_dhfr={primary!r} is not a bacterial DHFR "
                f"column in this run. Available: {', '.join(bact) or 'none'}"
            )

        wide["dhfr_selectivity"] = wide[f"selectivity_{primary}"]
        log.info("Added dhfr_selectivity vs %s (positive = favours bacterial); "
                 "per-receptor columns: %s",
                 primary, ", ".join(f"selectivity_{b}" for b in bact))
    wide.to_csv(dock_dir / "scores_wide.csv")

    redock_path = dock_dir / "redock_control.json"
    redock = json.loads(redock_path.read_text()) if redock_path.exists() else []

    summary = {
        "run_name": cfg.run_name,
        "n_ligand_states": int(len(df)),
        "n_molecules": int(best["ligand_base"].nunique()),
        "n_failed_states": int((df["status"] != "ok").sum()),
        "redock_control": redock,
        "per_receptor": {},
    }
    for recname, grp in best.groupby("receptor"):
        summary["per_receptor"][recname] = {
            "n_molecules": int(len(grp)),
            "score_mean": float(grp["best_score"].mean()),
            "score_min": float(grp["best_score"].min()),
            "score_q05": float(grp["best_score"].quantile(0.05)),
            "score_median": float(grp["best_score"].median()),
            "le_mean": float(grp["ligand_efficiency"].mean(skipna=True)),
            "le_max": float(grp["ligand_efficiency"].max(skipna=True)),
            "top10_by_le": grp.nlargest(10, "ligand_efficiency")[
                ["ligand_base", "best_score", "ligand_efficiency"]
            ].to_dict("records"),
        }

    (dock_dir / "summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    write_manifest(root / "manifests" / "merge.json",
                   {"stage": "merge", "summary": summary},
                   files=[merged, by_lig])

    print(f"\nMerged {len(df)} state rows -> {best['ligand_base'].nunique()} molecules")
    for name, s in summary["per_receptor"].items():
        print(f"  {name:<26} n={s['n_molecules']:<5} "
              f"best={s['score_min']:.2f}  median={s['score_median']:.2f}  "
              f"LEmax={s['le_max']:.3f}")


if __name__ == "__main__":
    main()
