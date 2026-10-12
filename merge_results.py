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

    shards = sorted(dock_dir.glob("scores*.csv"))
    shards = [s for s in shards if "merged" not in s.name
              and "by_ligand" not in s.name]
    if not shards:
        raise SystemExit(f"No shard CSVs found in {dock_dir}")
    log.info("Merging %d shard files", len(shards))

    df = pd.concat([pd.read_csv(s) for s in shards], ignore_index=True)
    df = df.drop_duplicates(subset=["ligand", "receptor"], keep="first")
    log.info("%d rows after dedup", len(df))

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
    bact = [c for c in wide.columns if "dhfr" in c and "human" not in c]
    human = [c for c in wide.columns if "dhfr_human" in c]
    if bact and human:
        wide["dhfr_selectivity"] = wide[human[0]] - wide[bact[0]]
        log.info("Added dhfr_selectivity (positive = favours bacterial)")
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
