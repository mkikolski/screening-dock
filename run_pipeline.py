#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

from dockpipe.common import PipelineConfig, setup_logging, write_manifest
from dockpipe import ligands as lig_mod
from dockpipe import receptors as rec_mod
from dockpipe import docking as dock_mod

log = logging.getLogger("run_pipeline")


def _paths(cfg: PipelineConfig) -> dict[str, Path]:
    root = Path(cfg.output_root) / cfg.run_name
    return {
        "root": root,
        "ligands": root / "ligands",
        "receptors": root / "receptors",
        "docking": root / "docking",
        "manifests": root / "manifests",
        "logs": root / "logs",
    }


def _load_receptor_index(p: dict[str, Path]) -> list[dict]:
    idx = p["receptors"] / "receptor_index.json"
    if not idx.exists():
        sys.exit(f"Missing {idx}. Run prep-receptors first.")
    recs = json.loads(idx.read_text())
    ok = [r for r in recs if r.get("status") != "failed"]
    if not ok:
        sys.exit("No successfully prepared receptors.")
    return ok


def _load_ligand_states(p: dict[str, Path]) -> list[str]:
    idx = p["ligands"] / "ligand_index.json"
    if not idx.exists():
        sys.exit(f"Missing {idx}. Run prep-ligands first.")
    states: list[str] = []
    for rec in json.loads(idx.read_text()):
        states.extend(rec.get("pdbqt_files", []))
    if not states:
        sys.exit("No ligand PDBQT states found.")
    return states


def cmd_prep_receptors(cfg: PipelineConfig, args) -> None:
    p = _paths(cfg)
    records = rec_mod.prepare_receptors(
        cfg.receptors, p["receptors"], use_pdb_redo=args.pdb_redo
    )
    write_manifest(
        p["manifests"] / "receptors.json",
        {"stage": "prep-receptors", "config": cfg.to_dict(), "receptors": records},
        files=[r["receptor_pdbqt"] for r in records if r.get("receptor_pdbqt")],
    )
    for r in records:
        status = "FAILED" if r.get("status") == "failed" else "ok"
        log.info("  %-24s %s", r["name"], status)


def cmd_prep_ligands(cfg: PipelineConfig, args) -> None:
    p = _paths(cfg)
    smi = args.ligands or cfg.ligands_smi
    records = lig_mod.process_ligands(
        smi,
        p["ligands"],
        ph=cfg.ligand_ph,
        skip_tautomers=cfg.skip_tautomers,
        n_attempts=cfg.n_scrub_attempts,
        n_workers=args.workers or cfg.n_workers,
    )
    write_manifest(
        p["manifests"] / "ligands.json",
        {
            "stage": "prep-ligands",
            "ligands_smi": str(smi),
            "n_input": len(records),
            "n_ok": sum(r["status"] == "ok" for r in records),
            "n_states": sum(r["n_states"] for r in records),
            "failures": [r for r in records if r["status"] != "ok"],
        },
        files=[smi],
    )


def cmd_redock(cfg: PipelineConfig, args) -> None:
    p = _paths(cfg)
    receptors = _load_receptor_index(p)
    results = []
    for r in receptors:
        if args.receptor and r["name"] != args.receptor:
            continue
        results.append(dock_mod.redock_control(r, p["docking"], cfg.vina))

    out = p["docking"] / "redock_control.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(results, indent=2), encoding="utf-8")
    write_manifest(
        p["manifests"] / "redock.json",
        {"stage": "redock", "vina": cfg.to_dict()["vina"], "results": results},
    )

    print("\n--- Redocking control (threshold 2.0 A) ---")
    for r in results:
        rmsd = r.get("rmsd_reported")
        print(
            f"  {r['receptor']:<24} "
            f"RMSD={rmsd if rmsd is None else f'{rmsd:.2f}':>6}  "
            f"{r.get('verdict') or r.get('status')}"
        )

    failed = [r for r in results if r.get("verdict") == "FAIL"]
    if failed and not args.allow_fail:
        sys.exit(
            f"\n{len(failed)} receptor(s) failed the redocking gate. "
            "Fix the setup or pass --allow-fail to proceed anyway."
        )


def cmd_dock(cfg: PipelineConfig, args) -> None:
    p = _paths(cfg)
    receptors = _load_receptor_index(p)
    states = _load_ligand_states(p)

    if args.receptor:
        receptors = [r for r in receptors if r["name"] == args.receptor]
        if not receptors:
            sys.exit(f"Receptor '{args.receptor}' not found in index.")

    if args.shard_count > 1:
        states = states[args.shard_index :: args.shard_count]
        log.info(
            "Shard %d/%d -> %d ligand states",
            args.shard_index,
            args.shard_count,
            len(states),
        )

    all_results = []
    for r in receptors:
        res = dock_mod.dock_many(
            states,
            r,
            p["docking"],
            cfg.vina,
            n_workers=args.workers or cfg.n_workers,
        )
        all_results.extend(res)

    suffix = f"_shard{args.shard_index:03d}" if args.shard_count > 1 else ""
    tag = f"_{args.receptor}" if args.receptor else ""
    csv_path = p["docking"] / f"scores{tag}{suffix}.csv"
    dock_mod.write_scores_csv(all_results, csv_path)

    write_manifest(
        p["manifests"] / f"dock{tag}{suffix}.json",
        {
            "stage": "dock",
            "vina": cfg.to_dict()["vina"],
            "receptors": [r["name"] for r in receptors],
            "n_ligand_states": len(states),
            "n_results": len(all_results),
            "n_failed": sum(r["status"] != "ok" for r in all_results),
            "shard": {"index": args.shard_index, "count": args.shard_count},
        },
        files=[csv_path],
    )


def cmd_all(cfg: PipelineConfig, args) -> None:
    cmd_prep_receptors(cfg, args)
    cmd_redock(cfg, args)
    cmd_prep_ligands(cfg, args)
    cmd_dock(cfg, args)


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "command", choices=["prep-ligands", "prep-receptors", "redock", "dock", "all"]
    )
    ap.add_argument("--config", required=True, help="YAML pipeline config")
    ap.add_argument("--ligands", help="override ligands_smi from config")
    ap.add_argument("--receptor", help="restrict to one receptor by name")
    ap.add_argument("--workers", type=int, help="override n_workers")
    ap.add_argument(
        "--pdb-redo",
        action="store_true",
        help="fetch PDB-REDO re-refined models instead of RCSB",
    )
    ap.add_argument(
        "--allow-fail",
        action="store_true",
        help="continue even if the redocking gate fails",
    )
    ap.add_argument("--shard-index", type=int, default=0)
    ap.add_argument("--shard-count", type=int, default=1)
    ap.add_argument("--log-level", default="INFO")
    args = ap.parse_args()

    cfg = PipelineConfig.from_yaml(args.config)
    p = _paths(cfg)
    p["logs"].mkdir(parents=True, exist_ok=True)
    tag = f"{args.command}{'_' + args.receptor if args.receptor else ''}"
    if args.shard_count > 1:
        tag += f"_shard{args.shard_index:03d}"
    setup_logging(args.log_level, p["logs"] / f"{tag}.log")

    log.info("=== %s | run=%s ===", args.command, cfg.run_name)
    {
        "prep-ligands": cmd_prep_ligands,
        "prep-receptors": cmd_prep_receptors,
        "redock": cmd_redock,
        "dock": cmd_dock,
        "all": cmd_all,
    }[args.command](cfg, args)
    log.info("=== %s done ===", args.command)


if __name__ == "__main__":
    main()
