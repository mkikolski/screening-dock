#!/usr/bin/env python3
"""Validate a config before submitting a job: fetch each PDB, split it, derive
the box. Catches the two silent setup errors -- a cofactor resname that is not
in the entry (dropped without an error) and a missing box reference -- in
seconds instead of after a queue wait.

Needs only pyyaml and numpy, so it runs on a login node.

    python preflight.py --config configs/panel500.yaml
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from dockpipe.common import PipelineConfig
from dockpipe.receptors import derive_box, fetch_pdb, split_pdb


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument(
        "--cache",
        help="PDB download dir (default: <output_root>/<run_name>/receptors/pdb_cache)",
    )
    args = ap.parse_args()

    cfg = PipelineConfig.from_yaml(args.config)
    cache = Path(
        args.cache
        or Path(cfg.output_root) / cfg.run_name / "receptors" / "pdb_cache"
    )

    problems: list[str] = []

    smi = Path(cfg.ligands_smi)
    if smi.is_file():
        n = sum(1 for line in smi.read_text().splitlines() if line.strip())
        print(f"ligands: {smi} ({n} lines)")
    else:
        print(f"ligands: {smi} -- NOT FOUND")
        problems.append(f"ligands_smi missing: {smi}")

    out_root = Path(cfg.output_root)
    if out_root.is_dir():
        print(f"output : {out_root}")
    else:
        print(f"output : {out_root} -- NOT FOUND")
        problems.append(f"output_root missing: {out_root}")
    print()

    for r in cfg.receptors:
        try:
            pdb = Path(r.local_pdb) if r.local_pdb else fetch_pdb(r.pdb_id, cache)
            rec, lig = split_pdb(
                pdb, r.keep_chains, r.keep_hetatms, r.redock_ligand_resname
            )
        except Exception as e:  # network, missing file, malformed entry
            print(f"{r.name:<22} FETCH/SPLIT FAILED: {e}")
            problems.append(f"{r.name}: {e}")
            continue

        het = sum(1 for line in rec if line.startswith("HETATM"))
        notes = []

        if not rec:
            notes.append("NO RECEPTOR ATOMS (check keep_chains)")
        # A cofactor named in the config but absent from the entry is dropped
        # silently by split_pdb, so the run would succeed with a broken site.
        if r.keep_hetatms and het == 0:
            notes.append(f"COFACTOR {','.join(r.keep_hetatms)} NOT IN ENTRY")

        try:
            center, size = derive_box(lig, r.box)
            box = (
                f"center=({center[0]:.2f},{center[1]:.2f},{center[2]:.2f}) "
                f"size={tuple(size)}"
            )
        except ValueError as e:
            box = "no box"
            notes.append(str(e))

        print(
            f"{r.name:<22} receptor={len(rec):5d} (het {het:3d})  "
            f"lig={len(lig):3d}  {box}"
        )
        for n in notes:
            print(f"{'':<22}   -> {n}")
            problems.append(f"{r.name}: {n}")

    print()
    if problems:
        print(f"{len(problems)} problem(s) found:")
        for p in problems:
            print(f"  - {p}")
        return 1
    print(f"OK: {len(cfg.receptors)} receptors ready")
    return 0


if __name__ == "__main__":
    sys.exit(main())
