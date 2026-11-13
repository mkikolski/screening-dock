#!/usr/bin/env python
"""Checks that an interrupted dock stage resumes correctly.

A full run is tens of hours, so this is the one piece of logic that must not
be wrong: if resume silently re-docks everything a restart costs days, and if
it silently skips work the scores table comes out short.

Run: python test_checkpoint.py
"""
from __future__ import annotations

import json
import tempfile
from pathlib import Path

from dockpipe.common import VinaSpec, read_checkpoint
from dockpipe.docking import dock_many

RECEPTOR = {
    "name": "r",
    "receptor_pdbqt": "x",
    "box_center": [0.0, 0.0, 0.0],
    "box_size": [20.0, 20.0, 20.0],
}


def main() -> None:
    d = Path(tempfile.mkdtemp())
    ck = d / ".ckpt.jsonl"

    # Two complete records plus a line cut off mid-write, which is what a
    # SIGKILL on a flushed append-only file actually leaves behind.
    ck.write_text(
        json.dumps({"ligand": "lig1", "seconds": 1800.0, "status": "ok"}) + "\n"
        + json.dumps({"ligand": "lig2", "seconds": 1700.0, "status": "ok"}) + "\n"
        + '{"ligand": "lig3", "seco'
    )

    done = read_checkpoint(ck)
    assert set(done) == {"lig1", "lig2"}, done
    assert "lig3" not in done, "truncated record must be redone, not trusted"

    # Everything recorded: must not build a Pool or import vina at all.
    ligs = [str(d / f"{n}.pdbqt") for n in ("lig1", "lig2")]
    res = dock_many(ligs, RECEPTOR, d, VinaSpec(), n_workers=2, checkpoint=ck)
    assert {r["ligand"] for r in res} == {"lig1", "lig2"}, res

    assert read_checkpoint(d / "absent.jsonl") == {}

    # The ligand stage keys on "name", not "ligand".
    lk = d / ".ckpt_ligands.jsonl"
    lk.write_text(
        json.dumps({"name": "m1", "smiles": "CCO", "pdbqt_files": ["m1_0.pdbqt"],
                    "n_states": 1, "status": "ok", "error": None}) + "\n"
    )
    assert set(read_checkpoint(lk, key="name")) == {"m1"}
    # ...and the docking key must not silently match it
    assert read_checkpoint(lk, key="ligand") == {}

    # retry_failed: a recorded failure must come back as work to redo, or a
    # rerun after fixing the cause silently does nothing. This is how the 8
    # scrub failures in panel500_v1 would otherwise stay lost forever.
    fk = d / ".ckpt_failed.jsonl"
    fk.write_text(
        json.dumps({"name": "ok1", "status": "ok", "n_states": 1}) + "\n"
        + json.dumps({"name": "bad1", "status": "scrub_failed",
                      "error": "ETKDG conformer generation failed"}) + "\n"
        + json.dumps({"name": "bad2", "status": "no_states_written"}) + "\n"
    )
    assert set(read_checkpoint(fk, key="name")) == {"ok1", "bad1", "bad2"}, \
        "default must still treat a failure as done"
    kept = read_checkpoint(fk, key="name", retry_failed=True)
    assert set(kept) == {"ok1"}, kept
    assert "bad1" not in kept and "bad2" not in kept, \
        "retry_failed must drop every non-ok status, not just one of them"

    print("checkpoint resume OK (both stages, retry_failed)")


if __name__ == "__main__":
    main()
