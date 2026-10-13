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

from dockpipe.common import VinaSpec
from dockpipe.docking import dock_many, read_checkpoint

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

    print("checkpoint resume OK")


if __name__ == "__main__":
    main()
