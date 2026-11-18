#!/usr/bin/env python3
"""Are two runs' scores comparable? Compare their receptor manifests.

Scores from two runs can only be read against each other if the receptors
were prepared identically. Two things have to match per receptor: the
SHA-256 of the receptor PDBQT, and the docking box. An identical PDBQT in a
shifted or resized box gives scores that look comparable and are not.

    python3 check_comparable.py a/receptors.json b/receptors.json

Exits non-zero if anything differs, so it can gate a pipeline.

The checksum keys are absolute paths that embed the run name, so they never
match across runs and are joined through each receptor's own record instead.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def load(path: str | Path) -> dict[str, dict]:
    """{receptor name: {sha256, box_center, box_size}} from one manifest."""
    doc = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(doc, dict) or "receptors" not in doc:
        raise SystemExit(
            f"{path}: not a prep-receptors manifest "
            f"(top-level keys: {sorted(doc)[:8] if isinstance(doc, dict) else type(doc).__name__})")
    sums = doc.get("checksums", {})
    out = {}
    for r in doc["receptors"]:
        pdbqt = r.get("receptor_pdbqt")
        out[r["name"]] = {
            "sha256": sums.get(str(pdbqt)) if pdbqt else None,
            "box_center": r.get("box_center"),
            "box_size": r.get("box_size"),
        }
    return out


def compare(a: dict[str, dict], b: dict[str, dict],
            names: tuple[str, str] = ("A", "B")) -> list[str]:
    """Human-readable lines; a leading '*' marks a blocking difference."""
    lines, shared = [], sorted(set(a) & set(b))
    for only, label in ((set(a) - set(b), names[0]), (set(b) - set(a), names[1])):
        for n in sorted(only):
            lines.append(f"  {n:<26} present only in {label} -- not compared")
    for n in shared:
        diffs = [k for k in ("sha256", "box_center", "box_size") if a[n][k] != b[n][k]]
        if not diffs:
            lines.append(f"  {n:<26} IDENTICAL")
        elif a[n]["sha256"] is None or b[n]["sha256"] is None:
            lines.append(f"* {n:<26} NO CHECKSUM RECORDED -- cannot verify")
        else:
            lines.append(f"* {n:<26} DIFFERS in {', '.join(diffs)}")
    if not shared:
        lines.append("* no receptor appears in both manifests")
    return lines


def selftest() -> None:
    base = {
        "receptors": [
            {"name": "gyrb_saureus_3u2d", "receptor_pdbqt": "/a/3u2d.pdbqt",
             "box_center": [1.0, 2.0, 3.0], "box_size": [22.0, 22.0, 22.0]},
            {"name": "gyrb_ecoli_4zvi", "receptor_pdbqt": "/a/4zvi.pdbqt",
             "box_center": [4.0, 5.0, 6.0], "box_size": [22.0, 22.0, 22.0]},
        ],
        "checksums": {"/a/3u2d.pdbqt": "aaa", "/a/4zvi.pdbqt": "bbb"},
    }
    # Same receptors, different run: paths differ, content does not.
    other = json.loads(json.dumps(base))
    for r in other["receptors"]:
        r["receptor_pdbqt"] = r["receptor_pdbqt"].replace("/a/", "/b/")
    other["checksums"] = {"/b/3u2d.pdbqt": "aaa", "/b/4zvi.pdbqt": "bbb"}

    import tempfile
    with tempfile.TemporaryDirectory() as d:
        p = Path(d)
        (p / "a.json").write_text(json.dumps(base))
        (p / "b.json").write_text(json.dumps(other))
        a, b = load(p / "a.json"), load(p / "b.json")
        out = compare(a, b)
        assert all("IDENTICAL" in l for l in out), out
        assert not any(l.startswith("*") for l in out), out

        # A shifted box with an identical file must still block.
        shifted = json.loads(json.dumps(other))
        shifted["receptors"][0]["box_center"] = [1.0, 2.0, 9.9]
        (p / "c.json").write_text(json.dumps(shifted))
        out = compare(a, load(p / "c.json"))
        assert any(l.startswith("*") and "box_center" in l for l in out), out

        # A changed file must block.
        changed = json.loads(json.dumps(other))
        changed["checksums"]["/b/3u2d.pdbqt"] = "zzz"
        (p / "d.json").write_text(json.dumps(changed))
        out = compare(a, load(p / "d.json"))
        assert any(l.startswith("*") and "sha256" in l for l in out), out
    print("selftest ok")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("manifests", nargs="*", help="two receptors.json paths")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()
    if args.selftest:
        selftest()
        return
    if len(args.manifests) != 2:
        ap.error("give exactly two receptors.json paths")
    a_path, b_path = args.manifests
    lines = compare(load(a_path), load(b_path), (a_path, b_path))
    print("\n".join(lines))
    if any(l.startswith("*") for l in lines):
        sys.exit("\nNOT comparable: scores from these runs must not be "
                 "read against each other until this is resolved.")
    print("\nComparable: receptors and boxes are identical.")


if __name__ == "__main__":
    main()
