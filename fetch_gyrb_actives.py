#!/usr/bin/env python3
"""Build a reference set of published GyrB ATP-site inhibitors from ChEMBL.

Why this exists: the panel's best GyrB scores (-8.75 S. aureus, -8.33 E. coli)
already beat both co-crystal ligands (08B -8.004, 4S4 -7.970). Either the
generated molecules outclass every pyrrolamide AstraZeneca made, or Vina is not
discriminating on this panel. Docking real inhibitors through the same pipeline
is how that question gets answered.

The filter that matters is mechanism, not potency. Fluoroquinolones and NBTIs
(gepotidacin) are annotated against DNA gyrase but bind the GyrA/DNA cleavage
interface, nowhere near the GyrB ATP pocket. Docking them into this box returns
numbers that mean nothing and would quietly corrupt the comparison, so they are
excluded structurally and the exclusion is recorded per compound.

Nothing unclassified ships by default. A compound matching neither the ATP-site
allowlist nor the off-site denylist goes to a review file. The chemotype column
in the output CSV is the auditable record of why each molecule is in the set.

Run on the cluster (needs rdkit; molscrub and meeko are not required):

    python fetch_gyrb_actives.py --selftest          # check the SMARTS first
    python fetch_gyrb_actives.py --out-dir /storage/users/kikolski/screening-dock

Writes gyrb_actives.smi (pipeline input), gyrb_actives.csv (provenance) and
gyrb_unclassified.csv (needs a human). configs/gyrb_knowns.yaml points at
the .smi.
"""
from __future__ import annotations

import argparse
import csv
import json
import logging
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

log = logging.getLogger("fetch_gyrb")

CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"

# Subunit B entries specifically, not the GyrA/B complex: a complex-level
# annotation cannot distinguish an ATP-site binder from a cleavage-complex
# poison, which is the only distinction this script exists to make.
TARGETS = {
    "P0AES6": "Escherichia coli",
    "P66936": "Staphylococcus aureus",
}

# Potency types worth keeping. MIC measures whole-cell growth and Inhibition %
# is not comparable across assays, so neither belongs in a score comparison.
POTENCY_TYPES = ("IC50", "Ki", "Kd")

# --- chemotype assignment -------------------------------------------------
#
# Allowlist = published GyrB ATP-site chemotypes. Denylist = binds gyrase
# somewhere else. The denylist is checked first, so anything matching both is
# excluded rather than included.
#
# These are deliberately loose. A chemotype tag claims "this is the series it
# belongs to", not a binding mode, and the CSV keeps the assay text so a human
# can disagree with any individual call.

ATP_SITE_SMARTS = {
    # Novobiocin, clorobiocin, coumermycin A1 -- the only registered GyrB
    # ATP-site drug is in this class.
    "aminocoumarin": "O=c1ccc2ccccc2o1",
    # AstraZeneca pyrrolamides including AZD5099, 3U2D's ligand 08B, and
    # Zidar's N-phenylpyrrolamides (4ZVI's 4S4).
    "pyrrolamide": "[nX3H1]1cccc1C(=O)[NX3]",
    "benzimidazole": "c1ccc2[nH]cnc2c1",
    "azaindole": "c1cc2cc[nH]c2nc1",
    "pyrazolopyridine": "c1cc2cn[nH]c2cn1",
    # Pyridyl and ethyl ureas. Broadest pattern here, so it is checked last.
    "urea": "[NX3][CX3](=O)[NX3]",
}

OFF_SITE_SMARTS = {
    # 4-oxo-1,4-dihydroquinoline-3-carboxylic acid: ciprofloxacin and every
    # other fluoroquinolone. Binds the cleavage complex, not GyrB.
    "fluoroquinolone": "O=c1c(cn(c2ccccc12))C(=O)[OX2H1,OX1-]",
    # NBTIs (gepotidacin and relatives): GyrA/DNA interface.
    "nbti_naphthyridone": "O=c1cc[nH]c2ncccc12",
}

# Compiled once, on first use. Recompiling per molecule would dominate runtime.
_PATTERNS: dict[str, list[tuple[str, object]]] = {}


def _compile(label: str, smarts_map: dict[str, str]) -> list[tuple[str, object]]:
    """Compile SMARTS once, failing loudly on a bad pattern.

    MolFromSmarts returns None rather than raising, and a None pattern matches
    nothing -- which would turn the denylist into a silent no-op and ship
    fluoroquinolones into a GyrB ATP box. Hence the explicit check.
    """
    if label not in _PATTERNS:
        from rdkit import Chem

        out = []
        for name, sma in smarts_map.items():
            patt = Chem.MolFromSmarts(sma)
            if patt is None:
                raise SystemExit(f"invalid SMARTS for {name!r}: {sma}")
            out.append((name, patt))
        _PATTERNS[label] = out
    return _PATTERNS[label]


def classify(mol) -> tuple[str, bool]:
    """Return (chemotype, keep). The denylist wins over the allowlist."""
    for name, patt in _compile("off", OFF_SITE_SMARTS):
        if mol.HasSubstructMatch(patt):
            return name, False
    for name, patt in _compile("atp", ATP_SITE_SMARTS):
        if mol.HasSubstructMatch(patt):
            return name, True
    return "unclassified", False


# --- ChEMBL ---------------------------------------------------------------

def _get(path: str, **params) -> dict:
    url = f"{CHEMBL}/{path}.json?{urllib.parse.urlencode(params)}"
    for attempt in range(1, 4):
        try:
            with urllib.request.urlopen(url, timeout=60) as fh:
                return json.loads(fh.read())
        except Exception as exc:  # noqa: BLE001 - network, retry then give up
            if attempt == 3:
                raise SystemExit(f"ChEMBL request failed: {url}\n{exc}") from exc
            log.warning("retry %d for %s (%s)", attempt, path, exc)
            time.sleep(2 * attempt)
    raise AssertionError("unreachable")


def target_ids(accession: str) -> list[str]:
    data = _get("target", target_components__accession=accession, limit=50)
    ids = [t["target_chembl_id"] for t in data.get("targets", [])
           if t.get("target_type") == "SINGLE PROTEIN"]
    if not ids:
        log.warning("no SINGLE PROTEIN target for %s", accession)
    return ids


def activities(target_chembl_id: str, limit: int) -> list[dict]:
    """Potency records carrying a structure, most potent first."""
    rows: list[dict] = []
    for stype in POTENCY_TYPES:
        offset = 0
        while len(rows) < limit:
            data = _get("activity", target_chembl_id=target_chembl_id,
                        standard_type=stype, standard_units="nM",
                        limit=1000, offset=offset)
            batch = data.get("activities", [])
            rows.extend(r for r in batch
                        if r.get("canonical_smiles") and r.get("standard_value"))
            if len(batch) < 1000:
                break
            offset += 1000
    rows.sort(key=lambda r: float(r["standard_value"]))
    return rows


# --- main -----------------------------------------------------------------

CSV_COLS = ["name", "chembl_id", "smiles", "inchikey", "heavy_atoms",
            "chemotype", "standard_type", "standard_value_nM",
            "target_uniprot", "target_organism", "assay_chembl_id",
            "assay_description"]


def build(out_dir: Path, max_per_target: int, include_unclassified: bool) -> None:
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")

    _compile("atp", ATP_SITE_SMARTS)    # fail before any network call
    _compile("off", OFF_SITE_SMARTS)

    kept: dict[str, dict] = {}          # inchikey -> row, dedup across targets
    rejected: list[dict] = []

    for accession, organism in TARGETS.items():
        for tid in target_ids(accession):
            acts = activities(tid, max_per_target)
            log.info("%s (%s, %s): %d potency records",
                     tid, accession, organism, len(acts))
            for a in acts:
                mol = Chem.MolFromSmiles(a["canonical_smiles"])
                if mol is None:
                    continue
                key = Chem.MolToInchiKey(mol)
                chemotype, keep = classify(mol)
                row = {
                    "name": "",
                    "chembl_id": a.get("molecule_chembl_id", ""),
                    "smiles": Chem.MolToSmiles(mol),
                    "inchikey": key,
                    "heavy_atoms": mol.GetNumHeavyAtoms(),
                    "chemotype": chemotype,
                    "standard_type": a.get("standard_type", ""),
                    "standard_value_nM": a.get("standard_value", ""),
                    "target_uniprot": accession,
                    "target_organism": organism,
                    "assay_chembl_id": a.get("assay_chembl_id", ""),
                    "assay_description": (a.get("assay_description") or "")[:200],
                }
                if keep or (include_unclassified and chemotype == "unclassified"):
                    prev = kept.get(key)
                    if prev is None or (float(row["standard_value_nM"])
                                        < float(prev["standard_value_nM"])):
                        kept[key] = row     # most potent record per structure
                else:
                    rejected.append(row)

    if not kept:
        raise SystemExit("no compounds survived the filter -- check the SMARTS")

    rows = sorted(kept.values(), key=lambda r: float(r["standard_value_nM"]))
    for i, r in enumerate(rows, 1):
        r["name"] = f"gyrb_known_{i:03d}"

    out_dir.mkdir(parents=True, exist_ok=True)
    smi_path = out_dir / "gyrb_actives.smi"
    with open(smi_path, "w", encoding="utf-8") as fh:
        for r in rows:
            fh.write(f"{r['smiles']} {r['name']}\n")

    with open(out_dir / "gyrb_actives.csv", "w", encoding="utf-8",
              newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=CSV_COLS)
        w.writeheader()
        w.writerows(rows)

    unclassified = [r for r in rejected if r["chemotype"] == "unclassified"]
    with open(out_dir / "gyrb_unclassified.csv", "w", encoding="utf-8",
              newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=CSV_COLS)
        w.writeheader()
        w.writerows(unclassified)

    by_type: dict[str, int] = {}
    for r in rows:
        by_type[r["chemotype"]] = by_type.get(r["chemotype"], 0) + 1
    log.info("kept %d compounds: %s", len(rows), by_type)
    log.info("excluded %d off-site, %d unclassified (see gyrb_unclassified.csv)",
             len(rejected) - len(unclassified), len(unclassified))
    log.info("wrote %s", smi_path)
    log.info("REVIEW gyrb_actives.csv before docking. The chemotype column is "
             "an automated guess, and a wrong mechanism call is worse than a "
             "missing compound.")


def selftest() -> None:
    """The one check that matters: does the mechanism filter actually work?

    A denylist that silently matches nothing is the failure mode here. It looks
    exactly like success and ships fluoroquinolones into a GyrB ATP box.
    """
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")

    cases = [
        # (SMILES, expected chemotype, expected keep)
        ("c1cc[nH]c1C(=O)NC", "pyrrolamide", True),
        ("O=c1ccc2ccccc2o1", "aminocoumarin", True),
        ("CCNC(=O)NCC", "urea", True),
        # 4-oxoquinoline-3-carboxylic acid core: must be EXCLUDED.
        ("OC(=O)c1cn(C)c2ccccc2c1=O", "fluoroquinolone", False),
        ("CCO", "unclassified", False),
    ]
    for smi, want_type, want_keep in cases:
        mol = Chem.MolFromSmiles(smi)
        assert mol is not None, f"test SMILES will not parse: {smi}"
        got_type, got_keep = classify(mol)
        assert got_type == want_type, f"{smi}: got {got_type}, want {want_type}"
        assert got_keep == want_keep, f"{smi}: keep={got_keep}, want {want_keep}"

    # A fluoroquinolone that also bears a urea must still be excluded.
    mol = Chem.MolFromSmiles("OC(=O)c1cn(CNC(=O)NC)c2ccccc2c1=O")
    assert mol is not None
    assert classify(mol) == ("fluoroquinolone", False), \
        "denylist must take precedence over the allowlist"

    print("selftest OK: chemotype filter excludes off-site binders")


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", default=".", type=Path)
    ap.add_argument("--max-per-target", type=int, default=2000,
                    help="cap on potency records pulled per ChEMBL target")
    ap.add_argument("--include-unclassified", action="store_true",
                    help="also keep compounds matching no known chemotype "
                         "(off by default: mechanism is the whole filter)")
    ap.add_argument("--selftest", action="store_true",
                    help="check the SMARTS filter and exit")
    ap.add_argument("--log-level", default="INFO")
    args = ap.parse_args()

    logging.basicConfig(level=args.log_level,
                        format="%(asctime)s %(levelname)s %(message)s")
    if args.selftest:
        selftest()
        return
    build(args.out_dir, args.max_per_target, args.include_unclassified)


if __name__ == "__main__":
    sys.exit(main())
