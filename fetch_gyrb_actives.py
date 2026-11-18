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
import hashlib
import json
import logging
import random
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

log = logging.getLogger("fetch_gyrb")

CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"

# Measured against the live service on 2026-11-17: it 500s in bursts lasting a
# minute or two, fast-failing in ~5s, then recovers and answers in ~0.3s. Its
# own /status.json endpoint 500s as well, so there is nothing to health-check
# against. Hence patient retries rather than few, and a cache so an outage
# part-way through does not discard the pages already fetched.
#
# Page size matters as much as the retries. limit=1000 (the documented maximum)
# times out; 100 and 50 both 500 even outside an outage window; 25 works.
PAGE = 25
TIMEOUT = 120
ATTEMPTS = 8
BACKOFF = (5, 10, 20, 30, 45, 60, 60)   # between attempts, ~4.5 min of patience

# Field selection. A full activity record carries large unused text, and asking
# only for what the output CSV needs keeps each page inside the timeout.
ACTIVITY_FIELDS = ("molecule_chembl_id,canonical_smiles,standard_type,"
                   "standard_relation,standard_value,standard_units,"
                   "assay_chembl_id,assay_description")
TARGET_FIELDS = "target_chembl_id,pref_name,target_type,organism"


class FetchError(RuntimeError):
    """A request that stayed broken through every retry.

    Deliberately not SystemExit. One dead (target, potency type) pair must not
    throw away the pages that already succeeded, which is what the first
    version did when Ki started 500ing after IC50 had been collected.
    """

# Subunit B entries specifically, not the GyrA/B complex: a complex-level
# annotation cannot distinguish an ATP-site binder from a cleavage-complex
# poison, which is the only distinction this script exists to make.
#
# Pinned as ChEMBL target IDs, not resolved from UniProt accessions. P66936 is
# the accession BindingDB shows for S. aureus GyrB, but it maps to ZERO ChEMBL
# targets -- an accession lookup therefore returned an empty S. aureus half
# while looking like it had worked. Verified 2026-11-17 and re-derivable with:
#
#   .../target.json?pref_name__icontains=gyrase%20subunit%20B&limit=25
#
# which lists six, of which these two are the organisms this panel docks.
# The accession is recorded for provenance only -- nothing is looked up by it.
TARGETS = {
    "CHEMBL1826": ("Escherichia coli", "P0AES6"),
    "CHEMBL1921666": ("Staphylococcus aureus", "P66936"),
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

# The strongest mechanism evidence is the assay, not the scaffold. 540 of the
# records pulled on 2026-11-17 are dominated by descriptions like "Inhibition
# of ATPase activity of Staphylococcus aureus GyrB" -- that measures the ATP
# site directly, and a fluoroquinolone cannot inhibit GyrB ATPase at all.
#
# So the assay gates inclusion and the chemotype tag is descriptive. An
# allowlist of hand-written scaffolds would reject most of this set: the modern
# GyrB series (pyrimidinoindoles, thiazolopyridines, pyrazolo-oxadiazolones)
# match none of the classical patterns above. --strict-chemotype restores
# allowlist-only behaviour.
# [\s-]* not \s*: ChEMBL writes "DNA gyrase-B" as well as "DNA gyrase B", and
# a bare \s* silently dropped the hyphenated ones.
ATP_SITE_ASSAY = re.compile(r"atpase|gyrase[\s-]*b|gyr[\s-]*b|subunit[\s-]*b",
                            re.IGNORECASE)

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
    """Return (chemotype, matches_atp_allowlist). The denylist wins."""
    for name, patt in _compile("off", OFF_SITE_SMARTS):
        if mol.HasSubstructMatch(patt):
            return name, False
    for name, patt in _compile("atp", ATP_SITE_SMARTS):
        if mol.HasSubstructMatch(patt):
            return name, True
    return "unclassified", False


def decide(mol, assay: str, strict: bool) -> tuple[str, str, bool]:
    """Return (chemotype, why, keep).

    Order of authority: an off-site scaffold excludes outright, a known
    ATP-site scaffold includes, and otherwise the assay decides. `why` goes
    into the CSV so every inclusion can be argued with individually.
    """
    chemotype, on_allowlist = classify(mol)
    if chemotype in OFF_SITE_SMARTS:
        return chemotype, "excluded: off-site chemotype", False
    if on_allowlist:
        return chemotype, "known ATP-site chemotype", True
    if strict:
        return chemotype, "excluded: no chemotype match (--strict-chemotype)", False
    if ATP_SITE_ASSAY.search(assay or ""):
        return chemotype, "GyrB/ATPase assay", True
    return chemotype, "excluded: no chemotype match, assay not ATP-specific", False


# --- ChEMBL ---------------------------------------------------------------

_CACHE: Path | None = None


def _get(path: str, **params) -> dict:
    """One ChEMBL page, cached on disk and retried through an outage window.

    Raises FetchError once patience runs out, so the caller can drop this one
    page and keep everything else.
    """
    query = urllib.parse.urlencode(params)
    url = f"{CHEMBL}/{path}.json?{query}"

    cached = None
    if _CACHE is not None:
        key = hashlib.sha256(f"{path}?{query}".encode()).hexdigest()[:32]
        cached = _CACHE / f"{path}_{key}.json"
        if cached.exists():
            try:
                return json.loads(cached.read_text(encoding="utf-8"))
            except json.JSONDecodeError:
                cached.unlink()      # truncated by an interrupt; refetch

    last = None
    for attempt in range(1, ATTEMPTS + 1):
        try:
            with urllib.request.urlopen(url, timeout=TIMEOUT) as fh:
                body = fh.read()
            data = json.loads(body)
            if cached is not None:
                cached.write_bytes(body)
            return data
        except Exception as exc:  # noqa: BLE001 - any network/server fault retries
            last = exc
            if attempt == ATTEMPTS:
                break
            # Jitter so a rerun does not line its retries up with the last one.
            pause = BACKOFF[min(attempt - 1, len(BACKOFF) - 1)]
            pause += random.uniform(0, pause * 0.25)
            log.warning("attempt %d/%d failed for %s (%s); waiting %.0fs",
                        attempt, ATTEMPTS, path, exc, pause)
            time.sleep(pause)

    raise FetchError(f"{url}\n  gave up after {ATTEMPTS} attempts: {last}")


def usable_records(batch: list[dict]) -> tuple[list[dict], int]:
    """(records with an exact value, count of censored ones dropped).

    A censored value -- "<0.05 nM", ">10 uM" -- is a bound, not a
    measurement, and must not be ranked against exact ones. It is typically a
    reference compound quoted below an assay's limit of quantitation, and
    because it sorts to the very top of an ascending potency ranking it would
    otherwise lead both the potency mode and the chemotype round-robin while
    carrying no actual number.
    """
    ok: list[dict] = []
    censored = 0
    for r in batch:
        if not (r.get("canonical_smiles") and r.get("standard_value")):
            continue
        if (r.get("standard_relation") or "=").strip() == "=":
            ok.append(r)
        else:
            censored += 1
    return ok, censored


def activities(target_chembl_id: str, limit: int, page: int = PAGE) -> list[dict]:
    """Potency records carrying a structure, most potent first.

    order_by=standard_value makes ChEMBL sort server-side, so the first pages
    are already the most potent and paging can stop at `limit` instead of
    pulling everything and sorting here.

    A potency type that stays broken is logged and skipped. Losing Ki is worth
    far less than losing the IC50 pages already in hand.
    """
    rows: list[dict] = []
    for stype in POTENCY_TYPES:
        got = 0
        offset = 0
        while got < limit:
            try:
                data = _get("activity", target_chembl_id=target_chembl_id,
                            standard_type=stype, standard_units="nM",
                            only=ACTIVITY_FIELDS, order_by="standard_value",
                            limit=page, offset=offset)
            except FetchError as exc:
                log.error("skipping %s/%s at offset %d: %s",
                          target_chembl_id, stype, offset, exc)
                break
            batch = data.get("activities", [])
            usable, censored = usable_records(batch)
            if censored:
                log.info("  %s/%s: %d censored record(s) skipped at offset %d",
                         target_chembl_id, stype, censored, offset)
            rows.extend(usable)
            got += len(usable)
            total = data.get("page_meta", {}).get("total_count")
            if offset == 0 and total is not None:
                log.info("  %s/%s: %s records available", target_chembl_id,
                         stype, total)
            if len(batch) < page:
                break
            offset += page
    rows.sort(key=lambda r: float(r["standard_value"]))
    return rows


# --- main -----------------------------------------------------------------

CSV_COLS = ["name", "chembl_id", "smiles", "inchikey", "heavy_atoms",
            "chemotype", "mechanism_evidence", "standard_type",
            "standard_relation", "standard_value_nM",
            "target_uniprot", "target_organism",
            "assay_chembl_id", "assay_description"]


def select(rows: list[dict], limit: int, mode: str) -> list[dict]:
    """Choose the reference set from the filtered pool.

    At a small limit this choice is the experiment, not a detail. ChEMBL's
    sub-nM end for GyrB is two or three patent series, so a plain potency
    ranking of the ~450-compound pool gives a top 10 that is five analogues of
    one ethylurea-thiazole-oxadiazolone plus two of one azaindole -- seven of
    ten compounds from two series, every one between 0.05 and 4 nM, while the
    pool spans to 1 mM. Such a set answers "do our molecules beat
    best-in-class" and nothing about chemotype or potency breadth.

    "diverse" (default) round-robins across chemotypes, most potent first
    within each, and skips a compound whose Murcko scaffold is already taken.
    "potency" is the naive ranking, kept for when the limit is large enough
    that selection stops mattering.
    """
    rows = sorted(rows, key=lambda r: float(r["standard_value_nM"]))
    if not limit or len(rows) <= limit:
        return rows
    if mode == "potency":
        return rows[:limit]

    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold

    def scaffold(smiles: str) -> str:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return smiles
        try:
            return MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
        except Exception:  # noqa: BLE001 - odd valences; fall back to the whole molecule
            return smiles

    buckets: dict[str, list[dict]] = {}
    for r in rows:
        buckets.setdefault(r["chemotype"], []).append(r)
    # Visit chemotypes best-first, so a limit smaller than the number of
    # chemotypes still picks up the strongest ones.
    order = sorted(buckets, key=lambda k: float(buckets[k][0]["standard_value_nM"]))

    out: list[dict] = []
    seen: set[str] = set()
    while len(out) < limit:
        progressed = False
        for key in order:
            if len(out) >= limit:
                break
            while buckets[key]:
                cand = buckets[key].pop(0)
                scaf = scaffold(cand["smiles"])
                if scaf in seen:
                    continue            # analogue of something already taken
                seen.add(scaf)
                out.append(cand)
                progressed = True
                break
        if not progressed:
            break

    if len(out) < limit:
        # Scaffold dedup exhausted the pool. Top up by potency rather than
        # silently returning a shorter set than asked for.
        taken = {r["inchikey"] for r in out}
        for r in rows:
            if len(out) >= limit:
                break
            if r["inchikey"] not in taken:
                out.append(r)
        log.info("scaffold dedup ran out of distinct series; topped up to %d "
                 "by potency", len(out))

    return sorted(out, key=lambda r: float(r["standard_value_nM"]))


def build(out_dir: Path, max_per_target: int, max_compounds: int,
          strict_chemotype: bool, select_mode: str = "diverse",
          page: int = PAGE) -> None:
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")

    _compile("atp", ATP_SITE_SMARTS)    # fail before any network call
    _compile("off", OFF_SITE_SMARTS)

    kept: dict[str, dict] = {}          # inchikey -> row, dedup across targets
    rejected: list[dict] = []

    failed_targets: list[str] = []
    for tid, (organism, accession) in TARGETS.items():
        log.info("%s -- %s GyrB (UniProt %s)", tid, organism, accession)
        acts = activities(tid, max_per_target, page)
        log.info("%s: %d usable potency records", tid, len(acts))
        if not acts:
            # Either an outage or a target that no longer carries data. Both
            # need saying: an empty half of the panel must not pass silently.
            failed_targets.append(f"{tid} ({organism})")
        for a in acts:
            mol = Chem.MolFromSmiles(a["canonical_smiles"])
            if mol is None:
                continue
            key = Chem.MolToInchiKey(mol)
            assay = a.get("assay_description") or ""
            chemotype, why, keep = decide(mol, assay, strict_chemotype)
            row = {
                "name": "",
                "chembl_id": a.get("molecule_chembl_id", ""),
                "smiles": Chem.MolToSmiles(mol),
                "inchikey": key,
                "heavy_atoms": mol.GetNumHeavyAtoms(),
                "chemotype": chemotype,
                "mechanism_evidence": why,
                "standard_type": a.get("standard_type", ""),
                "standard_relation": (a.get("standard_relation") or "=").strip(),
                "standard_value_nM": a.get("standard_value", ""),
                "target_uniprot": accession,
                "target_organism": organism,
                "assay_chembl_id": a.get("assay_chembl_id", ""),
                "assay_description": (a.get("assay_description") or "")[:200],
            }
            if keep:
                prev = kept.get(key)
                if prev is None or (float(row["standard_value_nM"])
                                    < float(prev["standard_value_nM"])):
                    kept[key] = row     # most potent record per structure
            else:
                rejected.append(row)

    if not kept:
        raise SystemExit(
            "no compounds survived the filter. If the log above shows ChEMBL "
            "errors this is an outage, not a filter problem -- rerun with the "
            "same --cache-dir and the fetched pages are reused. If the log is "
            "clean, check the SMARTS with --selftest."
        )
    if failed_targets:
        log.warning("INCOMPLETE: no data for %s. The set is usable but not "
                    "the full picture -- rerun with the same --cache-dir to "
                    "fill the gaps.", ", ".join(failed_targets))

    pool = sorted(kept.values(), key=lambda r: float(r["standard_value_nM"]))
    rows = select(pool, max_compounds, select_mode)
    if len(rows) < len(pool):
        log.info("selected %d of %d compounds by '%s'", len(rows), len(pool),
                 select_mode)

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

    # Everything rejected, with the reason, so an over-strict filter is
    # visible rather than invisible.
    with open(out_dir / "gyrb_excluded.csv", "w", encoding="utf-8",
              newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=CSV_COLS)
        w.writeheader()
        w.writerows(sorted(rejected,
                           key=lambda r: float(r["standard_value_nM"])))

    def tally(seq, field):
        out: dict[str, int] = {}
        for r in seq:
            out[r[field]] = out.get(r[field], 0) + 1
        return dict(sorted(out.items(), key=lambda kv: -kv[1]))

    log.info("kept %d compounds", len(rows))
    log.info("  by chemotype: %s", tally(rows, "chemotype"))
    log.info("  by evidence : %s", tally(rows, "mechanism_evidence"))
    log.info("excluded %d: %s", len(rejected),
             tally(rejected, "mechanism_evidence"))
    log.info("potency span: %.3g to %.3g nM",
             float(rows[0]["standard_value_nM"]),
             float(rows[-1]["standard_value_nM"]))
    log.info("wrote %s and gyrb_actives.csv (+ gyrb_excluded.csv)", smi_path)
    log.info("REVIEW gyrb_actives.csv before docking. mechanism_evidence says "
             "why each compound is in; most will read 'GyrB/ATPase assay', "
             "which is the assay measuring the ATP site directly. Skim "
             "gyrb_excluded.csv too -- a good compound dropped there means "
             "the filter is too strict.")


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

    # decide(): the assay now gates inclusion, so these are the cases that
    # determine what actually gets docked.
    fq = Chem.MolFromSmiles("OC(=O)c1cn(C)c2ccccc2c1=O")
    atpase = "Inhibition of ATPase activity of Staphylococcus aureus GyrB"
    unknown_scaffold = Chem.MolFromSmiles(
        "CNc1cc(F)cc2c1[nH]c1nc(Oc3cnc(C)nc3)nc(N3CCNCC3)c12")  # modern series

    # An ATPase assay must NOT rescue an off-site chemotype.
    assert decide(fq, atpase, False)[2] is False, \
        "a fluoroquinolone must stay excluded whatever the assay says"
    # A scaffold outside the allowlist is kept on assay evidence...
    assert decide(unknown_scaffold, atpase, False)[2] is True, \
        "GyrB ATPase assay evidence must admit an unlisted scaffold"
    # ...but not when the assay is silent about mechanism...
    assert decide(unknown_scaffold, "Antibacterial activity against E. coli",
                  False)[2] is False, "a whole-cell MIC assay is not evidence"
    # ...and not under --strict-chemotype.
    assert decide(unknown_scaffold, atpase, True)[2] is False, \
        "--strict-chemotype must require an allowlist match"

    # select(): at a small limit this is the experiment. Two analogues sharing
    # a Murcko scaffold must not both be chosen while a distinct series waits.
    pool = [
        {"smiles": "CCNC(=O)Nc1ccc(-c2ccccc2)cc1", "chemotype": "urea",
         "standard_value_nM": "1.0", "inchikey": "A"},
        {"smiles": "CCNC(=O)Nc1ccc(-c2ccccc2C)cc1", "chemotype": "urea",
         "standard_value_nM": "2.0", "inchikey": "B"},   # analogue of A
        {"smiles": "O=c1ccc2ccccc2o1", "chemotype": "aminocoumarin",
         "standard_value_nM": "500.0", "inchikey": "C"},
    ]
    batch = [
        {"canonical_smiles": "CC", "standard_value": "1.0", "standard_relation": "="},
        {"canonical_smiles": "CCC", "standard_value": "0.05", "standard_relation": "<"},
        {"canonical_smiles": "CCCC", "standard_value": "9.0", "standard_relation": ">"},
        {"canonical_smiles": "CCCCC", "standard_value": "2.0"},        # absent = exact
        {"canonical_smiles": "CCCCCC", "standard_value": None},        # no value
        {"canonical_smiles": None, "standard_value": "3.0", "standard_relation": "="},
    ]
    ok, censored = usable_records(batch)
    assert [r["standard_value"] for r in ok] == ["1.0", "2.0"], ok
    assert censored == 2, censored

    picked = select(pool, 2, "diverse")
    keys = {r["inchikey"] for r in picked}
    assert keys == {"A", "C"}, f"diverse must span series, got {keys}"
    assert {r["inchikey"] for r in select(pool, 2, "potency")} == {"A", "B"}, \
        "potency mode must keep the naive ranking"
    assert len(select(pool, 0, "diverse")) == 3, "limit 0 means no cap"

    print("selftest OK: off-site binders excluded, assay gates the rest, "
          "diverse selection spans series")


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", default=".", type=Path)
    ap.add_argument("--cache-dir", type=Path, default=None,
                    help="cache raw ChEMBL pages here (default: "
                         "<out-dir>/.chembl_cache). ChEMBL 500s in bursts, so "
                         "a rerun with the same cache resumes instead of "
                         "refetching. Delete it to force fresh data.")
    ap.add_argument("--page", type=int, default=PAGE,
                    help=f"records per request (default {PAGE}). Larger pages "
                         "time out or 500: limit=1000 times out, 50 and 100 "
                         "both 500 even when the service is healthy.")
    ap.add_argument("--max-per-target", type=int, default=2000,
                    help="cap on potency records pulled per ChEMBL target")
    ap.add_argument("--max-compounds", type=int, default=10,
                    help="size of the reference set (default 10 = 20 dockings "
                         "= ~1 h on two receptors). 0 for no cap, which is "
                         "~450 compounds and ~50 h. Cost is linear: ~2 "
                         "dockings per compound at ~18.4 dockings/h.")
    ap.add_argument("--select", choices=("diverse", "potency"),
                    default="diverse",
                    help="how to pick them. 'diverse' round-robins across "
                         "chemotypes and skips repeated Murcko scaffolds. "
                         "'potency' takes the strongest N, which at N=10 "
                         "returns 7 compounds from 2 patent series.")
    ap.add_argument("--strict-chemotype", action="store_true",
                    help="require a match against the ATP-site scaffold "
                         "allowlist instead of accepting GyrB/ATPase assay "
                         "evidence. Much stricter, and drops most modern "
                         "GyrB series.")
    ap.add_argument("--selftest", action="store_true",
                    help="check the SMARTS filter and exit")
    ap.add_argument("--log-level", default="INFO")
    args = ap.parse_args()

    logging.basicConfig(level=args.log_level,
                        format="%(asctime)s %(levelname)s %(message)s")
    if args.selftest:
        selftest()
        return

    global _CACHE
    _CACHE = args.cache_dir or (args.out_dir / ".chembl_cache")
    _CACHE.mkdir(parents=True, exist_ok=True)
    log.info("caching ChEMBL pages in %s", _CACHE)

    build(args.out_dir, args.max_per_target, args.max_compounds,
          args.strict_chemotype, args.select, args.page)


if __name__ == "__main__":
    sys.exit(main())
