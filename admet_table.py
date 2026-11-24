#!/usr/bin/env python3
"""Collect the screen's top molecules and the reference compounds into one
table of SMILES plus physicochemical properties, for an ADMET comparison.

Why this is worth doing carefully
---------------------------------
The panel was generated under an ADMET MPO objective -- every ligand name in
exp_c5_selected_500.smi reads exp_c5_admet_mpo_guided-*. So the question is not
"are these molecules drug-like", the generator was already rewarded for that.
The question is whether the molecules that score BEST against the targets are
still the ones the MPO liked, or whether docking score and ADMET pull apart at
the top of the ranking. That is the ADMET-degradation worry about conditioning
on Vina score, asked against this project's own data instead of in general.

Which is why the panel rows here are not simply "the top 15". A molecule is
included if it is top-N on ANY of the ranking axes the report uses -- raw
bacterial score, ligand efficiency, the bacterial-minus-human margin, and each
per-receptor score -- and the selected_by column records every axis that
picked it. If the raw-score winners and the ligand-efficiency winners are
different molecules with different property profiles, that is the answer, and a
single blended ranking would have hidden it.

Sign conventions, verified against report/data/per_molecule.csv rather than
assumed, because two of them run opposite ways:

  score columns        more negative is better (Vina kcal/mol)
  le_at_best_bacterial larger is better
  margin_vs_human      = score(human) - score(best bacterial), so POSITIVE
  dhfr_selectivity     = score(human) - score(bacterial DHFR), so POSITIVE
                         means bacterial-favouring. The panel's maximum is
                         -0.065: nothing crosses zero, which is the report's
                         headline. Ranking these ascending would therefore
                         return the most HUMAN-selective molecules and quietly
                         invert the comparison.

Charge is estimated, not read off the SMILES
--------------------------------------------
RDKit's formal charge on a neutral SMILES is 0, which is useless here: the
panel is full of molecules carrying three and four carboxylic acids. At pH 7.4
those are anions, and polyanions do not cross membranes. So acids and bases are
counted by substructure and net_charge_74 is their difference. It is an
estimate -- no pKa model -- and labelled as one, but the acid count alone
separates these molecules far more sharply than logP does.

Run (rdkit optional; without it the SMILES table and the selection are still
written, so the .smi can go straight to an external predictor):

    python admet_table.py --selftest
    python admet_table.py --top 15 --out-dir report/admet

Writes admet_input.csv (everything), admet_input.smi (for upload to an
external ADMET predictor) and prints a panel-vs-knowns property comparison.
"""
from __future__ import annotations

import argparse
import csv
import json
import logging
import re
import statistics
import sys
from pathlib import Path

log = logging.getLogger("admet_table")

STATE_RE = re.compile(r"^(?P<base>.+)_(?P<state>\d+)$")

# Columns in per_molecule.csv that are NOT receptor scores. Receptor columns are
# taken from summary.json["per_receptor"] instead of by pattern-matching names,
# the same way make_report.py does it -- the panel gained a sixth receptor once
# already and anything keyed on a hardcoded list silently went stale.
DERIVED = {"ligand_base", "best_bacterial", "best_bacterial_target",
           "margin_vs_human", "dhfr_selectivity", "heavy_atoms",
           "le_at_best_bacterial", "seed"}

# Substructures for the pH 7.4 charge estimate. Deliberately short: these are
# the groups that actually occur in this panel and in the antifolates.
ACIDS = {
    "carboxylic_acid": "[CX3](=O)[OX2H1]",
    "tetrazole": "[$(c1[nH]nnn1),$(c1n[nH]nn1)]",
    "sulfonic_acid": "[SX4](=O)(=O)[OX2H1]",
    "phosphonic_acid": "[PX4](=O)([OX2H1])[OX2H1]",
    "acyl_sulfonamide": "[CX3](=O)[NX3H1][SX4](=O)=O",
}
BASES = {
    # Aliphatic amine, excluding amides, anilines, sulfonamides, nitro and
    # anything already charged or multiply bonded.
    "aliphatic_amine":
        "[NX3;H0,H1,H2;!$(N[#6]=[O,N,S]);!$(N[S,P]=O);!$(N[O]);!$(N=*);"
        "!$(N#*);!$(Nc);!$([N+])]",
    "amidine": "[NX3][CX3]=[NX2]",
    "guanidine": "[NX3][CX3](=[NX2])[NX3]",
}


def base_name(ligand: str) -> str:
    """Drop the _<state> suffix the pipeline appends to each ligand state."""
    m = STATE_RE.match(ligand)
    return m.group("base") if m else ligand


# --- inputs ----------------------------------------------------------------

def read_smi(path: Path) -> dict[str, str]:
    """{name: smiles} from a <smiles> <name> file.

    Whitespace-split, not tab-split: exp_c5_selected_500.smi is space-separated
    while the files fetch_*_actives.py writes are tab-separated, and a hardcoded
    "\\t" reads the whole panel line as one field and silently matches nothing.
    """
    out: dict[str, str] = {}
    for n, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) < 2:
            log.warning("%s:%d has no name field, skipped: %.60s", path, n, line)
            continue
        out[parts[-1]] = parts[0]
    return out


def best_scores(scores_csv: Path) -> dict[str, dict[str, float]]:
    """{compound: {receptor: best score over its states}}.

    Keyed on the base name so the several protonation/tautomer states of one
    compound collapse to the one score the report uses.
    """
    out: dict[str, dict[str, float]] = {}
    for r in csv.DictReader(scores_csv.open(encoding="utf-8")):
        if r.get("status") != "ok" or not r.get("best_score"):
            continue
        b = base_name(r["ligand"])
        s = float(r["best_score"])
        rec = out.setdefault(b, {})
        if r["receptor"] not in rec or s < rec[r["receptor"]]:
            rec[r["receptor"]] = s
    return out


# --- panel selection -------------------------------------------------------

def criteria(receptors: list[str], columns: set[str]) -> list[tuple[str, bool, str]]:
    """(column, higher_is_better, human-readable label), in reporting order.

    Only columns actually present are returned, so this works on the 5-receptor
    per_molecule.csv and on the 6-receptor one without a flag.
    """
    out: list[tuple[str, bool, str]] = []
    if "best_bacterial" in columns:
        out.append(("best_bacterial", False, "best bacterial score"))
    if "le_at_best_bacterial" in columns:
        out.append(("le_at_best_bacterial", True, "ligand efficiency"))
    if "margin_vs_human" in columns:
        out.append(("margin_vs_human", True, "margin vs human"))
    for c in sorted(columns):
        if c == "dhfr_selectivity" or c.startswith("selectivity_"):
            out.append((c, True, c.replace("_", " ")))
    for rec in receptors:
        out.append((rec, False, f"score at {rec}"))
    return out


def pick_panel(rows: list[dict], receptors: list[str], top: int
               ) -> tuple[list[dict], list[tuple[str, bool, str]]]:
    """Union of the top-`top` molecules on every ranking axis.

    Union rather than one blended score on purpose: the point of the table is
    to see whether the axes disagree about which molecules are good, which a
    blend would average away.
    """
    cols = set(rows[0].keys())
    crits = criteria(receptors, cols)
    chosen: dict[str, dict] = {}
    for col, higher, label in crits:
        try:
            ranked = sorted(rows, key=lambda r: float(r[col]), reverse=higher)
        except (ValueError, KeyError) as exc:
            log.warning("skipping axis %s: %s", col, exc)
            continue
        for r in ranked[:top]:
            key = r["ligand_base"]
            entry = chosen.setdefault(key, {"row": r, "by": []})
            entry["by"].append(label)
    out = []
    for key, e in chosen.items():
        d = dict(e["row"])
        d["selected_by"] = "; ".join(e["by"])
        d["n_criteria"] = len(e["by"])
        out.append(d)
    out.sort(key=lambda d: (-d["n_criteria"], float(d.get("best_bacterial", 0))))
    return out, crits


# --- descriptors -----------------------------------------------------------

DESC_COLS = ["mw", "clogp", "tpsa", "hbd", "hba", "rotb", "aromatic_rings",
             "rings", "fraction_csp3", "n_acid", "n_base", "net_charge_74",
             "qed", "lipinski_violations", "veber_pass"]


def descriptors(rows: list[dict], smiles_key: str = "smiles") -> bool:
    """Fill DESC_COLS in place. Returns True if rdkit was available.

    Optional by design: the .smi for an external predictor and the whole
    selection do not need rdkit, so this degrades to blank columns rather than
    refusing to produce the table on a machine without it.
    """
    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem import Crippen, Descriptors, QED, rdMolDescriptors
    except ImportError:
        log.warning("rdkit absent -- property columns left blank. The SMILES "
                    "table and admet_input.smi are complete; run this on the "
                    "cluster env (or pip install rdkit) for the properties.")
        for r in rows:
            for c in DESC_COLS:
                r.setdefault(c, "")
        return False

    RDLogger.DisableLog("rdApp.*")
    acids = {k: Chem.MolFromSmarts(v) for k, v in ACIDS.items()}
    bases = {k: Chem.MolFromSmarts(v) for k, v in BASES.items()}
    for name, patt in list(acids.items()) + list(bases.items()):
        assert patt is not None, f"bad SMARTS for {name}"

    for r in rows:
        mol = Chem.MolFromSmiles(r[smiles_key] or "")
        if mol is None:
            log.error("rdkit cannot parse %s: %.70s", r.get("id"), r[smiles_key])
            for c in DESC_COLS:
                r[c] = ""
            continue
        # Count matches, not matching groups: a molecule with four carboxylic
        # acids must read as four, which is the whole point here.
        n_acid = sum(len(mol.GetSubstructMatches(p)) for p in acids.values())
        n_base = sum(len(mol.GetSubstructMatches(p)) for p in bases.values())
        mw = Descriptors.MolWt(mol)
        clogp = Crippen.MolLogP(mol)
        tpsa = rdMolDescriptors.CalcTPSA(mol)
        hbd = rdMolDescriptors.CalcNumHBD(mol)
        hba = rdMolDescriptors.CalcNumHBA(mol)
        rotb = rdMolDescriptors.CalcNumRotatableBonds(mol)
        r.update({
            "mw": f"{mw:.1f}", "clogp": f"{clogp:.2f}", "tpsa": f"{tpsa:.1f}",
            "hbd": hbd, "hba": hba, "rotb": rotb,
            "aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(mol),
            "rings": rdMolDescriptors.CalcNumRings(mol),
            "fraction_csp3": f"{rdMolDescriptors.CalcFractionCSP3(mol):.2f}",
            "n_acid": n_acid, "n_base": n_base,
            "net_charge_74": n_base - n_acid,
            "qed": f"{QED.qed(mol):.3f}",
            "lipinski_violations": sum((mw > 500, clogp > 5, hbd > 5, hba > 10)),
            "veber_pass": int(rotb <= 10 and tpsa <= 140),
        })
    return True


# --- report ----------------------------------------------------------------

NUMERIC_SUMMARY = ["mw", "clogp", "tpsa", "hbd", "hba", "rotb",
                   "aromatic_rings", "fraction_csp3", "n_acid", "n_base",
                   "net_charge_74", "qed", "heavy_atoms"]


def compare(rows: list[dict]) -> None:
    """Median of each property per set, printed side by side."""
    sets: dict[str, list[dict]] = {}
    for r in rows:
        sets.setdefault(r["set"], []).append(r)
    order = sorted(sets, key=lambda s: (s != "panel_top", s))
    width = max(len(s) for s in order) + 2

    print("\nMedian properties by set (n in brackets)")
    header = "  " + "property".ljust(20) + "".join(
        f"{s} [{len(sets[s])}]".rjust(width + 6) for s in order)
    print(header)
    print("  " + "-" * (len(header) - 2))
    for prop in NUMERIC_SUMMARY:
        cells = []
        any_val = False
        for s in order:
            vals = []
            for r in sets[s]:
                try:
                    vals.append(float(r[prop]))
                except (KeyError, TypeError, ValueError):
                    pass
            if vals:
                any_val = True
                cells.append(f"{statistics.median(vals):.2f}".rjust(width + 6))
            else:
                cells.append("--".rjust(width + 6))
        if any_val:
            print("  " + prop.ljust(20) + "".join(cells))


# --- build -----------------------------------------------------------------

def build(args) -> None:
    panel_rows = list(csv.DictReader(args.panel_csv.open(encoding="utf-8")))
    if not panel_rows:
        raise SystemExit(f"{args.panel_csv} has no rows")
    summary = json.loads(args.panel_summary.read_text(encoding="utf-8")) \
        if args.panel_summary.exists() else {}
    per_rec = summary.get("per_receptor", {})
    receptors = [c for c in panel_rows[0] if c in per_rec] or \
                [c for c in panel_rows[0] if c not in DERIVED
                 and not c.startswith("selectivity_")]
    if not per_rec:
        log.warning("%s absent or has no per_receptor block -- receptor columns "
                    "inferred as %s", args.panel_summary, receptors)

    panel_smi = read_smi(args.panel_smi)
    picked, crits = pick_panel(panel_rows, receptors, args.top)
    log.info("panel: %d molecules selected on %d axes (top %d each)",
             len(picked), len(crits), args.top)

    out_rows: list[dict] = []
    missing = 0
    for d in picked:
        smi = panel_smi.get(d["ligand_base"])
        if not smi:
            missing += 1
            continue
        out_rows.append({"id": d["ligand_base"], "set": "panel_top",
                         "smiles": smi, "label": "", "potency_nM": "",
                         "selected_by": d["selected_by"],
                         "n_criteria": d["n_criteria"],
                         **{k: v for k, v in d.items()
                            if k not in ("ligand_base", "selected_by",
                                         "n_criteria")}})
    if missing:
        # Loud, because a silent drop here would shrink the comparison set
        # without changing any number that looks wrong.
        log.error("%d selected panel molecule(s) had no SMILES in %s",
                  missing, args.panel_smi)

    # The rest of the panel, as the control that makes the comparison mean
    # anything. panel_top vs knowns alone cannot separate "docking score
    # selects against ADMET" from "the generator makes molecules that look like
    # this and the top-scoring ones are no different". Only panel_top vs
    # panel_rest isolates the score, because both come from the same generator
    # under the same ADMET MPO.
    if args.all_panel:
        chosen_ids = {d["ligand_base"] for d in picked}
        rest = 0
        for r in panel_rows:
            if r["ligand_base"] in chosen_ids:
                continue
            smi = panel_smi.get(r["ligand_base"])
            if not smi:
                continue
            out_rows.append({"id": r["ligand_base"], "set": "panel_rest",
                             "smiles": smi, "label": "", "potency_nM": "",
                             "selected_by": "", "n_criteria": "",
                             **{k: v for k, v in r.items()
                                if k != "ligand_base"}})
            rest += 1
        log.info("panel_rest: %d molecules included as the within-generator "
                 "control", rest)

    for d in args.knowns_dir:
        if not d.is_dir():
            log.warning("knowns dir %s absent, skipped", d)
            continue
        act = sorted(d.glob("*_actives.csv"))
        sc = d / "scores.csv"
        if not act or not sc.exists():
            log.warning("%s has no *_actives.csv and scores.csv pair, skipped", d)
            continue
        scores = best_scores(sc)
        for a in act:
            setname = f"knowns_{a.name.replace('_actives.csv', '')}"
            n = 0
            for r in csv.DictReader(a.open(encoding="utf-8")):
                row = {"id": r["name"], "set": setname,
                       "smiles": r.get("smiles", ""),
                       "label": r.get("pref_name") or r.get("chemotype") or
                                r.get("chembl_id", ""),
                       "potency_nM": r.get("standard_value_nM")
                                     or r.get("ecoli_nM", ""),
                       "selected_by": "reference compound", "n_criteria": "",
                       "heavy_atoms": r.get("heavy_atoms", "")}
                row.update(scores.get(r["name"], {}))
                # Censored potencies are bounds, not measurements -- the same
                # trap that once inverted the GyrB potency correlation.
                if (r.get("standard_relation") or "=").strip() != "=":
                    row["label"] = f"{row['label']} (potency censored)"
                out_rows.append(row)
                n += 1
            log.info("knowns: %d compounds from %s", n, a.name)

    had_rdkit = descriptors(out_rows)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    cols = ["id", "set", "label", "smiles", "selected_by", "n_criteria",
            "potency_nM", *DESC_COLS]
    cols += [c for c in dict.fromkeys(
        k for r in out_rows for k in r) if c not in cols]
    with (args.out_dir / "admet_input.csv").open("w", newline="",
                                                 encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(out_rows)
    (args.out_dir / "admet_input.smi").write_text(
        "".join(f"{r['smiles']}\t{r['id']}\n" for r in out_rows if r["smiles"]),
        encoding="utf-8")

    print(f"\n{len(out_rows)} compounds -> {args.out_dir}/admet_input.csv")
    print(f"{'':25}{args.out_dir}/admet_input.smi  (upload this)")
    counts: dict[str, int] = {}
    for r in out_rows:
        counts[r["set"]] = counts.get(r["set"], 0) + 1
    for s, n in sorted(counts.items()):
        print(f"  {s:<22} {n}")
    if had_rdkit:
        compare(out_rows)
    print("\nRanking axes used for the panel (top %d each):" % args.top)
    for col, higher, label in crits:
        print(f"  {label:<34} {'max' if higher else 'min'} of {col}")


# --- selftest --------------------------------------------------------------

def selftest() -> None:
    assert base_name("gyrb_known_001_0") == "gyrb_known_001"
    assert base_name("dhfr_known_018_12") == "dhfr_known_018"
    assert base_name("no_state_here") == "no_state_here"

    import tempfile
    with tempfile.TemporaryDirectory() as td:
        t = Path(td)
        # Space-separated (panel) and tab-separated (knowns) must both read.
        (t / "a.smi").write_text("CCO name-one\nCCN\tname-two\n\n", encoding="utf-8")
        got = read_smi(t / "a.smi")
        assert got == {"name-one": "CCO", "name-two": "CCN"}, got

        (t / "scores.csv").write_text(
            "ligand,receptor,best_score,status\n"
            "k_001_0,recA,-8.0,ok\n"
            "k_001_1,recA,-9.5,ok\n"       # better state wins
            "k_001_0,recB,-7.0,ok\n"
            "k_002_0,recA,-5.0,failed\n"   # non-ok ignored
            "k_002_1,recA,,ok\n",          # blank score ignored
            encoding="utf-8")
        bs = best_scores(t / "scores.csv")
        assert bs == {"k_001": {"recA": -9.5, "recB": -7.0}}, bs

    # Sign conventions: the two selectivity-like axes rank by MAXIMUM, the
    # score axes by MINIMUM. Getting this backwards returns the worst molecules.
    crits = criteria(["recA", "recB"],
                     {"best_bacterial", "le_at_best_bacterial",
                      "margin_vs_human", "dhfr_selectivity",
                      "selectivity_recA", "recA", "recB", "ligand_base"})
    d = {c: higher for c, higher, _ in crits}
    assert d["best_bacterial"] is False and d["recA"] is False
    assert d["le_at_best_bacterial"] is True
    assert d["margin_vs_human"] is True
    assert d["dhfr_selectivity"] is True and d["selectivity_recA"] is True

    # Absent columns must drop out rather than raise: the 5-receptor
    # per_molecule.csv has no selectivity_* columns at all.
    few = criteria(["recA"], {"recA", "ligand_base"})
    assert [c for c, _, _ in few] == ["recA"], few

    # The union keeps a molecule that wins on one axis only, and records every
    # axis that chose it.
    rows = [{"ligand_base": "m1", "best_bacterial": "-9.0",
             "le_at_best_bacterial": "0.10", "recA": "-9.0"},
            {"ligand_base": "m2", "best_bacterial": "-6.0",
             "le_at_best_bacterial": "0.40", "recA": "-6.0"},
            {"ligand_base": "m3", "best_bacterial": "-7.0",
             "le_at_best_bacterial": "0.20", "recA": "-7.0"}]
    picked, _ = pick_panel(rows, ["recA"], top=1)
    ids = {p["ligand_base"] for p in picked}
    assert ids == {"m1", "m2"}, ids
    m1 = next(p for p in picked if p["ligand_base"] == "m1")
    assert m1["n_criteria"] == 2, m1            # raw score and score at recA
    m2 = next(p for p in picked if p["ligand_base"] == "m2")
    assert m2["selected_by"] == "ligand efficiency", m2["selected_by"]

    # top >= len(rows) must return everything, not raise on the slice.
    allp, _ = pick_panel(rows, ["recA"], top=99)
    assert len(allp) == 3

    # rdkit optional: blanks, not an exception.
    rs = [{"id": "x", "smiles": "OC(=O)CN"}]
    had = descriptors(rs)
    assert all(c in rs[0] for c in DESC_COLS)
    if had:
        # Glycine: one acid, one aliphatic amine, net zero.
        assert rs[0]["n_acid"] == 1, rs[0]
        assert rs[0]["n_base"] == 1, rs[0]
        assert rs[0]["net_charge_74"] == 0, rs[0]
        poly = [{"id": "p", "smiles": "OC(=O)CC(C(=O)O)C(=O)O"}]
        descriptors(poly)
        assert poly[0]["n_acid"] == 3, poly[0]
        assert poly[0]["net_charge_74"] == -3, poly[0]
        # An amide nitrogen is not a base.
        am = [{"id": "a", "smiles": "CC(=O)NC"}]
        descriptors(am)
        assert am[0]["n_base"] == 0, am[0]
        bad = [{"id": "b", "smiles": "not a molecule"}]
        descriptors(bad)
        assert bad[0]["mw"] == "", bad[0]
    print("selftest OK" + ("" if had else " (rdkit absent: property checks skipped)"))


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--panel-csv", type=Path,
                   default=Path("report/data/per_molecule.csv"))
    p.add_argument("--panel-summary", type=Path,
                   default=Path("report/data/summary.json"),
                   help="for the authoritative receptor-column list")
    p.add_argument("--panel-smi", type=Path,
                   default=Path("exp_c5_selected_500.smi"))
    p.add_argument("--knowns-dir", type=Path, nargs="*",
                   default=[Path("report/data_knowns"), Path("report/data_dhfr")],
                   help="dirs holding <prefix>_actives.csv + scores.csv; "
                        "absent ones are skipped")
    p.add_argument("--top", type=int, default=15,
                   help="top N per ranking axis (default 15)")
    p.add_argument("--all-panel", action="store_true",
                   help="also include every non-selected panel molecule as the "
                        "'panel_rest' control. This is what separates a docking-"
                        "score effect from the generator's own chemistry, so use "
                        "it for the comparison; leave it off when the .smi is "
                        "going to an external predictor with an upload limit")
    p.add_argument("--out-dir", type=Path, default=Path("report/admet"))
    p.add_argument("--selftest", action="store_true")
    p.add_argument("-v", "--verbose", action="store_true")
    a = p.parse_args()
    logging.basicConfig(level=logging.DEBUG if a.verbose else logging.INFO,
                        format="%(levelname)s %(message)s")
    if a.selftest:
        selftest()
        return
    build(a)


if __name__ == "__main__":
    main()
