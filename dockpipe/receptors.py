from __future__ import annotations

import json
import logging
import shutil
import subprocess
import urllib.request
from pathlib import Path

import numpy as np

from .common import BoxSpec, ReceptorSpec

log = logging.getLogger(__name__)

RCSB_PDB = "https://files.rcsb.org/download/{pdb_id}.pdb"
PDB_REDO = "https://pdb-redo.eu/db/{pdb_id}/{pdb_id}_final.pdb"



def fetch_pdb(pdb_id: str, dest_dir: Path, use_pdb_redo: bool = False) -> Path:
    dest_dir.mkdir(parents=True, exist_ok=True)
    pdb_id = pdb_id.lower()
    suffix = "_redo" if use_pdb_redo else ""
    dest = dest_dir / f"{pdb_id}{suffix}.pdb"
    if dest.exists() and dest.stat().st_size > 0:
        log.info("Using cached %s", dest)
        return dest

    url = (PDB_REDO if use_pdb_redo else RCSB_PDB).format(pdb_id=pdb_id)
    log.info("Fetching %s", url)
    with urllib.request.urlopen(url, timeout=60) as resp:
        data = resp.read()
    dest.write_bytes(data)
    return dest



def _is_water(resname: str) -> bool:
    return resname.strip() in {"HOH", "WAT", "DOD"}


def split_pdb(
    pdb_path: Path,
    keep_chains: list[str] | None,
    keep_hetatms: list[str],
    ligand_resname: str | None,
) -> tuple[list[str], list[str]]:

    receptor: list[str] = []
    ligand: list[str] = []
    keep_het = {h.strip().upper() for h in keep_hetatms}
    lig_res = ligand_resname.strip().upper() if ligand_resname else None

    with open(pdb_path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            rec = line[:6]
            if rec not in ("ATOM  ", "HETATM"):
                continue
            altloc = line[16]
            if altloc not in (" ", "A"):
                continue
            chain = line[21]
            resname = line[17:20].strip().upper()

            if keep_chains and chain not in keep_chains:
                continue

            if rec == "ATOM  ":
                receptor.append(line)
            else:
                if lig_res and resname == lig_res:
                    ligand.append(line)
                elif resname in keep_het:
                    receptor.append(line)
    return receptor, ligand


def coords_from_lines(lines: list[str]) -> np.ndarray:
    return np.array(
        [[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in lines],
        dtype=float,
    )


def derive_box(ligand_lines: list[str], spec: BoxSpec) -> tuple[tuple, tuple]:
    if not ligand_lines:
        raise ValueError("No reference-ligand atoms found; cannot derive box.")
    xyz = coords_from_lines(ligand_lines)
    center = tuple(float(v) for v in xyz.mean(axis=0))
    if spec.size is not None:
        size = tuple(float(v) for v in spec.size)
    else:
        extent = xyz.max(axis=0) - xyz.min(axis=0)
        size = tuple(float(v) + 2 * spec.padding + 8.0 for v in extent)
    vol = size[0] * size[1] * size[2]
    if vol > 27000:
        log.warning(
            "Box volume %.0f A^3 exceeds ~27000; raise exhaustiveness or shrink.",
            vol,
        )
    return center, size



def protonate_with_pdbfixer(in_pdb: Path, out_pdb: Path, ph: float = 7.4) -> Path:
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile

    fixer = PDBFixer(filename=str(in_pdb))
    fixer.findMissingResidues()
    fixer.missingResidues = {}
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(ph)
    with open(out_pdb, "w", encoding="utf-8") as fh:
        PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)
    return out_pdb


def receptor_to_pdbqt(in_pdb: Path, out_pdbqt: Path) -> Path:
    out_pdbqt.parent.mkdir(parents=True, exist_ok=True)

    for cmd in (
        ["mk_prepare_receptor.py", "--read_pdb", str(in_pdb),
         "-o", str(out_pdbqt.with_suffix("")), "--write_pdbqt"],
        ["prepare_receptor", "-r", str(in_pdb), "-o", str(out_pdbqt),
         "-A", "checkhydrogens"],
    ):
        exe = shutil.which(cmd[0])
        if exe is None:
            continue
        log.info("Receptor PDBQT via %s", cmd[0])
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode == 0 and out_pdbqt.exists():
            return out_pdbqt
        log.warning("%s failed (rc=%d): %s", cmd[0], proc.returncode,
                    proc.stderr.strip()[:500])

    raise RuntimeError(
        "No receptor preparation tool succeeded. Install Meeko's "
        "mk_prepare_receptor.py or the ADFR suite's prepare_receptor."
    )


def ligand_lines_to_pdbqt(ligand_lines: list[str], out_pdbqt: Path,
                          template_smiles: str | None = None) -> Path | None:
    """Convert extracted co-crystal ligand HETATM lines to PDBQT for redocking."""
    from rdkit import Chem
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    block = "".join(ligand_lines) + "END\n"
    mol = Chem.MolFromPDBBlock(block, removeHs=False, sanitize=False)
    if mol is None:
        log.warning("Could not parse extracted ligand block")
        return None

    if template_smiles:
        try:
            ref = Chem.MolFromSmiles(template_smiles)
            mol = Chem.AllChem.AssignBondOrdersFromTemplate(ref, mol)
        except Exception as exc:  # noqa: BLE001
            log.warning("Bond-order assignment from template failed: %s", exc)

    try:
        Chem.SanitizeMol(mol)
        mol = Chem.AddHs(mol, addCoords=True)
    except Exception as exc:  # noqa: BLE001
        log.warning("Ligand sanitization failed: %s", exc)
        return None

    molprep = MoleculePreparation()
    for molsetup in molprep.prepare(mol):
        pdbqt, ok, err = PDBQTWriterLegacy.write_string(molsetup)
        if ok:
            out_pdbqt.parent.mkdir(parents=True, exist_ok=True)
            out_pdbqt.write_text(pdbqt, encoding="utf-8")
            return out_pdbqt
        log.warning("Reference ligand PDBQT write failed: %s", err)
    return None


# --------------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------------

def prepare_receptor(
    spec: ReceptorSpec,
    workdir: Path,
    use_pdb_redo: bool = False,
) -> dict:
    """Full receptor prep. Returns a record with paths and the docking box."""
    workdir = Path(workdir) / spec.name
    workdir.mkdir(parents=True, exist_ok=True)

    if spec.local_pdb:
        src = Path(spec.local_pdb)
        if not src.exists():
            raise FileNotFoundError(f"local_pdb not found: {src}")
    elif spec.pdb_id:
        src = fetch_pdb(spec.pdb_id, workdir / "raw", use_pdb_redo)
    else:
        raise ValueError(f"Receptor {spec.name}: need pdb_id or local_pdb")

    lig_res = spec.redock_ligand_resname or spec.box.reference_ligand_resname
    rec_lines, lig_lines = split_pdb(
        src, spec.keep_chains, spec.keep_hetatms, lig_res
    )
    if not rec_lines:
        raise ValueError(f"Receptor {spec.name}: no atoms retained after filtering")

    log.info("%s: %d receptor atoms, %d reference-ligand atoms",
             spec.name, len(rec_lines), len(lig_lines))

    clean_pdb = workdir / f"{spec.name}_clean.pdb"
    clean_pdb.write_text("".join(rec_lines) + "END\n", encoding="utf-8")

    prot_pdb = workdir / f"{spec.name}_H.pdb"
    try:
        protonate_with_pdbfixer(clean_pdb, prot_pdb, spec.ph)
    except Exception as exc:  # noqa: BLE001
        log.warning("%s: PDBFixer failed (%s); using unprotonated model", spec.name, exc)
        shutil.copy(clean_pdb, prot_pdb)

    rec_pdbqt = workdir / f"{spec.name}.pdbqt"
    receptor_to_pdbqt(prot_pdb, rec_pdbqt)

    # Box
    if spec.box.center is not None:
        center = tuple(float(v) for v in spec.box.center)
        size = tuple(float(v) for v in spec.box.size)
    else:
        center, size = derive_box(lig_lines, spec.box)
    log.info("%s: box center=%s size=%s", spec.name, center, size)

    # Reference ligand for the redocking control
    ref_pdbqt = None
    ref_sdf = None
    if lig_lines:
        ref_pdbqt = workdir / f"{spec.name}_ref.pdbqt"
        if ligand_lines_to_pdbqt(lig_lines, ref_pdbqt) is None:
            ref_pdbqt = None
        ref_sdf = workdir / f"{spec.name}_ref.pdb"
        ref_sdf.write_text("".join(lig_lines) + "END\n", encoding="utf-8")

    record = {
        "name": spec.name,
        "source_pdb": str(src),
        "pdb_id": spec.pdb_id,
        "pdb_redo": use_pdb_redo,
        "receptor_pdbqt": str(rec_pdbqt),
        "receptor_pdb": str(prot_pdb),
        "reference_ligand_pdbqt": str(ref_pdbqt) if ref_pdbqt else None,
        "reference_ligand_pdb": str(ref_sdf) if ref_sdf else None,
        "reference_ligand_resname": lig_res,
        "kept_hetatms": spec.keep_hetatms,
        "kept_chains": spec.keep_chains,
        "box_center": list(center),
        "box_size": list(size),
        "ph": spec.ph,
        "notes": spec.notes,
    }
    (workdir / "receptor_record.json").write_text(
        json.dumps(record, indent=2), encoding="utf-8"
    )
    return record


def prepare_receptors(
    specs: list[ReceptorSpec],
    workdir: str | Path,
    use_pdb_redo: bool = False,
) -> list[dict]:
    workdir = Path(workdir)
    out = []
    for spec in specs:
        try:
            out.append(prepare_receptor(spec, workdir, use_pdb_redo))
        except Exception as exc:  # noqa: BLE001
            log.error("Receptor %s failed: %s", spec.name, exc)
            out.append({"name": spec.name, "status": "failed", "error": str(exc)})
    (workdir / "receptor_index.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8"
    )
    return out
