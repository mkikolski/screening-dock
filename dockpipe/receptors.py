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
# The chemical component's idealised structure. Coordinates are useless here,
# but the bond orders are authoritative and cannot be recovered from the PDB
# entry, which stores no connectivity for ligands.
RCSB_LIGAND = "https://files.rcsb.org/ligands/download/{res}_ideal.sdf"



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



def fetch_ligand_template(resname: str, dest_dir: Path) -> Path | None:
    """Idealised SDF for a PDB chemical component, cached. None if unavailable."""
    dest_dir.mkdir(parents=True, exist_ok=True)
    res = resname.strip().upper()
    dest = dest_dir / f"{res}_ideal.sdf"
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    url = RCSB_LIGAND.format(res=res)
    try:
        with urllib.request.urlopen(url, timeout=60) as resp:
            dest.write_bytes(resp.read())
    except Exception as exc:  # noqa: BLE001
        log.warning("No ideal SDF for %s (%s): %s", res, url, exc)
        return None
    return dest


def _is_water(resname: str) -> bool:
    return resname.strip() in {"HOH", "WAT", "DOD"}


def split_pdb(
    pdb_path: Path,
    keep_chains: list[str] | None,
    keep_hetatms: list[str],
    ligand_resname: str | None,
) -> tuple[list[str], list[str], list[str]]:
    """Split into (protein, reference ligand, cofactors).

    Cofactors come out separately rather than merged into the protein: Meeko
    >=0.6 treats the receptor as a polymer and needs a residue template per
    residue, and it ships none for NAP/NDP, so a cofactor left in the protein
    file aborts the whole preparation.
    """
    receptor: list[str] = []
    ligand: list[str] = []
    cofactor: list[str] = []
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
            # Keeping altloc 'A' in column 17 makes Meeko stop and demand
            # --wanted_altloc for every such residue, even though we have
            # already picked one here. Blank it so the choice is final.
            line = line[:16] + " " + line[17:]
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
                    cofactor.append(line)
    return receptor, ligand, cofactor


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



def protonate_with_pdbfixer(in_pdb: Path, out_pdb: Path, ph: float = 7.4,
                            add_heavy: bool = True) -> Path:
    """Add hydrogens at `ph`, and optionally rebuild missing heavy atoms.

    add_heavy=False exists because addMissingAtoms() caps the C-terminus with
    OXT, and Meeko then cannot build an RDKit mol for that residue at all
    (3U2D GLU A:230, 4ZVI THR A:392 both die with a 5-valent backbone carbon).
    Without it the unresolved side chains stay unresolved and Meeko's
    --allow_bad_res drops them, which beats docking against invented atoms.
    """
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile

    fixer = PDBFixer(filename=str(in_pdb))
    fixer.findMissingResidues()
    fixer.missingResidues = {}
    if add_heavy:
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
    fixer.addMissingHydrogens(ph)
    with open(out_pdb, "w", encoding="utf-8") as fh:
        PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)
    return out_pdb


def receptor_to_pdbqt(in_pdb: Path, out_pdbqt: Path) -> Path:
    out_pdbqt.parent.mkdir(parents=True, exist_ok=True)

    for cmd in (
        # --allow_bad_res drops residues Meeko cannot template. In practice
        # those are chain termini and side chains the crystal never resolved
        # (3U2D: TYR A:14 missing 7 heavy atoms, GLU A:230 the C-terminus),
        # which are not part of any binding site. Cofactors are not in this
        # file at all -- they are prepared separately -- so nothing that
        # matters for docking can be silently deleted here.
        ["mk_prepare_receptor.py", "--read_pdb", str(in_pdb),
         "-o", str(out_pdbqt.with_suffix("")), "--write_pdbqt",
         "--allow_bad_res", "--default_altloc", "A"],
        ["prepare_receptor", "-r", str(in_pdb), "-o", str(out_pdbqt),
         "-A", "checkhydrogens"],
    ):
        exe = shutil.which(cmd[0])
        if exe is None:
            continue
        log.info("Receptor PDBQT via %s", cmd[0])
        proc = subprocess.run(cmd, capture_output=True, text=True)
        # Kept on success too: this is where Meeko names the residues
        # --allow_bad_res threw away, and that list wants reviewing.
        errlog = out_pdbqt.with_suffix(".prep.log")
        errlog.write_text(
            f"$ {' '.join(cmd)}\nrc={proc.returncode}\n"
            f"--- stdout ---\n{proc.stdout}\n--- stderr ---\n{proc.stderr}",
            encoding="utf-8",
        )
        if proc.returncode == 0 and out_pdbqt.exists():
            log.info("%s: deleted-residue report in %s", cmd[0], errlog)
            return out_pdbqt
        log.warning("%s failed (rc=%d); full output in %s\n%s",
                    cmd[0], proc.returncode, errlog,
                    proc.stderr.strip()[-800:])

    raise RuntimeError(
        "No receptor preparation tool succeeded. Install Meeko's "
        "mk_prepare_receptor.py or the ADFR suite's prepare_receptor."
    )


def add_cofactor_to_receptor(
    cofactor_lines: list[str],
    rec_pdbqt: Path,
    out_cofactor_pdbqt: Path,
    ph: float = 7.4,
) -> Path:
    """Type the cofactor with Open Babel and append it to the receptor PDBQT.

    Open Babel rather than RDKit+Meeko because bond orders cannot be recovered
    from PDB coordinates without a SMILES template, and Vina does not need
    them: its scoring reads the AutoDock atom type column only and ignores
    partial charges entirely. A rigid cofactor therefore only has to be typed,
    not charged correctly.
    """
    if shutil.which("obabel") is None:
        raise RuntimeError(
            "Cofactors declared in keep_hetatms but obabel is not installed; "
            "run: conda install -c conda-forge openbabel"
        )

    src_pdb = out_cofactor_pdbqt.with_suffix(".pdb")
    src_pdb.write_text("".join(cofactor_lines) + "END\n", encoding="utf-8")

    # -xr: no torsion tree, the cofactor stays rigid. -p: add H at this pH.
    cmd = ["obabel", str(src_pdb), "-O", str(out_cofactor_pdbqt),
           "-xr", "-p", str(ph)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if not out_cofactor_pdbqt.exists() or not out_cofactor_pdbqt.stat().st_size:
        raise RuntimeError(
            f"obabel failed to type the cofactor (rc={proc.returncode}): "
            f"{proc.stderr.strip()[-500:]}"
        )

    # Only the coordinate records: obabel may emit ROOT/ENDROOT/TORSDOF, which
    # a rigid Vina receptor must not contain. Serials are left clashing with
    # the protein's on purpose -- Vina ignores them for a rigid receptor.
    atoms = [l for l in out_cofactor_pdbqt.read_text(encoding="utf-8").splitlines(True)
             if l[:6] in ("ATOM  ", "HETATM")]
    if not atoms:
        raise RuntimeError(f"No atoms in {out_cofactor_pdbqt} after obabel")

    with open(rec_pdbqt, "a", encoding="utf-8") as fh:
        fh.writelines(atoms)
    log.info("Appended %d cofactor atoms to %s", len(atoms), rec_pdbqt.name)
    return out_cofactor_pdbqt


def ligand_lines_to_pdbqt(ligand_lines: list[str], out_pdbqt: Path,
                          template_sdf: Path | None = None
                          ) -> tuple[Path | None, str]:
    """Convert extracted co-crystal ligand HETATM lines to PDBQT for redocking.

    Returns (path, bond_order_source). A PDB entry stores no connectivity for
    ligands, so without `template_sdf` RDKit perceives bonds by proximity and
    makes every one of them single; sanitization then saturates the aromatic
    rings. The result is a hydrogenated analogue of the crystal ligand, and
    redocking it measures nothing useful -- hence the second return value.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    block = "".join(ligand_lines) + "END\n"
    mol = Chem.MolFromPDBBlock(block, removeHs=False, sanitize=False)
    if mol is None:
        log.warning("Could not parse extracted ligand block")
        return None, "unparsed"

    source = "proximity_guess"
    if template_sdf is not None:
        try:
            ref = Chem.MolFromMolFile(str(template_sdf), sanitize=True)
            if ref is None:
                raise ValueError(f"unreadable template {template_sdf}")
            mol = AllChem.AssignBondOrdersFromTemplate(Chem.RemoveHs(ref), mol)
            source = "ideal_sdf"
        except Exception as exc:  # noqa: BLE001
            # Loud: the redocking verdict is not interpretable without this.
            log.warning(
                "Bond-order assignment from %s failed (%s). Reference ligand "
                "will have GUESSED bond orders and its redocking RMSD is not "
                "a valid control.", template_sdf, exc)

    try:
        Chem.SanitizeMol(mol)
        mol = Chem.AddHs(mol, addCoords=True)
    except Exception as exc:  # noqa: BLE001
        log.warning("Ligand sanitization failed: %s", exc)
        return None, source

    molprep = MoleculePreparation()
    for molsetup in molprep.prepare(mol):
        pdbqt, ok, err = PDBQTWriterLegacy.write_string(molsetup)
        if ok:
            out_pdbqt.parent.mkdir(parents=True, exist_ok=True)
            out_pdbqt.write_text(pdbqt, encoding="utf-8")
            return out_pdbqt, source
        log.warning("Reference ligand PDBQT write failed: %s", err)
    return None, source


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
    rec_lines, lig_lines, cof_lines = split_pdb(
        src, spec.keep_chains, spec.keep_hetatms, lig_res
    )
    if not rec_lines:
        raise ValueError(f"Receptor {spec.name}: no atoms retained after filtering")
    if spec.keep_hetatms and not cof_lines:
        raise ValueError(
            f"Receptor {spec.name}: keep_hetatms={spec.keep_hetatms} matched "
            "nothing in the entry. Docking would run into a cofactor-less "
            "pocket -- check the residue names against the PDB entry."
        )

    log.info("%s: %d receptor atoms, %d reference-ligand atoms, %d cofactor atoms",
             spec.name, len(rec_lines), len(lig_lines), len(cof_lines))

    clean_pdb = workdir / f"{spec.name}_clean.pdb"
    clean_pdb.write_text("".join(rec_lines) + "END\n", encoding="utf-8")

    prot_pdb = workdir / f"{spec.name}_H.pdb"
    rec_pdbqt = workdir / f"{spec.name}.pdbqt"

    # Rebuilding missing heavy atoms helps some entries and makes Meeko reject
    # others outright, and which is which is a property of the crystal, not of
    # anything we can check up front. So try it, then retry without.
    for add_heavy in (True, False):
        try:
            protonate_with_pdbfixer(clean_pdb, prot_pdb, spec.ph, add_heavy)
        except Exception as exc:  # noqa: BLE001
            log.warning("%s: PDBFixer failed (%s); using unprotonated model",
                        spec.name, exc)
            shutil.copy(clean_pdb, prot_pdb)
        try:
            receptor_to_pdbqt(prot_pdb, rec_pdbqt)
            break
        except RuntimeError:
            if not add_heavy:
                raise
            log.warning("%s: receptor prep failed with rebuilt heavy atoms; "
                        "retrying without them", spec.name)

    cof_pdbqt = None
    if cof_lines:
        cof_pdbqt = add_cofactor_to_receptor(
            cof_lines, rec_pdbqt, workdir / f"{spec.name}_cofactor.pdbqt", spec.ph
        )

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
    bond_orders = None
    if lig_lines:
        template = fetch_ligand_template(lig_res, workdir / "raw") if lig_res else None
        ref_pdbqt = workdir / f"{spec.name}_ref.pdbqt"
        ref_pdbqt, bond_orders = ligand_lines_to_pdbqt(
            lig_lines, ref_pdbqt, template
        )
        ref_sdf = workdir / f"{spec.name}_ref.pdb"
        ref_sdf.write_text("".join(lig_lines) + "END\n", encoding="utf-8")
        log.info("%s: reference ligand %s, bond orders from %s",
                 spec.name, lig_res, bond_orders)

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
        "reference_bond_orders": bond_orders,
        "cofactor_pdbqt": str(cof_pdbqt) if cof_pdbqt else None,
        "cofactor_atoms": len(cof_lines),
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
