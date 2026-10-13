from __future__ import annotations

import csv
import json
import logging
import os
import time
from dataclasses import dataclass
from multiprocessing import Pool
from pathlib import Path

import numpy as np

from .common import VinaSpec

log = logging.getLogger(__name__)

_VINA = None
_VINA_KEY = None


@dataclass
class DockTask:
    ligand_pdbqt: str
    receptor_name: str
    receptor_pdbqt: str
    center: tuple[float, float, float]
    size: tuple[float, float, float]
    out_dir: str
    vina: VinaSpec


def _get_vina(task: DockTask):
    global _VINA, _VINA_KEY
    from vina import Vina

    key = (task.receptor_pdbqt, task.vina.scoring)
    if _VINA_KEY == key and _VINA is not None:
        return _VINA

    v = Vina(
        sf_name=task.vina.scoring, cpu=task.vina.cpu, seed=task.vina.seed, verbosity=0
    )
    v.set_receptor(task.receptor_pdbqt)
    v.compute_vina_maps(center=list(task.center), box_size=list(task.size))
    _VINA, _VINA_KEY = v, key
    return v


def dock_one(task: DockTask) -> dict:
    lig = Path(task.ligand_pdbqt)
    rec: dict = {
        "ligand": lig.stem,
        "ligand_file": str(lig),
        "receptor": task.receptor_name,
        "scoring": task.vina.scoring,
        "exhaustiveness": task.vina.exhaustiveness,
        "seed": task.vina.seed,
        "status": "ok",
        "error": None,
        "best_score": None,
        "pose_file": None,
        "n_poses": 0,
        "seconds": None,
    }
    t0 = time.perf_counter()
    try:
        v = _get_vina(task)
        v.set_ligand_from_file(str(lig))
        v.dock(
            exhaustiveness=task.vina.exhaustiveness,
            n_poses=task.vina.n_poses,
            min_rmsd=task.vina.min_rmsd,
            max_evals=task.vina.max_evals,
        )

        energies = np.asarray(v.energies(n_poses=task.vina.n_poses))
        if energies.size == 0:
            rec.update(status="no_poses", error="Vina returned no poses")
            return rec

        out_dir = Path(task.out_dir) / task.receptor_name
        out_dir.mkdir(parents=True, exist_ok=True)
        pose_file = out_dir / f"{lig.stem}_out.pdbqt"
        v.write_poses(str(pose_file), n_poses=task.vina.n_poses, overwrite=True)

        rec.update(
            best_score=float(energies[0][0]),
            all_scores=[float(e[0]) for e in energies],
            n_poses=int(energies.shape[0]),
            pose_file=str(pose_file),
        )
    except Exception as exc:  # noqa: BLE001
        rec.update(status="dock_failed", error=str(exc))
        log.warning(
            "Docking failed for %s on %s: %s", lig.stem, task.receptor_name, exc
        )
    finally:
        # Wall time per docking: the only honest basis for a --time estimate,
        # since Vina's cost is driven by the ligand's torsion count.
        rec["seconds"] = round(time.perf_counter() - t0, 2)
    return rec


def dock_many(
    ligand_pdbqts: list[str],
    receptor_record: dict,
    out_dir: str | Path,
    vina: VinaSpec,
    n_workers: int = 8,
) -> list[dict]:
    tasks = [
        DockTask(
            ligand_pdbqt=lp,
            receptor_name=receptor_record["name"],
            receptor_pdbqt=receptor_record["receptor_pdbqt"],
            center=tuple(receptor_record["box_center"]),
            size=tuple(receptor_record["box_size"]),
            out_dir=str(out_dir),
            vina=vina,
        )
        for lp in ligand_pdbqts
    ]

    n_workers = max(1, min(os.cpu_count() or 1, n_workers))
    log.info(
        "Docking %d ligand states into %s with %d workers "
        "(exhaustiveness=%d, scoring=%s)",
        len(tasks),
        receptor_record["name"],
        n_workers,
        vina.exhaustiveness,
        vina.scoring,
    )

    results: list[dict] = []
    # chunksize=1 so each worker reuses its cached Vina maps across many tasks
    with Pool(processes=n_workers) as pool:
        for rec in pool.imap_unordered(dock_one, tasks, chunksize=1):
            results.append(rec)
            if len(results) % 100 == 0:
                log.info("  docked %d/%d", len(results), len(tasks))
    return results


def _pose_blocks_pdbqt(path: str | Path) -> list[list[str]]:
    """Coordinate lines of each MODEL, in Vina's ranking order."""
    models: list[list[str]] = []
    cur: list[str] = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line[:6] in ("ATOM  ", "HETATM"):
                cur.append(line)
            elif line.startswith("ENDMDL") and cur:
                models.append(cur)
                cur = []
    if cur:
        models.append(cur)
    return models


def _heavy_atoms(lines: list[str], elem_col: int) -> np.ndarray:
    coords = [
        [float(l[30:38]), float(l[38:46]), float(l[46:54])]
        for l in lines
        if not (l[elem_col:elem_col + 2].strip()
                or l[12:14].strip()).upper().startswith("H")
    ]
    return np.asarray(coords, dtype=float)


def _heavy_atom_coords_pdb(path: str | Path) -> np.ndarray:
    with open(path, encoding="utf-8") as fh:
        lines = [l for l in fh if l[:6] in ("ATOM  ", "HETATM")]
    return _heavy_atoms(lines, 76)


def symmetry_naive_rmsd(a: np.ndarray, b: np.ndarray) -> float | None:
    if a.shape != b.shape or a.size == 0:
        return None
    return float(np.sqrt(((a - b) ** 2).sum(axis=1).mean()))


def redock_control(
    receptor_record: dict,
    out_dir: str | Path,
    vina: VinaSpec,
) -> dict:
    ref_pdbqt = receptor_record.get("reference_ligand_pdbqt")
    ref_pdb = receptor_record.get("reference_ligand_pdb")
    name = receptor_record["name"]

    if not ref_pdbqt:
        return {
            "receptor": name,
            "status": "skipped",
            "reason": "no reference ligand extracted",
        }

    out_dir = Path(out_dir) / "redock"
    out_dir.mkdir(parents=True, exist_ok=True)

    task = DockTask(
        ligand_pdbqt=ref_pdbqt,
        receptor_name=name,
        receptor_pdbqt=receptor_record["receptor_pdbqt"],
        center=tuple(receptor_record["box_center"]),
        size=tuple(receptor_record["box_size"]),
        out_dir=str(out_dir),
        vina=vina,
    )
    rec = dock_one(task)
    if rec["status"] != "ok":
        return {"receptor": name, "status": "failed", "detail": rec}

    crystal = _heavy_atom_coords_pdb(ref_pdb) if ref_pdb else None
    models = _pose_blocks_pdbqt(rec["pose_file"])

    # Every pose, not just Vina's top-ranked one. The literature convention is
    # that redocking succeeds when the BEST-SCORED pose is within 2 A, so that
    # stays the verdict -- but a native pose sitting at rank 3 means the box
    # and the receptor are right and only the ranking is off, which is a
    # completely different problem from a pose that is nowhere near the site.
    per_pose: list[float] = []
    for i, block in enumerate(models):
        naive = (
            symmetry_naive_rmsd(crystal, _heavy_atoms(block, 77))
            if crystal is not None
            else None
        )
        sym = None
        try:
            from spyrmsd import rmsd as sprmsd, molecule

            ref_mol = molecule.Molecule.from_rdkit(_rdkit_from_pdb(ref_pdb))
            pose_mol = molecule.Molecule.from_rdkit(_rdkit_from_pdb_block(block))
            ref_mol.strip()
            pose_mol.strip()
            sym = float(
                sprmsd.symmrmsd(
                    ref_mol.coordinates,
                    pose_mol.coordinates,
                    ref_mol.atomicnums,
                    pose_mol.atomicnums,
                    ref_mol.adjacency_matrix,
                    pose_mol.adjacency_matrix,
                )
            )
        except Exception as exc:
            log.debug("spyrmsd unavailable or failed on pose %d: %s", i + 1, exc)
        value = sym if sym is not None else naive
        if value is None:
            log.warning(
                "%s pose %d: RMSD not computed (crystal %s vs pose %d atoms)",
                name, i + 1, None if crystal is None else crystal.shape, len(block),
            )
        else:
            per_pose.append(value)

    final = per_pose[0] if per_pose else None
    best = min(per_pose) if per_pose else None
    best_rank = per_pose.index(best) + 1 if per_pose else None

    verdict = None
    if final is not None:
        verdict = "PASS" if final <= 2.0 else "FAIL"
        log.info("%s redocking RMSD = %.2f A -> %s "
                 "(best of %d poses: %.2f A at rank %d; bond orders: %s)",
                 name, final, verdict, len(per_pose), best, best_rank,
                 receptor_record.get("reference_bond_orders"))

    return {
        "receptor": name,
        "status": "ok",
        "best_score": rec["best_score"],
        "rmsd_reported": final,
        "rmsd_best": best,
        "best_pose_rank": best_rank,
        "n_poses_scored": len(per_pose),
        "rmsd_per_pose": per_pose,
        "bond_orders": receptor_record.get("reference_bond_orders"),
        "threshold": 2.0,
        "verdict": verdict,
        "pose_file": rec["pose_file"],
    }


def _rdkit_from_pdb(path):
    from rdkit import Chem

    return Chem.MolFromPDBFile(str(path), removeHs=True, sanitize=False)


def _rdkit_from_pdb_block(lines):
    """One pose's PDBQT coordinate lines as a mol.

    Truncated at column 66: past it a PDBQT carries the partial charge and the
    AutoDock type, which RDKit's PDB parser would read as occupancy/element.
    """
    from rdkit import Chem

    block = "".join(l[:66] + "\n" for l in lines) + "END\n"
    return Chem.MolFromPDBBlock(block, removeHs=True, sanitize=False)


def write_scores_csv(results: list[dict], path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    cols = [
        "ligand",
        "receptor",
        "scoring",
        "best_score",
        "n_poses",
        "seconds",
        "exhaustiveness",
        "seed",
        "status",
        "error",
        "pose_file",
        "ligand_file",
    ]
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in sorted(
            results,
            key=lambda x: (
                x["receptor"],
                x["best_score"] is None,
                x["best_score"] or 0.0,
            ),
        ):
            w.writerow(r)
    log.info("Wrote %d score rows to %s", len(results), path)
    return path
