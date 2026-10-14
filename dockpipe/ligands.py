"""Ligand preparation: SMILES -> protomers/conformers -> PDBQT.

Built on molscrub (protonation + 3D) and Meeko (PDBQT writing), keeping the
structure of the original process_ligand/process_ligands pair but fixing the
issues that make it fail at scale:

  * os.mkdirs            -> os.makedirs
  * unbound `scrubbed`   -> explicit success flag, clean skip on failure
  * rigid 3-field split  -> tolerant parsing of .smi lines
  * Scrub pickled per task -> one Scrub built per worker process
"""
from __future__ import annotations

import json
import logging
import os
from dataclasses import dataclass
from multiprocessing import Pool
from pathlib import Path

from rdkit import Chem, RDLogger

from .common import read_checkpoint

log = logging.getLogger(__name__)
RDLogger.DisableLog("rdApp.*")

# One Scrub per worker process, created lazily in the initializer.
_SCRUB = None


@dataclass
class LigandTask:
    smiles: str
    name: str
    outdir: str
    max_attempts: int = 5


def _init_worker(ph: float, skip_tautomers: bool) -> None:
    """Pool initializer: build the Scrub once per process, not once per task."""
    global _SCRUB
    from molscrub import Scrub
    _SCRUB = Scrub(ph_low=ph, ph_high=ph, skip_tautomers=skip_tautomers)


def read_smiles(path: str | Path) -> list[tuple[str, str]]:
    """Read a .smi file. Accepts 'SMILES', 'SMILES name', or extra columns.

    Blank lines and '#' comments are skipped. Unnamed molecules get lig_<n>.
    """
    out: list[tuple[str, str]] = []
    with open(path, encoding="utf-8") as fh:
        for idx, raw in enumerate(fh):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            smi = parts[0]
            name = parts[1] if len(parts) > 1 else f"lig_{idx:05d}"
            # Filesystem-safe name
            name = "".join(c if c.isalnum() or c in "-_." else "_" for c in name)
            out.append((smi, name))
    return out


def process_ligand(task: LigandTask) -> dict:
    """Prepare one ligand. Returns a record describing what was written."""
    rec: dict = {
        "name": task.name,
        "smiles": task.smiles,
        "pdbqt_files": [],
        "n_states": 0,
        "status": "ok",
        "error": None,
    }

    mol = Chem.MolFromSmiles(task.smiles)
    if mol is None:
        rec.update(status="parse_failed", error="RDKit could not parse SMILES")
        log.warning("[%s] SMILES parse failed", task.name)
        return rec

    # --- scrub (protonation states + 3D embedding), with retries -----------
    scrubbed = None
    last_err = None
    for attempt in range(1, task.max_attempts + 1):
        try:
            scrubbed = list(_SCRUB(mol))
            break
        except (RuntimeError, ValueError) as exc:
            last_err = str(exc)
            log.debug("[%s] scrub attempt %d/%d failed: %s",
                      task.name, attempt, task.max_attempts, exc)

    if not scrubbed:
        rec.update(status="scrub_failed", error=last_err or "no states produced")
        log.warning("[%s] scrub failed after %d attempts", task.name, task.max_attempts)
        return rec

    # --- Meeko preparation + PDBQT ----------------------------------------
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    molprep = MoleculePreparation()
    os.makedirs(task.outdir, exist_ok=True)
    count = 0

    for state in scrubbed:
        try:
            for molsetup in molprep.prepare(state):
                pdbqt, success, error = PDBQTWriterLegacy.write_string(molsetup)
                if not success:
                    log.warning("[%s] PDBQT write failed: %s", task.name, error)
                    continue
                path = os.path.join(task.outdir, f"{task.name}_{count}.pdbqt")
                with open(path, "w", encoding="utf-8") as fh:
                    fh.write(pdbqt)
                rec["pdbqt_files"].append(path)
                count += 1
        except Exception as exc:  # noqa: BLE001 - one bad state must not kill the run
            log.warning("[%s] Meeko prepare failed for a state: %s", task.name, exc)

    rec["n_states"] = count
    if count == 0:
        rec.update(status="no_states_written", error="all Meeko/PDBQT writes failed")
    return rec


def process_ligands(
    ligands_path: str | Path,
    output_path: str | Path,
    ph: float = 7.4,
    skip_tautomers: bool = True,
    n_attempts: int = 5,
    n_workers: int = 8,
    checkpoint: str | Path | None = None,
) -> list[dict]:
    """Prepare every ligand in a .smi file. Returns one record per input."""
    output_path = Path(output_path)
    output_path.mkdir(parents=True, exist_ok=True)

    entries = read_smiles(ligands_path)
    log.info("Read %d ligands from %s", len(entries), ligands_path)

    # ligand_index.json is only written at the end, so without a checkpoint an
    # interrupted prep leaves nothing behind and `dock` then refuses to start.
    done = read_checkpoint(checkpoint, key="name") if checkpoint else {}
    if done:
        log.info("Resuming: %d ligands already prepared", len(done))

    tasks = [
        LigandTask(smiles=smi, name=name, outdir=str(output_path),
                   max_attempts=n_attempts)
        for smi, name in entries
        if name not in done
    ]

    n_workers = max(1, min(os.cpu_count() or 1, n_workers))
    records: list[dict] = list(done.values())
    fh = None
    if checkpoint:
        Path(checkpoint).parent.mkdir(parents=True, exist_ok=True)
        fh = open(checkpoint, "a", encoding="utf-8")
    try:
        if tasks:
            with Pool(
                processes=n_workers,
                initializer=_init_worker,
                initargs=(ph, skip_tautomers),
            ) as pool:
                for i, rec in enumerate(
                    pool.imap_unordered(process_ligand, tasks, chunksize=4), 1
                ):
                    records.append(rec)
                    if fh is not None:
                        fh.write(json.dumps(rec) + "\n")
                        fh.flush()
                    if i % 50 == 0 or i == len(tasks):
                        log.info("prepared %d/%d", i, len(tasks))
    finally:
        if fh is not None:
            fh.close()

    ok = sum(r["status"] == "ok" for r in records)
    states = sum(r["n_states"] for r in records)
    log.info("Ligand prep done: %d/%d ok, %d PDBQT states written",
             ok, len(records), states)

    index = output_path / "ligand_index.json"
    index.write_text(json.dumps(records, indent=2), encoding="utf-8")
    return records
