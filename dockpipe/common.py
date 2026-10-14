from __future__ import annotations

import hashlib
import json
import logging
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass, field, asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import yaml

LOG_FORMAT = "%(asctime)s %(levelname)-7s [%(processName)s] %(name)s: %(message)s"


def setup_logging(level: str = "INFO", logfile: Path | None = None) -> None:
    handlers: list[logging.Handler] = [logging.StreamHandler(sys.stderr)]
    if logfile is not None:
        logfile.parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(logfile, encoding="utf-8"))
    logging.basicConfig(
        level=getattr(logging, level.upper()),
        format=LOG_FORMAT,
        handlers=handlers,
        force=True,
    )


def sha256(path: str | Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while block := fh.read(chunk):
            h.update(block)
    return h.hexdigest()


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def read_checkpoint(path: str | Path, key: str = "ligand") -> dict[str, dict]:
    """Records from an append-only JSONL checkpoint, indexed by `key`.

    Both long stages append a flushed line per finished item so an interrupted
    run resumes instead of starting over. A record counts as done whatever its
    status: every failure mode here is deterministic (unparseable SMILES, a
    Vina exception on one ligand), so retrying only burns the ~30 min the
    docking cost in the first place. Callers log the failure count instead.
    """
    done: dict[str, dict] = {}
    path = Path(path)
    if not path.exists():
        return done
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            try:
                rec = json.loads(line)
                done[rec[key]] = rec
            except (json.JSONDecodeError, KeyError):
                # A killed job can leave a half-written final line. Everything
                # before it is still good, so skip it rather than start over.
                logging.getLogger(__name__).warning(
                    "Ignoring truncated checkpoint line in %s", path
                )
    return done


@dataclass
class BoxSpec:
    center: tuple[float, float, float] | None = None
    size: tuple[float, float, float] = (22.0, 22.0, 22.0)
    reference_ligand_resname: str | None = None
    padding: float = 0.0


@dataclass
class ReceptorSpec:
    name: str
    pdb_id: str | None = None
    local_pdb: str | None = None
    box: BoxSpec = field(default_factory=BoxSpec)
    keep_chains: list[str] | None = None
    keep_hetatms: list[str] = field(default_factory=list)
    redock_ligand_resname: str | None = None
    ph: float = 7.4
    notes: str = ""


@dataclass
class VinaSpec:
    scoring: str = "vina"
    exhaustiveness: int = 32
    n_poses: int = 9
    seed: int = 42
    cpu: int = 1
    min_rmsd: float = 1.0
    max_evals: int = 0


@dataclass
class PipelineConfig:
    run_name: str
    ligands_smi: str
    output_root: str
    receptors: list[ReceptorSpec]
    vina: VinaSpec = field(default_factory=VinaSpec)
    ligand_ph: float = 7.4
    skip_tautomers: bool = True
    n_scrub_attempts: int = 5
    n_workers: int = 8

    @staticmethod
    def from_yaml(path: str | Path) -> "PipelineConfig":
        raw = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
        receptors = []
        for r in raw.pop("receptors"):
            box_raw = r.pop("box", {}) or {}
            if box_raw.get("center") is not None:
                box_raw["center"] = tuple(float(x) for x in box_raw["center"])
            if box_raw.get("size") is not None:
                box_raw["size"] = tuple(float(x) for x in box_raw["size"])
            receptors.append(ReceptorSpec(box=BoxSpec(**box_raw), **r))
        vina = VinaSpec(**(raw.pop("vina", {}) or {}))
        return PipelineConfig(receptors=receptors, vina=vina, **raw)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

CODE_ROOT = Path(__file__).resolve().parent.parent


def _safe(cmd: list[str]) -> str | None:
    """Run cmd in the source tree, not the cwd, so git provenance is the code's."""
    try:
        return subprocess.run(
            cmd, capture_output=True, text=True, timeout=10, check=True,
            cwd=CODE_ROOT,
        ).stdout.strip()
    except Exception:
        return None


PACKAGES = (
    "rdkit", "meeko", "molscrub", "vina", "numpy", "scipy", "pandas",
    # these three change the numbers but are easy to forget:
    # openmm/pdbfixer do protonation, spyrmsd decides whether the redocking
    # RMSD is symmetry-corrected (it falls back silently if absent).
    "openmm", "pdbfixer", "spyrmsd",
    # meeko's own dependency; a missing gemmi surfaces as an import error from
    # meeko, which is confusing, so record it explicitly.
    "gemmi",
)

# Receptor PDBQT conversion shells out to whichever of these is on PATH, so the
# result depends on a binary that no package version captures.
EXTERNAL_TOOLS = ("mk_prepare_receptor.py", "prepare_receptor", "vina", "obabel")


def package_versions() -> dict[str, str]:
    out: dict[str, str] = {"python": sys.version.split()[0]}
    for mod in PACKAGES:
        try:
            m = __import__(mod)
            out[mod] = getattr(m, "__version__", "unknown")
        except Exception:
            out[mod] = "not-installed"
    return out


def external_tools() -> dict[str, str | None]:
    """Resolved paths of the CLI tools the pipeline shells out to."""
    return {t: shutil.which(t) for t in EXTERNAL_TOOLS}


def environment_provenance() -> dict[str, Any]:
    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "hostname": os.uname().nodename,
        "packages": package_versions(),
        "tools": external_tools(),
        # Module systems export these, so the loaded modules land in the
        # manifest without the script having to report them.
        "modules": os.environ.get("LOADEDMODULES") or None,
        "conda_prefix": os.environ.get("CONDA_PREFIX"),
        "git_commit": _safe(["git", "rev-parse", "HEAD"]),
        "git_dirty": bool(_safe(["git", "status", "--porcelain"])),
        "slurm": {
            k: v for k, v in os.environ.items()
            if k.startswith("SLURM_") and k in {
                "SLURM_JOB_ID", "SLURM_JOB_NAME", "SLURM_ARRAY_TASK_ID",
                "SLURM_ARRAY_JOB_ID", "SLURM_CPUS_PER_TASK", "SLURM_NNODES",
                "SLURM_JOB_PARTITION", "SLURM_SUBMIT_DIR",
            }
        },
    }


def write_manifest(
    path: str | Path,
    payload: dict[str, Any],
    files: Iterable[str | Path] = (),
) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    checksums = {}
    for f in files:
        f = Path(f)
        if f.is_file():
            checksums[str(f)] = sha256(f)
    doc = {
        "provenance": environment_provenance(),
        **payload,
        "checksums": checksums,
    }
    path.write_text(json.dumps(doc, indent=2, default=str), encoding="utf-8")
    return path
