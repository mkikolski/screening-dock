help([[
AutoDock Vina screening pipeline (dockpipe).

Provides the dockpipe Python environment: RDKit, Meeko, molscrub, the Vina
Python bindings, OpenMM/PDBFixer and spyrmsd.

Docking uses the Vina *Python API* (`from vina import Vina`), not the vina
executable, so a CLI-only vina module is not a substitute for this one.

Sets DOCKPIPE_PROJECT, which run_all.slurm picks up as its default project
path; override it per run if the checkout lives elsewhere.
]])

whatis("Name: dockpipe")
whatis("Version: phase1")
whatis("Category: Python / Molecular docking")
whatis("Description: AutoDock Vina panel screening with Meeko ligand prep")

-- Any two conda-env modules fight over PATH and CONDA_PREFIX, so they must not
-- be loaded together. For this guard to be mutual the other conda-env modules
-- need the same family() line; otherwise only an explicit conflict() applies.
family("conda_env")
conflict("python_md")

-- conda_root is only the miniforge install; the env itself does not live under
-- it, since /opt/modules is admin-owned and `conda env create` there is
-- refused. Change env_root below if the env moves.
local conda_root = "/opt/modules/miniforge3"
-- Resolved from $HOME so this does not go stale when the env moves, and so
-- the file is not tied to one account. Override with DOCKPIPE_ENV and use
-- run_all.slurm instead if the env lives somewhere else entirely; the slurm
-- script probes several locations, a modulefile cannot.
local env_name = "dockpipe"
local env_root = pathJoin(os.getenv("HOME") or "", ".conda/envs")
local env_path = pathJoin(env_root, env_name)
local project_path = "/storage/users/kikolski/screening-dock"

prepend_path("PATH", pathJoin(env_path, "bin"))
prepend_path("PATH", pathJoin(conda_root, "condabin"))

setenv("CONDA_PREFIX", env_path)
setenv("CONDA_DEFAULT_ENV", env_name)
setenv("CONDA_EXE", pathJoin(conda_root, "bin", "conda"))
setenv("CONDA_ENVS_PATH", env_root)

setenv("DOCKPIPE_PROJECT", project_path)
setenv("PYTHONNOUSERSITE", "1")
setenv("PYTHONUNBUFFERED", "1")

-- Thread limits are deliberately NOT set here: Vina runs single-threaded with
-- parallelism coming from the worker Pool, which is a property of the job, not
-- of the environment. run_all.slurm exports OMP_NUM_THREADS=1 itself.
