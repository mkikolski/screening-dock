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

-- conda_root is only the miniforge install; the env itself need not live under
-- it. /opt/modules is usually admin-owned, so if `conda env create` there is
-- refused, build the env in your own storage and swap env_root for the second
-- line below. Nothing else in this file changes.
local conda_root = "/opt/modules/miniforge3"
-- The env lives on shared storage, not under miniforge3/envs or in $HOME.
local env_root = "/storage/users/kikolski/envs"
local env_name = "dockpipe"
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
