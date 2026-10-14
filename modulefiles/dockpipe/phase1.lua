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

-- The env lives inside the shared miniforge install, which every compute node
-- mounts. $HOME does not travel to the node a job lands on, so an env under
-- ~/.conda is unusable from a batch job even though it works interactively.
local conda_root = "/opt/modules/miniforge3"
local env_root = pathJoin(conda_root, "envs")
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
