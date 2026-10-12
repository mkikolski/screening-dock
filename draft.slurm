#!/usr/bin/env bash
#SBATCH --export=ALL
#SBATCH --job-name=docking_c5_selected
#SBATCH --output=docking_c5_v1_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=30G
#SBATCH --gres=gpu:2
#SBATCH --time=48:00:00
#SBATCH --partition=phd
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mikolaj.kikolski@pwr.edu.pl

set -euo pipefail

PROJECT="/storage/kikolski/selected_docking"
ENV_DEF="$PROJECT/environment.yml"
ENV_NAME="selected-docking"
SMI_FILE="$PROJECT/input_files/exp_c5_selected_500.smi"
RECEPTOR_PDB="$PROJECT/input_files/receptor.pdb"
RECEPTOR_LIG="$PROJECT/input_files/receptor.lig"
PDBQT_DIR="$PROECT/pdbqt_mols"
TEMP_DIR="$PROJECT/intermediate_files"

export CONDA_EXE=/opt/modules/miniforge3/bin/conda
export CONDA_ENVS_PATH=/opt/modules/miniforge3/envs

cd "$PROJECT"
mkdir $TEMP_DIR

conda env create --name "$ENV_NAME" --file="$ENV_DEF"
conda activate "$ENV_NAME"

scrub.py $SMI_FILE -o $TEMP_DIR/mols.sdf
mk_prepare_ligand.py -i $TEMP_DIR/mols.sdf --multimol_outdir $PDBQT_DIR

