from rdkit import Chem
from molscrub import Scrub
from meeko import MoleculePreparation, PDBQTWriterLegacy
import os
from time import sleep
import multiprocessing
from multiprocessing import Pool

def process_ligand(args: tuple[str, str, str, Scrub, int]) -> list[str]:
    smi, lname, outdir, scrubber_instance, max_scrub = args
    outdirs = []

    mol = Chem.MolFromSmiles(smi)
    if mol is None: 
        print(f"[MolFromSmiles] Can't parse SMILES for {lname}.")
        return []
    
    for i in range(max_scrub):
        try:
            scrubbed = list(scrubber_instance(mol))
            print(f"[Scrub] Scrub finished for {lname}.")
            break
        except RuntimeError as e:
            print(f"[Scrub] Scrub failed. Attempt {i + 1}/{max_scrub} for {lname}. Error: {e}.")
            sleep(.001)

    molprep = MoleculePreparation()
    count = 0

    for state in scrubbed:
        molsetup_list = molprep.prepare(state)
        for molsetup in molsetup_list:
            pdbqt, success, error = PDBQTWriterLegacy.write_string(molsetup)
            if success:
                path = os.path.join(outdir, f"{lname}_{count}.pdbqt")
                with open(path, "w") as f:
                    f.write(pdbqt)
                outdirs.append(path)
                count += 1
            else:
                print(f"[PDBQTWrite] Writing PDBQT failed for {lname}. Error: {error}.")

    return outdirs


def process_ligands(ligands_path: str, output_path: str, scrub_ph: float = 7.4, n_scrub_attempts: int = 5, n_threads: int = 8):
    n_threads = min(multiprocessing.cpu_count(), n_threads)

    scrub = Scrub(ph_low=scrub_ph, ph_high=scrub_ph, skip_tautomers=True)
    ligands_list = []

    os.mkdirs(output_path, exist_ok=True)

    with open(ligands_path) as f:
        for l in f.readlines():
            if len(l) > 2:
                smi, name, _ = l.split(" ")
                ligands_list.append(
                    (smi, name, output_path, scrub, n_scrub_attempts)
                )

    with Pool(processes=n_threads) as pool:
        for results in pool.imap_unordered(process_ligand, ligands_list):
            pass