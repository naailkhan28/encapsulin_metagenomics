import pandas as pd
import os
from glob import glob
import subprocess
from Bio.PDB import PDBParser

def calculate_average_b_factor(pdb_file):
    parser = PDBParser()
    structure = parser.get_structure('protein', pdb_file)
    model = structure[0]
    atoms = []
    for chain in model:
        for residue in chain:
            for atom in residue:
                atoms.append(atom)
    b_factors = [atom.get_bfactor() for atom in atoms if atom.get_bfactor() != 0.0]
    average_b_factor = sum(b_factors) / len(b_factors)
    return average_b_factor

df_records = []

for mgyp_accession in os.listdir("AF2_predictions"):

    af2_filename = f"AF2_predictions/{mgyp_accession}/ranked_0.pdb"
    esmfold_filename = glob(f"structures/confident/{mgyp_accession}*.pdb")[0]

    try:
        command_output = subprocess.check_output(f"./TMalign '{esmfold_filename}' '{af2_filename}'", shell=True).decode().split("\n")
        tm_line = command_output[13]
        tm_score = float(tm_line.split("=")[1].split()[0])

    except IndexError:
        continue

    esmfold_plddt = calculate_average_b_factor(esmfold_filename)
    af2_plddt = calculate_average_b_factor(af2_filename)

    df_records.append({
        "Accession": mgyp_accession,
        "TM-Score": tm_score,
        "ESMFold pLDDT": esmfold_plddt,
        "AF2 pLDDT": af2_plddt
    })

outfile = pd.DataFrame(df_records).to_csv("metadata/AF2_vs_ESMFold_predictions.csv", index=False)