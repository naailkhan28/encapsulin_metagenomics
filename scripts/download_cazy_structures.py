import urllib.request
from Bio import PDB
import pandas as pd
import os
import warnings
warnings.filterwarnings("ignore")


#Load CAZY DataFrame and get PDB codes
cazy_df = pd.read_csv("metadata/cazy_structures.csv")
cazy_df["PDB Code"] = cazy_df["PDB Code"].str.split(",").str[0]
cazy_df["PDB Code"] = cazy_df["PDB Code"].str.replace("[", "_").str.replace("]", "")
pdb_codes = cazy_df.sort_values(by="Resolution", ascending=True).groupby("Protein Name").first()["PDB Code"].unique()

#Check for output directory and create it if it doesn't exist
if not os.path.exists("structures/CAZY"):
    os.mkdir("structures/CAZY")

#Download each PDB file to the output folder
base_url = "https://files.rcsb.org/download/"
missing_structures = []

for pdb_code in pdb_codes:
    pdb_id = pdb_code.split("_")[0]
    chain_id = pdb_code.split("_")[1]
    
    print(f"Downloading {pdb_id}")
    try:
        urllib.request.urlretrieve(f"{base_url}{pdb_id}.pdb", f"structures/CAZY/{pdb_id}.pdb")
    except urllib.error.HTTPError:
        missing_structures.append(f"{pdb_code}")
        continue

    #Extract only the chain we need
    parser = PDB.PDBParser()
    structure = parser.get_structure("structure", f"structures/CAZY/{pdb_id}.pdb")

    extracted_structure = PDB.Structure.Structure("extracted_structure")
    extracted_model = PDB.Model.Model("new_model")
    
    try:
        for model in structure:
            for input_chain in model:
                if input_chain.id == chain_id:
                    extracted_model.add(input_chain)
                    break
    except PDB.PDBExceptions.PDBConstructionException:
        missing_structures.append(f"{pdb_code}")
        continue

    extracted_structure.add(extracted_model)
    io = PDB.PDBIO()
    io.set_structure(extracted_structure)
    io.save(f"structures/CAZY/{pdb_id}.pdb")

print(f"_________________Missing Structures:_________________")
print(missing_structures)