import urllib.request
from Bio import PDB
import os
import warnings
warnings.filterwarnings("ignore")



#Check for output directory and create it if it doesn't exist
if not os.path.exists("structures/CAZY"):
    os.mkdir("structures/CAZY")

parser = PDB.MMCIFParser()

#Download each PDB file to the output folder
base_url = "https://files.rcsb.org/download/"

with open("structures/CAZY/missing_structures.txt", "r") as infile:
    pdb_codes = [line.rstrip() for line in infile]

for pdb_code in pdb_codes:
    pdb_id = pdb_code.split("_")[0]
    chain_id = pdb_code.split("_")[1]
    
    print(f"Downloading {pdb_id}")
    urllib.request.urlretrieve(f"{base_url}{pdb_id}.cif", f"structures/CAZY/{pdb_id}.cif")

    #Extract only the chain we need
    structure = parser.get_structure("structure", f"structures/CAZY/{pdb_id}.cif")

    extracted_structure = PDB.Structure.Structure("extracted_structure")
    extracted_model = PDB.Model.Model("new_model")
    
    for model in structure:
        for input_chain in model:
            if input_chain.id == chain_id:
                extracted_model.add(input_chain)
                break
        else:
            continue
        break

    extracted_structure.add(extracted_model)
    io = PDB.PDBIO()
    io.set_structure(extracted_structure)
    io.save(f"structures/CAZY/{pdb_id}.pdb")
    os.remove(f"structures/CAZY/{pdb_id}.cif")