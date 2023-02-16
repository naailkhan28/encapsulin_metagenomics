import subprocess
import os

all_pdb_files = os.listdir("structures/confident")
all_pdb_files = [file for file in all_pdb_files if file.endswith(".pdb")]

ids = list(range(len(all_pdb_files)))
ids = [str(id).rjust(4, "0") for id in ids]

id_mapping_dict = {}

for id, pdb_file in zip(ids, all_pdb_files):
   id_mapping_dict[id] = pdb_file
   subprocess.run(["perl", "../DALI/DaliLite.v5/bin/import.pl", "--pdbfile", f"structures/confident/{pdb_file}", "--pdbid", id, "--dat", "DALI/data/", "--clean"])


with open("DALI/chain_id_mapping.txt", "w") as outfile:
    for id, pdb_file in id_mapping_dict.items():
        outfile.write(f"{id}A\t{pdb_file}\n")