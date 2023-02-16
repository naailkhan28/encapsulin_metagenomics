import os
import numpy as np

confident = []

for pdb_file in os.listdir("structures/"):
    if pdb_file.endswith(".pdb"):
        with open(f"structures/{pdb_file}", "r") as structure_file:
            all_plddts = []

            for line in structure_file:
                line = line.rstrip().split()

                if line[0] != "ATOM":
                    continue
                
                plddt = float(line[-2])
                if plddt >= 1:
                    plddt = plddt / 100
                all_plddts.append(plddt)
        
        if np.mean(all_plddts) >= 0.7:
            confident.append(pdb_file)

print(len(confident))

with open("structures/confident_predicted_structures.txt", "w") as outfile:
    outfile.write("\n".join(confident))