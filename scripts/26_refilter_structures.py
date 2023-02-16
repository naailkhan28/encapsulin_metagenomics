import os
import glob
import re

with open("final_filtered_mgyp_list.txt", "r") as mgypfile:
    mgyp_list = [line.rstrip() for line in mgypfile]

mgyp_pattern = re.compile(r"MGYP\d+")
structures = []

for file in os.listdir("structures/"):
    if file.endswith(".pdb"):
        structures.append(re.match(mgyp_pattern, file).group(0))

mgyp_list = set(mgyp_list)
structure_list = set(structures)
unneeded_structures = structure_list.difference(mgyp_list)


for structure in unneeded_structures:
    os.remove(glob.glob(f"structures/{structure}*.pdb")[0])
    try:
        os.remove(glob.glob(f"structures/pae/{structure}*.txt")[0])
    except IndexError:
        continue