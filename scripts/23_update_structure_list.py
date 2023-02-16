from Bio import SeqIO
import os
import re

with open("final_filtered_mgyp_list.txt", "r") as mgyp_file:
    all_mgyps = [line.rstrip() for line in mgyp_file]

mgyp_pattern = re.compile(r"MGYP\d+")
mgyps_with_structure = []

for file in os.listdir("structures/"):
    if file.endswith(".pdb"):
        mgyps_with_structure.append(re.match(mgyp_pattern, file).group(0))

mgyps_with_structure = set(mgyps_with_structure)
all_mgyps = set(all_mgyps)

missing_mgyps = all_mgyps.difference(mgyps_with_structure)

all_records = []

for record in SeqIO.parse("seqs/encapsulin_hits_filtered.fasta", "fasta"):
    if str(record.id) in missing_mgyps:
        all_records.append(record)

outfile = SeqIO.write(all_records, "seqs/encapsulin_hits_for_structure_prediction.fasta", "fasta")