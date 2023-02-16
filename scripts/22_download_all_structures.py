from urllib.request import urlretrieve
from urllib.error import HTTPError
from Bio import SeqIO

with open("final_filtered_mgyp_list.txt", "r") as mgyp_file:
    mgyp_list = [line.rstrip() for line in mgyp_file]

missing_mgyps = []

for i, mgyp in enumerate(mgyp_list):
    if i == 0 or (i + 1) % 500 == 0:
        print(f"Downloading structure {i+1}/{len(mgyp_list)}")
    try:
        request = urlretrieve(f"https://api.esmatlas.com/fetchPredictedStructure/{mgyp}", f"structures/{mgyp}.pdb")
    except HTTPError:
        missing_mgyps.append(mgyp)
        continue

with open("structures/missing_mgyps.txt", "w") as missing_mgyps_file:
    missing_mgyps_file.write("\n".join(missing_mgyps))

all_records = []

for record in SeqIO.parse("seqs/encapsulin_hits_filtered.fasta", "fasta"):
    if str(record.id) in missing_mgyps:
        all_records.append(record)

outfile = SeqIO.write(all_records, "seqs/encapsulin_hits_for_structure_prediction.fasta", "fasta")