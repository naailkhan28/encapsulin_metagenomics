from urllib.request import urlretrieve
from urllib.error import HTTPError
from Bio import SeqIO
import ssl
import os

#We need to disable SSL due to a certificate error with the ESM Atlas API
#See https://github.com/facebookresearch/esm/discussions/627
ssl._create_default_https_context = ssl._create_unverified_context

#Check for output directory and create it if it doesn't exist
if not os.path.exists("structures/glycosyltransferases"):
    os.mkdir("structures/glycosyltransferases")

with open("metadata/saccharide_BGC_glycosyltransferases.txt", "r") as mgyp_file:
    mgyp_list = [line.rstrip() for line in mgyp_file]

missing_mgyps = []

for i, mgyp in enumerate(mgyp_list):
    if i == 0 or (i + 1) % 5 == 0:
        print(f"Downloading structure {i+1}/{len(mgyp_list)}")
    try:
        request = urlretrieve(f"https://api.esmatlas.com/fetchPredictedStructure/{mgyp}", f"structures/glycosyltransferases/{mgyp}.pdb")
    except HTTPError:
        missing_mgyps.append(mgyp)
        continue

with open("structures/CAZY/missing_mgyps.txt", "w") as missing_mgyps_file:
    missing_mgyps_file.write("\n".join(missing_mgyps))

print(f"Writing missing sequences to a file")
all_records = []

for record in SeqIO.parse("seqs/filtered_cargo_proteins.fasta", "fasta"):
    if str(record.id) in missing_mgyps:
        all_records.append(record)

outfile = SeqIO.write(all_records, "seqs/glycosyltransferases_for_structure_prediction.fasta", "fasta")