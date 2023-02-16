import pandas as pd
import gzip
import json
from collections import defaultdict
from Bio import SeqIO

print("Loading DataFrames")
operon_df = pd.read_csv("metadata/operon_df.csv")
pfam_df = pd.read_csv("pfams/cargo_pfams.tsv", sep="\t", names=["MGYP", "Pfam", "Start", "End"])

print("Mapping Pfams...")
with open("DBs/label_descriptions.json.gz", 'rb') as f:
    with gzip.GzipFile(fileobj=f, mode='rb') as gzip_file:
      labels_dict = json.load(gzip_file)

def get_label(pfam):
    try:
        return(labels_dict[pfam])
    except KeyError:
        return(None)

def get_pfams(mgyp):
    try:
        return(pfam_mapping_dict[mgyp])
    except KeyError:
        return None

pfam_df["Description"] = pfam_df["Pfam"].apply(get_label)

pfam_mapping_dict = defaultdict(list)
pfam_df_records = pfam_df.to_dict(orient="records")

for record in pfam_df_records:
    pfam_mapping_dict[record["MGYP"]].append(record["Pfam"])

for i in [-10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]:
    operon_df[f"Pfam {i}"] = operon_df[str(i)].apply(get_pfams)

print("Filtering out phage operons...")
with open("pfams/phage_pfams.txt", "r") as phage_pfams_file:
    phage_pfams = [line.rstrip().split()[0] for line in phage_pfams_file]
pfam_records = operon_df.to_dict(orient="records")

phage_mgyps = []
for record in pfam_records:
    all_pfams = []
    for i in [-10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]:
        all_pfams.extend(record[f"Pfam {i}"])

    for pfam in all_pfams:
        if pfam in phage_pfams:
            phage_mgyps.append(record["Encapsulin MGYP"])
            break

phage_mgyps = list(set(phage_mgyps))
operon_df = operon_df.query("`Encapsulin MGYP` not in @phage_mgyps")

non_phage_mgyps = list(operon_df["Encapsulin MGYP"].unique())

print(len(non_phage_mgyps))

print("Writing output DataFrame")
operon_outfile = operon_df.to_csv("metadata/operon_df_filtered.csv", index=False)

print("Writing MGYP list")
with open("metadata/final_filtered_mgyp_list.txt", "w") as mgyp_file:
    mgyp_file.write("\n".join(non_phage_mgyps))

print("Re-writing output sequences")
all_records = []

for record in SeqIO.parse("seqs/all_encapsulin_hits.fasta", "fasta"):
    if str(record.id) in non_phage_mgyps:
        all_records.append(record)
    
sequences_outfile = SeqIO.write(all_records, "seqs/encapsulin_hits_filtered.fasta", "fasta")