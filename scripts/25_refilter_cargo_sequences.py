import pandas as pd
from Bio import SeqIO


operon_df = pd.read_csv("operon_df_filtered.csv")

with open("final_filtered_mgyp_list.txt", "r") as mgypfile:
    mgyp_list = [line.rstrip() for line in mgypfile]


operon_df = operon_df.query("`Encapsulin MGYP` in @mgyp_list")
outfile = operon_df.to_csv("operon_df_filtered.csv", index=False)

all_cargos = []
for i in [-10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]:
    all_cargos.extend(operon_df[f"{i}"].to_list())

all_cargos = set(all_cargos)


all_records = []
for record in SeqIO.parse("seqs/all_putative_cargo_proteins.fasta", "fasta"):
    if str(record.id) in all_cargos:
        all_records.append(record)
    

outfile = SeqIO.write(all_records, "seqs/filtered_cargo_proteins.fasta", "fasta")