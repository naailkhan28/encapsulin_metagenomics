import pandas as pd
from collections import defaultdict
import re

metadata_df = pd.read_csv("metadata/final_mgy_seq_metadata_table.csv")

with open("metadata/final_mgyp_list_curated.txt", "r") as mgypfile:
    mgyp_list = [line.rstrip() for line in mgypfile]

metadata_df = metadata_df.query("MGYP in @mgyp_list")


contig_pattern = re.compile(r"ERZ\d+\.[\w-]+\d\.\d+")
cds_counts = defaultdict(int)

with open("seqs/all_search_contigs.fasta", "r") as CDSes:
    text = "".join(CDSes.readlines())
    data = re.findall(contig_pattern, text)

for contig in data:
    cds_counts[contig] += 1

cds_count_records = []
for contig, count in cds_counts.items():
    cds_count_records.append({"Contig": contig, "Number of CDSes": count})

cds_count_records = pd.DataFrame(cds_count_records)
metadata_df["Contig"] = metadata_df["ERZ"] + "." + metadata_df["Contig"]
metadata_df = metadata_df.merge(cds_count_records, on="Contig")
metadata_df = metadata_df[metadata_df["Number of CDSes"] > 10]

all_mgycs = set(metadata_df["MGYC"].to_list())

with open("metadata/cargo_MGYCs.txt", "w") as outfile:
    outfile.write("\n".join([mgyc for mgyc in all_mgycs]))