import pandas as pd
import os

all_records = []

for folder in os.listdir("deepbgc"):
    try:
        df = pd.read_csv(f"deepbgc/{folder}/{folder}.bgc.tsv", sep="\t")
    except FileNotFoundError:
        continue
    
    all_records.extend(df.to_dict(orient="records"))

cluster_df = pd.DataFrame(all_records).rename(columns={"sequence_id": "Contig", "nucl_start": "Cluster Start", "nucl_end": "Cluster End", "nucl_length": "Cluster Length"})

#Load Encapsulin Metadata Table and filter
print("Loading Metadata Table")
metadata_df = pd.read_csv("metadata/mgy_seq_metadata_with_contigs_mgya_filtered_nonans.csv").rename(columns={"Start": "Encapsulin Start", "End": "Encapsulin End"})

print("Filtering Metadata DataFrame")
with open("final_filtered_mgyp_list.txt", "r") as cargo_file:
    mgyp_list = [line.rstrip() for line in cargo_file]

metadata_df = metadata_df.query("MGYP in @mgyp_list")

#Add encapsulin metadata to deepBGC cluster DataFrame
metadata_df["Contig"] = metadata_df[["ERZ", "Contig"]].agg(".".join, axis=1)

cluster_df = pd.merge(left=cluster_df, right=metadata_df, on="Contig")
cluster_df["Hit?"] = (cluster_df["Cluster Start"] < cluster_df["Encapsulin Start"]) & (cluster_df["Cluster End"] > cluster_df["Encapsulin End"])
cluster_df = cluster_df[cluster_df["Hit?"] == True]

#Remove any already annotated encapsulins
annotated_encapsulins = list(pd.read_csv("encapsulin_families.csv")["Encapsulin MGYP"].unique())
cluster_df = cluster_df.query("MGYP not in @annotated_encapsulins")

outfile = cluster_df.to_csv("metadata/deepbgc_hits.csv", index=False)