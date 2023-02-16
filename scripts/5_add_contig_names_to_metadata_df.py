import pandas as pd
import numpy as np

print("Loading DataFrames")
metadata_df = pd.read_csv("metadata/mgy_seq_metadata_filtered.csv")
contig_df = pd.read_csv("metadata/mgy_contigs_filtered.tsv", sep="\t", names=["MGYC", "Contig"])

print("Processing Contigs")
contig_df["Contig"] = contig_df["Contig"].str.extract(r"ERZ\d+\.(.+)", expand=True)

print("Merging Tables")
metadata_df = metadata_df.merge(contig_df, on="MGYC", how="left")

print(metadata_df[metadata_df["Contig"].isna()])

print("Writing Output File")
outfile = metadata_df.to_csv("metadata/mgy_seq_metadata_with_contigs_filtered.csv")