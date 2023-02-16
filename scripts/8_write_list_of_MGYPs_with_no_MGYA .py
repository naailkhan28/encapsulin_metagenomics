import pandas as pd
import numpy as np

metadata_df = pd.read_csv("metadata/mgy_seq_metadata_with_contigs_filtered_nonans.csv")
mgya_df = pd.read_csv("metadata/mgya_mappings.csv", names=["ERZ", "MGYA"])

print("Merging Tables")
metadata_df = metadata_df.merge(mgya_df, on="ERZ", how="left")
metadata_df = metadata_df.dropna(subset="MGYA")

metadata_outfile = metadata_df.to_csv("metadata/mgy_seq_metadata_with_contigs_mgya_filtered_nonans.csv", index=False)