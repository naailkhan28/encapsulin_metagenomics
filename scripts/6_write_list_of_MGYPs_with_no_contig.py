import pandas as pd
import numpy as np

metadata_df = pd.read_csv("metadata/mgy_seq_metadata_with_contigs_filtered.csv")
metadata_df = metadata_df.iloc[:, 1:]

orphan_df = metadata_df[metadata_df["Contig"].isna()]
metadata_df = metadata_df.dropna(subset="Contig")

print(len(orphan_df))

metadata_outfile = metadata_df.to_csv("metadata/mgy_seq_metadata_with_contigs_filtered_nonans.csv", index=False)
orphan_outfile = orphan_df.to_csv("metadata/orphan_sequences.csv", index=False)