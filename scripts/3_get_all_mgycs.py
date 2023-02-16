import pandas as pd
import numpy as np


metadata_df = pd.read_csv("metadata/mgy_seq_metadata_filtered.csv")
all_mgycs = set(metadata_df["MGYC"].to_list())

with open("metadata/all_MGYCs.txt", "w") as outfile:
    outfile.write("\n".join([mgyc for mgyc in all_mgycs]))