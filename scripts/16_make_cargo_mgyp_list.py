import pandas as pd

metadata_df = pd.read_csv("metadata/cargo_seq_metadata_filtered.csv")
all_mgyps = set(metadata_df["MGYP"].to_list())

with open("metadata/cargo_MGYPs.txt", "w") as outfile:
    outfile.write("\n".join([mgyp for mgyp in all_mgyps]))