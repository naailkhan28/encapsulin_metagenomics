import pandas as pd

operon_df = pd.read_csv("operon_df_filtered.csv")

all_mgyps = []

for i in [-10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]:
    all_mgyps.extend(operon_df[str(i)].dropna().to_list())

all_mgyps = set(all_mgyps)

with open("metadata/cargo_MGYPs_updated.txt", "w") as outfile:
    outfile.write("\n".join(all_mgyps))