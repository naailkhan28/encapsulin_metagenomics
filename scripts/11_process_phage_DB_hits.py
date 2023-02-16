import pandas as pd
import re

hits_df = pd.read_csv("metadata/contig_phage_hits.m8", sep="\t", names=["Query", "Target", "% Identity", "Alignment Length", "Mismatches",
                        "Gap Openings", "Query Start", "Query End", "Target Start", "Target End", "E-Value", "Bit Score"])

metadata_df = pd.read_csv("metadata/final_mgy_seq_metadata_table.csv")

erz_regex = re.compile(r"(ERZ\d+)\.")
contig_regex = re.compile(r"\.(.+)_")

hits_df["ERZ"] = hits_df["Query"].str.extract(erz_regex)
hits_df["Contig"] = hits_df["Query"].str.extract(contig_regex)
hits_df["Phage"] = True

subset_hits_df = hits_df.loc[:, ["ERZ", "Contig", "Phage"]]
subset_hits_df = subset_hits_df.drop_duplicates()

phage_df = metadata_df.merge(subset_hits_df, how="left", on="Contig")
phage_df["Phage"] = phage_df["Phage"].fillna(False)

phage_mgyps = set(phage_df[phage_df["Phage"] == True]["MGYP"].to_list())
non_phage_mgyps = set(phage_df[phage_df["Phage"] == False]["MGYP"].to_list())

non_phage_mgyps = non_phage_mgyps.difference(phage_mgyps)

print(f"Phage-associated MGYPs: {len(phage_mgyps)}")
print(f"Encapsulin MGYPs: {len(non_phage_mgyps)}")

outfile = phage_df.to_csv("metadata/phage_hits_df.csv", index=False)