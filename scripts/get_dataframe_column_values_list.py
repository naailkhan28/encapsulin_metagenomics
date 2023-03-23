import pandas as pd

df = pd.read_csv("encapsulin_UniRef90_hits.m8", sep="\t", names=["Query", "Target", "Identity", "Alignment Length", "Mismatches", "Gap Openings",
                                                              "Query Start", "Query End", "Target Start", "Target End", "E-Value", "Bitscore"])

IDs = df["Target"].unique()

with open("encapsulin_UniRef90_hits_IDs.txt", "w") as outfile:
    outfile.write("\n".join(list(IDs)))