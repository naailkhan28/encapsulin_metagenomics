import pandas as pd
import numpy as np

all_records = []

columns_dict = {
    0: "Biome Category",
    1: "Biome Type",
    2: "Biome Subtype",
    3: "Ecosystem Type",
    4: "Specific Ecosystem"
}

with open("biomes/encapsulin_hit_biomes_filtered.tsv", "r") as infile:
    for line in infile:
        line = line.rstrip().split("\t")

        mgyp = line[0]
        biomes = line[2].split(":")[1:]

        record = {}

        record["MGYP"] = mgyp
        for i, biome in enumerate(biomes):
            record[columns_dict[i]] = biome
        
        all_records.append(record)

biome_df = pd.DataFrame(all_records)

outfile = biome_df.to_csv("biomes/processed_biome_dataframe.csv", index=False)