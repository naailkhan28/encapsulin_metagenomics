import pandas as pd
import numpy as np

print("Processing Input File")

records = []

with open("metadata/cargo_metadata_filtered.tsv", "r") as metadatafile:
    for i, line in enumerate(metadatafile):
        if i == 0 or (i + 1) % 500 == 0:
            print(f"Reading line {i+1}")
        
        line = line.rstrip().split()

        mgyp = line[0]
        metadata_records = line[1].split(";")

        for record in metadata_records:
            record = record.split(".")

            erz = record[0]

            non_erz_data = record[1].split(":")

            mgyc = non_erz_data[0]
            start = non_erz_data[1].split("-")[0]
            end = non_erz_data[1].split("-")[1]
            strand = non_erz_data[2]

            records.append({"MGYP": mgyp, 
                            "ERZ": erz, 
                            "MGYC": mgyc, 
                            "Start": start, 
                            "End": end, 
                            "Strand": strand})

print("Creating DataFrame")
metadata_df = pd.DataFrame(records)

with open("metadata/final_mgyp_list_curated.txt", "r") as mgypfile:
    mgyp_list = [line.rstrip() for line in mgypfile]

metadata_df = metadata_df.query("MGYP not in @mgyp_list")

print("Writing Output File")
outfile = metadata_df.to_csv("metadata/cargo_seq_metadata_filtered.csv", index=False)