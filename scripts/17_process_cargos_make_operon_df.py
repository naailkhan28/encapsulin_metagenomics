import pandas as pd

print("Loading DataFrames")
encapsulin_df = pd.read_csv("metadata/encapsulin_hits_contig_metadata.csv").drop(["ERZ_x", "Strand", "Contig", "MGYA", "Phage", "Contig Length"],axis=1).rename(columns={"MGYP": "Encapsulin MGYP", "Start": "Encapsulin Start", "End": "Encapsulin End"})
cargo_df = pd.read_csv("metadata/cargo_seq_metadata_filtered.csv").drop(["ERZ", "Strand"], axis=1).rename(columns={"MGYP": "Cargo MGYP", "Start": "Cargo Start", "End": "Cargo End"})

print("Merging DataFrames")
merged_df = encapsulin_df.merge(cargo_df, on="MGYC", how="right")
merged_df = merged_df.dropna(how="any")

all_records = []
grouped = merged_df.groupby(by=["Encapsulin MGYP", "MGYC"])

print("Retrieving Operons")
for name, group in grouped:    
    encapsulin_mgyp = group.iloc[0, 0]
    encapsulin_start = group.iloc[0, 2]

    dummy_record = [{"Cargo MGYP": encapsulin_mgyp, "Cargo Start": encapsulin_start}]
    new_group = pd.concat([group, pd.DataFrame(dummy_record)]).sort_values(by="Cargo Start")
    
    cargo_mgyps = new_group["Cargo MGYP"].to_list()
    
    for i, mgyp in enumerate(cargo_mgyps):
        if mgyp == encapsulin_mgyp:
            mgyps_before = cargo_mgyps[max(0, i-10):i]
            mgyps_after = cargo_mgyps[i+1:min(i+11, len(cargo_mgyps))]

            desired_mgyps = mgyps_before + mgyps_after
            indices = [-1 * x for x in range(len(mgyps_before), 0, -1)] + list(range(1, len(mgyps_after) + 1))

            record = {index: mgyp for index, mgyp in zip(indices, desired_mgyps)}
            record["Encapsuling MGYP"] = encapsulin_mgyp
            all_records.append(record)

operon_df = pd.DataFrame(all_records)

print("Writing output CSV...")
outfile = operon_df.to_csv("metadata/operon_df.csv", index=False)