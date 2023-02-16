import pandas as pd

print("Processing empty files list")
empty_files_records = []
with open("metadata/processed_empty_file_contigs.txt", "r") as empty_files_listfile:
    for line in empty_files_listfile:
        line = line.split("-")
        empty_files_records.append({"MGYA": line[0], "ERZ": line[1][:-11]})

print("Loading tables")
metadata_df = pd.read_csv("metadata/mgy_seq_metadata_with_contigs_mgya_filtered_nonans.csv")
orphan_df = pd.read_csv("metadata/orphan_sequences.csv")
missing_accessions = pd.read_csv("metadata/all_unknown_accessions.txt", sep="\t", names=["MGYA", "ERZ"])
missing_accessions = pd.concat([missing_accessions, pd.DataFrame(empty_files_records)])
missing_accessions = missing_accessions[missing_accessions["MGYA"] != "None"].drop_duplicates()

print("Creating new orphan table")
new_missing_values = metadata_df.merge(missing_accessions, on=["MGYA", "ERZ"])

orphan_df = pd.concat([orphan_df, new_missing_values]).drop_duplicates()

print("Filtering original metadata table")
orphan_mgyas = orphan_df["MGYA"].unique()
orphan_erzs = orphan_df["ERZ"].unique()
metadata_df = metadata_df.query("MGYA not in @orphan_mgyas")
metadata_df = metadata_df.query("ERZ not in @orphan_erzs")

print("Writing output tables")
orphan_df_outfile = orphan_df.to_csv("metadata/final_orphan_sequences.csv", index=False)
metadata_df_outfile = metadata_df.to_csv("metadata/final_mgy_seq_metadata_table.csv", index=False)

mgyps = len(metadata_df["MGYP"].unique())
orphans = len(orphan_df["MGYP"].unique())

print(f"MGYPs with contig: {mgyps}")
print(f"Orphan MGYPs with no contig data: {orphans}")