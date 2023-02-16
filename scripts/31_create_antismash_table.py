import pandas as pd
import re
import os
from bs4 import BeautifulSoup

contig_pattern = re.compile(r"original\sname\swas:\s(.+)\.\.\.")
region_pattern = re.compile(r"Region&nbsp(.+)")
contigs = os.listdir("antiSMASH/")

cluster_records = []

#Parse all antiSMASH HTML file outputs
print("Parsing antiSMASH HTML Files...")
for contig in contigs:
    try:
        with open(f"antiSMASH/{contig}/index.html", "r") as infile:
            soup = BeautifulSoup(infile, "html.parser")
    except FileNotFoundError:
        continue
    
    if soup.find_all(class_="overview-layout"):
        region_names = soup.find_all("div", class_="record-overview-header")
        tables = soup.find_all("table", class_="region-table")

        for region, table in zip(region_names, tables):

            name = re.search(contig_pattern, region.text)
            try:
                name = name.group(1)
            except AttributeError:
                continue

            table_body = table.find("tbody")
            rows = table_body.find_all("tr")
            for row in rows:
                cols = row.find_all("td")
                cols = [ele.text.strip() for ele in cols]
                try:
                    cluster_records.append({"Partial Contig": name, "Region": re.findall(region_pattern, cols[0])[0], 
                                            "Cluster Type": cols[1], 
                                            "Cluster Start": int(cols[2].replace(",", "")), 
                                            "Cluster End": int(cols[3].replace(",", "")), 
                                            "Closest Match": cols[4], "Identity": int(cols[6][:-1]) / 100})
                except IndexError:
                    cluster_records.append({"Partial Contig": name, "Region": re.findall(region_pattern, cols[0])[0], 
                                            "Cluster Type": cols[1], 
                                            "Cluster Start": int(cols[2].replace(",", "")), 
                                            "Cluster End": int(cols[3].replace(",", "")), 
                                            "Closest Match": cols[4], "Identity": ""})

cluster_df = pd.DataFrame(cluster_records)

#Load Encapsulin Metadata Table and filter
print("Loading Metadata Table")
metadata_df = pd.read_csv("metadata/mgy_seq_metadata_with_contigs_mgya_filtered_nonans.csv").rename(columns={"Start": "Encapsulin Start", "End": "Encapsulin End"})

print("Filtering Metadata DataFrame")
with open("final_filtered_mgyp_list.txt", "r") as cargo_file:
    mgyp_list = [line.rstrip() for line in cargo_file]

metadata_df = metadata_df.query("MGYP in @mgyp_list")

#Add encapsulin metadata to antiSMASH cluster DataFrame
metadata_df["Contig"] = metadata_df[["ERZ", "Contig"]].agg(".".join, axis=1)
partial_contigs = list(cluster_df["Partial Contig"].unique())
full_contigs = list(metadata_df["Contig"].unique())

contigs_mapping_dict = {}

for partial_contig in partial_contigs:
    for full_contig in full_contigs:
        if full_contig.startswith(partial_contig):
            contigs_mapping_dict[partial_contig] = full_contig
            break

def get_full_contig(partial_contig):
    return(contigs_mapping_dict[partial_contig])

cluster_df["Contig"] = cluster_df["Partial Contig"].apply(get_full_contig)
cluster_df = cluster_df.drop("Partial Contig", axis=1)
cluster_df = pd.merge(left=cluster_df, right=metadata_df, on="Contig")

cluster_df["Hit?"] = (cluster_df["Cluster Start"] < cluster_df["Encapsulin Start"]) & (cluster_df["Cluster End"] > cluster_df["Encapsulin End"])
cluster_df = cluster_df[cluster_df["Hit?"] == True]
cluster_df = cluster_df.loc[:, ["MGYP", "MGYA", "ERZ", "Region", "Encapsulin Start", "Encapsulin End", "Cluster Start", "Cluster End", "Cluster Type", "Closest Match", "Identity"]]

#Remove any already annotated encapsulins
annotated_encapsulins = list(pd.read_csv("encapsulin_families.csv")["Encapsulin MGYP"].unique())
cluster_df = cluster_df.query("MGYP not in @annotated_encapsulins")

outfile = cluster_df.to_csv("metadata/encapsulin_antiSMASH_predictions.csv", index=False)