import pandas as pd
import requests

metadata_df = pd.read_csv("metadata/mgy_seq_metadata_with_contigs_filtered_nonans.csv")

all_erzs = list(metadata_df["ERZ"].unique())
mgya_dict = {}

def get_mgya_accession(erz):
    url = f"https://www.ebi.ac.uk/metagenomics/api/v1/assemblies/{erz}/analyses"

    api_request = requests.get(url=url)
    api_request = api_request.json()

    try:
        mgya = api_request["data"][0]["id"]
    except IndexError:
        return(None)     

    return(mgya)

for i, erz in enumerate(all_erzs):
    if i == 0 or (i + 1) % 50 == 0:
        print(f"Requesting MGYA {i} / {len(all_erzs)}")
    
    try:
        mgya_dict[erz] = get_mgya_accession(erz)
    except ValueError:
        continue

with open("metadata/mgya_mappings.csv", "w") as outfile:
    outfile.write("\n".join([f"{erz},{mgya}" for erz, mgya in mgya_dict.items()]))
