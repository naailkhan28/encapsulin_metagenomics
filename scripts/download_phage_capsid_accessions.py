import requests
from collections import defaultdict
import re

def get_uniprot_accessions(url, prev_accessions=defaultdict(list), count=0):
    print(f"Making request #{count+1} at url {url}")
    print(f"Currently have {len(prev_accessions)} accessions")
    response_data = requests.get(url).json()

    if response_data["next"]: 
        for item in response_data["results"]:
            accession = item["metadata"]["accession"]
            taxid = item["metadata"]["source_organism"]["taxId"]
            prev_accessions[taxid].append(accession)

        return(get_uniprot_accessions(response_data["next"], prev_accessions, count+1))
    else:
        return(prev_accessions)
    

phage_coat_url = "http://www.ebi.ac.uk/interpro/api/protein/UniProt/set/pfam/CL0373/?page_size=100"
all_cl0373_accessions = get_uniprot_accessions(phage_coat_url)

taxid_pattern = re.compile(r"\(taxid=(\d+),")
with open("metadata/phage_capsid_taxids.txt", "r") as phage_taxids_file:
    phage_taxids = [re.findall(taxid_pattern, line)[0] for line in phage_taxids_file]

print(phage_taxids)

phage_accessions = []
for taxid in set(phage_taxids):
    phage_accessions.extend(all_cl0373_accessions[taxid])

with open("metadata/phage_capsid_uniprot_accessions.txt", "w") as outfile:
    outfile.write("\n".join(phage_accessions))
    outfile.write("\n")