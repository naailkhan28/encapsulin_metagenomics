import urllib.request

with open("metadata/phage_capsid_uniprot_accessions.txt", "r") as accessions_file:
    accessions = [line.rstrip() for line in accessions_file]

base_url="https://www.uniprot.org/uniprot/"

for i, accession in enumerate(accessions):
    if i == 0 or (i+1) % 100 == 0:
        print(f"Downloading accession {i+1} / {len(accessions)}")

    url = f"{base_url}{accession}.fasta"
    filename = f"seqs/phage_capsids/{accession}.fasta"
    urllib.request.urlretrieve(url, filename)