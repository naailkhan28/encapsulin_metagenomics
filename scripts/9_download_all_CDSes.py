import pandas as pd
import gzip
from time import sleep
from Bio.SeqIO import parse, write
from collections import defaultdict
from urllib import request, error
import os

print("Loading metadata table")
metadata_df = pd.read_csv("metadata/mgy_seq_metadata_with_contigs_mgya_filtered_nonans.csv")
grouped_df = metadata_df.groupby(["ERZ", "MGYA", "Contig"])

print("Creating mapping dictionary")
contigs_dict = defaultdict(list)
for name, _ in grouped_df:
    contigs_dict[f"{name[1]}-{name[0]}"].append(name[2])

unresolved_contigs = []
empty_contigs = []
for i, (name, contig_list) in enumerate(contigs_dict.items()):
    all_records = []

    if i == 0 or (i+1) % 100 == 0:
        print(f"Processing CDS set {i+1} / {len(contigs_dict)}")

    
    name = name.split("-")
    mgya = name[0]
    erz = name[1]

    url = f"https://www.ebi.ac.uk/metagenomics/api/v1/analyses/{mgya}/file/{erz}_FASTA_CDS_unannotated.faa.gz"
    filename = f"contigs/{mgya}-{erz}.faa.gz"

    try:
        request.urlretrieve(url, filename)
    except error.ContentTooShortError:
        sleep(30)
        request.urlretrieve(url, filename)
    except error.HTTPError:
        url = f"https://www.ebi.ac.uk/metagenomics/api/v1/analyses/{mgya}/file/{erz}_FASTA_predicted_cds.faa.gz"
        try:
            request.urlretrieve(url, filename)
        except error.ContentTooShortError:
            sleep(30)
            request.urlretrieve(url, filename)
        except error.HTTPError:
            print(f"Cannot resolve CDSes for entry {mgya} - {erz}")
            unresolved_contigs.append(f"{mgya}\t{erz}")
            continue
    
    with gzip.open(filename, "rt") as cdsfile:
        for record in parse(cdsfile, "fasta"):
            for contig in contig_list:
                if contig in str(record.id):
                    all_records.append(record)
    
    if not all_records:
        empty_contigs.append(f"{mgya}\t{erz}")
    
    outfile = write(all_records, f"contigs/{mgya}-{erz}_CDS.fasta", "fasta")
    os.remove(filename)


with open("metadata/unresolved_contigs.txt", "w") as unresolved_file:
    unresolved_file.write("\n".join([missing_contig for missing_contig in unresolved_contigs]))

with open("metadata/empty_contigs.txt", "w") as empty_file:
    empty_file.write("\n".join([empty_contig for empty_contig in empty_contigs]))
