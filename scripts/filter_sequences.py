from Bio import SeqIO
import gzip
import argparse

parser = argparse.ArgumentParser(prog = 'filter_mgnify_sequences', description = 'Filters a downloaded list of ')
parser.add_argument("input_file",  type=str, help="An input compressed fasta.gz file to filter")
parser.add_argument("output_path",  type=str, help="Filename of the output processed file to write")

args = parser.parse_args()
input_path = args.input_file
output_path = args.output_path

with open("metadata/cargo_MGYPs_updated.txt", "r") as mgyp_file:
    mgyp_list = [line.rstrip() for line in mgyp_file]

print("Creating Dictionary")
with gzip.open(input_path, "rt") as seqs_file:
    records_dict = {str(record.id).split()[0]: record for record in SeqIO.parse(seqs_file, "fasta")}

all_records = []

for i, mgyp in enumerate(mgyp_list):
    if i == 0 or (i + 1) % 5000 == 0:
        print(f"Processing sequence {i + 1} / {len(mgyp_list)}")
    
    try:
        all_records.append(records_dict[mgyp])
    except KeyError:
        continue

outfile = SeqIO.write(all_records, output_path, "fasta")
