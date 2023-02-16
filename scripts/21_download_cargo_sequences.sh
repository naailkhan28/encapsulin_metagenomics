#!/bin/bash

for i in {1..25}
do
	echo "Downloading file ${i}"
	wget "http://ftp.ebi.ac.uk/pub/databases/metagenomics/peptide_database/current_release/mgy_proteins_${i}.fa.gz"
	echo "Processing file ${i}"
	python scripts/filter_sequences.py "mgy_proteins_${i}.fa.gz" "seqs/cargo_proteins_${i}.fasta"
	echo "Cleaning up..."
	rm "mgy_proteins_${i}.fa.gz"
done
