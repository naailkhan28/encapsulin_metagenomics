#!/bin/bash

for i in {1..19}
do
	echo "Downloading file ${i}"
	wget "http://ftp.ebi.ac.uk/pub/databases/metagenomics/peptide_database/current_release/mgy_seq_metadata_${i}.tsv.gz"
	echo "Processing file ${i}"
	zgrep -F -f hits/all_non_redundant_hits.txt "mgy_seq_metadata_${i}.tsv.gz" > metadata/mgy_seq_metadata_${i}.tsv
	echo "Cleaning up..."
	rm "mgy_seq_metadata_${i}.tsv.gz"
done
