#!/bin/bash

for i in {1..9}
do
	echo "Downloading file ${i}"
	wget "http://ftp.ebi.ac.uk/pub/databases/metagenomics/peptide_database/current_release/mgy_contig_map_${i}.tsv.gz"
	echo "Processing file ${i}"
	zgrep -F -f metadata/all_MGYCs.txt "mgy_contig_map_${i}.tsv.gz" > metadata/mgy_contig_map_${i}.tsv
	echo "Cleaning up..."
	rm "mgy_contig_map_${i}.tsv.gz"
done
