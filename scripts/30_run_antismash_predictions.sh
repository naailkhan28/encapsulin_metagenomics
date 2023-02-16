#!/bin/bash

FILES="../contigs/*.fasta"
for f in $FILES
do
	echo "Processing file ${f}"
	antismash --cb-general --cb-knownclusters --cb-subclusters --asf --pfam2go --clusterhmmer --pfam2go --genefinding-tool prodigal ${f}
done