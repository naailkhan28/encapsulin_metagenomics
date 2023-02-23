#!/bin/bash

FILES="seqs/DALI_structure_clusters/*.fasta"
for f in $FILES
do
	echo "Processing file ${f}"
    mmseqs easy-cluster ${f} ${f}_clustered tmp --min-seq-id 0.8
done