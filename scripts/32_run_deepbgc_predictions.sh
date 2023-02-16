#!/bin/bash

FILES="../contigs/*.fasta"
for f in $FILES
do
	echo "Processing file ${f}"
	deepbgc pipeline --prodigal-meta-mode --detector deepbgc --classifier product_class --classifier product_activity ${f}
done