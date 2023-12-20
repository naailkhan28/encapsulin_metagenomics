#!/usr/bin/env bash

#Make sure you run this from the AlphaFold2 directory!

for FASTA_FILE in /blockstorage/encapsulin_metagenomics/seqs/AF2_vs_ESMFold/*.fasta; do
    ./run_alphafold.sh -d /blockstorage/sequence_DBs -o /blockstorage/encapsulin_metagenomics/AF2_predictions -f $FASTA_FILE -t 2023-12-01 -a 0
done