#!/bin/bash
wget "http://ftp.ebi.ac.uk/pub/databases/metagenomics/peptide_database/current_release/mgy_biomes.tsv.gz";
zgrep -F -f final_filtered_mgyp_list.txt mgy_biomes.tsv.gz > metadata/encapsulin_hit_biomes_filtered.tsv;
rm mgy_biomes.tsv.gz;