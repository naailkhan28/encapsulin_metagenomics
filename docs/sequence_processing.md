# Sequence Processing and Filtering

We have a list of initial encapsulin hits in the MGnify protein database - these are found in `all_non_redundant_hits.txt`. We now need to filter these to remove any phage sequences and leave only our true positive encapsulin hits.

**Note:** All shell scripts and python scripts are in `scripts/`

To do this, we do the following steps:

## Download the ERZ and MGYC data for each of our MGYP hits

    ./1_download_seq_metadata.sh

## Load all of this data into a DataFrame and write a CSV

    cd metadata/
    cat mgy_seq_metadata_* > mgy_metadata_filtered.tsv
    python 2_create_metadata_table.py
    rm mgy_metadata_filtered.tsv

## Download the contig name mapping data for all of these sequences

    python 3_get_all_mgycs.py
    ./4_download_contig_mapping.sh

## Add this contig name information to the original CSV

    python 5_add_contig_names_to_metadata_df.py

## Make a list of sequences with no known contig information

We'll call these "orphan" proteins:

    python 6_write_list_of_MGYPs_with_no_contig.py

## Get an MGYA accession for each ERZ accession

    python 7_get_MGYA_accessions_from_ERZs.py

## Merge these MGYA accessions into the original metadata table

    python 8_write_list_of_MGYPs_with_no_MGYA.py

**Note:** This script originally also wrote an output table with a list of all MGYPs that had no matching MGYA accession, however this list turned out to be empty.

## Download all protein sequences found on contigs in the metadata table

    mkdir contigs
    python 9_download_all_CDSes.py

This script will also make two text file - `empty_contigs.txt` and `unresolved_contigs.txt` - these contain the MGYA and ERZ accessions which either had no contigs matching those in the table, or whose contigs could not be downloaded from the MGnify database respectively.

## Remove all missing contigs from the metadata table and make a final table of "orphan" contigs

    cat metadata/empty_contigs.txt metadata/unresolved_contigs.txt > metadata/all_unknown_accessions.txt
    cd contigs
    find -empty > empty_contigs.txt
    mv empty_contigs.txt ../metadata
    cd ../metadata
    cat empty_contigs.txt | cut -c 3- > ../processed_empty_file_contigs.txt
    python 10_make_final_orphan_contig_table.py

MGYPs with contig: 372241
Orphan MGYPs with no contig data: 271734

## Search our contigs against phage proteomes using mmseqs2

To make the target database, you'll need the phage proteome data. The two proteomes are [UP000391682](https://www.uniprot.org/proteomes/UP000391682) and [UP000002576](https://www.uniprot.org/proteomes/UP000002576) which can be found at UniProt - download all proteins as FASTA and cat into a single file. 

    mkdir seqs
    cd contigs
    cat MGYA*.fasta > all_search_contigs.fasta
    mv all_search_contigs.fasta ../seqs
    mkdir tmp
    mkdir DBs
    cd DBs/

    conda activate mmseqs2
    mmseqs createdb ../seqs/all_search_contigs.fasta query_contigs
    mmseqs createdb ../seqs/all_phage_proteins.fasta target_phage_proteins
    mmseqs createindex target_phage_proteins ../tmp

    mmseqs search query_contigs target_phage_proteins contig_phage_hits tmp -a --start-sens 4 --sens-steps 5 -s 7
    mmseqs convertalis  query_contigs target_phage_proteins contig_phage_hits ../contig_phage_hits.m8
    

### Merge phage protein hit data into metadata table and write to a new file

    python 11_process_phage_DB_hits.py

Phage-associated MGYPs: 44340
Encapsulin MGYPs: 327901

## Explore phage hits, contig metadata, and write a final curated list of MGYPs

See notebook `phage_contigs_EDA.ipynb` - some figures + writes an outfile of 15,736 MGYP accessions

## Download curated sequences from ESM Atlas

    ./12_download_encapsulin_sequences.sh