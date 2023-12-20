# Analyzing Metagenomic Encapsulin Hits

## Search Encapsulin Hit Sequences against Existing Encapsulin Sequences

    conda activate mmseqs2
    mkdir DBs/ #If not already present - put all DB files here
    mkdir tmp/ #If not already present
    cd DBs/

    mmseqs createdb ../all_natural_encapsulins_family_labelled.fasta target_natural_encapsulins
    mmseqs createdb ../curated_sequences.fasta query_all_encapsulin_metagenomic_hits
    mmseqs createindex target_natural_encapsulins ../tmp/
    mmseqs search query_all_encapsulin_metagenomic_hits target_natural_encapsulins natural_encapsulin_hits ../tmp -a --start-sens 4 --sens-steps 5 -s 7
    mmseqs convertalis query_all_encapsulin_metagenomic_hits target_natural_encapsulins natural_encapsulin_hits natural_encapsulin_hits.m8

## Search Encapsulin Hit Sequences against Uniref90

Requires download of UniRef90 (can be downloaded using AF2 DB download script)

    conda activate mmseqs2
    mkdir /blockstorage/sequence_DBs/tmp
    mkdir /blockstorage/sequence_DBs/encapsulin_hits
    mmseqs createdb encapsulin_metagenomics/seqs/encapsulin_hits_filtered.fasta /blockstorage/sequence_DBs/encapsulin_hits/v0_encapsulin_hits
    mmseqs createdb /blockstorage/sequence_DBs/uniref90/UniRef90.fasta /blockstorage/sequence_DBs/uniref90/UniRef90
    mmseqs createindex /blockstorage/sequence_DBs/uniref90/UniRef90 /blockstorage/sequence_DBs/tmp

    mmseqs search /blockstorage/sequence_DBs/encapsulin_hits/v0_encapsulin_hits blockstorage/sequence_DBs/uniref90/UniRef90 /blockstorage/sequence_DBs/encapsulin_hits/encapsulin_UniRef90_hits tmp -a  --start-sens 1 --sens-steps 3 -s 7 --max-accept 1 

    mmseqs convertalis /blockstorage/sequence_DBs/encapsulin_hits/v0_encapsulin_hits /blockstorage/sequence_DBs/uniref90/UniRef90 \
    /blockstorage/sequence_DBs/encapsulin_hits/encapsulin_UniRef90_hits  \
    encapsulin_UniRef90_hits.tsv \
    --format-output "query,theader,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"

Note here that `v0_encapsulin_hits` refers to the initial encapsulin hit dataset before any large scale filtering or searches using Oracle Cloud.

Analysis of these hits is in `notebooks/encapsulin_uniref90_hits.ipynb` and use of this hit data to remove viral sequences is in `notebooks/encapsulin_filtering_by_sequence_identity.ipynb`

## Search Existing Natural Encapsulins against Phage HK97 Proteins

This requires the file `family_1_2_3_natural_encapsulins.fasta` which is included in the repo. This contains all family 1, 2, and 3 encapsulins from Andreas and Giessen, 2021.

We'll also need to download a set of phage capsid protein sequences. We can use the InterPro API to get all proteins in [Pfam clan CL0373](http://www.ebi.ac.uk/interpro/set/pfam/CL0373/entry/pfam/#table) that are viral in origin. First, we query the API and get all taxids that have proteins in this clan that are viruses:

    python scripts/download_CL0373_tax_ids.py

Then, we can query the InterPro API again and get all UniProt accessions in CL0373 and filter it to only include sequences from viruses:

    python scripts/download_phage_capsid_accessions.py

And finally, we can download these sequences from UniProt directly:

    mkdir seqs/phage_capsids
    python scripts/download_phage_capsid_sequences.py
    cat seqs/phage_capsids/*.fasta > seqs/phage_capsids.fasta

Now let's build mmseqs databases and search:

    conda activate mmseqs2
    mkdir tmp
    mmseqs createdb family_1_2_3_natural_encapsulins.fasta DBs/natural_encapsulin_sequences
    mmseqs createdb seqs/phage_capsid_sequences.fasta DBs/phage_capsids
    mmseqs search DBs/natural_encapsulin_sequences DBs/phage_capsids DBs/natural_encapsulin_phage_hits tmp -a  --start-sens 1 --sens-steps 3 -s 7 --max-accept 1
    mmseqs convertalis DBs/natural_encapsulin_sequences DBs/phage_capsids DBs/natural_encapsulin_phage_hits natural_encapsulin_phage_hits.tsv --format-output "query,theader,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"

Analysis of this data is in `notebooks/encapsulin_phage_hits.ipynb` 