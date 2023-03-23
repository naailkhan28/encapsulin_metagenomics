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

## Analyze Operon Types and Cargo Proteins

Initial exploration is in `encapsulin_hits_contig_EDA.ipynb`

### Make Cargo Metadata Table

    python 13_get_cargo_mgycs.py
    ./14/download_cargo_metadata.sh

    cat metadata/cargo_seq_metadata_* > metadata/cargo_metadata_filtered.tsv
    python 15_make_cargo_metadata_table.py
    rm metadata/cargo_metadata_filtered.tsv

### Get Cargo Pfam Annotations

    python 16_make_cargo_mgyp_list.py
    wget http://ftp.ebi.ac.uk/pub/databases/metagenomics/peptide_database/current_release/mgy_ProtENN_pfam.tsv.gz
    setenv LC_ALL C
    zgrep -F -f metadata/cargo_MGYPs.txt mgy_ProtENN_pfam.tsv.gz > pfams/cargo_pfams.tsv
    rm mgy_ProtENN_pfam.tsv.gz

### Make Operon Table

This table has a row for each encapsulin MGYP containing potential cargo MGYPs, and their position either side of the encapsulin (from -10 to +10)

    python 17_process_cargos_make_operon_df.py

### Re-filter operon table, sequences, and MGYP list to remove final phage-associated sequences

    python 18_refilter_sequences_write_new_files.py

### Download all cargo sequences

    python 19_make_filtered_cargo_list.py
    ./21_download_cargo_sequences.sh

### Search cargo sequences for targeting peptide sequences

You'll need a file `cargo_targeting_peptides.fasta` and the substitution matrix `PAM30.out` which you can get from [here](https://github.com/soedinglab/MMseqs2/blob/master/data/PAM30.out)

    mkdir tmp
    cd DBs/

    conda activate mmseqs2
    mmseqs createdb ../seqs/cargo_targeting_peptides.fasta query_targeting_peptides
    mmseqs createdb ../seqs/all_putative_cargo_proteins.fasta target_cargo_proteins
    mmseqs createindex target_cargo_proteins ../tmp --spaced-kmer-mode 0

    mmseqs search query_targeting_peptides target_cargo_proteins cargo_hits tmp -a \
    -s 7 --spaced-kmer-mode 0 --max-seqs 2000 -e 200000 --sub-mat ../PAM30.out

    mmseqs convertalis  query_targeting_peptides target_cargo_proteins cargo_hits ../cargo_hits.m8

### Search cargo sequences with Family 2 HMMs

There are 4 types of family 2 cargo protein:

 * Cysteine Desulfurase

 * Xylulose Kinase

 * Terpene Cyclase

 * Polyprenyl Transferase

I've previously downloaded all the sequences for these types as listed in the Excel spreadsheet for the supplementary data of Andreas and Giessen, 2021 from UniProt. Make alignments using ClustalO:

    clustalo -i cysteine desulfurase.fasta -o alignments/cysteine_desulfurase.fasta
    #Repeat for other 3 types

Install hmmer:

    conda create -n hmmer python=3.10
    conda activate hmmer
    conda install -c bioconda hmmer

Build HMMs:

    mkdir hmms/

    hmmbuild HMMs/cysteine_desulfurase.hmm alignments/cysteine_desulfurase.fasta
    hmmbuild HMMs/xylulose_kinase.hmm alignments/xylulose_kinase.fasta
    hmmbuild HMMs/polyprenyl_transferase.hmm alignments/polyprenyl_transferase.fasta
    hmmbuild HMMs/terpene_cyclase.hmm alignments/terpene_cyclase.fasta

Search HMMs against cargo proteins:

    hmmsearch HMMs/cysteine_desulfurase.hmm seqs/all_putative_cargo_proteins.fasta > HMMs/cysteine_desulfurase.out
    hmmsearch HMMs/xylulose_kinase.hmm seqs/all_putative_cargo_proteins.fasta > HMMs/xylulose_kinase.out
    hmmsearch HMMs/polyprenyl_transferase.hmm seqs/all_putative_cargo_proteins.fasta > HMMs/polyprenyl_transferase.out
    hmmsearch HMMs/terpene_cyclase.hmm seqs/all_putative_cargo_proteins.fasta > HMMs/terpene_cyclase.out

### Search cargo sequences against BLAST non-redundant database

(Requires BLAST nr database download - FASTA can be downloaded from their FTP server)

    conda activate mmseqs2
    mmseqs createdb encapsulin_metagenomics/seqs/filtered_cargo_proteins.fasta /blockstorage/sequence_DBs/encapsulin_hits/v0_all_cargo_proteins
    mmseqs createdb /blockstorage/sequence_DBs/blast_nr/nr.fasta /blockstorage/sequence_DBs/blast_nr/blast_nr
    mmseqs createindex /blockstorage/sequence_DBs/blast_nr/blast_nr /blockstorage/sequence_DBs/tmp

    mmseqs search /blockstorage/sequence_DBs/encapsulin_hits/v0_all_cargo_proteins blockstorage/sequence_DBs/blast_nr/blast_nr /blockstorage/sequence_DBs/encapsulin_hits/cargo_blast_nr_hits tmp -a  --start-sens 1 --sens-steps 3 -s 7 --max-accept 30 
    mmseqs convertalis /blockstorage/sequence_DBs/encapsulin_hits/v0_all_cargo_proteins /blockstorage/sequence_DBs/blast_nr/blast_nr \
    /blockstorage/sequence_DBs/encapsulin_hits/cargo_blast_nr_hits \
    cargo_blast_nr_hits.tsv \
    --format-output "query,theader,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"

### Run disorder predictions using ADOPT

Install ADOPT as per instructions


    mkdir disorder
    mkdir disorder/reprs
    mkdir disorder/seqs
    python scripts/split_sequences.py

    conda activate adopt
    cd ../ADOPT
    ./scripts/adopt_inference_loop.sh

### Predict Biosynthetic Gene Clusters (BGCs)

Download and install [antiSMASH](https://academic.oup.com/nar/article/49/W1/W29/6274535?login=true):

    conda create -n antismash python=3.8
    conda activate antismash
    conda install pip numpy hmmer2 hmmer diamond fasttree prodigal blast muscle=3.8.1551 glimmerhmm meme java-jdk

    cd ..
    mkdir antismash
    cd antismash
    wget https://dl.secondarymetabolites.org/releases/6.1.1/antismash-6.1.1.tar.gz
    tar -zxf antismash-6.1.1.tar.gz
    rm antismash-6.1.1.tar.gz

    ../anaconda3/envs/antismash/bin/pip install ./antismash-6.1.1
    download-antismash-databases
    antismash --check-prereqs

Download all contig nucleotide sequences for BGC prediction with antiSMASH:

    python scripts/29_download_all_contigs_nucleotide.py

Run antiSMASH on all contigs:

    mkdir antiSMASH
    cd antiSMASH
    find ../contigs/ -maxdepth 1 -size 0 -delete
    sh ../scripts/30_run_antismash_predictions.sh

Process antiSMASH HTML files to make an output DataFrame

    python 31_create_antismash_table.py

Explore antiSMASH results in notebook `antiSMASH_EDA.ipynb`

Download and install [DeepBGC](https://github.com/Merck/deepbgc):

    conda create -n deepbgc python=3.7 hmmer prodigal
    conda activate deepbgc
    ../anaconda3/envs/deepbgc/bin/pip install deepbgc protobuf=3.20 hmmlearn

Run DeepBGC predictions for all contigs:

    mkdir deepbgc
    cd deepbgc
    sh ../scripts/32_run_deepbgc_predictions.sh

Compile all hits into a DataFrame:

    python 37_create_deepBGC_table.py

Explore antiSMASH results in notebook `deepBGC_EDA.ipynb`

## Analyze ESMFold Predicted Structures

### Download All Available Structures From ESM Atlas

    mkdir structures
    python scripts/22_download_all_structures.py

### Download PAE data for these structures

    mkdir structures/pae
    python scripts/24_download_structure_PAEs.py

### Predict all remaining structures (this will probably be most of them)

You'll need to use Google Colab for this - there's a notebook [here](https://colab.research.google.com/drive/1LGX6xflEwnUDDDOq2DCOVUM7KZfgiXqx#scrollTo=7nEDhqwESdrq).

Upload the file `seqs/encapsulin_hits_for_structure_prediction.fasta` and run all the cells. This will predict all sequences under 900aa in length. 

If you're using the free tier of Colab then you'll run into usage limits and the scripts will eventually time out - to avoid this, I recommend interrupting the structure prediction and downloading the results as a zip file every 1-2 hours. The notebook has code to write these zip files and download them. Put all the resulting PDBs and PAE files in the respective folders inside this repository, and then run `23_update_structure_list.py`. This will write a new, updated version of `seqs/encapsulin_hits_for_structure_prediction.fasta` with the already predicted structures removed. You can then wait a day or so and then re-upload and predict in Google Colab as before.

Rinse and repeat until all structure data is present!

### Filter out any structures with average pLDDT <= 0.7:

    python scripts/33_filter_predicted_structures.py
    mkdir structures/confident
    cd structures
    sh ../scripts_34_move_confident_structure_predictions.sh

### Install DALI

TODO

### Import all predicted structures into DALI

For this we need to generate unique four-character identifiers for each input structure - the Python script will handle this and write an output file `DALI/chain_id_mapping.txt` containing a tab-separated list of each ID and its corresponding predicted structure filename.

We also download and import the 4x experimental encapsulin structures that we used for the structure search against ESM Atlas (these use their PDB IDs as identifiers).

    mkdir DALI
    mkdir DALI/data
    python scripts/35_dali_import_structures.py
    sh scripts/36_dali_import_natural_encapsulin_structures.sh

You'll also need to add the following to the end of the `chain_id_mapping.txt` file (copy paste into a text editor or cat using terminal):

    7MU1A   T_maritima_T1
    6X8MA   S_elongatus_T1
    7S2TA   M_xanthus_T3
    6NJ8A   Q_thermotolerans_T4

### Run DALI all-against-all comparison

    perl ../DALI/DaliLite.v5/bin/dali.pl --matrix --query DALI/chain_id_mapping.txt --dat1 DALI/data --clean

Be patient - this takes 5+ days!!

### Process DALI output data and explore

See notebook `notebooks/DALI_data_analysis.ipynb` for details and processing steps. To make output DataFrame:

    python scripts/40_write_DALI_dataframe.py

Let's also write sequence files containing the sequences of each cluster:

    mkdir seqs/DALI_structure_clusters
    python 41_make_structure_cluster_files.py

And we can also cluster these sets of sequences at 80% identity to get a reduced list of structures in each cluster:

    sh scripts/42_cluster_seqs_from_structure_clusters.sh

Lastly, let's create a set of folders and put the reduced structure sets in each one so we can download them and interpret them:

    python scripts/43_move_reduced_structure_clusters.py

Check `structures/clusters` - there'll be a folder for each cluster with a set of reduced structures in each one for download and manual inspection.

### Foldseek clustering

Install Foldseek:

    conda create -n foldseek python=3.9
    conda activate foldseek
    conda install -c conda-forge -c bioconda foldseek

Cluster all structures using Foldseek with a minumum alignment length cutoff of 0.8:

    mkdir foldseek
    cd foldseek
    mkdir DBs
    mkdir tmp
    foldseek createdb ../structures/confident DBs/confident_esmfold_predictions_foldseek_db
    foldseek search DBs/confident_esmfold_predictions_foldseek_db DBs/confident_esmfold_predictions_foldseek_db DBs/all_against_all tmp -c 0.8 
    foldseek clust DBs/confident_esmfold_predictions_foldseek_db DBs/all_against_all DBs/clusters
    foldseek createtsv DBs/confident_esmfold_predictions_foldseek_db DBs/confident_esmfold_predictions_foldseek_db DBs/clusters clusters.tsv

## Analyzing Biome Distribution Data

Download biome metadata and create dataframe:

    ./scripts/27_download_biome_data.sh
    python scripts/28_create_biome_table.py

Explore biome data and plots in file `notebooks/biomes_EDA.ipynb`

### Extract ESM-IF Representations

Install ESM-IF as per instructions [here](https://github.com/facebookresearch/esm/tree/main/examples/inverse_folding).

    mkdir structures/representations
    python scripts/38_get_esm_if_representations.py

## Phylogenetic Analysis

### Prepare Input Data

Cluster all known family 2 encapsulins at 60% identity (requires file `fixed_family2.fasta`):

    mmseqs easy-cluster fixed_family2.fasta family2_clustered_60 tmp --min-seq-id 0.6

Cluster all new unannotated encapsulins at 30% identity:

    python scripts/create_family_sequence_files.py
    cd seqs
    mkdir tmp
    mmseqs easy-cluster unannotated.fasta unannotated.fasta tmp --min-seq-id 0.3