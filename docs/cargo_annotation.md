# Annotating Cargo Proteins for Encapsulin Hits

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