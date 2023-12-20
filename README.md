# Mining Metagenomics Data for Novel Bacterial Nanocompartments

This repository contains scripts, code, and data from our upcoming manuscript *"Mining Metagenomics Data for Novel Bacterial Nanocompartments"*.

## Data

This repo contains a lot of raw data and intermediate data files generated throughout the analysis that was done for this paper, however the following files are the most important ones:

 * `seqs/encapsulin_hits_filtered.fasta` contains the final dataset of 1548 putative encapsulin sequences in FASTA format. The FASTA headers contain each sequence's MGnify Protein Database accession, as well as whether the sequence is a full-length or partial sequence (`FL=1` or `FL=0`) and whether the sequence is a cluster representative (`CR=1` or `CR=0`). See the [MGnify Protein DB readme](http://ftp.ebi.ac.uk/pub/databases/metagenomics/peptide_database/current_release/README.txt) for more info.

 * `encapsulin_families.csv` contains putative family assignments for ≈200 of our new encapsulin sequences. This file is in comma-separated format with column headings `Encapsulin MGYP`, `Cargo Description`, and `Cargo Search Method` which describes how the cargo type was assigned.

 * `operon_df_filtered.csv` is a comma-separated file containing data on the surrounding genes for each putative encapsulin. The column headings in this file are `Encapsulin MGYP` and `[-10, -9, -8, -7... -1]` `[1, 2, 3, 4...10]` which correspond to the MGYP accessions for up to 10 genes either side of the putative encapsulin. This file also contains columns `[Pfam -10, Pfam -9, Pfam -8, Pfam -7... Pfam -1]` `[Pfam 1, Pfam 2, Pfam 3, Pfam 4...Pfam 10]` which contain the Pfam annotations for each gene surrounding the encapsulin.


## Documentation And Scripts

The `docs/` directory contains Markdown files documenting the scripts and analyses for this study. The folder contains the following files:

1. `sequence_processing.md` documents the process of filtering the initial non-redundant hits from structure and annotation search, to remove potentially phage-related proteins.
2. `cargo_annotation.md` describes how the putative encapsulin dataset was annotated with hypothetical cargo protein functions.
3. `encapsulin_sequence_searches.md` details the sequence searches with `mmseqs2` that were carried out in this study - searching the new encapsulin dataset against various databases, as well as searching the pre-existing encapsulin dataset against UniRef90, among others. 
4. `predicted_structure_analysis.md` goes through the steps to filter and cluster the ESMFold predicted structures using either DALI (slow) or Foldseek (fast).
5. `misc_analysis.md` contains details around a few other miscellaneous analyses, most of which weren't included in the manuscript. This includes extracting ESM-IF representations for the ESMFold predicted structures, and searching the full-size ESM Atlas using the experimentally solved encapsulin structures as queries.

The `scripts/` folder contains all bash and Python scripts necessary to reproduce the analyses described above - scripts are numbered in the order that they were written and can be run sequentially to reproduce this study - although note that some scripts will download large amounts of data so make sure you have enough storage space and compute time for these.

## Notebooks

The `notebooks/` directory contains Jupyter Notebooks with exploratory data analysis and plots for the sequence searches, predicted structure clustering, phage sequence filtering, BGC prediction with antiSMASH and DeepBGC, and a few others. The code in these notebooks may not run without acquiring and processing the required data, which will be detailed in the documentation outlined above.