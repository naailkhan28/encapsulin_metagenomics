# Metagenomic Encapsulin Predicted Structures

## Predict Structures

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

## Cluster Structures

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

## Compare AF2 and ESMFold Predictions

Write a list of cluster representatives:

    mkdir seqs/AF2_vs_ESMfold
    python scripts/write_structure_cluster_centres.py

See `notebooks/AlphaFold2_vs_ESMFold_comparison.ipynb` for more information about these cluster "centre" sequences.

Predict the structure of these sequences using AF2:

    mkdir AF2_predictions
    cd ~/alphafold
    sh /blockstorage/encapsulin_metagenomics/scripts/af2_predictions.sh