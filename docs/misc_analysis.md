# Miscellaneous Analysis

## Extract ESM-IF Representations

Install ESM-IF as per instructions [here](https://github.com/facebookresearch/esm/tree/main/examples/inverse_folding).

    mkdir structures/representations
    python scripts/38_get_esm_if_representations.py

## Analyzing Biome Distribution Data

Download biome metadata and create dataframe:

    ./scripts/27_download_biome_data.sh
    python scripts/28_create_biome_table.py

Explore biome data and plots in file `notebooks/biomes_EDA.ipynb`

## Large-scale structure searches against ESM Atlas

    foldseek easy-search natural_structures/*.pdb /confident_structure_DB ./esm_atlas_hits.m8 tmp --max-seqs 1000000 -c 0.5 

This required downloading all confident structures from ESM Atlas and making a foldseek DB, which was done on Oracle Cloud.

This analysis didn't give any extra hits compared to the previous strategy of searching the Atlas using the web API and searching the MGnify database using Pfam annotations.