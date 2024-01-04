# Analyzing Saccharide BGC Glycosyltransferases
## Background
### Saccharide BGCs
In the manuscript we presented a dataset of 29 predicted Saccharide BGCs - these are predicted BGCs containing a putative encapsulin as well as predicted carbohydrate-related enzymes. All of these BGCs contained at least one glycosyltransferase enzyme (but sometimes multiple)

Doing some more research and background reading into glycosyltransferase enzymes was fascinating - it turns out that bacteria encode lots of them for various different reasons, but also bacteriophages also have their own too! This is an interesting but also slightly nefarious sign - what if some of the 29 Saccharide BGCs contain phage-associated glycosyltransferases?

In a recent addition to the notebook file `notebooks/cargo_BLAST_nr_hits.ipynb`, I explored this possibility and showed that quite a few of our putative encapsulin hits are associated with cargo proteins which show >90% identity to proteins from viruses! And some of these are our Saccharide BGC-associated encapsulins.

So, my idea now is to look deeper into the Saccharide BGC-associated glycosyltransferase proteins and try to analyze them, understand what they do and where they come from (and indeed, if they come from phages or not).

### CAZY

In order to try to better understand these glycosyltransferases in our BGCs, we'll look to [CAZY](http://www.cazy.org/), the Carbohydrate-Active enZYmes database. CAZY covers all enzymes related to carbohydrates but one of their main focuses is glycosyltransferases which is what we're interested in here. CAZY has over [100 different glycosyltransferase families](http://www.cazy.org/GlycosylTransferases). Each of these families contains a set of representative protein structures, along with details on substrate and ligands and species of origin.

The basic premise is to download all of these representative structures, and also download or predict structures for our Saccharide BGC-associated glycosyltransferase enzymes. We can then use Foldseek to search these BGC enzymes against CAZY and try to discover information about each one based on its nearest hit.

## Searching Saccharide BGC Glycosyltransferases against CAZY representatives

Download a table of all CAZZY representative glycosyltransferase structures from each family:

    python scripts/get_cazy_dataframe.py

Download all of these representative structures from the PDB:

    python scripts/download_cazy_structures.py

This script misses a few structures, so download these missing CAZY structures:

    python scripts/download_missing_cazy_structures.py

Make a list of all Saccharide BGC glycosyltransferases (MGYP accessions):

    python scripts/make_saccharide_glycosyltransferase_list.py

A notebook containing code snippets and more background on the above scripts can be found at `notebooks/glycosyltransferase_classification.ipynb`

Download all available Saccharide BGC glycosyltransferase predicted structures from ESM Atlas:

    python scripts/download_saccharide_BGC_structures.py

Predict missing ESM Atlas structures:

    python scripts/esmfold_predict_glycosyltransferases.py

Foldseek search:

    conda activate foldseek
    foldseek easy-search structures/glycosyltransferases/ structures/CAZY/ metadata/cazy_foldseek_search.tsv tmp \
    --alignment-type 1 \
    --format-mode 4 \
    --format-output "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,prob"

NOTE: I'm not 100% sure but I'm fairly certain that the `bits` field here corresponds to the query-normalized TM-Score * 100. The explanation [here](https://github.com/steineggerlab/foldseek?tab=readme-ov-file#alignment-mode) is a bit cryptic but the way I interpret this is that since we're using `--alignment-type 1` then everything is scored using TM-Score instead of 3Di. I assume that the `"score"` header referenced in the above linked GitHub repo is actually the `bits` heading we see here.

Exploration of the Foldseek hits can be found in `notebooks/saccharide_BGC_CAZY_foldseek_hits_EDA.ipynb`