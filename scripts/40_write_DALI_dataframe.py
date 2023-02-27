#This script contains functions and processing steps explained in more detail in the notebook ../notebooks/DALI_data_analysis.ipynb
#I know this is ugly spaghetti code but it does the job alright :)

import pandas as pd
import numpy as np
from scipy.cluster import hierarchy
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import PDB
from glob import glob
import warnings
from collections import defaultdict
import matplotlib as mpl
from matplotlib.pyplot import cm

#We need to skip the first column and the first row of the input file since this is where the labels are
print("Loading DALI matrix and labels")
columns = tuple(range(1, 952))
similarity_matrix = np.loadtxt("DALI/ordered", usecols=columns, skiprows=1)

#Let's also load the labels into a list for reference
with open("DALI/ordered", "r") as matrixfile:
    matrixfile.readline()
    identifiers = [line.rstrip().split()[0] for line in matrixfile]

#We cluster a submatrix of the full similarity matrix since the first 130 columns/rows contain structures highly dissimilar to the others
similar_submatrix = similarity_matrix[130:, 130:]

#Carry out hierarchical clustering with complete linkage
print("Clustering data")
Z = hierarchy.linkage(similar_submatrix, "complete")
cmap = cm.tab20b(np.linspace(0, 1, 16))
hierarchy.set_link_color_palette([mpl.colors.rgb2hex(rgb[:3]) for rgb in cmap])
tree = hierarchy.dendrogram(Z, p=5, truncate_mode="none", color_threshold=150)


#Get a dictionary mapping each structure to a cluster using the dendrogram
print("Loading clustered data into DataFrame")
cluster_dict = defaultdict(list)
for cluster, item in zip(tree["leaves_color_list"], tree["leaves"]):
    #Each cluster will be named after its colour
    cluster_dict[cluster].append(item)

#Load our file mapping each DALI four-letter identifier to a PDB filename
with open("DALI/chain_id_mapping.txt", "r") as mapping_file:
    chain_mapping_dict = {line.split()[0]: line.split()[1].split("_ptm0")[0].replace(".pdb", "") for line in mapping_file}

def get_sequence_name_from_matrix_index(index):
    return(chain_mapping_dict[identifiers[index]])

cluster_records = []

for i, (cluster, indices) in enumerate(cluster_dict.items()):
    for index in indices:
        #We add 130 to the matrix index since the submatrix we've used here is missing the first 130 records
        cluster_records.append({"Cluster": str(i), "MGYP": get_sequence_name_from_matrix_index(index + 130)})

#Let's also add one more cluster containing the dissimilar structures not included in the dendrogram:
cluster_records.extend([{"Cluster": "Dissimilar", "MGYP": get_sequence_name_from_matrix_index(x)} for x in range(130)])
cluster_df = pd.DataFrame(cluster_records)

#Load all of our encapsulin hit sequences into a  dictionary for access
print("Getting length, mW, and pI data for sequences")
seq_record_dict = {str(record.id).split()[0]: str(record.seq) for record in SeqIO.parse("seqs/encapsulin_hits_filtered.fasta", "fasta")}

#Add the experimentally solved structures to this dictionary
for record in SeqIO.parse("seqs/experimental.fasta", "fasta"):
    seq_record_dict[str(record.id)] = str(record.seq)

def get_sequence_length(mgyp):
    return(len(seq_record_dict[mgyp]))

def get_sequence_pI(mgyp):
    return(ProteinAnalysis(seq_record_dict[mgyp].replace("X", "")).isoelectric_point())

def get_sequence_mW(mgyp):
    return(ProteinAnalysis(seq_record_dict[mgyp].replace("X", "I")).molecular_weight() / 1000)

cluster_df["Length"] = cluster_df["MGYP"].apply(get_sequence_length)
cluster_df["mW"] = cluster_df["MGYP"].apply(get_sequence_mW)
cluster_df["pI"] = cluster_df["MGYP"].apply(get_sequence_pI)

#We need to make a dictionary mapping MGYP to our original DALI four-letter identifiers
#This is the inverse of the chain_mapping_dict we made earlier from our text file
print("Adding similarity values for experimentally solved encapsulin structures")
identifier_mapping_dict = {value: key for key, value in chain_mapping_dict.items()}

def get_experimental_similarity(mgyp, experimental_id):
    row = identifiers.index(identifier_mapping_dict[experimental_id])
    column = identifiers.index(identifier_mapping_dict[mgyp])

    return(similarity_matrix[row, column] / similarity_matrix[row, row])

experimental_ids = ["T_maritima_T1", "M_xanthus_T3", "S_elongatus_T1", "Q_thermotolerans_T4"]

for idx in experimental_ids:
    cluster_df[idx] = cluster_df["MGYP"].apply(get_experimental_similarity, args=(idx, ))

cluster_df["Closest Match"] = cluster_df.loc[:, experimental_ids].idxmax(axis="columns")

print("Getting pLDDT Disorder Data")
def get_longest_unconfident_stretch(mgyp):
    if len(mgyp) != 16:
        return float("NaN")

    path = glob(f"structures/confident/{mgyp}*.pdb")[0]

    #Catch any warnings to avoid cluttering the output
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        structure = PDB.PDBParser().get_structure("Input", path)

    pLDDTs = [residue["CA"].get_bfactor() for residue in structure[0]["A"]]
    pLDDTs = [pLDDT / 100 for pLDDT in pLDDTs if pLDDT > 1]

    unconfident_stretch = []
    lengths = []

    for pLDDT in pLDDTs:
        if pLDDT < 0.7:
            unconfident_stretch.append(pLDDT)

        else:
            if len(unconfident_stretch) > 0:
                lengths.append(len(unconfident_stretch))
                unconfident_stretch = []

    if len(unconfident_stretch) > 0:
                lengths.append(len(unconfident_stretch))
                unconfident_stretch = []
    
    return(max(lengths) if lengths else 0)

cluster_df["Longest Disordered Region (pLDDT)"] = cluster_df["MGYP"].apply(get_longest_unconfident_stretch)
outfile = cluster_df.to_csv("metadata/DALI_cluster_table_2.csv", index=False)