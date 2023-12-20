from Bio import SeqIO
import numpy as np
from scipy.cluster.hierarchy import linkage, dendrogram, set_link_color_palette
import matplotlib as mpl
from matplotlib.pyplot import cm
from collections import defaultdict

sequence_ids = []

#Get all sequence IDs from any clusters which have a small enough number of sequences
for i in range(17):
    infile = f"seqs/DALI_structure_clusters/cluster_{i}.fasta_clustered_rep_seq.fasta"

    if len(list(SeqIO.parse(infile, "fasta"))) < 15:
        sequence_ids.extend([str(record.id).split()[0] for record in SeqIO.parse(infile, "fasta")])


#Load and cluster the similarity matrix
columns = tuple(range(1, 952))
similarity_matrix = np.loadtxt("DALI/ordered", usecols=columns, skiprows=1)

#The first 130 columns contain structures with low similarity to all the rest so I've manually removed these into their own cluster
similar_submatrix = similarity_matrix[130:, 130:]

Z = linkage(similar_submatrix, "complete")

cmap = cm.tab20b(np.linspace(0, 1, 16))
set_link_color_palette([mpl.colors.rgb2hex(rgb[:3]) for rgb in cmap])
tree = dendrogram(Z, p=5, truncate_mode="none", color_threshold=150)

cluster_dict = defaultdict(list)
for cluster, item in zip(tree["leaves_color_list"], tree["leaves"]):
    #Each cluster will be named after its colour
    cluster_dict[cluster].append(item)

clusters_similarity = {}

for cluster_name, cluster_members in cluster_dict.items():
    clusters_similarity[cluster_name] = similar_submatrix[np.ix_(cluster_members, cluster_members)]

#Load the labels into a list for reference
with open("DALI/ordered", "r") as matrixfile:
    matrixfile.readline()
    identifiers = [line.rstrip().split()[0] for line in matrixfile]

#Load our file mapping each DALI four-letter identifier to a PDB filename
with open("DALI/chain_id_mapping.txt", "r") as mapping_file:
    chain_mapping_dict = {line.split()[0]: line.split()[1].split("_ptm0")[0].replace(".pdb", "") for line in mapping_file}

def get_sequence_name_from_matrix_index(index):
    return(chain_mapping_dict[identifiers[index]])

for cluster_name, similarities in clusters_similarity.items():
    mean_distances = np.mean(similarities, axis=1)

    centre_index = np.argmin(mean_distances)
    centre_structure_index = cluster_dict[cluster_name][centre_index]
    sequence_ids.append(get_sequence_name_from_matrix_index(centre_structure_index + 130)) #We need to add 130 to the index since we're using a subset of the full DALI similarity matrix

sequence_ids = set(sequence_ids)

sequence_records = []

for record in SeqIO.parse("seqs/all_encapsulin_hits.fasta", "fasta"):
    if str(record.id).split()[0] in sequence_ids:
        outfile = SeqIO.write(record, f"seqs/AF2_vs_ESMFold/{str(record.id).split()[0]}.fasta", "fasta")