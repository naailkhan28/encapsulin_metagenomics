import pandas as pd
from Bio import SeqIO

all_records = SeqIO.parse("seqs/encapsulin_seqs_for_disorder_prediction.fasta", "fasta")
sequence_dict = {str(record.id).split()[0]: record for record in all_records}

cluster_df = pd.read_csv("metadata/DALI_cluster_table_2.csv")
clusters = list(cluster_df["Cluster"].unique())

for cluster in clusters:
    mgyps = cluster_df[cluster_df["Cluster"] == cluster]["MGYP"].unique()
    seq_records = [sequence_dict[mgyp] for mgyp in mgyps]

    outfile = SeqIO.write(seq_records, f"seqs/DALI_structure_clusters/cluster_{cluster}.fasta", "fasta")