import os
from Bio import SeqIO
from glob import glob
import shutil

base_path = "seqs/DALI_structure_clusters"
cluster_sequence_files = [f"{base_path}/{file}" for file in os.listdir(base_path) if file.endswith("clustered_rep_seq.fasta")]

os.mkdir("structures/clusters")

for file in cluster_sequence_files:
    seq_records = SeqIO.parse(file, "fasta")

    mgyps = [str(record.id).split()[0] for record in seq_records]

    cluster = file.split("cluster_")[1].split(".")[0]
    os.mkdir(f"structures/clusters/{cluster}")

    for mgyp in mgyps:
        try:
            path = glob(f"structures/confident/{mgyp}*.pdb")[0]
        except IndexError:
            continue


        shutil.copy2(path, f"structures/clusters/{cluster}/{mgyp}.pdb")