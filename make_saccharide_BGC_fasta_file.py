import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#This is the dirtiest, jankiest script I possibly could have written for this - sorry if you're having to read or use this

encapsulin_sequences = {str(record.id): str(record.seq) for record in SeqIO.parse("seqs/all_encapsulin_hits.fasta", "fasta")}
cargo_sequences = {str(record.id): str(record.seq) for record in SeqIO.parse("seqs/all_putative_cargo_proteins.fasta", "fasta")}

#Load operon DataFrame
operon_df = pd.read_csv("operon_df_filtered.csv").iloc[:, :21]

#Filter operon DataFrame to only include Saccharide BGCs
family_df = pd.read_csv("encapsulin_families.csv")
saccharide_BGC_encapsulins = family_df[family_df["Cargo Description"].str.contains("Saccharide BGC")]["Encapsulin MGYP"].unique()

operon_df = operon_df.query("`Encapsulin MGYP` in @saccharide_BGC_encapsulins")
accessions = list(map(set,operon_df.values.T))

accessions = [accession for column in accessions for accession in column]
accessions = [accession for accession in accessions if accession != float("nan") and accession != "nan" and accession]

seq_records = []

for accession in accessions:
    try:
        seq_records.append(SeqRecord(
            Seq(encapsulin_sequences[accession]),
            id=accession, description=""
        ))
    except KeyError:
        try:
            seq_records.append(SeqRecord(
                Seq(cargo_sequences[accession]),
                id=accession, description=""
            ))
        except KeyError:
            continue

outfile = SeqIO.write(seq_records, "seqs/saccharide_BGC_genes.fasta", "fasta")