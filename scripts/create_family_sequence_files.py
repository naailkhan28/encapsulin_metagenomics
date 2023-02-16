import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

family_df = pd.read_csv("encapsulin_families.csv")

cargo_dict = {"DyP Peroxidase": 1, "Alternative oxidase": 1, "Amine oxidoreductase": 1, "Ferritin": 1,
 "Ubiquinone biosynthesis protein COQ7": 1, "Rubrerythrin": 1,
 "Manganese containing catalase": 1, "Cysteine Desulfurase": 2,
 "Polyprenyl Transferase": 2, "Xylulose Kinase": 2, "Saccharide BGC (Redox)": 3,
 "Saccharide BGC": 3, "DeoC": 4, "OsmC": 4}

def get_family(cargo):
    return(cargo_dict[cargo])

family_df["Family"] = family_df["Cargo Description"].apply(get_family)

for family in family_df["Family"].unique():
    mgyps = family_df[family_df["Family"] == family]["Encapsulin MGYP"].unique()
    records = []

    for record in SeqIO.parse("seqs/encapsulin_hits_filtered.fasta", "fasta"):
        if str(record.id).split()[0] in mgyps:
            id = f"family_{family}_{str(record.id)}"
            sequence = record.seq
            records.append(SeqRecord(record.seq, id=id, description=""))
    
    outfile = SeqIO.write(records, f"seqs/annotated_family_{family}.fasta", "fasta")

annotated_mgyps = list(family_df["Encapsulin MGYP"].unique())
unannotated_mgyps = list(pd.read_csv("metadata/mgy_seq_metadata_with_contigs_filtered_nonans.csv").query("MGYP not in @annotated_mgyps")["MGYP"].unique())

unannotated_records = []
for record in SeqIO.parse("seqs/encapsulin_hits_filtered.fasta", "fasta"):
    if str(record.id).split()[0] in unannotated_mgyps:
        id = f"unannotated_{str(record.id)}"
        unannotated_records.append(SeqRecord(record.seq, id=id, description=""))
    
outfile = SeqIO.write(unannotated_records, f"seqs/unannotated.fasta", "fasta")