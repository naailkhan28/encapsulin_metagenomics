import pandas as pd
import gzip
import json
import re
from Bio import SeqIO
from glob import glob

#We can drop a bunch of unnecessary columns from the CSV file
deepbgc_df = pd.read_csv("metadata/deepbgc_hits.csv").drop(
    ["detector", "detector_version", "detector_label", "bgc_candidate_id", "antibacterial", "cytotoxic",
     "inhibitor", "antifungal","Alkaloid", "NRP", "Other", "Polyketide", "RiPP", "Saccharide", "Terpene", "Hit?", "Contig"],
axis=1)

# The below code will load the Pfam labels dictionary downloaded from the GoogleResearch/Proteinfer GitHub repo
# I've also manually added some code below to fill in some important missing labels that I manually curated



with open("DBs/label_descriptions.json.gz", 'rb') as f:
    with gzip.GzipFile(fileobj=f, mode='rb') as gzip_file:
      labels_dict = json.load(gzip_file)

labels_dict["PF19821"] = "Phage capsid protein"
labels_dict["PF19307"] = "Phage capsid-like protein"
labels_dict['PF19289'] = "PmbA/TldA metallopeptidase C-terminal domain"
labels_dict['PF19290'] = "PmbA/TldA metallopeptidase central domain"
labels_dict['PF20211'] = "Family of unknown function (DUF6571)"
labels_dict['PF19782'] = "Family of unknown function (DUF6267)"
labels_dict['PF19343'] = "Family of unknown function (DUF5923)"
labels_dict['PF20036'] = "Major capsid protein 13-like"
labels_dict['PF18960'] = "Family of unknown function (DUF5702)"
labels_dict['PF18906'] = "Phage tail tube protein"
labels_dict['PF19753'] = "Family of unknown function (DUF6240)"

def get_label(pfam):
    try:
        return(labels_dict[pfam])
    except KeyError:
        return(None)
    
pattern = r"Glycosyl\s*transferase"

#This function will access the Pfam family lists for each DataFrame entry and convert it to a free text description
def get_pfam_list_labels(pfam_list):
    pfam_list = pfam_list.split(";")

    try:
        return(" | ".join([get_label(pfam) for pfam in pfam_list]))
    except TypeError:
        return("None")

deepbgc_df["bio_pfam_ids"] = deepbgc_df["bio_pfam_ids"].fillna("None")
deepbgc_df["Descriptions"] = deepbgc_df["bio_pfam_ids"].apply(get_pfam_list_labels)

#Make a list of all glycosyltransferase proteins
#This list contains the CDS name from the contig, in the format ERZ.contigname_CDSnumber
proteins = []

saccharide_df = deepbgc_df[deepbgc_df["product_class"] == "Saccharide"].sort_values(by="MGYP")

for row in saccharide_df.to_dict(orient="records"):
    mgya = row["MGYA"]
    erz = row["ERZ"]
    mgyp  = row["MGYP"]
    protein_ids = row["protein_ids"].split(";")
    pfam_df = pd.read_csv(f"deepbgc/{mgya}-{erz}/{mgya}-{erz}.pfam.tsv", sep="\t").query("protein_id in @protein_ids").query("deepbgc_score > 0.5")

    for i, protein in enumerate(protein_ids):
        pfams = pfam_df[pfam_df["protein_id"] == protein]["pfam_id"].values
        if len(pfams) > 0:  
            pfam_labels = " | ".join([get_label(pfam) for pfam in pfams])
        
        if re.match(pattern, pfam_labels):
            proteins.append(protein)


#Load the metadata mapping MGYPs to their ERZ, start, end, and strand data
#We'll also filter this DataFrame to make it easier to process
metadata_df = pd.read_csv("metadata/cargo_seq_metadata_filtered.csv")

erzs = [protein.split(".")[0] for protein in proteins]
metadata_df = metadata_df.query("ERZ in @erzs")

#Now, we need to map each CDS name to its ERZ, start, end, and strand
#The ERZ is in the CDS name itself but the start, end, and strand data needs to be loaded from the CDS FASTA header
#We can load these FASTA headers from the FASTA files containing all CDSes from each contig
protein_names = [f"{protein.split('_')[1]}_{protein.split('_')[2]}" for protein in proteins]
erz_dict = {}

for protein_name in protein_names:
    erz = protein_name.split(".")[0]
    seq_records = SeqIO.parse(glob(f"contigs/CDS/*-{erz}_CDS.fasta")[0], "fasta")

    fasta_header = [str(record.description) for record in seq_records if str(record.description).split()[0] == protein_name][0]

    start = fasta_header.split(" # ")[1]
    end = fasta_header.split(" # ")[2]
    strand = fasta_header.split(" # ")[3]

    erz_dict["_".join((erz, start, end, strand))] = protein_name


#And now FINALLY we can make a dictionary mapping the CDS names to MGYPs
mgyp_dict = {}

for _, row in metadata_df.iterrows():
    try:
        mgyp_dict[erz_dict["_".join((row["ERZ"], str(row["Start"]), str(row["End"]), str(row["Strand"])))]] = row["MGYP"]
    except KeyError:
        continue


with open("metadata/saccharide_BGC_glycosyltransferases.txt", "w") as outfile:
    outfile.write("\n".join(mgyp for mgyp in set(mgyp_dict.values())))