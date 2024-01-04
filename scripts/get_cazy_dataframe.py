import pandas as pd

pattern = r"(\w\w\w\w\[[\w,]+\])\s([A-Za-z1-9<>()\-_Α-Ωα-ω\s]+\s+)?(\d\.\d\d)?"
urls = [f"http://www.cazy.org/GT{i}_structure.html" for i in range(118)]

dfs = []

for i, url in enumerate(urls):
    output = pd.read_html(url)
    
    enzyme_table = output[1]
    try:
        enzyme_table.rename(columns=enzyme_table.iloc[2], inplace = True)
        enzyme_table = enzyme_table.iloc[3:-1, :-2].reset_index(drop=True)
    except IndexError:
        continue

    enzyme_table = enzyme_table[enzyme_table["Protein Name"] != "Eukaryota"]
    enzyme_table = enzyme_table[enzyme_table["Protein Name"] != "Archaea"]
    enzyme_table["Family"] = i

    regex_matches = enzyme_table["PDB/3D Carbohydrate Ligands  Resolution (Å)"].str.extractall(pattern).rename(
    columns={0: "PDB Code", 1: "Ligand", 2: "Resolution"})

    regex_matches.index.names=["index", "Match"]
    enzyme_table.index.name="index"

    df = regex_matches.merge(enzyme_table, left_index=True, right_index=True).reset_index(drop=True)
    dfs.append(df)

all_enzymes = pd.concat(dfs)
all_enzymes = all_enzymes.replace(0, "Unclassified").iloc[:, [0, 1, 2, 3, 4, 5, 6, 7, 9]]

outfile = all_enzymes.to_csv("metadata/cazy_structures.csv", index=False)