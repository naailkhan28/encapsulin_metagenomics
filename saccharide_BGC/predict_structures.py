import torch
import esm
import numpy as np
from Bio import SeqIO
import os

#Load ESMFold model and send it to the GPU (Change device before you start)
device = "cuda:0"
model = esm.pretrained.esmfold_v1().eval().to(device)

#Load input FASTA file
records = SeqIO.parse("saccharide_BGC/bgc_genes.fasta", "fasta")

#Initialize a list to keep track of timings
#We want to keep track of how long each sequence takes to predict

for record in records:

    name = str(record.id).split("_")[0]
    sequence = str(record.seq)
    print(f"Processing structure {name}")

    #Reduce chunk sizes to improve memory requirements
    if len(sequence) > 700:
        model.set_chunk_size(1)

    with torch.no_grad():
        output = model.infer(sequence)

    #Retrieve PAE values from the model outputs
    pae = (output["aligned_confidence_probs"][0].cpu().numpy() * np.arange(64)).mean(-1) * 31
    mask = output["atom37_atom_exists"][0,:,1] == 1
    mask = mask.cpu()
    pae = pae[mask,:][:,mask]

    #Retrieve PDB string to write from model outputs
    pdb_str = model.output_to_pdb(output)[0]


    #Write files
    with open(f"saccharide_BGC/predicted_structures/{name}.pdb", "w") as f:
        f.write(pdb_str)

    np.savetxt(f"saccharide_BGC/predicted_structures/{name}_pae.txt", pae, "%.3f")

    #Clean up
    torch.cuda.empty_cache()