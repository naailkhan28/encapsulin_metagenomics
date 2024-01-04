import torch
import esm
import numpy as np
from Bio import SeqIO
import os


#Check for output directory and create it if it doesn't exist
if not os.path.exists("structures/glycosyltransferases/pae"):
    os.mkdir("structures/glycosyltransferases/pae")

#Load ESMFold model and send it to the GPU (Change device before you start)
device = "cuda:0"
model = esm.pretrained.esmfold_v1().eval().to(device)

#Load input FASTA file
records = SeqIO.parse("seqs/glycosyltransferases_for_structure_prediction.fasta", "fasta")

for record in records:

    name = str(record.id)
    sequence = str(record.seq)
    print(f"Processing structure {name}")

    #Skip any sequences longer than 900
    if len(sequence) > 900:
        continue

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
    with open(f"structures/glycosyltransferases/{name}.pdb", "w") as f:
        f.write(pdb_str)

    np.savetxt(f"structures/glycosyltransferases/pae/{name}_pae.txt", pae, "%.3f")

    #Clean up
    torch.cuda.empty_cache()