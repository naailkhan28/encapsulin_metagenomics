import torch
import esm
import esm.inverse_folding.util as util
import os
import numpy as np

#Set up cuda device
device = "cuda:0"

#Initialize model
print("Loading model")
model, alphabet = esm.pretrained.esm_if1_gvp4_t16_142M_UR50()
model.eval()

pdbs = [file for file in os.listdir("../encapsulin_generation/ESM-IF/generated_outputs/discriminator_structures/") if file.endswith(".pdb")]
print("Extracting representations")
for i, pdb in enumerate(pdbs):
    with torch.no_grad():
        model = model.to(device)
        if i == 0 or (i+1) % 100 == 0:
            print(f"Processing structure {i+1} / {len(pdbs)}")

        #Import protein structure
        seed_model = f"../encapsulin_generation/ESM-IF/generated_outputs/discriminator_structures/{pdb}"
        seed_chain = "A"

        structure = util.load_structure(seed_model, seed_chain)
        coords, seq = util.extract_coords_from_structure(structure)

        print(len(seq))

        #Exract representations and calculate mean to get vector
        representation = util.get_encoder_output(model, alphabet, coords, device=device)
        representation = torch.mean(representation, dim=1).to("cpu").detach().numpy()
        
        #Save representation to file
        filename = pdb.split("_")[0]
        outfile = np.savetxt(f"structures/representations/{filename}.txt", representation)