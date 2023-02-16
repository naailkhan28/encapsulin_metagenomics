import requests
import json
import os
import numpy as np

mgyp_list = []

for file in os.listdir("structures/"):
    if len(file) == 20:
        mgyp_list.append(file[:-4])

for i, mgyp in enumerate(mgyp_list):
    if i == 0 or (i + 1) % 10 == 0:
        print(f"Downloading structure {i+1}/{len(mgyp_list)}")
        request = requests.get(f"https://api.esmatlas.com/fetchConfidencePrediction/{mgyp}")

        request = json.loads(request.content)
        pae = np.array(request["pae"])
        ptm = request["ptm"]

        outfile = np.savetxt(f"structures/pae/{mgyp}_ptm{ptm}.txt", pae)