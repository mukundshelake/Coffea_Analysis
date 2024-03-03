from coffea.util import load
import json
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import os


outputDir = 'outputs'
coffeaFile = "Output.coffea"


# out = load("Output.coffea")
out = load(os.path.join(outputDir,coffeaFile))

ymax = 2.5j
Fu = 0.6034
Fd = 0.3966
Du = 0.6264
Dd = 0.4098

for era in ['UL2016preVFP', 'UL2016postVFP', 'UL2017', 'UL2018']:
    print(f"Working on {era}")
    df = pd.DataFrame(columns = ['y0', 'A', 'A_', 'B', 'D', 'C', 'C_', 'inFactor', 'outFactor', 'C1', 'C2', 'C3', 'C4', 'Ain', 'Aout'])
    ytHigher = out[f"{era}"][f"{era}"]["yMatrix_higherYt"]
    ytbarHigher = out[f"{era}"][f"{era}"]["yMatrix_higherYtbar"]

    ytHigher_yt = ytHigher.project('y_t')
    ytHigher_ytbar = ytHigher.project('y_tbar')

    ytbarHigher_yt = ytbarHigher.project('y_t')
    ytbarHigher_ytbar = ytbarHigher.project('y_tbar')
    for i in range(24):
        y0 = (i+1)*0.1
        slice_y0 = y0*1j
        # print(y0, slice_y0)
        A = ytHigher[0j:slice_y0,0:slice_y0].sum()
        B = ytHigher[slice_y0:ymax,0:slice_y0].sum()
        C = ytHigher[slice_y0:ymax,slice_y0:ymax].sum()


        A_ = ytbarHigher[0j:slice_y0,0:slice_y0].sum()
        D =  ytbarHigher[0:slice_y0, slice_y0:ymax].sum()
        C_ = ytbarHigher[slice_y0:ymax,slice_y0:ymax].sum()

        inFactor = -(1/(1+2*((A+A_)/(B+D))))
        outFactor = 1/(1+2*((C+C_)/(B+D)))

        C1 = inFactor*Du*Fu
        C2 = inFactor*Dd*Fd
        C3 = outFactor*Du*Fu
        C4 = outFactor*Dd*Fd

        Ain = (D-B)/(D+B+A+A_)
        Aout = (B-D)/(D+B+C+C_)

        lst = [y0, A, A_, B, D, C, C_, inFactor, outFactor, C1, C2, C3, C4, Ain, Aout]
        df.loc[i] = lst
    # print(f"{era}")
    jsonFile = f"{era}_results.json"
    jsonPath = os.path.join(outputDir, jsonFile)
    df.to_json(jsonPath)
    with open(jsonPath, 'r') as file:
        data = json.load(file)

    # Writing it back with indentation for better readability
    with open(jsonPath, 'w') as file:
        json.dump(data, file, indent=4) 