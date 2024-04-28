from coffea.util import load, save
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import os

outputDir = 'outputs'
coffeaFile = "Fq_results.coffea"


# out = load("Output.coffea")
out = load(os.path.join(outputDir,coffeaFile))

plotDictionary = {
    "Fu Distribution": {'xAxis': 'y0', 'yAxis' : 'Fu', 'figName' : 'Fu'},
    "Fd Distribution": {'xAxis': 'y0', 'yAxis' : 'Fd', 'figName' : 'Fd'},
    "Fu_ud Distribution" : {'xAxis' : 'y0', 'yAxis' : 'Fu_ud', 'figName' : 'Fu_ud'},
    "Fd_ud Distribution" : {'xAxis' : 'y0', 'yAxis' : 'Fd_ud', 'figName' : 'Fd_ud'},
    "Du Distribution" : {'xAxis' : 'y0', 'yAxis' : 'Du', 'figName' : 'Du'},
    "Dd Distribution" : {'xAxis' : 'y0', 'yAxis' : 'Dd', 'figName' : 'Dd'},
    "Fu vs Fd Distribution" : {'xAxis' : 'Fu', 'yAxis' : 'Fd', 'figName' : 'Fu_vs_Fd'},
    "Fu_ud vs Fd_ud Distribution" : {'xAxis' : 'Fu_ud', 'yAxis' : 'Fd_ud', 'figName' : 'Fu_ud_vs_Fd_ud'},
    "Du vs Dd Distribution" : {'xAxis' : 'Du', 'yAxis' : 'Dd' , 'figName' : 'Du_vs_Dd'},
    "FuDu Distribution" : {'xAxis' : 'y0', 'yAxis' : 'FuDu', 'figName' : 'FuDu'},
    "FdDd Distribution" : {'xAxis' : 'y0', 'yAxis' : 'FdDd', 'figName' : 'FdDd'},
    "FuDu vs FdDd Distribution" : {'xAxis' : 'FuDu', 'yAxis' : 'FdDd', 'figName' : 'FuDu_vs_FdDd'},
}

regions = ['A', 'B', 'C', 'D']

# Create the main "figures" directory if it doesn't exist
figures_folder = os.path.join("plots", "paramDistributions")
os.makedirs(figures_folder, exist_ok=True)

for era in out:
    era_folder = os.path.join(figures_folder, era)
    os.makedirs(era_folder, exist_ok=True)  # Create era sub-directory
    results = out[era]['Y0_variations']
    for region in regions:
        region_folder = os.path.join(era_folder, f"region_{region}")
        os.makedirs(region_folder, exist_ok=True)  # Create region sub-sub-directory
        
        for plot in plotDictionary:
            xAxis = plotDictionary[plot]['xAxis']
            yAxis = plotDictionary[plot]['yAxis']
            x_array = np.zeros(24)
            y_array = np.zeros(24)
            for index, (key, value) in enumerate(results.items()):
                x_array[index] = key
                if xAxis != 'y0':
                    x_array[index] = value[region][xAxis]
                y_array[index] = value[region][yAxis]
            plt.figure()
            plt.title(f"{plot} ({region}) in {era}")
            plt.plot(x_array, y_array)
            plt.xlabel(xAxis)
            plt.ylabel(yAxis)
            
            # Save the figure in the appropriate sub-folder
            plt.savefig(os.path.join(region_folder, f"{plotDictionary[plot]['figName']}.png"))
            plt.close()

print("Figures saved in separate folders based on era and regions.")