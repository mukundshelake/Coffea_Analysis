from coffea.util import load
import os

outputDir = '../outputs'
coffeaFile = "extractorOutput.coffea"


# out = load("Output.coffea")
out = load(os.path.join(outputDir,coffeaFile))


Data = {"Ac" : out['Ac'],
        "Ain" : "../plots/Ain.png",
        "Aout" : "../plots/Aout.png"
        }