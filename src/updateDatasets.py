### This script generates/updates file lists to be used for the analysis; in the form compatible with coffea. That is
### {'Dataset_label':[file1, file2, ..], ...}

### Input: The Data Path json files for the eras to be updated
import argparse
from Scripts.lib.helpers import getFiles
import json
import os

# Create the parser
parser = argparse.ArgumentParser(description="Update the datasets")

# Add arguments to the parser
parser.add_argument('-s', '--sample', action='store_true', help="update only the sample datasets")
parser.add_argument('--skimmed', action='store_true', help="update the skimmed datasets")

# Parse the arguments
args = parser.parse_args()

eras = ["UL2016preVFP", "UL2016postVFP", "UL2017", "UL2018"]

# Choose flag to update
## For sample from each dataset:'sample', Only TTbar MC files: "TTbar"; any other flag == Entire dataset

flag = 'all'

if args.sample:
    if args.skimmed:
        print(f"Updating the sample skimmed datasets for eras {eras}.")
    else:
        print(f"Updating the sample unskimmed datasets for eras {eras}.")
    flag = 'sample'
else:
    if args.skimmed:
        print(f"Updating all skimmed datasets for eras {eras}.")
    else:
        print(f"Updating all unskimmed datasets for eras {eras}.")

# Determine the input file prefix based on whether skimmed datasets are to be updated
file_prefix = "skimmed_filePaths_" if args.skimmed else "filePaths_"

# Determine the output file prefix based on whether skimmed datasets are to be updated
output_prefix = "skimmed_" if args.skimmed else ""

inoutFolder = "Datasets"

for era in eras:
    print(f"Updating the filepaths for the {era}")
    pathFile = f"{file_prefix}{era}.json"
    print(pathFile)
    with open(os.path.join(inoutFolder, pathFile), 'r') as json_file:
        DatasetPaths = json.load(json_file) 
    DataFiles = {}
    for DataMC in DatasetPaths:
        DataFiles[DataMC] = {}
        for process in DatasetPaths[DataMC]:
            print(process)
            DataFiles[DataMC][process] = getFiles(DatasetPaths[DataMC][process], flag)  
    json_file_path = f"{output_prefix}dataFiles_{era}.json"
    if flag == 'sample':
        json_file_path = f"{output_prefix}sampleFiles_{era}.json"
    with open(os.path.join(inoutFolder, json_file_path), 'w') as json_file:
        json.dump(DataFiles, json_file, indent=4)