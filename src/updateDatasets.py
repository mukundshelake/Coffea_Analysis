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
parser.add_argument('--selected', action='store_true', help="update the selected datasets")
parser.add_argument('-e', '--eras', nargs='+', default=["UL2016preVFP", "UL2016postVFP", "UL2017", "UL2018"], choices=["UL2016preVFP", "UL2016postVFP", "UL2017", "UL2018"], help="Specify the era(s) to process. Defaults to all eras if not specified.")

# Parse the arguments
args = parser.parse_args()

# Eras are now handled by args.eras with a default value

# Choose flag to update
## For sample from each dataset:'sample', Only TTbar MC files: "TTbar"; any other flag == Entire dataset

flag = 'all'

if args.sample:
    if args.skimmed:
        print(f"Updating the sample skimmed datasets for eras {args.eras}.")
    elif args.selected:
        print(f"Updating the sample selected datasets for eras {args.eras}.")
    else:
        print(f"Updating the sample unskimmed datasets for eras {args.eras}.")
    flag = 'sample'
else:
    if args.skimmed:
        print(f"Updating all skimmed datasets for eras {args.eras}.")
    elif args.selected:
        print(f"Updating all selected datasets for eras {args.eras}.")
    else:
        print(f"Updating all unskimmed datasets for eras {args.eras}.")

# Determine the input file prefix based on whether skimmed datasets are to be updated
file_prefix = "skimmed_filePaths_" if args.skimmed else "filePaths_"
if args.selected:
    file_prefix = "selected_filePaths_"

# Determine the output file prefix based on whether skimmed datasets are to be updated
output_prefix = "skimmed_" if args.skimmed else ""
if args.selected:
    output_prefix = "selected_"

inoutFolder = "Datasets"

for era in args.eras:
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