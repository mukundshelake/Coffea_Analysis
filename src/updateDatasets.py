### This script generates/updates file lists to be used for the analysis; in the form compatible with coffea. That is
### {'Dataset_label':[file1, file2, ..], ...}

### Input: The Data Path json files for the eras to be updated

from Scripts.lib.helpers import getFiles
import json
import os

eras = ["UL2016preVFP", "UL2016postVFP", "UL2017", "UL2018"]

# Choose flag to update
## For sample from each dataset:'sample', Only TTbar MC files: "TTbar"; any other falg == Entire dataset

flag = 'all'

#### If you want to update specific eras
# eras = ["UL2016preVFP"]

inoutFolder = "Datasets"


for era in eras:
    print(f"Updating the filepaths for the {era}")
    pathFile = f"filePaths_{era}.json"
    with open(os.path.join(inoutFolder,pathFile), 'r') as json_file:
        DatasetPaths = json.load(json_file) 
    DataFiles = {}
    for DataMC in DatasetPaths:
        DataFiles[DataMC] = {}
        for process in DatasetPaths[DataMC]:
            # print(process)
            DataFiles[DataMC][process] = getFiles(DatasetPaths[DataMC][process], 4, flag)  
    json_file_path = f"dataFiles_{era}.json"
    if flag=='sample':
        json_file_path = f'sampleFiles_{era}.json'
    with open(os.path.join(inoutFolder,json_file_path), 'w') as json_file:
        json.dump(DataFiles, json_file, indent=4) 


# for era in DatasetPaths:
#     DataFiles[era] = {}
#     for DataMC in DatasetPaths[era]:
#         DataFiles[era][DataMC] = {}
#         for process in DatasetPaths[era][DataMC]:
#             # print(process)
#             DataFiles[era][DataMC][process] = getFiles(DatasetPaths[era][DataMC][process])

# json_file_path = f"dataFiles_{era}.json"

# with open(json_file_path, 'w') as json_file:
#     json.dump(DataFiles[era], json_file, indent=4) 