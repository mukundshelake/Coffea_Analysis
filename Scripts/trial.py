import glob, os, json
from lib.helpers import getFiles


with open('Datasets/sampleFiles_UL2016preVFP.json', 'r') as json_file:
    fileset = json.load(json_file) 
print({**fileset['Data_el'], **fileset['MC_el']})