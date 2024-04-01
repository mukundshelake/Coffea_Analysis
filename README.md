# Coffea_Analysis
Top asymmetry analysis using the Coffea framework 


## To-do
* Remove the sensitive info from the repo and public the repo
* Add snakemake
* Add a twiki


## Project structure
```
Coffea_Analysis/
│
├── src/
│   ├── __init__.py
│   ├── getSFs.py
│   ├── updateDatasets.py
│   │ 
│   ├── Datasets/
│   │   ├── filePaths_{era}.json files for eras = ['UL2016preVFP', 'UL2016postVFP', 'UL2017', 'UL2018']
│   │   ├── dataFiles_{era}.json files for eras = ['UL2016preVFP', 'UL2016postVFP', 'UL2017', 'UL2018']
│   │   └── sampleFiles_{era}.json files for eras = ['UL2016preVFP', 'UL2016postVFP', 'UL2017', 'UL2018']
│   │ 
│   │ 
│   └── SFs/
│       ├── __init__.py
│       └── sfPaths.py
│
├── docs
├── Logs
├── tests
├── Makefile
├── LICENSE
├── README.md
├── requirements.txt
├── .gitignore
└── binder.ipynb

```

  
