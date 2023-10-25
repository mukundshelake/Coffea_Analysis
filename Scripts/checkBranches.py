import awkward as ak
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
from lib.helpers import getfileset, getFiles


era = 'SIXTEEN_preVFP'
lep = 'el'
# DataDir = f'/nfs/home/common/RUN2_UL/Tree_crab/{era}/Data_{lep}'
MCDir = f'/nfs/home/common/RUN2_UL/Tree_crab/{era}/MC'
inputDir = [MCDir]
fileset = getfileset(inputDir)

for s in fileset:
    for file in (fileset[s]):
        # print(file)
        events = NanoEventsFactory.from_root(
            file,
            schemaclass=NanoAODSchema.v6,
            metadata={"dataset": "DYJets"},
        ).events()
        if not hasattr(events, "LHEWeight"):
            print(f"{file} has no LHEWeight")
