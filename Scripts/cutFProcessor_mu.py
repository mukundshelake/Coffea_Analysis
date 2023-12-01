import numpy as np
import hist
from coffea.analysis_tools import Weights
from coffea import processor
import awkward as ak
from coffea.nanoevents import NanoAODSchema
from coffea.lookup_tools import extractor
from coffea.analysis_tools import PackedSelection
import warnings
warnings.filterwarnings("ignore")
import correctionlib

class NanoProcessor(processor.ProcessorABC):
    def __init__(self):
        pass
        # Define histograms with hist class- https://hist.readthedocs.io/

    @property
    ## accumulating processed results
    def accumulator(self):
        return self._accumulator

    ## Place to process the Nanoevents, selections, weights, fill histogram is done here, considered as a loop with operating on columnar dataset
    def process(self, events):
        ## create the dictionary contains 
        # print(f"Working with {events.metadata['dataset']}")
        output = {
            "sumw" : processor.defaultdict_accumulator(float),
            "NoSel":{
                "selEvents" : processor.defaultdict_accumulator(float),
                "wtEvents" : processor.defaultdict_accumulator(float)
                },
            "HLT":{
                "selEvents" : processor.defaultdict_accumulator(float),
                "wtEvents" : processor.defaultdict_accumulator(float),
                },
            "HLTandGoodLep":{
                "selEvents" : processor.defaultdict_accumulator(float),
                "wtEvents" : processor.defaultdict_accumulator(float)
                },
            "HLTGoodLandThreeJ":{
                "selEvents" : processor.defaultdict_accumulator(float),
                "wtEvents" : processor.defaultdict_accumulator(float),
                },
            "HLTGoodLandGood3J":{
                "selEvents" : processor.defaultdict_accumulator(float),
                "wtEvents" : processor.defaultdict_accumulator(float),
                },
            "Total":{
                "selEvents" : processor.defaultdict_accumulator(float),
                "wtEvents" : processor.defaultdict_accumulator(float),
                "muon_pt":  hist.Hist(
                        hist.axis.Regular(60, 0, 300, name="pt", label="$p_T$ [GeV]"),
                        hist.storage.Weight(),
                        ),
                "muon_eta": hist.Hist(
                        hist.axis.Regular(60, -3.0, 3.0, name="eta", label="$ \eta $ "),
                        hist.storage.Weight(),
                        ),
                "jet_pt":  hist.Hist(
                        hist.axis.Regular(60, 0, 300, name="pt", label="$p_T$ [GeV]"),
                        hist.storage.Weight(),
                        ),
                "jet_eta": hist.Hist(
                        hist.axis.Regular(60, -3.0, 3.0, name="eta", label="$ \eta $ "),
                        hist.storage.Weight(),
                        )                    
                }
        }
        isRealData = not hasattr(events, "Generator")
        dataset = events.metadata["dataset"]
        if isRealData:
            output["sumw"] = len(events)
        else:
            output["sumw"] = ak.sum(events.Generator.weight)
        ####################
        #    Selections    #
        ####################
        selection = PackedSelection()
        selection.add("atleastOnelep", ak.num(events.Muon) > 0)
        selection.add("leadPtandEta",ak.sum((events.Muon.pt >= 26.0) & (events.Muon.pt <= 200.0) & (abs(events.Muon.eta) <= 2.4) & (events.Muon.tightId), axis=1) > 0)
        selection.add("atleastThreeJ", ak.num(events.Jet) > 2)
        selection.add("JetPtandEta", ak.sum((events.Jet.pt >= 30.0) & (abs(events.Jet.eta) < 3.0), axis = 1) > 2)
        selection.add("HLTTrigger", events.HLT.IsoTkMu24 | events.HLT.IsoMu24),
        selection.add("bTagMedium", ak.sum(events.Jet.btagDeepFlavB > 0.2598, axis=1) > 1)
        selectionList = {
        "NoSel":{},
        "HLT":{"HLTTrigger": True},
        "HLTandGoodLep":{"HLTTrigger": True,'leadPtandEta': True},
        "HLTGoodLandThreeJ":{'atleastThreeJ': True, "HLTTrigger": True,'leadPtandEta': True},
        "HLTGoodLandGood3J":{'atleastOnelep': True, 'leadPtandEta': True, "atleastThreeJ": True, "JetPtandEta": True, "HLTTrigger": True},
        "Total": {'atleastOnelep': True, 'leadPtandEta': True, "atleastThreeJ": True, "JetPtandEta": True, "HLTTrigger": True, "bTagMedium": True}
        }
        for region, cuts in selectionList.items():
            event_level = selection.require(**cuts)
            output[region]["selEvents"] = float(sum(event_level))
            if sum(event_level) == 0:
                print("No event selected")
                output[region]["wtEvents"] = float(sum(event_level))
                print(output)
                if region == 'Total':
                    return {dataset: output}
                continue
    
            ####################
            # Selected objects #
            ####################
            if "muon_pt" in output[region]:
                smu = events.Muon[(events.Muon.pt >= 26.0) & (events.Muon.pt <= 200.0) & (abs(events.Muon.eta) <= 2.4) & (events.Muon.tightId)][event_level][:, 0]
            if 'jet_pt' in output[region]:
                sjet = events.Jet[(events.Jet.pt >= 30.0) & (abs(events.Jet.eta) < 3.0)][event_level][:, 0]
            ####################
            # Weight & Geninfo #
            ####################
            weights = Weights(sum(event_level), storeIndividual=True)
            if isRealData:
                try:
                    weights.add("L1preFireWt", weight=events[event_level].L1PreFiringWeight.Nom, weightUp = events[event_level].L1PreFiringWeight.Up, weightDown = events[event_level].L1PreFiringWeight.Dn)
                    # print("Added L1preFireW")
                except:
                    print(f"L1preFireW is not there for dataset {dataset}; adding 1s as L1preFireW")
                    weights.add("L1preFireW", weight = np.ones(sum(event_level), dtype = float))
            else:
                try:
                    weights.add("PUWt", weight = events[event_level].puWeight, weightUp = events[event_level].puWeightUp, weightDown = events[event_level].puWeightDown)
                    # print("Added PUWt")
                except:
                    print(f"puWeight is not there for dataset {dataset}; adding 1s as puWeights")
                    weights.add("PUWt", weight = np.ones(sum(event_level), dtype = float))
                
                try:
                    weights.add("L1preFireWt", weight=events[event_level].L1PreFiringWeight.Nom, weightUp = events[event_level].L1PreFiringWeight.Up, weightDown = events[event_level].L1PreFiringWeight.Dn)
                    # print("Added L1preFireW")
                except:
                    print(f"L1preFireW is not there for dataset {dataset}; adding 1s as L1preFireW")
                    weights.add("L1preFireW", weight = np.ones(sum(event_level), dtype = float))
                    
                try:
                    weights.add("LHEWeightSign", weight = events[event_level].LHEWeight.originalXWGTUP/abs(events[event_level].LHEWeight.originalXWGTUP))
                    # print("Added LHEWeightSign")
                except:
                    print(f"LHEWeight is not there for dataset {dataset}; adding +1s as LHEWeightSign")
                    weights.add("LHEWeightSign", weight = np.ones(sum(event_level), dtype = float))
                
                if region == 'Total':
                    try: 
                        ceval = correctionlib.CorrectionSet.from_file("SFs/Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_SingleMuonTriggers_schemaV2.json")
    #                     ptflat, counts = ak.flatten(smu.pt), ak.num(smu.pt)
    #                     etaflat, counts = ak.flatten(smu.eta), ak.num(smu.eta)
    #                     muSF = ak.unflatten(ceval["NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight"].evaluate(ptflat, etaflat, "nominal"),counts=counts,)
                        muSF = ceval["NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight"].evaluate(abs(smu.eta), smu.pt, 'nominal')
                        weights.add("HLT_SF",weight=muSF)
                        print("Added HLT weights")
                    except:
                        print(f"HLT Scale factors is not working for dataset {dataset}; adding +1s as Scale factor weight")
                        weights.add("HLT_SF", weight = np.ones(sum(event_level), dtype = float))
                
            ####################
            #  Fill histogram  #
            ####################
            ## Filling the output dictionary

            output[region]["wtEvents"] = float(sum(weights.weight()))
            if "muon_pt" in output[region]:
                output[region]["muon_pt"].fill(ak.flatten(smu.pt, axis=-1), weight=weights.weight())
                output[region]["muon_eta"].fill(ak.flatten(smu.eta, axis=-1), weight=weights.weight())
        
            if "jet_pt" in output[region]:
                output[region]["jet_pt"].fill(ak.flatten(sjet.pt, axis=-1), weight=weights.weight())
                output[region]["jet_eta"].fill(ak.flatten(sjet.eta, axis=-1), weight=weights.weight())


        return {dataset: output}

    def postprocess(self, accumulator):
        return accumulator