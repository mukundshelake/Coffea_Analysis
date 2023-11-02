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
        print(f"Working with {events.metadata['dataset']}")
        output = {
            "sumw" : processor.defaultdict_accumulator(float),
            "NoSel":{
                "selEvents" : processor.defaultdict_accumulator(float),
                "wtEvents" : processor.defaultdict_accumulator(float)
            },
            "AtleastOnelep":{
                "selEvents" : processor.defaultdict_accumulator(float),
                "wtEvents" : processor.defaultdict_accumulator(float)
            },
            "LeadLep":{
                "selEvents" : processor.defaultdict_accumulator(float),
                "wtEvents" : processor.defaultdict_accumulator(float)
            },
            "AtleastThreeJ":{
                "selEvents" : processor.defaultdict_accumulator(float),
                "wtEvents" : processor.defaultdict_accumulator(float)
            },
            "JetPtEta":{
                "selEvents" : processor.defaultdict_accumulator(float),
                "wtEvents" : processor.defaultdict_accumulator(float)
            },
            "HLTcut":{
                "selEvents" : processor.defaultdict_accumulator(float),
                "wtEvents" : processor.defaultdict_accumulator(float)
            },
            "LeadLandHLT":{
                "selEvents" : processor.defaultdict_accumulator(float),
                "wtEvents" : processor.defaultdict_accumulator(float)
            },
            "LeadLandthreeJandHLT":{
                "selEvents" : processor.defaultdict_accumulator(float),
                "wtEvents" : processor.defaultdict_accumulator(float)
            },
            "LeadLandthreeJ":{
                "selEvents" : processor.defaultdict_accumulator(float),
                "wtEvents" : processor.defaultdict_accumulator(float)
            },
            "LeadLandgoodJ":{
                "selEvents" : processor.defaultdict_accumulator(float),
                "wtEvents" : processor.defaultdict_accumulator(float)
            },
            "Total":{
                "selEvents" : processor.defaultdict_accumulator(float),
                "wtEvents" : processor.defaultdict_accumulator(float)
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
        selection.add(
            "leadPtandEta",
            ak.sum((events.Muon.pt >= 35.0) & (abs(events.Muon.eta) <= 3.0), axis=1) > 0
        )

        selection.add("atleastThreeJ", ak.num(events.Jet) > 2)
        selection.add(
            "JetPtandEta",
            ak.sum((events.Jet.pt >= 20.0) & (abs(events.Jet.eta) < 3.0), axis = 1) > 2
        )
        selection.add("HLTIsoTk24", events.HLT.IsoTkMu24)
        selectionList = {
        "NoSel":{},
        "AtleastOnelep":{'atleastOnelep': True},
        "LeadLep":{'leadPtandEta': True},
        "AtleastThreeJ":{'atleastThreeJ': True},
        "JetPtEta":{"JetPtandEta": True},
        "HLTcut":{"HLTIsoTk24": True},
        "LeadLandHLT":{"HLTIsoTk24": True, 'leadPtandEta': True},
        "LeadLandthreeJandHLT":{'atleastThreeJ': True, 'leadPtandEta': True, "HLTIsoTk24": True},
        "LeadLandthreeJ":{'atleastThreeJ': True, 'leadPtandEta': True},
        "LeadLandgoodJ":{'atleastThreeJ': True, 'leadPtandEta': True, "JetPtandEta": True},
        "Total":{'atleastOnelep': True, 'leadPtandEta': True, "atleastThreeJ": True, "JetPtandEta": True, "HLTIsoTk24": True}
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

            ####################
            # Weight & Geninfo #
            ####################
            weights = Weights(sum(event_level), storeIndividual=True)
            if isRealData:
                try:
                    weights.add("L1preFireWt", weight=events[event_level].L1PreFiringWeight.Nom, weightUp = events[event_level].L1PreFiringWeight.Up, weightDown = events[event_level].L1PreFiringWeight.Dn)
                    print("Added L1preFireW")
                except:
                    print(f"L1preFireW is not there for dataset {dataset}; adding 1s as L1preFireW")
                    weights.add("L1preFireW", weight = np.ones(sum(event_level), dtype = float))
            else:
                try:
                    weights.add("PUWt", weight = events[event_level].puWeight, weightUp = events[event_level].puWeightUp, weightDown = events[event_level].puWeightDown)
                    print("Added PUWt")
                except:
                    print(f"puWeight is not there for dataset {dataset}; adding 1s as puWeights")
                    weights.add("PUWt", weight = np.ones(sum(event_level), dtype = float))
                
                try:
                    weights.add("L1preFireWt", weight=events[event_level].L1PreFiringWeight.Nom, weightUp = events[event_level].L1PreFiringWeight.Up, weightDown = events[event_level].L1PreFiringWeight.Dn)
                    print("Added L1preFireW")
                except:
                    print(f"L1preFireW is not there for dataset {dataset}; adding 1s as L1preFireW")
                    weights.add("L1preFireW", weight = np.ones(sum(event_level), dtype = float))
                    
                try:
                    weights.add("LHEWeightSign", weight = events[event_level].LHEWeight.originalXWGTUP/abs(events[event_level].LHEWeight.originalXWGTUP))
                    print("Added LHEWeightSign")
                except:
                    print(f"LHEWeight is not there for dataset {dataset}; adding +1s as LHEWeightSign")
                    weights.add("LHEWeightSign", weight = np.ones(sum(event_level), dtype = float))
                
                # if "HLTIsoTk24" in cuts:
                #     try:
                #         ext = extractor()
                #         ext.add_weight_sets(["* * SFs/UL2016_postVFP_Tight.root"])
                #         ext.finalize()
                #         evaluator = ext.make_evaluator()
                #         eleSF = evaluator["EGamma_SF2D"](events.Muon[event_level][:, 0].eta, events.Muon[event_level][:, 0].pt)
                #         eleSFerror = evaluator["EGamma_SF2D_error"](events.Muon[event_level][:, 0].eta, events.Muon[event_level][:, 0].pt)
                #         weights.add("eleSF",weight=eleSF,weightUp=eleSF + eleSFerror,weightDown = eleSF - eleSFerror)
                #         print("Added HLT SF")
                #     except:
                #         print(f"HLT Scale factors is not working for dataset {dataset}; adding +1s as Scale factor weight")
                #         weights.add("eleSF", weight = np.ones(sum(event_level), dtype = float))
            ####################
            #  Fill histogram  #
            ####################
            ## Filling the output dictionary
            output[region]["wtEvents"] = float(sum(weights.weight()))
        return {dataset: output}

    def postprocess(self, accumulator):
        return accumulator