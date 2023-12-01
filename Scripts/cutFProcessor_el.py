import numpy as np
import hist
from coffea.analysis_tools import Weights
from coffea import processor
import awkward as ak
from coffea.nanoevents import NanoAODSchema
from coffea.lookup_tools import extractor
from coffea.analysis_tools import PackedSelection
import warnings
import correctionlib
warnings.filterwarnings("ignore")


class NanoProcessor(processor.ProcessorABC):
    def __init__(self):
        
        self.nLep = 1 # greater than equal to
        self.nJet = 3 # greater than equal to
        self.nbJet = 2 # greater than equal to
        # Lepton pT thresholds
        self.lep_pt_min = 35.0
        self.lep_pt_max = 200.0
        self.lep_ptHist_nBins = 60
        self.lep_ptHist_min = 0
        self.lep_ptHist_max = 300
        
        # Lepton eta thresholds
        self.lep_abseta_max = 3.0
        self.lep_etaHist_nBins = 60
        self.lep_etaHist_min = -3.0
        self.lep_etaHist_max = 3.0
        
        # AK4 jets pt kinematics
        
        self.jet_pt_min = 30.0
        self.jet_pt_max = 300.0
        self.jet_ptHist_nBins = 60
        self.jet_ptHist_min = 0.0
        self.jet_ptHist_max = 300.0
        
        # AK4 jets eta kinematics
        self.jet_abseta_max = 3.0
        self.jet_etaHist_nBins = 60
        self.jet_etaHist_min = -3.0
        self.jet_etaHist_max = 3.0
        
        # Tight ID configs
        self.eleTightID_idx = 4
        self.eleTightID_file = 'SFs/UL2016preVFP_ele_ID_SFs.json'
        
        # B-Tagging configs
        self.btagThreshold = 0.2598
        
        # SFs
        ##For "Medium" requires (sel.eta, sel.pt) ranges (-inf --> inf, 10 --> inf)
        
        
        
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
                "wtEvents" : processor.defaultdict_accumulator(float)
                },
            "HLTGoodLandGood3J":{
                "selEvents" : processor.defaultdict_accumulator(float),
                "wtEvents" : processor.defaultdict_accumulator(float)
                },
            "Total":{
                "selEvents" : processor.defaultdict_accumulator(float),
                "wtEvents" : processor.defaultdict_accumulator(float),
                "electron_pt":  hist.Hist(
                        hist.axis.Regular(self.lep_ptHist_nBins, self.lep_ptHist_min, self.lep_ptHist_max, name="pt", label="$p_T$ [GeV]"),
                        hist.storage.Weight(),
                        ),
                "electron_eta": hist.Hist(
                        hist.axis.Regular(self.lep_etaHist_nBins, self.lep_etaHist_min, self.lep_etaHist_max, name="eta", label="$ \eta $ "),
                        hist.storage.Weight(),
                        ),
                "jet_pt":  hist.Hist(
                        hist.axis.Regular(self.jet_ptHist_nBins, self.jet_ptHist_min, self.jet_ptHist_max, name="pt", label="$p_T$ [GeV]"),
                        hist.storage.Weight(),
                        ),
                "jet_eta": hist.Hist(
                        hist.axis.Regular(self.jet_etaHist_nBins, self.jet_etaHist_min, self.jet_etaHist_max, name="eta", label="$ \eta $ "),
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
        selection.add("atleastOnelep", ak.num(events.Electron) >= self.nLep)
        selection.add("TightLpt_eta", ak.sum((events.Electron.pt >= self.lep_pt_min) & (abs(events.Electron.eta) <= self.lep_abseta_max) & (events.Electron.cutBased == self.eleTightID_idx), axis=1) >= self.nLep)
        selection.add("atleastThreeJ", ak.num(events.Jet) >= self.nJet)
        selection.add("JetPtandEta", ak.sum((events.Jet.pt >= self.jet_pt_min) & (abs(events.Jet.eta) < self.jet_abseta_max), axis = 1) >= self.nJet)
        selection.add("HLT", events.HLT.Ele32_eta2p1_WPTight_Gsf)
        selection.add("mediumBTagGoodJets", ak.sum((events.Jet.pt >= self.jet_pt_min) & (abs(events.Jet.eta) < self.jet_abseta_max) & (events.Jet.btagDeepFlavB > self.btagThreshold), axis = 1) >= self.nbJet)
        selectionList = {
        "NoSel":{},
        "HLT":{"HLT": True},
        "HLTandGoodLep":{"HLT": True,'TightLpt_eta': True},
        "HLTGoodLandThreeJ":{'atleastThreeJ': True, "HLT": True,'TightLpt_eta': True},
        "HLTGoodLandGood3J":{'atleastOnelep': True, 'TightLpt_eta': True, "atleastThreeJ": True, "JetPtandEta": True, "HLT": True},
        "Total":{'atleastOnelep': True, 'TightLpt_eta': True, "atleastThreeJ": True, "JetPtandEta": True, "HLT": True, "mediumBTagGoodJets": True}
        }
        for region, cuts in selectionList.items():
            event_level = selection.require(**cuts)
            output[region]["selEvents"] = float(sum(event_level))
            if sum(event_level) == 0:
                # print("No event selected")
                output[region]["wtEvents"] = float(sum(event_level))
                # print(output)
                if region == 'Total':
                    return {dataset: output}
                continue
    
            ####################
            # Selected objects #
            ####################
            if "electron_pt" in output[region]:
                sel = events.Electron[(events.Electron.pt >= self.lep_pt_min) & (abs(events.Electron.eta) <= self.lep_abseta_max) & (events.Electron.cutBased == self.eleTightID_idx)][event_level][:, 0]
            if 'jet_pt' in output[region]:
#                 sjet = events.Jet[(events.Jet.pt >= 30.0) & (abs(events.Jet.eta) < 3.0)][event_level][:, 0]
                sjet = events.Jet[(events.Jet.pt >= self.jet_pt_min) & (abs(events.Jet.eta) < self.jet_abseta_max) & (events.Jet.btagDeepFlavB > self.btagThreshold)][event_level][:, 0]
                
                

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
                
                if region=="Total":
                    # print("%^%^%^%^%^*&^%^%^$#%#%")
                    try:
                        ext = extractor()
                        ext.add_weight_sets(["* * SFs/UL2016_preVFP_Tight.root"])
                        ext.finalize()
                        evaluator = ext.make_evaluator()
                        HLTSF = evaluator["EGamma_SF2D"](sel.eta, sel.pt)
                        HLTSFerror = evaluator["EGamma_SF2D_error"](sel.eta, sel.pt)
                        weights.add("TriggerSF",weight=HLTSF,weightUp=HLTSF + HLTSFerror,weightDown = HLTSF - HLTSFerror)
                        # print("Added HLT SF")
                    except Exception as e:
                        print(f"error is {e}")
                        print(f"HLT Scale factors is not working for dataset {dataset}; adding +1s as Scale factor weight")
                        weights.add("TriggerSF", weight = np.ones(sum(event_level), dtype = float))
                        
                if region=="Total":
                        print("gfgdgdfgdfg")
                        elIDeval = correctionlib.CorrectionSet.from_file(self.eleTightID_file)
                        IDsf = ak.unflatten(
                            elIDeval["UL-Electron-ID-SF"].evaluate("2016preVFP","sf","Medium",sel.eta, sel.pt),
                            counts= ak.num(sel, axis=-1))[0,:]
                        IDsfup = ak.unflatten(
                            elIDeval["UL-Electron-ID-SF"].evaluate("2016preVFP","sfup","Medium",sel.eta, sel.pt),
                            counts= ak.num(sel, axis=-1))[0,:]
                        IDsfdown = ak.unflatten(
                            elIDeval["UL-Electron-ID-SF"].evaluate("2016preVFP","sfdown","Medium",sel.eta, sel.pt),
                            counts= ak.num(sel, axis=-1))[0,:]  
                        weights.add("IDSF",weight= IDsf, weightUp= IDsfup, weightDown = IDsfdown)
            ####################
            #  Fill histogram  #
            ####################
            ## Filling the output dictionary

            output[region]["wtEvents"] = float(sum(weights.weight()))
            if "electron_pt" in output[region]:
                output[region]["electron_pt"].fill(ak.flatten(sel.pt, axis=-1), weight=weights.weight())
                output[region]["electron_eta"].fill(ak.flatten(sel.eta, axis=-1), weight=weights.weight())
        
            if "jet_pt" in output[region]:
                output[region]["jet_pt"].fill(ak.flatten(sjet.pt, axis=-1), weight=weights.weight())
                output[region]["jet_eta"].fill(ak.flatten(sjet.eta, axis=-1), weight=weights.weight())


        return {dataset: output}

    def postprocess(self, accumulator):
        return accumulator