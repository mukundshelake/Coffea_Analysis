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
        self.lep_abseta_max = 2.4 # HLT SF limit
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
        self.jet_abseta_max = 2.5  # HLT SF limit
        self.jet_etaHist_nBins = 60
        self.jet_etaHist_min = -2.5
        self.jet_etaHist_max = 2.5 
        
        # Tight ID configs
        self.muTightID_file = 'SFs/UL2016postVFP_mu_ID_SFs.json'
        
        # HLT configs
        self.muHLT_file = 'SFs/UL2016postVFP_muon_triggers.json'
        
        # B-Tagging configs
#         self.btagThreshold = 0.2598
        self.btagThreshold = 0.2489
        self.btagSF_file = 'SFs/UL2016postVFP_btagging.json'
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
                "muon_pt":  hist.Hist(
                        hist.axis.Regular(self.lep_ptHist_nBins, self.lep_ptHist_min, self.lep_ptHist_max, name="pt", label="$p_T$ [GeV]"),
                        hist.storage.Weight(),
                        ),
                "muon_eta": hist.Hist(
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
        selection.add("atleastOnelep", ak.num(events.Muon) >= self.nLep)
        selection.add("TightLpt_eta", ak.sum((events.Muon.pt >= self.lep_pt_min) & (events.Muon.pt <= self.lep_pt_max) & (abs(events.Muon.eta) <= self.lep_abseta_max) & (events.Muon.tightId), axis=1) >= self.nLep)
        selection.add("atleastThreeJ", ak.num(events.Jet) >= self.nJet)
        selection.add("JetPtandEta", ak.sum((events.Jet.pt >= self.jet_pt_min) & (abs(events.Jet.eta) < self.jet_abseta_max), axis = 1) >= self.nJet)
        selection.add("HLT", events.HLT.IsoTkMu24 | events.HLT.IsoMu24)
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
                smu = events.Muon[(events.Muon.pt >= self.lep_pt_min) & (events.Muon.pt <= self.lep_pt_max) & (abs(events.Muon.eta) <= self.lep_abseta_max) & (events.Muon.tightId)][event_level][:, 0]
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
#                     print("%^%^%^%^%^*&^%^%^$#%#%")
                        TReval = correctionlib.CorrectionSet.from_file(self.muHLT_file)
                        TRsf = ak.unflatten(TReval["NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight"].evaluate(smu.eta, smu.pt, "nominal"),counts = ak.num(smu, axis=0))[0,:]
                        weights.add("TrSF",weight= TRsf)
                        
                if region=="Total":
                        print("gfgdgdfgdfg")
                        IDeval = correctionlib.CorrectionSet.from_file(self.muTightID_file)
                        IDsf = ak.unflatten(IDeval["NUM_TightID_DEN_genTracks"].evaluate("2016postVFP_UL", abs(smu.eta), smu.pt, 'sf'),counts = ak.num(smu, axis=0))[0,:]
                        IDsfup = ak.unflatten(IDeval["NUM_TightID_DEN_genTracks"].evaluate("2016postVFP_UL", abs(smu.eta), smu.pt, 'systup'),counts = ak.num(smu, axis=0))[0,:] 
                        IDsfdown = ak.unflatten(IDeval["NUM_TightID_DEN_genTracks"].evaluate("2016postVFP_UL", abs(smu.eta), smu.pt, 'systdown'),counts = ak.num(smu, axis=0))[0,:] 
                        weights.add("IDSF",weight= IDsf, weightUp= IDsfup, weightDown = IDsfdown)
                if region=='Total':
                        Btageval = correctionlib.CorrectionSet.from_file(self.btagSF_file) 
                        Btagsf = ak.unflatten(Btageval["deepJet_mujets"].evaluate("central","M",5,abs(sjet.eta), sjet.pt),counts= ak.num(sjet, axis=0))[0,:]
                        Btagsfdown = ak.unflatten(Btageval["deepJet_mujets"].evaluate("down","M",5,abs(sjet.eta), sjet.pt),counts= ak.num(sjet, axis=0))[0,:]
                        Btagsfup = ak.unflatten(Btageval["deepJet_mujets"].evaluate("up","M",5,abs(sjet.eta), sjet.pt),counts= ak.num(sjet, axis=0))[0,:]
                        weights.add("BtagSF",weight= Btagsf, weightUp= Btagsfup, weightDown = Btagsfdown)
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