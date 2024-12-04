import json, os, argparse
import hist
import dask
import awkward as ak
import hist.dask as hda
import dask_awkward as dak
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
from coffea import processor
from coffea.nanoevents.methods import candidate
from coffea.dataset_tools import (
    apply_to_fileset,
    max_chunks,
    preprocess,
)
from distributed import Client
from coffea.analysis_tools import PackedSelection, Weights
from coffea.util import save
import numpy as np
from coffea.lookup_tools import extractor
import correctionlib

class MyProcessor(processor.ProcessorABC):
    def __init__(self):
        pass

    def process(self, events):
        dataset = events.metadata['dataset']
        print(dataset)

        parts = dataset.split('_', 1)

        # Extract the era and channel
        era = parts[0]
        channel = parts[1]

        isData = False
        if 'Run' in channel:
            isData = True

        maps = {
            'btagThreshold': {
                'UL2016preVFP': 0.2598,
                'UL2016postVFP': 0.2489,
                'UL2017': 0.3040,
                'UL2018': 0.2783,
            }
        }

        selection = PackedSelection()
        # weight = Weights.weight()

        selection.add_multiple(
            {
                "atleastOneLep": ak.num(events.Muon) > 0,
                "atleastThreeJ": ak.num(events.Jet) > 2,
                "goodLeps" : ak.sum((events.Muon.pt >= 35.0) & (abs(events.Muon.eta) <= 2.4) & (events.Muon.tightId), axis=1) >= 1,
                "goodJets" : ak.sum((events.Jet.pt >= 30.0) & (abs(events.Jet.eta) < 2.4), axis = 1) >= 3,
                "BTag"  : ak.sum((events.Jet.pt >= 30.0) & (abs(events.Jet.eta) < 2.4) & (events.Jet.btagDeepFlavB > maps['btagThreshold'][era]), axis = 1) >= 2
            }
        )
        if 'UL2016' in dataset:
            selection.add("HLT", events.HLT.IsoTkMu24 | events.HLT.IsoMu24)
        elif 'UL2017' in dataset:
            selection.add("HLT", events.HLT.IsoMu27)
        else:
            selection.add("HLT", events.HLT.IsoMu24)
        # mask = selection.all("atleastOneLep", "atleastThreeJ")
        cutflow = selection.cutflow("atleastOneLep", "atleastThreeJ", "goodLeps", "goodJets", "BTag", "HLT")

        # honecut, hcutflow, labels = cutflow.yieldhist()

        results = cutflow.result()

        finalMask = results[-1][-1]

        nSelected = np.sum(finalMask.compute())

        # muon_hist = (
        #     hda.Hist.new
        #     .Reg(20, 0, 600.0, label="$Muon_pt$", name="muonPt")
        #     .Reg(24, -2.4, 2.4, label="$Muon_eta$", name="muonEta")
        #     .Double()
        # )

        # leading_muon_hist = (
        #     hda.Hist.new
        #     .Reg(20, 0, 600.0, label="$leadingMuon_pt$", name="leadingmuonPt")
        #     .Reg(24, -2.4, 2.4, label="$leadingMuon_eta$", name="leadingmuonEta")
        #     .Double()
        # )

        # jet_hist = (
        #     hda.Hist.new
        #     .Reg(20, 0, 1000.0, label="$jet_pt$", name="jetPt")
        #     .Reg(24, -2.4, 2.4, label="$jet_eta$", name="jetEta")
        #     .Double()
        # )
        if nSelected == 0:
            print(f"No events selected for {dataset}")


        if nSelected == 0:
            return {
                dataset: {
                    "entries": ak.num(events, axis=0),
                    "nSelected": nSelected,
                    "nWeighted": 0,
                    # "weightStats": weights.weightStatistics
                }
            }
        # muons = events.Muon[(events.Muon.pt >= 35.0) & (abs(events.Muon.eta) <= 2.4) & (events.Muon.tightId)][finalMask]
        # spt_flat = ak.flatten(muons.pt)
        # seta_flat = ak.flatten(muons.eta)
        # leading_muon = muons[ak.argmax(muons.pt, axis=1, keepdims=True)]
        # leading_muon_pt = ak.flatten(leading_muon.pt)
        # leading_muon_eta = ak.flatten(leading_muon.eta)
        

        # sjets = events.Jet[(events.Jet.pt >= 30.0) & (abs(events.Jet.eta) < 2.4) & (events.Jet.btagDeepFlavB > maps['btagThreshold'][era])][finalMask]
        # leadingJet = sjets[ak.argmax(sjets.pt, axis=1, keepdims=True)]
        # leadingJet_pt = ak.flatten(leadingJet.pt)
        # leadingJet_eta = ak.flatten(leadingJet.eta)


        # met_hist = (
        #     hda.Hist.new
        #     .Reg(20, 0, 500.0, label="$MET_pt$", name="metPt")
        #     .Reg(24, -3.2, 3.2, label="$MET_phi$", name="metPhi")
        #     .Double()
        # )

        # met_pt = events.MET.pt[finalMask]
        # met_phi = events.MET.phi[finalMask]

        # met_hist.fill(
        #     metPt=met_pt,
        #     metPhi=met_phi
        # )



        weights = Weights(nSelected, storeIndividual=True)
        

        if isData:
            weights.add("L1preFireWt", weight=events[finalMask].L1PreFiringWeight.Nom.compute(), weightUp = events[finalMask].L1PreFiringWeight.Up.compute(), weightDown = events[finalMask].L1PreFiringWeight.Dn.compute())
        else:
            weights.add("PUWt", weight = events[finalMask].puWeight.compute(), weightUp = events[finalMask].puWeightUp.compute(), weightDown = events[finalMask].puWeightDown.compute())
            weights.add("L1preFireWt", weight=events[finalMask].L1PreFiringWeight.Nom.compute(), weightUp = events[finalMask].L1PreFiringWeight.Up.compute(), weightDown = events[finalMask].L1PreFiringWeight.Dn.compute())

            if hasattr(events, "LHEWeight"):
                weights.add("LHEWeightSign", weight = events[finalMask].LHEWeight.originalXWGTUP.compute()/abs(events[finalMask].LHEWeight.originalXWGTUP.compute()))
            spt = events.Muon[(events.Muon.pt >= 35.0)  & (abs(events.Muon.eta) <= 2.4) & (events.Muon.tightId)].pt[finalMask]
            seta = events.Muon[(events.Muon.pt >= 35.0) & (abs(events.Muon.eta) <= 2.4) & (events.Muon.tightId)].eta[finalMask]

            # leadingMuon_pt = ak.max(spt, axis=1)
            # leadingMuon_eta = ak.max(seta, axis=1)


            IDFile = f"/nfs/home/mukund/Projects/updatedCoffea/Coffea_Analysis/src/SFs/{era}_mu_ID.json"
            IDeval = correctionlib.CorrectionSet.from_file(IDFile)
            IDSF = IDeval["NUM_TightID_DEN_TrackerMuons"].evaluate(f"{era[2:]}_UL", abs(seta), spt, 'sf')
            IDup = IDeval["NUM_TightID_DEN_TrackerMuons"].evaluate(f"{era[2:]}_UL", abs(seta), spt, 'systup')
            IDdown = IDeval["NUM_TightID_DEN_TrackerMuons"].evaluate(f"{era[2:]}_UL", abs(seta), spt, 'systdown')
            IDweight = ak.prod(IDSF, axis=1)
            IDweightUp = ak.prod(IDup, axis=1)
            IDweightDown = ak.prod(IDdown, axis=1)


            weights.add("ID",weight=IDweight.compute(),weightUp=IDweightUp.compute(),weightDown = IDweightDown.compute())

            # try:
            #     HLTFile = f"/nfs/home/mukund/Projects/updatedCoffea/Coffea_Analysis/src/SFs/{era}_mu_HLT.json"
            #     HLTeval = correctionlib.CorrectionSet.from_file(HLTFile)
            #     if era == "UL2016preVFP" or era == "UL2016postVFP":
            #         sfString = "NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight"
            #     elif era == "UL2017":
            #         sfString = "NUM_IsoMu27_DEN_CutBasedIdTight_and_PFIsoTight"
            #     else:
            #         sfString = "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight"
            #     HLTSF = HLTeval[sfString].evaluate(abs(seta), spt, 'nominal')
            #     HLTweight = ak.prod(HLTSF, axis=1)
            #     weights.add("HLT",weight=HLTweight.compute())
            # except Exception as e:
            #     print(f"HLT SF not found for {era}_{channel}; skipping")
            #     print(e)

            bjets = events.Jet[(events.Jet.pt >= 30.0) & (abs(events.Jet.eta) < 2.4) & (events.Jet.btagDeepFlavB > maps['btagThreshold'][era])][finalMask]
            nonbjets = events.Jet[(events.Jet.pt >= 30.0) & (abs(events.Jet.eta) < 2.4) & (events.Jet.btagDeepFlavB <= maps['btagThreshold'][era])][finalMask]


            leadingJet_pt = ak.max(bjets.pt, axis=1)
            leadingJet_eta = ak.max(bjets.eta, axis=1)

            bjets_BC = bjets[(bjets.hadronFlavour == 4) | (bjets.hadronFlavour == 5)]
            bjets_L = bjets[(bjets.hadronFlavour == 0)]
            bjets_BC_pt = bjets_BC.pt.compute()
            bjets_L_pt = bjets_L.pt.compute()
            bjets_BC_eta = bjets_BC.eta.compute()
            bjets_L_eta = bjets_L.eta.compute()
            bjets_BC_fl = bjets_BC.hadronFlavour.compute()
            bjets_L_fl = bjets_L.hadronFlavour.compute()

            nonbjets_B = nonbjets[nonbjets.hadronFlavour == 4]
            nonbjets_C = nonbjets[nonbjets.hadronFlavour == 5]
            nonbjets_L = nonbjets[nonbjets.hadronFlavour == 0]
            nonbjets_B_pt = nonbjets_B.pt.compute()
            nonbjets_C_pt = nonbjets_C.pt.compute()
            nonbjets_L_pt = nonbjets_L.pt.compute()
            nonbjets_B_eta = nonbjets_B.eta.compute()
            nonbjets_C_eta = nonbjets_C.eta.compute()
            nonbjets_L_eta = nonbjets_L.eta.compute()
            nonbjets_B_fl = nonbjets_B.hadronFlavour.compute()
            nonbjets_C_fl = nonbjets_C.hadronFlavour.compute()
            nonbjets_L_fl = nonbjets_L.hadronFlavour.compute()

            btagging_evaluator = correctionlib.CorrectionSet.from_file(f'/nfs/home/mukund/Projects/updatedCoffea/Coffea_Analysis/src/SFs/{era}_btagging.json.gz')

            effiFile = f"/nfs/home/mukund/Projects/updatedCoffea/Coffea_Analysis/src/SFs/Efficiency/{era}/{channel}.root"
            b_eff_ext = extractor()
            b_eff_ext.add_weight_sets(["* * "+effiFile])
            b_eff_ext.finalize()
            b_eff_evaluator = b_eff_ext.make_evaluator()

            btag_effi_Bpass = b_eff_evaluator['Efficiency/FlavourB_Wp_pass_BM'](abs(nonbjets_B_eta), nonbjets_B_pt)
            btag_effi_Ball = b_eff_evaluator['Efficiency/FlavourB_Wp_pass_No'](abs(nonbjets_B_eta), nonbjets_B_pt) 
            btag_effi_B = btag_effi_Bpass/btag_effi_Ball

            btag_effi_Cpass = b_eff_evaluator['Efficiency/FlavourC_Wp_pass_BM'](abs(nonbjets_C_eta), nonbjets_C_pt)
            btag_effi_Call = b_eff_evaluator['Efficiency/FlavourC_Wp_pass_No'](abs(nonbjets_C_eta), nonbjets_C_pt)
            btag_effi_C = btag_effi_Cpass/btag_effi_Call

            btag_effi_Lpass = b_eff_evaluator['Efficiency/FlavourL_Wp_pass_BM'](abs(nonbjets_L_eta), nonbjets_L_pt)
            btag_effi_Lall = b_eff_evaluator['Efficiency/FlavourL_Wp_pass_No'](abs(nonbjets_L_eta), nonbjets_L_pt)
            btag_effi_L = btag_effi_Lpass/btag_effi_Lall


            Sel_btag_jets_SFBC = btagging_evaluator['deepJet_mujets'].evaluate('central', 'M', bjets_BC_fl, abs(bjets_BC_eta), bjets_BC_pt)
            Sel_btag_jets_SFL = btagging_evaluator['deepJet_incl'].evaluate('central', 'M', bjets_L_fl, abs(bjets_L_eta), bjets_L_pt)
            Rej_btag_jets_SFB = btagging_evaluator['deepJet_mujets'].evaluate('central', 'M', nonbjets_B_fl, abs(nonbjets_B_eta), nonbjets_B_pt) 
            Rej_btag_B = (1 - Rej_btag_jets_SFB*btag_effi_B)/(1 - btag_effi_B)
            Rej_btag_jets_SFC = btagging_evaluator['deepJet_mujets'].evaluate('central', 'M', nonbjets_C_fl, abs(nonbjets_C_eta), nonbjets_C_pt) 
            Rej_btag_C = (1 - Rej_btag_jets_SFC*btag_effi_C)/(1 - btag_effi_C)
            Rej_btag_jets_SFL = btagging_evaluator['deepJet_incl'].evaluate('central', 'M', nonbjets_L_fl, abs(nonbjets_L_eta), nonbjets_L_pt) 
            Rej_btag_L = (1 - Rej_btag_jets_SFL*btag_effi_L)/(1 - btag_effi_L)

            btag_SF = ak.prod(Sel_btag_jets_SFBC, axis=1)*ak.prod(Sel_btag_jets_SFL, axis=1)*ak.prod(Rej_btag_B, axis=1)*ak.prod(Rej_btag_C, axis=1)*ak.prod(Rej_btag_L, axis=1)

            Sel_btag_jets_SFBC_up = btagging_evaluator['deepJet_mujets'].evaluate('up', 'M', bjets_BC_fl, abs(bjets_BC_eta), bjets_BC_pt)
            Sel_btag_jets_SFL_up = btagging_evaluator['deepJet_incl'].evaluate('up', 'M', bjets_L_fl, abs(bjets_L_eta), bjets_L_pt)
            Rej_btag_jets_SFB_up = btagging_evaluator['deepJet_mujets'].evaluate('up', 'M', nonbjets_B_fl, abs(nonbjets_B_eta), nonbjets_B_pt) 
            Rej_btag_B_up = (1 - Rej_btag_jets_SFB_up*btag_effi_B)/(1 - btag_effi_B)
            Rej_btag_jets_SFC_up = btagging_evaluator['deepJet_mujets'].evaluate('up', 'M', nonbjets_C_fl, abs(nonbjets_C_eta), nonbjets_C_pt) 
            Rej_btag_C_up = (1 - Rej_btag_jets_SFC_up*btag_effi_C)/(1 - btag_effi_C)
            Rej_btag_jets_SFL_up = btagging_evaluator['deepJet_incl'].evaluate('up', 'M', nonbjets_L_fl, abs(nonbjets_L_eta), nonbjets_L_pt) 
            Rej_btag_L_up = (1 - Rej_btag_jets_SFL_up*btag_effi_L)/(1 - btag_effi_L)

            btag_SF_up = ak.prod(Sel_btag_jets_SFBC_up, axis=1)*ak.prod(Sel_btag_jets_SFL_up, axis=1)*ak.prod(Rej_btag_B_up, axis=1)*ak.prod(Rej_btag_C_up, axis=1)*ak.prod(Rej_btag_L_up, axis=1)

            Sel_btag_jets_SFBC_down = btagging_evaluator['deepJet_mujets'].evaluate('down', 'M', bjets_BC_fl, abs(bjets_BC_eta), bjets_BC_pt)
            Sel_btag_jets_SFL_down = btagging_evaluator['deepJet_incl'].evaluate('down', 'M', bjets_L_fl, abs(bjets_L_eta), bjets_L_pt)
            Rej_btag_jets_SFB_down = btagging_evaluator['deepJet_mujets'].evaluate('down', 'M', nonbjets_B_fl, abs(nonbjets_B_eta), nonbjets_B_pt) 
            Rej_btag_B_down = (1 - Rej_btag_jets_SFB_down*btag_effi_B)/(1 - btag_effi_B)
            Rej_btag_jets_SFC_down = btagging_evaluator['deepJet_mujets'].evaluate('down', 'M', nonbjets_C_fl, abs(nonbjets_C_eta), nonbjets_C_pt) 
            Rej_btag_C_down = (1 - Rej_btag_jets_SFC_down*btag_effi_C)/(1 - btag_effi_C)
            Rej_btag_jets_SFL_down = btagging_evaluator['deepJet_incl'].evaluate('down', 'M', nonbjets_L_fl, abs(nonbjets_L_eta), nonbjets_L_pt) 
            Rej_btag_L_down = (1 - Rej_btag_jets_SFL_down*btag_effi_L)/(1 - btag_effi_L)

            btag_SF_down = ak.prod(Sel_btag_jets_SFBC_down, axis=1)*ak.prod(Sel_btag_jets_SFL_down, axis=1)*ak.prod(Rej_btag_B_down, axis=1)*ak.prod(Rej_btag_C_down, axis=1)*ak.prod(Rej_btag_L_down, axis=1)

            weights.add("btag", weight=btag_SF, weightUp = btag_SF_up, weightDown = btag_SF_down)

            weights.add("btag", weight=btag_SF)


        # muon_hist.fill(
        #     muonPt = spt_flat,
        #     muonEta = seta_flat
        # )

        # leading_muon_hist.fill(
        #     leadingmuonPt = leading_muon_pt,
        #     leadingmuonEta = leading_muon_eta
        # )

        # jet_hist.fill(
        #     jetPt = leadingJet_pt,
        #     jetEta = leadingJet_eta
        # )


        return {
            dataset: {
                "entries": ak.num(events, axis=0),
                "nSelected": nSelected,
                "nWeighted": weights.weight().sum()
                # "muonHist": muon_hist,
                # "leadingMuonHist": leading_muon_hist,
                # "jetHist": jet_hist,
                # "metHist": met_hist
                # "weightStats": weights.weightStatistics
            }
        }

    def postprocess(self, accumulator):
        pass

def main():
    parser = argparse.ArgumentParser(description="Process some eras.")
    
    # Define the allowed choices
    allowed_eras = ['UL2016preVFP', 'UL2016postVFP', 'UL2017', 'UL2018']

    # Add the --eras argument with choices and default value
    parser.add_argument(
        '-e', '--eras', 
        choices=allowed_eras, 
        nargs='*', 
        default=allowed_eras,
        help="Specify one or more eras. Allowed values are: 'UL2016preVFP', 'UL2016postVFP', 'UL2017', 'UL2018'. If not provided, all eras will be used by default."
    )

    parser.add_argument(
        '-c', '--channels',
        nargs='+',  # Accept one or more values
        type=str,
        default=[],  # Default to an empty list if no arguments are provided
        help="Specify one or more channels to process. By default, all channels"
    )


    # Add the --sample flag argument
    parser.add_argument(
        '-s', '--sample',
        action='store_true',
        help="If provided, the sample mode will be enabled."
    )

    parser.add_argument(
        '--onlyData',
        action='store_true',
        help="If provided, only the data will be processed."
    )

    parser.add_argument(
        '--onlyMC',
        action='store_true',
        help="If provided, only the MC will be processed."
    )

    parser.add_argument(
        '-t','--timestamp',
        type=str,
        default='timestamp',
        help="Specify the timestamp'."
    )
    args = parser.parse_args()

    # Display the parsed arguments
    print(f"Selected eras: {args.eras}")
    print(f"Sample mode: {args.sample}")
    print(f"Output file: Output_{args.timestamp}.coffea")
    if len(args.channels) > 0:
        print(f"Channels: {args.channels}")
    else:
        print("Channels: All")

    outputDir = "outputs"
    datasetFlag = 'data'

    if args.sample:
        datasetFlag = 'sample'


    excludeChannels = ['WWTolnulnu']

    fileset = {}
    for era in args.eras:
        with open(f'../../Datasets/{datasetFlag}Files_{era}.json', 'r') as json_file:
            dicti = json.load(json_file)
            for pr in dicti['Data_mu']:
                if args.onlyMC:
                    print("Working on only the MC, ignoring data")
                    continue
                if len(args.channels) > 0:
                    if pr not in args.channels:
                        # print(f'skipping {era}_{pr} as not in list')
                        continue
                datasetName = f'{era}_{pr}'
                if dicti['Data_mu'][pr] == {}:
                    print(f"Skipping {datasetName} as no files found")
                    continue
                fileset[datasetName] = {"files": dicti['Data_mu'][pr]}
            for pr in dicti['MC_mu']:
                if args.onlyData:
                    print("Working on only the data, ignoring MC")
                    continue
                if len(args.channels) > 0:
                    if pr not in args.channels:
                        # print(f'skipping {era}_{pr} as not in list')
                        continue
                if pr in excludeChannels:
                    print(f"Skipping {era}_{pr} as it is in the exclude list")
                    continue
                datasetName = f'{era}_{pr}'
                if dicti['MC_mu'][pr] == {}:
                    print(f"Skipping {datasetName} as no files found")
                    continue
                fileset[datasetName] = {"files": dicti['MC_mu'][pr]}

    dataset_runnable, dataset_updated = preprocess(
        fileset,
        align_clusters=False,
        files_per_batch=1,
        skip_bad_files=True,
        save_form=False,
    )

    to_compute = apply_to_fileset(
        MyProcessor(),
        max_chunks(dataset_runnable, 300),
        schemaclass=NanoAODSchema,
    )

    (out,) = dask.compute(to_compute, scheduler='threads')
    
    outputFile = f"Output_{args.timestamp}.coffea"
    save(out, os.path.join(outputDir, outputFile))
    print(f"Output file is stored in {os.path.join(outputDir, outputFile)}")

if __name__ == '__main__':
    from multiprocessing import freeze_support
    freeze_support()
    main()