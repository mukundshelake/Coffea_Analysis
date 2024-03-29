from coffea.util import load
import matplotlib.pyplot as plt
import mplhep as hep
import hist
import os
### Parameters for Lumi factor calculations
#### Obtaned from "https://github.com/mintutifr/nanoAOD-tools/blob/Mintu_dev/crab/minitree/Xsec_Nevent.py"
UL2016preVFP = {
		'Tchannel'      : ['52437432','134.2'],
		'Tbarchannel'   : ['29205915','80.0'],
		'Schannel'      : ['3592772','2.215836'],
		'tw_top'        : ['3294485','39.65'],
		'tw_antitop'    : ['3176335','39.65'],
		'ttbar_SemiLeptonic'  : ['131106831','366.3'],
		'ttbar_FullyLeptonic' : ['37202073','88.5'],
		'WJetsToLNu_0J' : ['121208493','52780.0'],
		'WJetsToLNu_1J' : ['84198168','8832.0'],
		'WJetsToLNu_2J' : ['27463756','3276.0'],
		'DYJetsToLL'    : ['95170552','6424.0'],
		'WWTo2L2Nu'     : ['3006596','11.09'],
		'WZTo2Q2L'      : ['9780392','6.565'],
		'ZZTo2Q2L'      : ['10406942','3.676'],
		'QCD_Pt-15To20_MuEnriched'   : ['4749558','2881500'],
		'QCD_Pt-20To30_MuEnriched'   : ['31926643','3657550'],
		'QCD_Pt-30To50_MuEnriched'   : ['28664840','1154520'],
		'QCD_Pt-50To80_MuEnriched'   : ['19724658','347191'],
		'QCD_Pt-80To120_MuEnriched'  : ['21978558','90131'],
		'QCD_Pt-120To170_MuEnriched' : ['19136907','21509.488'],
		'QCD_Pt-170To300_MuEnriched' : ['37080544','7224.3'],
		'QCD_Pt-300To470_MuEnriched' : ['28465385','617.1102'],
		'QCD_Pt-470To600_MuEnriched' : ['19695298','57.84961'],
		'QCD_Pt-600To800_MuEnriched' : ['19598815','22.21919375'],
		'QCD_Pt-800To1000_MuEnriched': ['38933940','3.86242284'],
		'QCD_Pt-1000_MuEnriched'     : ['13660999','1.34449145'],
		'QCD_Pt-30to50_EMEnriched'   : ['4361931','6485600'],
		'QCD_Pt-50to80_EMEnriched'   : ['5440758','1991360'],
		'QCD_Pt-80to120_EMEnriched'  : ['4847354','380538'],
		'QCD_Pt-120to170_EMEnriched' : ['4852573','72499.4'],
		'QCD_Pt-170to300_EMEnriched' : ['35640892','21349.8'],
		'QCD_Pt-300toInf_EMEnriched' : ['1142775','1350'],

		'Tchannel_TuneCP5CR2' : ['20230062','134.2'],
		'Tchannel_TuneCP5CR1' : ['18671220','134.2'],
		'Tchannel_TuneCP5up' : ['20545577','134.2'],
		'Tchannel_TuneCP5down' : ['20681075','134.2'],
		#'Tchannel_hdampup' : ['20779593','134.2'],
		#'Tchannel_hdampdown' : ['20784203','134.2'],
		'Tchannel_erdON' : ['21487232','134.2'],
		'Tbarchannel_TuneCP5CR2' : ['10300876','80.0'],
		'Tbarchannel_TuneCP5CR1' : ['9739578','80.0'],
		'Tbarchannel_TuneCP5up' : ['10086078','80.0'],
		'Tbarchannel_TuneCP5down' : ['10216908','80.0'],
		#'Tbarchannel_hdampup' : ['9947894','80.0'],
		#'Tbarchannel_hdampdown' : ['10124656','80.0'],
		'Tbarchannel_erdON' : ['10410358','80.0']
}


### Lumi info and generating the scale factor for each dataset mentioned above. 
###### Note that: we may not be using all the factors.
intLumi = 19521.2283 
factor = {}
# print(f"The lumi scale factors for the datasets")
for key in UL2016preVFP:
    nMC = float(UL2016preVFP[key][0])
    XMC = float(UL2016preVFP[key][1])
    x = intLumi*XMC/nMC
    factor[key] = x
    # print(f"{key}:{x} by using {XMC} and {nMC}")

### load .coffea file generated by the processor
output = load("../coffeaOutputs/ElHLT32eta2p1_pt35_20230926_105202_SIXTEEN_preVFP_el.coffea")
plt.style.use(hep.style.ROOT)



###### Bunch the datasets into proper buckets
eventtypes = {
    'Data':{
        'UL2016preVFP_el':{'Datasets':[s for s in output.keys() if ("Run2016" in s) and ("el" in s) and ("ver1" not in s)]}
    },
    'MC':{
    'singleT':{'Datasets':["tw_top", "tw_antitop", "Tchannel", "Tbarchannel", "Schannel"]},
    'QCD_EMEnriched':{'Datasets':[s for s in output.keys() if ("QCD" in s) and ("EM" in s)]},
    'QCD_MuEnriched':{'Datasets':[s for s in output.keys() if ("QCD" in s) and ("Mu" in s)]},
    'WJets':{'Datasets':[s for s in output.keys() if ("WJet" in s)]},
    'Diboson':{'Datasets':[s for s in output.keys() if ("WW" in s) or ("ZZ" in s) or ("WZ" in s)]},
    'DYJets':{'Datasets':["DYJetsToLL"]},
    'TTbar_semileptonic':{'Datasets':["ttbar_SemiLeptonic"]},
    'TTbar_fullyleptonic':{'Datasets':["ttbar_FullyLeptonic"]},
    }
}


### List the features you want to add to the eventtypes dictionary
histFeatures = ['electron_eta', 'electron_pt']



###### Following loop populates the eventtypes dictionary
for DataMC in eventtypes:
#     print(key)
    for key in eventtypes[DataMC]:
        for feature in histFeatures:
            noEvents = 0
            ## To get the empty histogram of same format as the output[first_dataset in the eventtypes-->Data/MC --> key][feature]
            ## This will make sure that the histograms can be properly added.
            hist_ = output[eventtypes[DataMC][key]['Datasets'][0]][feature]*1.0 
            hist_.reset()

            for s in eventtypes[DataMC][key]['Datasets']:
                if s == 'ZZTo2L2Nu' or s == 'WWTolnulnu':
                    print(f"Ignoring the {s} dataset as we don't have the numbers for them")
                    continue
                if DataMC == 'Data':
                    scale = 1.0
                else:
                    scale = factor[s]
                    #scale = 1.0
                noEvents += sum(output[s][feature].values())
                hist_ = hist_ + output[s][feature]*scale
                # print(output[s][feature]*scale)
            #print(hist_)
            # print(f"{key}: {noEvents}")
            eventtypes[DataMC][key][feature] = hist_
            eventtypes[DataMC][key]['noEvents'] = noEvents

            
#### To print the populated eventtypes dictionary
print(eventtypes)
################################################################
## THE PLOTTER

### for plot naming
era = 'UL2016preVFP'
lep = 'el'
feature = 'electron_pt'
plotDir = '../plots'

### Fig1: Only Data plot 
fig1, ax = plt.subplots()
hep.histplot(eventtypes['Data']['UL2016preVFP_el'][feature] , ax=ax)
fig1_name = era + '_' + lep + '_' + feature + '.pdf'
fig1_file = os.path.join(plotDir, fig1_name)
plt.savefig(fig1_file)


### Fig2: stack the MC and overlay data with errorbar
fig2, ax2 = plt.subplots()
hep.histplot(
    [eventtypes['MC'][s][feature] for s in eventtypes['MC'].keys()],
    histtype="fill",
    stack=True,
    label=[s for s in eventtypes['MC'].keys()],
    ax=ax2,
)

hep.histplot(
#     [output[s][feature] for s in output.keys() if "Run" in s],
    eventtypes['Data']['UL2016preVFP_el'][feature],
    histtype="errorbar",
    label='Data',
    ax=ax2,
    color="k",
)
plt.legend()
fig2_name = era + '_' + lep + '_' + feature + '_stacked'+'.pdf'
fig2_file = os.path.join(plotDir, fig2_name)
plt.savefig(fig2_file)

#################################################
### Counting the total Data and MC events
DataEvents = sum(eventtypes['Data'][s]['noEvents'] for s in eventtypes['Data'])
MCEvents = sum(eventtypes['MC'][s]['noEvents'] for s in eventtypes['MC'])
print(f"The Data/MC ratio: {DataEvents/MCEvents}")

#################################################
