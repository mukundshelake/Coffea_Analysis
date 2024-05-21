import os
import dask
import numpy as np
import awkward as ak 
import hist.dask as hda
from coffea import processor
# import warnings
# warnings.filterwarnings("error", module="coffea.*")
from coffea.nanoevents.methods import candidate
from coffea.dataset_tools import (
    apply_to_fileset,
    max_chunks,
    preprocess,
)
# from distributed import Client
from coffea.nanoevents import NanoEventsFactory, BaseSchema
import matplotlib.pyplot as plt
from coffea.util import save, rich_bar
import json, argparse

class MyProcessor(processor.ProcessorABC):
    def __init__(self):
        pass

    def process(self, events):
        dataset = events.metadata['dataset']
        # print(dataset)
        tops = ak.zip(
            {
                "pt" : events.GenPart_pt[:, 2],
                "eta": events.GenPart_eta[:, 2],
                "mass": events.GenPart_mass[:, 2],
                "phi" : events.GenPart_phi[:, 2]
            }
        )
        antitops = ak.zip(
            {
                "pt" : events.GenPart_pt[:, 3],
                "eta": events.GenPart_eta[:, 3],
                "mass": events.GenPart_mass[:, 3],
                "phi" : events.GenPart_phi[:, 3]
            }
        )
        tpt = tops["pt"]
        teta = tops["eta"]
        tmass = tops["mass"]
        tphi = tops["phi"]
        yt = abs(teta - 0.50*np.tanh(teta)*np.square(tmass/tpt))
        tbarpt = antitops["pt"]
        tbareta = antitops["eta"]
        tbarmass = antitops["mass"]
        tbarphi = antitops["phi"]
        ytbar = abs(tbareta - 0.50*np.tanh(tbareta)*np.square(tbarmass/tbarpt))


        tpx = tpt*np.cos(tphi)
        tpy = tpt*np.sin(tphi)
        tpz = tpt*np.sinh(teta)
        tE = np.sqrt(tpt*tpt*np.cosh(teta)*np.cosh(teta) + tmass*tmass)

        tbarpx = tbarpt*np.cos(tbarphi)
        tbarpy = tbarpt*np.sin(tbarphi)
        tbarpz = tbarpt*np.sinh(tbareta)
        tbarE = np.sqrt(tbarpt*tbarpt*np.cosh(tbareta)*np.cosh(tbareta) + tbarmass*tbarmass)

        ttbarpx = tpx + tbarpx
        ttbarpy = tpy + tbarpy
        ttbarpz = tpz + tbarpz
        ttbarE = tE + tbarE

        mtt = np.sqrt(ttbarE*ttbarE - (ttbarpx*ttbarpx + ttbarpy*ttbarpy + ttbarpz*ttbarpz))

        betatt = abs(ttbarpz)/ttbarE

        yt2D = (
            hda.Hist.new
            .Bool(name="is_Yt_Higher")
            .Reg(20, 0, 900, label = "$m_tt$", name = "m_tt")
            .Reg(25, 0, 1.5, label = "$beta_tt$", name = "beta_tt")
            .Double()
        )
        Yt_cut = yt > ytbar

        yt2D.fill(
            is_Yt_Higher = Yt_cut,
            m_tt = mtt,
            beta_tt = betatt
        )

        return {
                "entries": ak.num(events, axis=0),
                "yMatrix": yt2D,
            }
        rich_bar()
    def postprocess(self, accumulator):
        pass


def main():
    parser = argparse.ArgumentParser(description='Do you want to run it on sample')
    parser.add_argument('-s', '--sample', action='store_true', help='run on sample')

    args = parser.parse_args()


    outputDir = "outputs"
    to_analyze = 'fullRun2'

    if args.sample:
        to_analyze = 'ttbarSample_UL2016preVFP'

    print(f"\n\nWorking on {to_analyze}")


    if to_analyze=='fullRun2':
        fileset = {}
        for era in ['UL2016preVFP', 'UL2016postVFP' , 'UL2017', 'UL2018']:
            with open(f'../Datasets/dataFiles_{era}.json', 'r') as json_file:
                fileset[era]= {'files':json.load(json_file)['MC_el']['ttbar_SemiLeptonic']}
    elif to_analyze=='ttbarSample_UL2016preVFP':
        fileset = {
            'ttbarSample_UL2016preVFP': {
                "files": {
                    'file://../../tests/UL2016_preVFP_ttbarSemileptonic.root': "Events",
                }
            }
        }
    else:
        fileset = {
            'ttbarSample_UL2016preVFP': {
                "files": {
                    'file://../../tests/UL2016_preVFP_ttbarSemileptonic.root': "Events",
                }
            }
        }

    # Your code that starts new processes goes here.
    dataset_runnable, dataset_updated = preprocess(
        fileset,
        align_clusters=False,
        # maybe_step_size=100_000,
        files_per_batch=1,
        skip_bad_files=True,
        save_form=False,)

    to_compute = apply_to_fileset(
                    MyProcessor(),
                    max_chunks(dataset_runnable, 300),
                    schemaclass=BaseSchema,
                )
    (out,) = dask.compute(to_compute, scheduler='threads')
    # print(out)
    outputFile = "skimmerOutput.coffea"
    save(out, os.path.join(outputDir, outputFile))
    print(f"output file is stored in {os.path.join(outputDir, outputFile)}")
    pass


if __name__ == '__main__':
    # On Windows, this prevents the fork bomb issue when starting new processes.
    from multiprocessing import freeze_support
    freeze_support()
    main()
