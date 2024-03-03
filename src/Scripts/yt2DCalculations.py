import hist
import dask
import numpy as np
import awkward as ak
import hist.dask as hda
import dask_awkward as dak

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
import json


class MyProcessor(processor.ProcessorABC):
    def __init__(self):
        pass

    def process(self, events):
        dataset = events.metadata['dataset']
        print(dataset)
        tops = ak.zip(
            {
                "pt" : events.GenPart_pt[:, 2],
                "eta": events.GenPart_eta[:, 2],
                "mass": events.GenPart_mass[:, 2]
            }
        )
        antitops = ak.zip(
            {
                "pt" : events.GenPart_pt[:, 3],
                "eta": events.GenPart_eta[:, 3],
                "mass": events.GenPart_mass[:, 3]
            }
        )
        tpt = tops["pt"]
        teta = tops["eta"]
        tmass = tops["mass"]
        yt = abs(teta - 0.50*np.tanh(teta)*np.square(tmass/tpt))
        tbarpt = antitops["pt"]
        tbareta = antitops["eta"]
        tbarmass = antitops["mass"]
        ytbar = abs(tbareta - 0.50*np.tanh(tbareta)*np.square(tbarmass/tbarpt)) 

        yt2D_higherYt = (
            hda.Hist.new
            # .StrCat(["higher_yt", "higher_ytbar"], name="sign")
            .Reg(80, 0, 4., label="$y_t$", name = "y_t")
            .Reg(80, 0, 4., label="$y_tbar$", name = "y_tbar")
            .Double()
        )
        yt2D_higherYtbar = (
            hda.Hist.new
            # .StrCat(["higher_yt", "higher_ytbar"], name="sign")
            .Reg(80, 0, 4., label="$y_t$", name = "y_t")
            .Reg(80, 0, 4., label="$y_tbar$", name = "y_tbar")
            .Double()
        )
        cut = yt > ytbar
        yt2D_higherYt.fill(y_t = yt[cut], y_tbar = ytbar[cut])
        cut = ytbar > yt
        yt2D_higherYtbar.fill(y_t = yt[cut], y_tbar = ytbar[cut])
        return {
                "entries": ak.num(events, axis=0),
                "yMatrix_higherYt": yt2D_higherYt,
                "yMatrix_higherYtbar": yt2D_higherYtbar,
            }
        rich_bar()
    def postprocess(self, accumulator):
        pass

def main():
    # to_analyze = 'fullRun2'
    to_analyze = 'ttbarSample_UL2016preVFP'


    if to_analyze=='fullRun2':
        fileset = {}
        for era in ['UL2016preVFP', 'UL2016postVFP', 'UL2017', 'UL2018']:
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
    print(out)
    save(out, f"Output.coffea")
    pass

if __name__ == '__main__':
    # On Windows, this prevents the fork bomb issue when starting new processes.
    from multiprocessing import freeze_support
    freeze_support()
    main()
