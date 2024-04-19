import os
import dask
import numpy as np
import awkward as ak 
import hist.dask as hda
import dask_awkward as dak
import argparse
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

        yt2D = (
            hda.Hist.new
            .Bool(name="is_Yt_Higher")
            .Bool(name = "is_x1_Higher")
            .Bool(name = "is_uubar")
            .Bool(name = "is_ddbar")
            .Bool(name = "is_ubaru")
            .Bool(name = "is_dbard")
            .Reg(250, 0, 2.5, label="$y_t$", name = "y_t")
            .Reg(250, 0, 2.5, label="$y_tbar$", name = "y_tbar")
            .Double()
        )
        Yt_cut = yt > ytbar
        x_cut = events.Generator_x1 > events.Generator_x2
        uubar_cut = (events.Generator_id1 == 2) & (events.Generator_id2 == -2)
        ubaru_cut = (events.Generator_id1 == -2) & (events.Generator_id2 == 2)
        ddbar_cut = (events.Generator_id1 == 1) & (events.Generator_id2 == -1)
        dbard_cut = (events.Generator_id1 == -1) & (events.Generator_id2 == 1)

        yt2D.fill(
            is_Yt_Higher = Yt_cut,
            is_x1_Higher = x_cut,
            is_uubar = uubar_cut,
            is_ddbar = ddbar_cut,
            is_ubaru = ubaru_cut,
            is_dbard = dbard_cut,
            y_t = yt,
            y_tbar = ytbar
        )

        return {
                "entries": ak.num(events, axis=0),
                "yMatrix": yt2D,
            }
        rich_bar()
    def postprocess(self, accumulator):
        pass

def main():
    outputDir = "outputs"
    to_analyze = 'fullRun2'
    # to_analyze = 'ttbarSample_UL2016preVFP'


    if to_analyze=='fullRun2':
        fileset = {}
        for era in ['UL2016preVFP']:
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
    outputFile = "Output_new.coffea"
    save(out, os.path.join(outputDir, outputFile))
    pass

if __name__ == '__main__':
    # On Windows, this prevents the fork bomb issue when starting new processes.
    from multiprocessing import freeze_support
    freeze_support()
    main()
