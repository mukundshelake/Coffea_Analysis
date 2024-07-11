import os
import dask
import numpy as np
import awkward as ak
import hist.dask as hda
from coffea import processor
from coffea.nanoevents import NanoEventsFactory, BaseSchema
from coffea.dataset_tools import apply_to_fileset, max_chunks, preprocess
import matplotlib.pyplot as plt
import json, argparse
from coffea.util import save, rich_bar


class MyProcessor(processor.ProcessorABC):
    def __init__(self):
        pass

    def lorentz_transform(self, E, p, beta):
        beta_mag = dask.array.linalg.norm(beta)
        gamma = 1.0 / dask.array.sqrt(1 - beta_mag**2)
        bp = dask.array.dot(beta, p)
        E_prime = gamma * (E - bp)
        p_prime = p + ((gamma - 1) * bp / beta_mag**2 - gamma * E) * beta
        return E_prime, p_prime

    def process(self, events):
        dataset = events.metadata['dataset']
        
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

        deltay = yt - ytbar

        tbarpx = tbarpt*np.cos(tbarphi)
        tbarpy = tbarpt*np.sin(tbarphi)
        tbarpz = tbarpt*np.sinh(tbareta)
        tbarE = np.sqrt(tbarpt*tbarpt*np.cosh(tbareta)*np.cosh(tbareta) + tbarmass*tbarmass)

        # Four-momentum of the tbar system in lab frame
        ttbarpx = tpx + tbarpx
        ttbarpy = tpy + tbarpy
        ttbarpz = tpz + tbarpz
        ttbarE = tE + tbarE


        # Boost velocity of the tbar system
        # beta_ttbar_x = ttbarpx/ttbarE
        # beta_ttbar_y = ttbarpy/ttbarE
        # beta_ttbar_z = ttbarpz/ttbarE

        # beta = np.sqrt(beta_ttbar_x*beta_ttbar_x + beta_ttbar_y*beta_ttbar_y + beta_ttbar_z*beta_ttbar_z)

        # gamma = 1.0 / np.sqrt(1 - beta**2)

        # bp = beta_ttbar_x*tpx + beta_ttbar_y*tpy + beta_ttbar_z*tpz
        # E_prime = gamma * (tE - bp)
        # p_prime_x = tpx + ((gamma - 1) * bp / beta**2 - gamma * tE) * beta_ttbar_x
        # p_prime_y = tpy + ((gamma - 1) * bp / beta**2 - gamma * tE) * beta_ttbar_y
        # p_prime_z = tpz + ((gamma - 1) * bp / beta**2 - gamma * tE) * beta_ttbar_z

        # Calculate the angle between the top quark and the z-axis in the tbar rest frame
        # cos_theta = p_prime_z/ np.sqrt(p_prime_x*p_prime_x + p_prime_y*p_prime_y + p_prime_z*p_prime_z)

        mtt = np.sqrt(ttbarE*ttbarE - (ttbarpx*ttbarpx + ttbarpy*ttbarpy + ttbarpz*ttbarpz))

        betattz = abs(ttbarpz)/ttbarE

        yt2D = (
            hda.Hist.new
            # .Reg(20, -1.0, 1.0, label = "$c*$", name = "c")
            .Bool(name="is_Yt_Higher")
            .Reg(20, 250, 1250, label = "$m_tt$", name = "m_tt")
            .Reg(17, 0, 1.02, label = "$beta_ttz$", name = "beta_ttz")
            .Reg(8, 0, 2.4, label="$y_t$", name = "y_t")
            .Reg(8, 0, 2.4, label="$y_tbar$", name = "y_tbar")
            .Reg(2, -2.4, 2.4, label="$deltay$", name = "deltay")
            .Double()
        )
        Yt_cut = yt > ytbar
        yt2D.fill(
            # c = cos_theta,
            is_Yt_Higher = Yt_cut,
            m_tt = mtt,
            beta_ttz = betattz,
            y_t = yt,
            y_tbar = ytbar,
            deltay = deltay
        )

        return {
                "entries": ak.num(events, axis=0),
                "yMatrix": yt2D,
            }

    def postprocess(self, accumulator):
        pass

def main():
    parser = argparse.ArgumentParser(description='Run processor on sample')
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
        schemaclass=BaseSchema,
    )

    (out,) = dask.compute(to_compute, scheduler='threads')
    
    outputFile = "skimmerOutput.coffea"
    save(out, os.path.join(outputDir, outputFile))
    print(f"Output file is stored in {os.path.join(outputDir, outputFile)}")

if __name__ == '__main__':
    from multiprocessing import freeze_support
    freeze_support()
    main()














# import os
# import dask
# import numpy as np
# import awkward as ak 
# import hist.dask as hda
# from coffea import processor
# # import warnings
# # warnings.filterwarnings("error", module="coffea.*")
# from coffea.nanoevents.methods import candidate
# from coffea.dataset_tools import (
#     apply_to_fileset,
#     max_chunks,
#     preprocess,
# )
# # from distributed import Client
# from coffea.nanoevents import NanoEventsFactory, BaseSchema
# import matplotlib.pyplot as plt
# from coffea.util import save, rich_bar
# import json, argparse

# class MyProcessor(processor.ProcessorABC):
#     def __init__(self):
#         pass

#     def lorentz_transform(self, E, p, beta):
#         beta_mag = np.linalg.norm(beta)
#         gamma = 1.0 / np.sqrt(1 - beta_mag**2)
#         bp = np.dot(beta, p)
#         E_prime = gamma * (E - bp)
#         p_prime = p + ((gamma - 1) * bp / beta_mag**2 - gamma * E) * beta
#         return E_prime, p_prime

#     def four_momentum(self, pt, eta, phi, mass):
#         px = pt * np.cos(phi)
#         py = pt * np.sin(phi)
#         pz = pt * np.sinh(eta)
#         E = np.sqrt(px**2 + py**2 + pz**2 + mass**2)
#         return E, np.array([px, py, pz])

#     def process(self, events):
#         dataset = events.metadata['dataset']
#         # print(dataset)
#         tops = ak.zip(
#             {
#                 "pt" : events.GenPart_pt[:, 2],
#                 "eta": events.GenPart_eta[:, 2],
#                 "mass": events.GenPart_mass[:, 2],
#                 "phi" : events.GenPart_phi[:, 2]
#             }
#         )
#         antitops = ak.zip(
#             {
#                 "pt" : events.GenPart_pt[:, 3],
#                 "eta": events.GenPart_eta[:, 3],
#                 "mass": events.GenPart_mass[:, 3],
#                 "phi" : events.GenPart_phi[:, 3]
#             }
#         )
#         tpt = tops["pt"]
#         teta = tops["eta"]
#         tmass = tops["mass"]
#         tphi = tops["phi"]
#         tbarpt = antitops["pt"]
#         tbareta = antitops["eta"]
#         tbarmass = antitops["mass"]
#         tbarphi = antitops["phi"]

#         E_t, p_t = self.four_momentum(tpt, teta, tphi, tmass)
#         E_tbar, p_tbar = self.four_momentum(tbarpt, tbareta, tbarphi, tbarmass)

#         # Four-momentum of the tbar system in lab frame
#         E_ttbar = E_t + E_tbar
#         p_ttbar = p_t + p_tbar


#         # Boost velocity of the tbar system
#         beta_ttbar = p_ttbar / E_ttbar

#         # Boost the top quark to the tbar rest frame
#         E_t_prime, p_t_prime = self.lorentz_transform(E_t, p_t, beta_ttbar)

#         # Calculate the angle between the top quark and the z-axis in the tbar rest frame
#         cos_theta = p_t_prime[2] / np.linalg.norm(p_t_prime)

#         yt2D = (
#             hda.Hist.new
#             .Reg(20, -1.0, 1.0, label = "$c*$", name = "c")
#             .Double()
#         )

#         yt2D.fill(
#             c = cos_theta
#         )

#         return {
#                 "entries": ak.num(events, axis=0),
#                 "yMatrix": yt2D,
#             }
#         rich_bar()
#     def postprocess(self, accumulator):
#         pass


# def main():
#     parser = argparse.ArgumentParser(description='Do you want to run it on sample')
#     parser.add_argument('-s', '--sample', action='store_true', help='run on sample')

#     args = parser.parse_args()


#     outputDir = "outputs"
#     to_analyze = 'fullRun2'

#     if args.sample:
#         to_analyze = 'ttbarSample_UL2016preVFP'

#     print(f"\n\nWorking on {to_analyze}")


#     if to_analyze=='fullRun2':
#         fileset = {}
#         for era in ['UL2016preVFP', 'UL2016postVFP' , 'UL2017', 'UL2018']:
#             with open(f'../Datasets/dataFiles_{era}.json', 'r') as json_file:
#                 fileset[era]= {'files':json.load(json_file)['MC_el']['ttbar_SemiLeptonic']}
#     elif to_analyze=='ttbarSample_UL2016preVFP':
#         fileset = {
#             'ttbarSample_UL2016preVFP': {
#                 "files": {
#                     'file://../../tests/UL2016_preVFP_ttbarSemileptonic.root': "Events",
#                 }
#             }
#         }
#     else:
#         fileset = {
#             'ttbarSample_UL2016preVFP': {
#                 "files": {
#                     'file://../../tests/UL2016_preVFP_ttbarSemileptonic.root': "Events",
#                 }
#             }
#         }

#     # Your code that starts new processes goes here.
#     dataset_runnable, dataset_updated = preprocess(
#         fileset,
#         align_clusters=False,
#         # maybe_step_size=100_000,
#         files_per_batch=1,
#         skip_bad_files=True,
#         save_form=False,)

#     to_compute = apply_to_fileset(
#                     MyProcessor(),
#                     max_chunks(dataset_runnable, 300),
#                     schemaclass=BaseSchema,
#                 )
#     (out,) = dask.compute(to_compute, scheduler='threads')
#     # print(out)
#     outputFile = "skimmerOutput.coffea"
#     save(out, os.path.join(outputDir, outputFile))
#     print(f"output file is stored in {os.path.join(outputDir, outputFile)}")
#     pass


# if __name__ == '__main__':
#     # On Windows, this prevents the fork bomb issue when starting new processes.
#     from multiprocessing import freeze_support
#     freeze_support()
#     main()
