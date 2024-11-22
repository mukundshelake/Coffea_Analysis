from coffea.util import load, save
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import awkward as ak
from coffea.nanoevents.methods import vector
from scipy.special import legendre
import os
from coffea.util import load, save
from scipy.special import eval_legendre

timestamp = 'Nov17'

outputDir = f"outputs/{timestamp}"
era = 'UL2016preVFP'
era_channel = 'UL2016preVFP_ttbar_SemiLeptonic'
for i in range(2):
    coffeaFile = f"BDTSkimmerOutput_{timestamp}_chunk{i}.coffea"
    parquetFile = f"events_{era}_chunk{i}.parquet"
    out = load(os.path.join(outputDir,coffeaFile))
    params = out['UL2016preVFP_ttbar_SemiLeptonic']['params']
    params_ak = ak.from_iter(params)
    params_ak = params_ak[:100]
    parquetFile = f"events_{era}_chunk{i}_sample.parquet"
    sum_p1_p2 = params_ak['p1_id'] + params_ak['p2_id']
    condition1 = (sum_p1_p2 == 0)
    condition2 = (sum_p1_p2 == 42)
    new_array = ak.Array({'y': ak.where(condition1, 0, ak.where(condition2, 1, -1))})


    new_array['ttbarpz'] = params_ak['ttbarpz']
    jet_HT = ak.sum(params_ak['jet_details']['pt'], axis=1)
    nJets = ak.num(params_ak['jet_details']['pt'])
    muon_HT = ak.sum(params_ak['muon_details']['pt'], axis=1)
    pT_sum = muon_HT + params_ak['MET'] 
    new_array['Jet_HT'] = jet_HT
    new_array['nJet'] = nJets
    new_array['pT_sum'] = pT_sum

    jet_pt = params_ak['jet_details']['pt']
    jet_eta = params_ak['jet_details']['eta']
    jet_phi = params_ak['jet_details']['phi']

    jet_px = jet_pt * np.cos(jet_phi)
    jet_py = jet_pt * np.sin(jet_phi)
    jet_pz = jet_pt * np.sinh(jet_eta)
    jet_p = np.sqrt(jet_px**2 + jet_py**2 + jet_pz**2)

    # Calculate the dot product of jet momenta
    dot_product = jet_px * jet_px[:, np.newaxis] + jet_py * jet_py[:, np.newaxis] + jet_pz * jet_pz[:, np.newaxis]

    # Calculate the magnitude of jet momenta
    jet_p_magnitude = ak.sum(jet_p, axis=1)

    # Calculate the cosine of the angle between jet momenta
    cos_theta = dot_product / (jet_p * jet_p[:, np.newaxis])

    for l in [1]:
        Pl_cos_theta = eval_legendre(l, cos_theta)
        hl = ak.sum(jet_p * jet_p[:, np.newaxis] * Pl_cos_theta, axis=2)
        Hl = ak.sum(hl, axis=1)/jet_p_magnitude**2
        new_array[f'FW{l}'] = Hl

    # // Calculate sphericity matrix elements
    Sxx = ak.sum(jet_px * jet_px, axis=1)/ak.sum(jet_p * jet_p, axis=1)
    Syy = ak.sum(jet_py * jet_py, axis=1)/ak.sum(jet_p * jet_p, axis=1)
    Szz = ak.sum(jet_pz * jet_pz, axis=1)/ak.sum(jet_p * jet_p, axis=1)
    Sxy = ak.sum(jet_px * jet_py, axis=1)/ak.sum(jet_p * jet_p, axis=1)
    Sxz = ak.sum(jet_px * jet_pz, axis=1)/ak.sum(jet_p * jet_p, axis=1)
    Syz = ak.sum(jet_py * jet_pz, axis=1)/ak.sum(jet_p * jet_p, axis=1)

    new_array['Sxx'] = Sxx       


    new_array['Syy'] = Syy
    new_array['Szz'] = Szz
    new_array['Sxy'] = Sxy
    new_array['Sxz'] = Sxz
    new_array['Syz'] = Syz

    # Construct the sphericity tensor
    sphericity_tensor = np.array([[Sxx, Sxy, Sxz],
                                [Sxy, Syy, Syz],
                                [Sxz, Syz, Szz]])
    
    # Initialize an array to store the eigenvalues
    eigenvalues = np.zeros((3, sphericity_tensor.shape[2]))

    # Iterate over the third dimension and compute eigenvalues for each 3x3 matrix
    for i in range(sphericity_tensor.shape[2]):
        eigenvalues[:, i] = np.sort(np.linalg.eigvals(sphericity_tensor[:, :, i]))[::-1]

    lambda1, lambda2, lambda3 = eigenvalues[0], eigenvalues[1], eigenvalues[2]

    sphericity = 1.5 * (lambda2 + lambda3)
    planarity = lambda3/lambda2
    alignment = lambda2/lambda1
    p2in = lambda2/nJets
    p2out = lambda3/nJets

    new_array['S'] = sphericity
    new_array['P'] = planarity
    new_array['A'] = alignment
    new_array['p2in'] = p2in
    new_array['p2out'] = p2out

    ak.to_parquet(new_array, os.path.join(outputDir, parquetFile))
    print(f"Parquet file saved to {os.path.join(outputDir, parquetFile)}")
