from coffea.util import load, save
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import awkward as ak
from coffea.nanoevents.methods import vector
from scipy.special import legendre

timestamp = 'Nov17'

outputDir = f'outputs/{timestamp}'
coffeaFile = f"BDTSkimmerOutput_{timestamp}.coffea"


# out = load("Output.coffea")
out = load(os.path.join(outputDir,coffeaFile))


def calculate_angle(px1, py1, pz1, px2, py2, pz2):
    vector1 = ak.zip({"x": px1, "y": py1, "z": pz1}, with_name="Vector3D")
    vector2 = ak.zip({"x": px2, "y": py2, "z": pz2}, with_name="Vector3D")
    
    cos_theta = (vector1.x * vector2.x + 
                 vector1.y * vector2.y + 
                 vector1.z * vector2.z) / (
                 np.sqrt(vector1.x**2 + vector1.y**2 + vector1.z**2) * 
                 np.sqrt(vector2.x**2 + vector2.y**2 + vector2.z**2))
    
    angle = np.arccos(cos_theta)
    return angle


for key in out:
    era = key.split('_')[0]
    outputFile = f'events_{era}.parquet'
    awkward_array = ak.Array(out[key]['params'])
    print(out[key]['nSelected'])
    sum_p1_p2 = awkward_array['p1_id'] + awkward_array['p2_id']
    condition1 = (sum_p1_p2 == 0)
    condition2 = (sum_p1_p2 == 42)
    awkward_array['y'] = ak.where(condition1, 0, ak.where(condition2, 1, -1))
    leadingMuon_px = awkward_array['leading_muon_px']
    leadingMuon_py = awkward_array['leading_muon_py']
    leadingMuon_pz = awkward_array['leading_muon_pz']
    leadingMuon_pt = np.sqrt(leadingMuon_px**2 + leadingMuon_py**2)
    awkward_array['leading_muon_pt'] = leadingMuon_pt
    leadingMuon_p = np.sqrt(leadingMuon_px**2 + leadingMuon_py**2 + leadingMuon_pz**2)



    leadingJet_px = awkward_array['leadingJet_px']
    leadingJet_py = awkward_array['leadingJet_py']
    leadingJet_pz = awkward_array['leadingJet_pz']
    leadingJet_pt = np.sqrt(leadingJet_px**2 + leadingJet_py**2)
    awkward_array['leadingJet_pt'] = leadingJet_pt
    leadingJet_p = np.sqrt(leadingJet_px**2 + leadingJet_py**2 + leadingJet_pz**2)

    subleadingJet_px = awkward_array['subleadingJet_px']
    subleadingJet_py = awkward_array['subleadingJet_py']
    subleadingJet_pz = awkward_array['subleadingJet_pz']
    subleadingJet_pt = np.sqrt(subleadingJet_px**2 + subleadingJet_py**2)
    awkward_array['subleadingJet_pt'] = subleadingJet_pt
    subleadingJet_p = np.sqrt(subleadingJet_px**2 + subleadingJet_py**2 + subleadingJet_pz**2)

    subsubleadingJet_px = awkward_array['subsubleadingJet_px']
    subsubleadingJet_py = awkward_array['subsubleadingJet_py']
    subsubleadingJet_pz = awkward_array['subsubleadingJet_pz']
    subsubleadingJet_pt = np.sqrt(subsubleadingJet_px**2 + subsubleadingJet_py**2)
    awkward_array['subsubleadingJet_pt'] = subsubleadingJet_pt
    subsubleadingJet_p = np.sqrt(subsubleadingJet_px**2 + subsubleadingJet_py**2 + subsubleadingJet_pz**2)
    # pTprod_lj1 = leadingMuon_pt * leadingJet_pt
    # pTprod_lj2 = leadingMuon_pt * subleadingJet_pt
    # pTprod_lj3 = leadingMuon_pt * subsubleadingJet_pt
    pTprod_j1j2 = leadingJet_pt * subleadingJet_pt
    pTprod_j1j3 = leadingJet_pt * subsubleadingJet_pt
    pTprod_j2j3 = subleadingJet_pt * subsubleadingJet_pt

    # angle_lj1 = calculate_angle(leadingMuon_px, leadingMuon_py, leadingMuon_pz, leadingJet_px, leadingJet_py, leadingJet_pz)
    # angle_lj2 = calculate_angle(leadingMuon_px, leadingMuon_py, leadingMuon_pz, subleadingJet_px, subleadingJet_py, subleadingJet_pz)
    # angle_lj3 = calculate_angle(leadingMuon_px, leadingMuon_py, leadingMuon_pz, subsubleadingJet_px, subsubleadingJet_py, subsubleadingJet_pz)
    angle_j1j2 = calculate_angle(leadingJet_px, leadingJet_py, leadingJet_pz, subleadingJet_px, subleadingJet_py, subleadingJet_pz)
    angle_j2j3 = calculate_angle(subleadingJet_px, subleadingJet_py, subleadingJet_pz, subsubleadingJet_px, subsubleadingJet_py, subsubleadingJet_pz)
    angle_j1j3 = calculate_angle(leadingJet_px, leadingJet_py, leadingJet_pz, subsubleadingJet_px, subsubleadingJet_py, subsubleadingJet_pz)

    # FW3 = pTprod_lj1 * legendre(3)(np.cos(angle_lj1)) + pTprod_lj2 * legendre(3)(np.cos(angle_lj2)) + pTprod_lj3 * legendre(3)(np.cos(angle_lj3)) + pTprod_j1j2 * legendre(3)(np.cos(angle_j1j2)) + pTprod_j1j3 * legendre(3)(np.cos(angle_j1j3)) + pTprod_j2j3 * legendre(3)(np.cos(angle_j2j3))
    # FW1 = pTprod_lj1 * legendre(1)(np.cos(angle_lj1)) + pTprod_lj2 * legendre(1)(np.cos(angle_lj2)) + pTprod_lj3 * legendre(1)(np.cos(angle_lj3)) + pTprod_j1j2 * legendre(1)(np.cos(angle_j1j2)) + pTprod_j1j3 * legendre(1)(np.cos(angle_j1j3)) + pTprod_j2j3 * legendre(1)(np.cos(angle_j2j3))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
    # FW2 = pTprod_lj1 * legendre(2)(np.cos(angle_lj1)) + pTprod_lj2 * legendre(2)(np.cos(angle_lj2)) + pTprod_lj3 * legendre(2)(np.cos(angle_lj3)) + pTprod_j1j2 * legendre(2)(np.cos(angle_j1j2)) + pTprod_j1j3 * legendre(2)(np.cos(angle_j1j3)) + pTprod_j2j3 * legendre(2)(np.cos(angle_j2j3))

    # FW4 = pTprod_j1j2 * legendre(4)(np.cos(angle_j1j2)) + pTprod_j1j3 * legendre(4)(np.cos(angle_j1j3)) + pTprod_j2j3 * legendre(4)(np.cos(angle_j2j3))
    # FW3 = pTprod_j1j2 * legendre(3)(np.cos(angle_j1j2)) + pTprod_j1j3 * legendre(3)(np.cos(angle_j1j3)) + pTprod_j2j3 * legendre(3)(np.cos(angle_j2j3))
    FW1 = pTprod_j1j2 * legendre(1)(np.cos(angle_j1j2)) + pTprod_j1j3 * legendre(1)(np.cos(angle_j1j3)) + pTprod_j2j3 * legendre(1)(np.cos(angle_j2j3))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
    # FW2 = pTprod_j1j2 * legendre(2)(np.cos(angle_j1j2)) + pTprod_j1j3 * legendre(2)(np.cos(angle_j1j3)) + pTprod_j2j3 * legendre(2)(np.cos(angle_j2j3))
    FW0 = pTprod_j1j2 * legendre(0)(np.cos(angle_j1j2)) + pTprod_j1j3 * legendre(0)(np.cos(angle_j1j3)) + pTprod_j2j3 * legendre(0)(np.cos(angle_j2j3))

    # awkward_array['FW3'] = FW3/FW0
    awkward_array['FW1'] = FW1/FW0
    # awkward_array['FW2'] = FW2/FW0
    # awkward_array['FW4'] = FW4/FW0


    # Calculate elements of sphericity tensor
    total_p2 = leadingJet_p*leadingJet_p + subleadingJet_p*subleadingJet_p + subsubleadingJet_p*subsubleadingJet_p
    Sxx = (leadingJet_px**2 + subleadingJet_px**2 + subsubleadingJet_px**2)/total_p2
    Syy = (leadingJet_py**2 + subleadingJet_py**2 + subsubleadingJet_py**2)/total_p2
    Szz = (leadingJet_pz**2 + subleadingJet_pz**2 + subsubleadingJet_pz**2)/total_p2
    Sxy = (leadingJet_px * leadingJet_py + subleadingJet_px * subleadingJet_py + subsubleadingJet_px * subsubleadingJet_py)/total_p2
    Sxz = (leadingJet_px * leadingJet_pz + subleadingJet_px * subleadingJet_pz + subsubleadingJet_px * subsubleadingJet_pz)/total_p2
    Syz = (leadingJet_py * leadingJet_pz + subleadingJet_py * subleadingJet_pz + subsubleadingJet_py * subsubleadingJet_pz)/total_p2


    # Fill the sphericity tensor elements
    awkward_array['Sxx'] = Sxx
    awkward_array['Syy'] = Syy
    awkward_array['Szz'] = Szz
    awkward_array['Sxy'] = Sxy
    awkward_array['Sxz'] = Sxz
    awkward_array['Syz'] = Syz


    ak.to_parquet(awkward_array, os.path.join(outputDir,outputFile))
    print(f"Parquet file saved to {outputFile}")