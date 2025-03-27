import uproot
import awkward as ak
import numpy as np
from scipy.optimize import minimize
from concurrent.futures import ProcessPoolExecutor

def optimize_event(event_data):
    event_idx, permutations, resolutions = event_data
    best_chi2 = np.inf
    best_result = None

    for perm_idx, (perm, res) in enumerate(zip(permutations, resolutions)):
        fitted_momenta = perm[event_idx]
        measured_momenta = perm[event_idx]
        resolution = res[event_idx]

        if not np.all(resolution > 0):
            continue

        result = minimize(
            chi2,
            fitted_momenta,
            args=(measured_momenta, resolution),
            method="Powell",
            options={"disp": False},
        )

        if result.success and result.fun < best_chi2:
            best_chi2 = result.fun
            best_result = {
                "event": event_idx,
                "permutation": perm_idx + 1,
                "success": result.success,
                "optimized_momenta": result.x,
                "chi_squared": result.fun,
                "message": result.message,
            }

    return best_result


def solve_npz(mpx, mpy, mpz, mpe, metx, mety):
    # Calculate the neutrino pz assuming a W-boson mass constraint
    alpha = MW**2 - mpe**2 + mpx**2 + mpy**2 + mpz**2 + 2*(mpx*metx + mpy*mety)
    a = mpz**2 - mpe**2
    b = alpha*mpz 
    c = alpha**2/4 - mpe**2*(metx**2 + mety**2)

    # Calculate the discriminant
    discriminant = b**2 - 4*a*c
    # print(discriminant)
    # print(a)
    # print(b)
    # print(c)

    pz1 = ak.where(discriminant < 0, 0, (-b + np.sqrt(discriminant))/(2*a))
    pz2 = ak.where(discriminant < 0, 0, (-b - np.sqrt(discriminant))/(2*a))

    return pz1, pz2

def invariant_mass(candidate):
    mass2 = candidate[0]**2 - candidate[1]**2 - candidate[2]**2 - candidate[3]**2
    if mass2 < 0:
        # print(f"Warning: Negative mass squared encountered. Setting mass to 0. Candidate: {candidate}, Mass2: {mass2}")
        return 1000000  # Return 0 or some default handling
    return np.sqrt(mass2)

def chi2(fitted_momenta, measured_momenta, resolutions):
    chi = (fitted_momenta - measured_momenta) / resolutions
    chi_sq = np.sum(chi**2)

    w_had_candidate = fitted_momenta[12:16] + fitted_momenta[16:20]
    w_had_mass = invariant_mass(w_had_candidate)

    top_had_candidate = w_had_candidate + fitted_momenta[4:8] # First b is hadronic b
    top_had_mass = invariant_mass(top_had_candidate)

    top_lep_candidate = fitted_momenta[0:4] + fitted_momenta[20:24] + fitted_momenta[8:12] # Second b is leptonic b
    top_lep_mass = invariant_mass(top_lep_candidate)

    chi_sq += (w_had_mass - MW)**2 / delta_MW**2
    chi_sq += (top_had_mass - MT)**2 / delta_MT**2
    chi_sq += (top_lep_mass - MT)**2 / delta_MT**2

    # print(w_had_mass, top_had_mass, top_lep_mass, chi_sq)

    return chi_sq

filePath = '/nfs/home/mukund/Projects/Coffea_Analysis/src/Scripts/Reco/tree_14_Skim.root'
file = uproot.open(filePath)
tree = file['Events']
jets = ak.zip({
    'pt': tree['Jet_pt'].array(),
    'eta': tree['Jet_eta'].array(),
    'phi': tree['Jet_phi'].array(),
    'mass': tree['Jet_mass'].array(),
    'btag': tree['Jet_btagDeepFlavB'].array()
})

sorted_jets = jets[ak.argsort(jets.btag, axis=1, ascending=False)]
leading_bjets = sorted_jets[:, :2]
non_bjets = sorted_jets[:, 2:]

# Sort non_bjets by pt in descending order
sorted_non_bjets = non_bjets[ak.argsort(non_bjets.pt, axis=1, ascending=False)]

# Extract 2 leading non_bjets with highest pt
leading_non_bjets = sorted_non_bjets[:, :2]

bjet1 = leading_bjets[:, :1]
bjet2 = leading_bjets[:, 1:2]
ljet1 = leading_non_bjets[:, :1]
ljet2 = leading_non_bjets[:, 1:2]

bjet1_px = bjet1.pt * np.cos(bjet1.phi)
bjet1_py = bjet1.pt * np.sin(bjet1.phi)
bjet1_pz = bjet1.pt * np.sinh(bjet1.eta)
bjet1_e = np.sqrt(bjet1_px**2 + bjet1_py**2 + bjet1_pz**2 + bjet1.mass**2)

bjet2_px = bjet2.pt * np.cos(bjet2.phi)
bjet2_py = bjet2.pt * np.sin(bjet2.phi)
bjet2_pz = bjet2.pt * np.sinh(bjet2.eta)
bjet2_e = np.sqrt(bjet2_px**2 + bjet2_py**2 + bjet2_pz**2 + bjet2.mass**2)

ljet1_px = ljet1.pt * np.cos(ljet1.phi)
ljet1_py = ljet1.pt * np.sin(ljet1.phi)
ljet1_pz = ljet1.pt * np.sinh(ljet1.eta)
ljet1_e = np.sqrt(ljet1_px**2 + ljet1_py**2 + ljet1_pz**2 + ljet1.mass**2)

ljet2_px = ljet2.pt * np.cos(ljet2.phi)
ljet2_py = ljet2.pt * np.sin(ljet2.phi)
ljet2_pz = ljet2.pt * np.sinh(ljet2.eta)
ljet2_e = np.sqrt(ljet2_px**2 + ljet2_py**2 + ljet2_pz**2 + ljet2.mass**2)


jet_resolution = 0.05

delta_bjet1_px = np.sqrt((jet_resolution*bjet1.pt * np.cos(bjet1.phi))**2)
delta_bjet1_py = np.sqrt((jet_resolution*bjet1.pt * np.sin(bjet1.phi))**2)
delta_bjet1_pz = np.sqrt((jet_resolution*bjet1.pt * np.sinh(bjet1.eta))**2)
delta_bjet1_e = np.sqrt((bjet1.pt*bjet1.pt * jet_resolution*np.cosh(bjet1.eta)*np.cosh(bjet1.eta)/bjet1_e)**2)

delta_bjet2_px = np.sqrt((jet_resolution*bjet2.pt * np.cos(bjet2.phi))**2)
delta_bjet2_py = np.sqrt((jet_resolution*bjet2.pt * np.sin(bjet2.phi))**2)
delta_bjet2_pz = np.sqrt((jet_resolution*bjet2.pt * np.sinh(bjet2.eta))**2)
delta_bjet2_e = np.sqrt((bjet2.pt*bjet2.pt * jet_resolution*np.cosh(bjet2.eta)*np.cosh(bjet2.eta)/bjet2_e)**2)

delta_ljet1_px = np.sqrt((jet_resolution*ljet1.pt * np.cos(ljet1.phi))**2)
delta_ljet1_py = np.sqrt((jet_resolution*ljet1.pt * np.sin(ljet1.phi))**2)
delta_ljet1_pz = np.sqrt((jet_resolution*ljet1.pt * np.sinh(ljet1.eta))**2)
delta_ljet1_e = np.sqrt((ljet1.pt*ljet1.pt * jet_resolution*np.cosh(ljet1.eta)*np.cosh(ljet1.eta)/ljet1_e)**2)

delta_ljet2_px = np.sqrt((jet_resolution*ljet2.pt * np.cos(ljet2.phi))**2)
delta_ljet2_py = np.sqrt((jet_resolution*ljet2.pt * np.sin(ljet2.phi))**2)
delta_ljet2_pz = np.sqrt((jet_resolution*ljet2.pt * np.sinh(ljet2.eta))**2)
delta_ljet2_e = np.sqrt((ljet2.pt*ljet2.pt * jet_resolution*np.cosh(ljet2.eta)*np.cosh(ljet2.eta)/ljet2_e)**2)

muons = ak.zip({
    'pt': tree['Muon_pt'].array(),
    'eta': tree['Muon_eta'].array(),
    'phi': tree['Muon_phi'].array(),
    'mass': tree['Muon_mass'].array(),
    'charge': tree['Muon_charge'].array(),
    'tightid': tree['Muon_tightId'].array(),
})

muon = muons[ak.argmax(muons.pt, axis=1, keepdims=True)]

muon_px = muon.pt * np.cos(muon.phi)
muon_py = muon.pt * np.sin(muon.phi)
muon_pz = muon.pt * np.sinh(muon.eta)
muon_e = np.sqrt(muon_px**2 + muon_py**2 + muon_pz**2 + muon.mass**2)


muon_resolution = 0.05

delta_muon_px = np.sqrt((muon_resolution*muon.pt * np.cos(muon.phi))**2)
delta_muon_py = np.sqrt((muon_resolution*muon.pt * np.sin(muon.phi))**2)
delta_muon_pz = np.sqrt((muon_resolution*muon.pt * np.sinh(muon.eta))**2)
delta_muon_e = np.sqrt((muon.pt*muon.pt * muon_resolution*np.cosh(muon.eta)*np.cosh(muon.eta)/muon_e)**2)

mets = ak.zip({
    'pt': tree['MET_pt'].array(),
    'phi': tree['MET_phi'].array(),
})

met_px = mets.pt * np.cos(mets.phi)
met_py = mets.pt * np.sin(mets.phi)

# Define physical constants
MW = 80.4  # W-boson mass (GeV)
MT = 172.5  # Top quark mass (GeV)
delta_MW = 2.085  # W-boson mass resolution (GeV)
delta_MT = 1.43  # Top quark mass resolution (GeV)




met_pz1 , met_pz2 = solve_npz(muon_px, muon_py, muon_pz, muon_e, met_px, met_py)

met_e1 = np.sqrt(met_px**2 + met_py**2 + met_pz1**2)
met_e2 = np.sqrt(met_px**2 + met_py**2 + met_pz2**2)

met_resolution = 0.05
delta_met_px = np.sqrt((met_resolution*mets.pt * np.cos(mets.phi))**2)
delta_met_py = np.sqrt((met_resolution*mets.pt * np.sin(mets.phi))**2)
delta_met_pz1 = np.sqrt((met_resolution*mets.pt)**2)
delta_met_pz2 = np.sqrt((met_resolution*mets.pt)**2)
delta_met_e1 = np.sqrt((met_resolution*mets.pt)**2)
delta_met_e2 = np.sqrt((met_resolution*mets.pt)**2)

met_px_2d = ak.to_numpy(met_px).reshape(-1, 1)
met_py_2d = ak.to_numpy(met_py).reshape(-1, 1)

delta_met_e1_2d = ak.to_numpy(delta_met_e1).reshape(-1, 1)
delta_met_e2_2d = ak.to_numpy(delta_met_e2).reshape(-1, 1)
delta_met_px_2d = ak.to_numpy(delta_met_px).reshape(-1, 1)
delta_met_py_2d = ak.to_numpy(delta_met_py).reshape(-1, 1)
delta_met_pz1_2d = ak.to_numpy(delta_met_pz1).reshape(-1, 1)
delta_met_pz2_2d = ak.to_numpy(delta_met_pz2).reshape(-1, 1)


perm1 = np.concatenate([muon_e, muon_px, muon_py, muon_pz, bjet1_e, bjet1_px, bjet1_py, bjet1_pz, 
                        bjet2_e, bjet2_px, bjet2_py, bjet2_pz, ljet1_e, ljet1_px, ljet1_py, ljet1_pz, 
                        ljet2_e, ljet2_px, ljet2_py, ljet2_pz, met_e1, met_px_2d, met_py_2d, met_pz1], axis=1)

perm2 = np.concatenate([muon_e, muon_px, muon_py, muon_pz, bjet1_e, bjet1_px, bjet1_py, bjet1_pz,
                        bjet2_e, bjet2_px, bjet2_py, bjet2_pz, ljet1_e, ljet1_px, ljet1_py, ljet1_pz,
                        ljet2_e, ljet2_px, ljet2_py, ljet2_pz, met_e2, met_px_2d, met_py_2d, met_pz2], axis=1)

perm3 = np.concatenate([muon_e, muon_px, muon_py, muon_pz, bjet2_e, bjet2_px, bjet2_py, bjet2_pz,
                        bjet1_e, bjet1_px, bjet1_py, bjet1_pz, ljet1_e, ljet1_px, ljet1_py, ljet1_pz,
                        ljet2_e, ljet2_px, ljet2_py, ljet2_pz, met_e1, met_px_2d, met_py_2d, met_pz1], axis=1)

perm4 = np.concatenate([muon_e, muon_px, muon_py, muon_pz, bjet2_e, bjet2_px, bjet2_py, bjet2_pz,
                        bjet1_e, bjet1_px, bjet1_py, bjet1_pz, ljet1_e, ljet1_px, ljet1_py, ljet1_pz,
                        ljet2_e, ljet2_px, ljet2_py, ljet2_pz, met_e2, met_px_2d, met_py_2d, met_pz2], axis=1)


res1 = np.concatenate([delta_muon_e, delta_muon_px, delta_muon_py, delta_muon_pz, delta_bjet1_e, delta_bjet1_px, delta_bjet1_py, delta_bjet1_pz,
                        delta_bjet2_e, delta_bjet2_px, delta_bjet2_py, delta_bjet2_pz, delta_ljet1_e, delta_ljet1_px, delta_ljet1_py, delta_ljet1_pz,
                        delta_ljet2_e, delta_ljet2_px, delta_ljet2_py, delta_ljet2_pz, delta_met_e1_2d, delta_met_px_2d, delta_met_py_2d, delta_met_pz1_2d], axis=1)

res2 = np.concatenate([delta_muon_e, delta_muon_px, delta_muon_py, delta_muon_pz, delta_bjet1_e, delta_bjet1_px, delta_bjet1_py, delta_bjet1_pz,
                        delta_bjet2_e, delta_bjet2_px, delta_bjet2_py, delta_bjet2_pz, delta_ljet1_e, delta_ljet1_px, delta_ljet1_py, delta_ljet1_pz,
                        delta_ljet2_e, delta_ljet2_px, delta_ljet2_py, delta_ljet2_pz, delta_met_e2_2d, delta_met_px_2d, delta_met_py_2d, delta_met_pz2_2d], axis=1)

res3 = np.concatenate([delta_muon_e, delta_muon_px, delta_muon_py, delta_muon_pz, delta_bjet2_e, delta_bjet2_px, delta_bjet2_py, delta_bjet2_pz,
                        delta_bjet1_e, delta_bjet1_px, delta_bjet1_py, delta_bjet1_pz, delta_ljet1_e, delta_ljet1_px, delta_ljet1_py, delta_ljet1_pz,
                        delta_ljet2_e, delta_ljet2_px, delta_ljet2_py, delta_ljet2_pz, delta_met_e1_2d, delta_met_px_2d, delta_met_py_2d, delta_met_pz1_2d], axis=1)

res4 = np.concatenate([delta_muon_e, delta_muon_px, delta_muon_py, delta_muon_pz, delta_bjet2_e, delta_bjet2_px, delta_bjet2_py, delta_bjet2_pz,
                        delta_bjet1_e, delta_bjet1_px, delta_bjet1_py, delta_bjet1_pz, delta_ljet1_e, delta_ljet1_px, delta_ljet1_py, delta_ljet1_pz,
                        delta_ljet2_e, delta_ljet2_px, delta_ljet2_py, delta_ljet2_pz, delta_met_e2_2d, delta_met_px_2d, delta_met_py_2d, delta_met_pz2_2d], axis=1)


np.savez('output_arrays.npz', perm1=perm1, perm2=perm2, perm3=perm3, perm4=perm4, res1=res1, res2=res2, res3=res3, res4=res4)
print("Saving permutations and resolutions to 'output_arrays.npz'.")
# Assume perm1, perm2, perm3, perm4 and their resolutions (res1, res2, res3, res4) are given
permutations = [perm1, perm2, perm3, perm4]
resolutions = [res1, res2, res3, res4]

n_events = len(perm1)
# n_events = 100

# Prepare data for all events
event_data = [(i, permutations, resolutions) for i in range(n_events)]

# Use all available CPU cores
with ProcessPoolExecutor() as executor:
    results = list(executor.map(optimize_event, event_data))

np.savez('optimization_results.npz', results=results)
print("Optimization process completed and results saved to 'optimization_results.npz'.")