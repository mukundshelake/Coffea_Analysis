import os
import h5py
import numpy as np
import awkward as ak
from scipy.optimize import minimize
from concurrent.futures import ProcessPoolExecutor
import argparse
from datetime import datetime

# Define physical constants
MW = 80.4  # W-boson mass (GeV)
MT = 172.5  # Top quark mass (GeV)
delta_MW = 2.085  # W-boson mass resolution (GeV)
delta_MT = 1.43  # Top quark mass resolution (GeV)

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
    alpha = MW**2 - mpe**2 + mpx**2 + mpy**2 + mpz**2 + 2*(mpx*metx + mpy*mety)
    a = mpz**2 - mpe**2
    b = alpha*mpz 
    c = alpha**2/4 - mpe**2*(metx**2 + mety**2)

    discriminant = b**2 - 4*a*c

    pz1 = ak.where(discriminant < 0, 0, (-b + np.sqrt(discriminant))/(2*a))
    pz2 = ak.where(discriminant < 0, 0, (-b - np.sqrt(discriminant))/(2*a))

    return pz1, pz2

def invariant_mass(candidate):
    mass2 = candidate[0]**2 - candidate[1]**2 - candidate[2]**2 - candidate[3]**2
    if mass2 < 0:
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

    return chi_sq

def process_file(file_path):
    print(f"Processing file: {file_path}")
    with h5py.File(file_path, 'r') as f:
        perm1 = f['perm1'][:]
        perm2 = f['perm2'][:]
        perm3 = f['perm3'][:]
        perm4 = f['perm4'][:]
        res1 = f['res1'][:]
        res2 = f['res2'][:]
        res3 = f['res3'][:]
        res4 = f['res4'][:]

    permutations = [perm1, perm2, perm3, perm4]
    resolutions = [res1, res2, res3, res4]

    n_events = len(perm1)
    event_data = [(i, permutations, resolutions) for i in range(n_events)]

    print(f"Starting optimization for {n_events} events...")
    with ProcessPoolExecutor() as executor:
        results = list(executor.map(optimize_event, event_data))

    output_file = file_path.replace('.h5', '_results.npz')
    np.savez(output_file, results=results)
    print(f"Optimization process completed and results saved to '{output_file}'.")

def main(input_dir):

    h5_files = []
    for root, dirs, files in os.walk(input_dir):
        for file in files:
            if file.endswith('.h5'):
                h5_files.append(os.path.join(root, file))

    print(f"Found {len(h5_files)} HDF5 files in '{input_dir}'")
    for h5_file in h5_files:
        process_file(h5_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process HDF5 files and save optimization results.")
    parser.add_argument('-t', '--timestamp', type=str, required=True, help="Timestamp to locate the input and output folders.")
    args = parser.parse_args()

    timestamp = args.timestamp
    input_dir = f'outputs/{timestamp}'

    print(f"Input directory: {input_dir}")

    main(input_dir)