import uproot
import awkward as ak
import json
import argparse
import os
from tqdm import tqdm

def get_puid_counts(file_path, tree_name="Events"):
    """
    Counts jets based on puId from a single ROOT file, after applying event selection.
    Event selection: Sum$(Jet_pt > 25 && abs(Jet_eta) < 2.4) > 3

    Args:
        file_path (str): Path to the ROOT file.
        tree_name (str): Name of the TTree containing events.

    Returns:
        tuple: (nTight, nMedium, nLoose, nTotal) counts for jets in selected events.
               Returns (0, 0, 0, 0) if the file/tree cannot be opened or required branches are missing.
    """
    try:
        with uproot.open(f"{file_path}:{tree_name}") as tree:
            # Check for required branches
            required_branches = ["Jet_pt", "Jet_eta", "Jet_puId"]
            missing_branches = [b for b in required_branches if b not in tree]
            if missing_branches:
                print(f"Warning: Missing branches {missing_branches} in {file_path}. Skipping.")
                return 0, 0, 0, 0

            # Read necessary arrays
            events = tree.arrays(required_branches, library="ak")

            # Apply jet selection within each event
            jet_selection_mask = (events.Jet_pt > 25) & (abs(events.Jet_eta) < 2.4)

            # Count selected jets per event
            n_selected_jets_per_event = ak.sum(jet_selection_mask, axis=1)

            # Create event selection mask
            event_selection_mask = (n_selected_jets_per_event > 3)

            # Filter events based on the event selection mask
            selected_events = events[event_selection_mask]

            # Also filter the jet selection mask to match the selected events
            jet_mask_in_selected_events = jet_selection_mask[event_selection_mask]

            # Apply the jet selection mask *within* the selected events to get the puIds to count
            puids_to_count = selected_events.Jet_puId[jet_mask_in_selected_events]

            # Flatten the puId array from selected jets in selected events
            flat_puid = ak.flatten(puids_to_count, axis=None)

            if len(flat_puid) == 0: # Handle case where no jets pass selection in selected events
                 return 0, 0, 0, 0

            # Define masks for specific puId values for the selected jets
            is_tight = (flat_puid == 7)
            is_medium = (flat_puid == 3)
            is_loose = (flat_puid == 1)

            # Count jets inclusively for each category within selected events
            n_tight = ak.sum(is_tight)
            # Medium includes Medium OR Tight
            n_medium = ak.sum(is_medium | is_tight)
            # Loose includes Loose OR Medium OR Tight
            n_loose = ak.sum(is_loose | is_medium | is_tight)
            # Total is the count of all jets passing kinematic cuts in selected events
            n_total = len(flat_puid)

            return n_tight, n_medium, n_loose, n_total

    except Exception as e:
        print(f"Error processing file {file_path}: {e}")
        return 0, 0, 0, 0

def main(input_json_path, output_json_path, data_key="MC_mu"):
    """
    Main function to process datasets and generate efficiency JSON.

    Args:
        input_json_path (str): Path to the input JSON file containing dataset paths.
        output_json_path (str): Path to save the output JSON file.
        data_key (str): The top-level key in the input JSON to process (e.g., "MC_mu").
    """
    try:
        with open(input_json_path, 'r') as f:
            all_datasets = json.load(f)
    except FileNotFoundError:
        print(f"Error: Input JSON file not found at {input_json_path}")
        return
    except json.JSONDecodeError:
        print(f"Error: Could not decode JSON from {input_json_path}")
        return

    if data_key not in all_datasets:
        print(f"Error: Key '{data_key}' not found in {input_json_path}")
        return

    mc_datasets = all_datasets[data_key]
    results = {}

    print(f"Processing datasets under key '{data_key}'...")
    for dataset_name, files_dict in tqdm(mc_datasets.items(), desc="Datasets"):
        print(f"\nProcessing dataset: {dataset_name}")
        total_tight = 0
        total_medium = 0
        total_loose = 0
        total_jets = 0

        file_paths = list(files_dict.keys())
        for file_path in tqdm(file_paths, desc=f"  Files in {dataset_name}", leave=False):
            # Assuming the JSON provides the tree name as value, but we use default "Events"
            # tree_name = files_dict[file_path] # Use if tree name varies
            tree_name = "Events" # Hardcoding for now based on example
            n_tight, n_medium, n_loose, n_total = get_puid_counts(file_path, tree_name)
            total_tight += n_tight
            total_medium += n_medium
            total_loose += n_loose
            total_jets += n_total

        results[dataset_name] = {
            "nTight": float(total_tight),
            "nMedium": float(total_medium),
            "nLoose": float(total_loose),
            "nTotal": float(total_jets)
        }
        print(f"  Finished {dataset_name}: Tight={total_tight}, Medium={total_medium}, Loose={total_loose}, Total={total_jets}")


    # Ensure the output directory exists
    output_dir = os.path.dirname(output_json_path)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")

    # Write results to output JSON
    try:
        with open(output_json_path, 'w') as f:
            json.dump(results, f, indent=4)
        print(f"\nSuccessfully wrote results to {output_json_path}")
    except IOError as e:
        print(f"Error writing output JSON to {output_json_path}: {e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate Jet puId efficiencies from ROOT files listed in a JSON.")
    parser.add_argument("input_json", help="Path to the input JSON file containing dataset file paths.")
    parser.add_argument("output_json", help="Path to the output JSON file to save results.")
    parser.add_argument("--data_key", default="MC_mu", help="The top-level key in the input JSON to process (default: MC_mu).")
    # parser.add_argument("--era", required=True, help="Specify the era (e.g., UL2016preVFP, UL2016postVFP) - used for output filename logic if needed.") # Example if era needed for output path

    args = parser.parse_args()

    main(args.input_json, args.output_json, args.data_key)