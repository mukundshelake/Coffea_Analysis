import uproot
import glob
import json
import os
from multiprocessing import Pool, cpu_count

era = "SIXTEEN_postVFP"

# Root directory where your NanoAOD files are stored
root_dir = f"/nfs/home/common/RUN2_UL/Tree_crab/{era}/Data_mu"  # Change this to your actual directory

# Recursively find all ROOT files
files = glob.glob(os.path.join(root_dir, "**/*.root"), recursive=True)
# files = [files[0], files[1]]

def process_file(file):
    """Process a single ROOT file and extract run and lumisection information."""
    print(f"Processing: {file}")
    lumi_dict = {}
    try:
        with uproot.open(file) as f:
            tree = f["Events"]

            # Read 'run' and 'luminosityBlock' branches
            runs = tree["run"].array(library="np")
            lumis = tree["luminosityBlock"].array(library="np")

            # Store in dictionary
            for run, lumi in zip(runs, lumis):
                if run not in lumi_dict:
                    lumi_dict[run] = set()
                lumi_dict[run].add(lumi)
    except Exception as e:
        print(f"Error processing {file}: {e}")
    return lumi_dict

def convert_to_ranges(lst):
    if not lst:
        return []
    ranges = []
    start = lst[0]
    for i in range(1, len(lst)):
        if lst[i] != lst[i - 1] + 1:
            ranges.append(f"{start}-{lst[i - 1]}" if start != lst[i - 1] else f"{start}")
            start = lst[i]
    ranges.append(f"{start}-{lst[-1]}" if start != lst[-1] else f"{start}")
    return ranges

def merge_dicts(dicts):
    """Merge a list of dictionaries into a single dictionary."""
    merged = {}
    for d in dicts:
        for run, lumis in d.items():
            if run not in merged:
                merged[run] = set()
            merged[run].update(lumis)
    return merged

if __name__ == "__main__":
    # Use multiprocessing to process files in parallel
    with Pool(cpu_count()) as pool:
        results = pool.map(process_file, files)

    # Merge results from all processes
    lumi_dict = merge_dicts(results)

    # Convert sets to sorted lists for JSON output
    for run in lumi_dict:
        lumi_dict[run] = sorted(list(lumi_dict[run]))

    print(lumi_dict)


    output_file = f"lumi_info_{era}.json"
    with open(output_file, "w") as outfile:
        lumi_dict = {int(k): [int(x) for x in v] for k, v in lumi_dict.items()}
        json.dump(lumi_dict, outfile, indent=4)

    print(f"Lumisection extraction complete! Saved as {output_file}")

    converted_data = {key: convert_to_ranges(value) for key, value in lumi_dict.items()}

    with open(f'lumi_range_{era}.json', 'w') as file:
        json.dump(converted_data, file, indent=4)