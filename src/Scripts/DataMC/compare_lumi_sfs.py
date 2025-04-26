import json
from collections import defaultdict

def normalize_key(key):
    """Remove year-specific prefix from a process name"""
    for prefix in ["UL2016preVFP_", "UL2017_"]:
        if key.startswith(prefix):
            return key.replace(prefix, "")
    return key  # fallback


def calculate_sfs(config):
    """Calculate scale factors for all processes in a config"""
    sfs = {}
    luminosity = config["Luminosity"]
    for process, xsec in config["cross_sections"].items():
        if process in config["generated_events"]:
            gen_events = config["generated_events"][process]
            sfs[normalize_key(process)] = luminosity * xsec / gen_events
    return sfs

def group_by_category(config, sfs):
    """Group scale factors by their category"""
    grouped = defaultdict(dict)
    for category, processes in config["category_map"].items():
        for process in processes:
            norm_process = normalize_key(process)
            if norm_process in sfs:
                grouped[category][norm_process] = sfs[norm_process]
    return grouped

def print_comparison(sfs_2016, sfs_2017, grouped_2016, grouped_2017):
    """Print comparison table of scale factors"""
    print("\n{:<30} {:<20} {:<20} {:<15}".format(
        "Process", "2016preVFP SF", "2017 SF", "Ratio (2017/2016)"))
    print("-" * 85)
    
    # Print all processes
    for process in sorted(set(sfs_2016.keys()).union(sfs_2017.keys())):
        sf16 = sfs_2016.get(process, "N/A")
        sf17 = sfs_2017.get(process, "N/A")
        ratio = float(sf17)/float(sf16) if (process in sfs_2016 and process in sfs_2017) else "N/A"
        print("{:<30} {:<20} {:<20} {:<15}".format(
            process,
            f"{sf16:.4g}" if isinstance(sf16, (int, float)) else sf16,
            f"{sf17:.4g}" if isinstance(sf17, (int, float)) else sf17,
            f"{ratio:.4g}" if isinstance(ratio, (int, float)) else ratio
        ))
    
    # Print grouped by category
    print("\n\n{:<30} {:<20} {:<20} {:<15}".format(
        "Category", "2016preVFP Avg SF", "2017 Avg SF", "Ratio (2017/2016)"))
    print("-" * 85)
    
    for category in sorted(set(grouped_2016.keys()).union(grouped_2017.keys())):
        avg16 = sum(grouped_2016[category].values())/len(grouped_2016[category]) if category in grouped_2016 else "N/A"
        avg17 = sum(grouped_2017[category].values())/len(grouped_2017[category]) if category in grouped_2017 else "N/A"
        ratio = float(avg17)/float(avg16) if (category in grouped_2016 and category in grouped_2017) else "N/A"
        
        print("{:<30} {:<20.4g} {:<20.4g} {:<15.4g}".format(
            category, avg16, avg17, ratio))

def main():
    # Load configuration files
    with open('UL2016preVFP_configs.json') as f:
        config_2016 = json.load(f)
    with open('UL2017_configs.json') as f:
        config_2017 = json.load(f)
    print("Loaded configurations for 2016 and 2017.")
    # Calculate scale factors
    sfs_2016 = calculate_sfs(config_2016)
    sfs_2017 = calculate_sfs(config_2017)
    
    # Group by category
    grouped_2016 = group_by_category(config_2016, sfs_2016)
    grouped_2017 = group_by_category(config_2017, sfs_2017)
    
    # Print comparison
    print_comparison(sfs_2016, sfs_2017, grouped_2016, grouped_2017)

if __name__ == "__main__":
    main()