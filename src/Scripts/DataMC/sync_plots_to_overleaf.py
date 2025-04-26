#!/usr/bin/env python3
import os
import shutil
import yaml
from pathlib import Path

# Configuration
CONFIG_FILE = "/home/mukund/Projects/updatedCoffea/Coffea_Analysis/src/Scripts/DataMC/plot_tags.yaml"
SOURCE_DIR = "/home/mukund/Projects/updatedCoffea/Coffea_Analysis/src/Scripts/DataMC/plots"
TARGET_DIR = "/home/mukund/Projects/AN_overleaf/6734a241cdf53075ce1278fe/figures/DataMC"

def load_config():
    """Load YAML config file with target directories for each era"""
    with open(CONFIG_FILE) as f:
        return yaml.safe_load(f) or {}

def get_clean_filename(filename, tag):
    """Remove tag from filename"""
    return filename.replace(f"_{tag}", "")

def sync_plots():
    """Main function to sync plots to target directories"""
    config = load_config()
    
    for era, tag in config.items():
        print(f"Syncing plots for era: {era} with tag: {tag}")
        target_dir = os.path.join(TARGET_DIR, era)
            
        Path(target_dir).mkdir(parents=True, exist_ok=True)
        
        for filename in os.listdir(SOURCE_DIR):
            if tag not in filename:
                continue

            if era not in filename:
                continue
                
            clean_name = get_clean_filename(filename, tag)
            source_path = os.path.join(SOURCE_DIR, filename)
            target_path = os.path.join(target_dir, clean_name)
            
            if os.path.exists(target_path):
                os.remove(target_path)
                
            shutil.copy2(source_path, target_path)
            print(f"Copied {filename} to {target_path}")

if __name__ == "__main__":
    sync_plots()