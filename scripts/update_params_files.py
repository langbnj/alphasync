#!/usr/bin/env python3
"""
One-time update script to add e.g. alphafold_version to params JSON files that didn't contain this yet
"""

# Initialize
# import pandas as pd
# import numpy as np
# import networkx as nx
# import seaborn as sns
# import matplotlib.pyplot as mp
# import json
# import os
import json
import sys
from pathlib import Path
from blang_mysql import *
from blang import *

def update_json_keys(infile):
    """Update JSON files with version information."""
    # Define the desired key order
    key_order = [
        "alphafold_version",
        "hhblits_version",
        "hhblits_limits",
        "uniref30_version",
        "bfd_version",
        "uniref90_version",
        "mgnify_version",
        "pdb_type",
        "pdb_version",
        "max_template_date"
    ]
    
    with open(infile, 'r') as f:
        data = json.load(f)
    
    modified = False
    new_data = {}

    if Switch('debug'):
        print(f" >> {infile}")
    Log("total files", infile)

    # Log number of keys
    if Switch('debug'):
        print(f"   >> {len(data)} input keys")
    Log(f"input number of keys is {len(data)} for file", infile)
    
    # Add alphafold_version if missing
    if 'alphafold_version' not in data:
        if Switch('debug'):
            print(f"   >> added alphafold_version")
        Log("update: added alphafold_version", infile)
        new_data['alphafold_version'] = '2.3.2'
        modified = True
    
    # Convert pdb70_version if present
    if 'pdb70_version' in data:
        if Switch('debug'):
            print(f"   >> converted pdb70_version to pdb_type and pdb_version")
        Log("update: converted pdb70_version to pdb_type and pdb_version", infile)
        new_data['pdb_type'] = 'pdb70'
        new_data['pdb_version'] = 'pdb70_from_mmcif_220313'
        data.pop('pdb70_version')
        modified = True

    # Remove hhblits_limits from hhblits_version, if it contains it
    if 'hhblits_version' in data:
        if data['hhblits_version'].endswith('_high_limits') or data['hhblits_version'].endswith('_high_cpus') or data['hhblits_version'].endswith('_low_limits'):
            if Switch('debug'):
                print(f"   >> removed hhblits_limits from hhblits_version")
            Log("update: removed hhblits_limits from hhblits_version", infile)
            new_data['hhblits_version'] = data['hhblits_version'].replace('_high_limits', '').replace('_high_cpus', '').replace('_low_limits', '')
            data.pop('hhblits_version')
            modified = True

    # Rename "hhblits_limits": "alphafold_default" to "default"
    if 'hhblits_limits' in data and data['hhblits_limits'] == 'alphafold_default':
        if Switch('debug'):
            print(f"   >> renamed hhblits_limits from 'alphafold_default' to 'default'")
        Log("update: renamed hhblits_limits from 'alphafold_default' to 'default'", infile)
        new_data['hhblits_limits'] = 'default'
        data.pop('hhblits_limits')
        modified = True

    # Check if there's a final newline in infile, and mark it modified if not (so it'll get written out with one)
    with open(infile, 'rb') as f:
        if not f.read().endswith(b'\n'):
            modified = True
            if Switch('debug'):
                print(f"   >> added final newline")
            Log("update: added final newline", infile)
    
    # If any modifications were made, update and write back
    if modified:

        # Add remaining keys from original data
        new_data.update(data)

        # Create final ordered dictionary
        final_data = {}
        for key in key_order:
            if key in new_data or key in data:
                # Use new_data if available, otherwise use original data
                final_data[key] = new_data.get(key, data.get(key))
        
        # Log number of output keys
        if Switch('debug'):
            print(f"   >> {len(final_data)} output keys")
        Log(f"output number of keys is {len(final_data)} for file", infile)

        if not Switch('debug'):
            # print(f"     >> updated")
            Log("updated files", infile)

            # Write back with matching formatting and ordered keys (unless -debug is active)
            with open(infile, 'w') as f:
                json.dump(final_data, f, indent=4, separators=(',', ': '))
                # Final newline
                f.write('\n')
        else:
            print(f"     >> would have been updated (but -debug is active)")
            Log("would have updated files (but -debug is active)", infile)

def process_folder(inpath):
    folder = Path(inpath)
    if not folder.is_dir():
        print(f"Error: {inpath} is not a directory")
        return
        
    # Process all JSON files in the directory
    print(f"\nProcessing all JSON files from '{inpath}':")
    infiles = list(folder.glob('*.json'))
    for infile in tq(infiles):
        try:
            update_json_keys(infile)
        except Exception as e:
            print(f"Error processing {infile}: {e}")

if __name__ == "__main__":
    # if len(sys.argv) != 2:
    #     print("Usage: python test_update_params_files.py input/alphasync/_backup_modelonly_expanded_iso_nofrag_peptides/params3")
    #     sys.exit(1)
    
    inpath = Args(1, "Path to the folder containing JSON files", ".")

    Starttime()
    process_folder(inpath)
    Stoptime()

    Show(sort=1)
    print("\nDone!")

# Run script on test folder & show differences:
# cd input/alphasync/_backup_modelonly_expanded_iso_nofrag_peptides/params3 ; \cp -p ../*.json .; scripts/tmp_update_params_files.py . ; find . -type f|xa sh -c 'echo -e "\n\n{}"; diff ../{} {}'

# Run script on final folder:
# scripts/tmp_update_params_files.py input/alphasync/params

# Compare new output folder to the updated backup:
# input/alphasync/params >> rsync -navic --delete ./ ../_backup_modelonly_expanded_iso_nofrag_peptides/params/ | cut -f1 -d' ' | suq

# Show unique combinations
# 1|xa sh -c "cat {}|perl -npe 's/\n/\t/g'"|suq
