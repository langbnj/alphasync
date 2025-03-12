#!/usr/bin/env python3
"""
alphamap_uniprot_cleanup.py:
Removes AlphaSync prediction (CIF, PAE and params files) that are no longer necessary since there are better structures available for their sequences.
"""

# Initialize
import pandas as pd
# import numpy as np
# import networkx as nx
# import seaborn as sns
# import matplotlib.pyplot as mp
import tarfile
import gzip
import io
import requests
import json
# from Bio import SeqIO
from blang_mysql import *
from blang import *

# type = "alphauniprot"
type = "uniprot"

# alphamap = "alphasync_compact.alphamap"
alphamap = "alphamap"
alphafrag = "alphafrag"
alphaseq = "alphaseq"
alphasa = "alphasa"
alphacon = "alphacon"

# Paths for deleting outdated AlphaSync structure predictions (where an accession has since been made obsolete, or a sequence was changed, in UniProt)
cifdir = "input/alphasync/cif"
paedir = "input/alphasync/pae"
paramdir = "input/alphasync/params"

(version) = Args(1, "[Current UniProt version]\n -debug: Don't actually make any changes, just simulate", "2025_01")
print()

# Get AFDB=0 accessions in alphaseq
print(f" >> alphaseq afdb=0 accessions:")
alphaseq_afdb0 = FetchSet(Query(f"SELECT DISTINCT acc FROM {alphaseq} WHERE afdb=0"))
print(f"   >> {len(alphaseq_afdb0):,}")

# Get AFDB=0 accessions that are actually mapped to in alphamap
print(f" >> alphamap afdb=0 accessions that are actually mapped to:")
alphamap_afdb0 = FetchSet(Query(f"SELECT DISTINCT map FROM {alphamap} WHERE type='{type}' AND version='{version}' AND afdb=0 AND best=1"))
print(f"   >> {len(alphamap_afdb0):,}")

print(f" >> unnecessary structures (alphaseq afdb=0 accessions that are not mapped to in alphamap):")
unnecessary_structures = alphaseq_afdb0 - alphamap_afdb0
print(f"   >> {len(unnecessary_structures):,}")

print()
print("Deleting unnecessary structure CIF/PAE/params output files:")
for acc in tq(unnecessary_structures):
    if not Switch('debug'):

        # Delete output files (CIF/PAE/params)
        Run(f"Delete all output files for obsolete acc '{acc}' (if they exist)", f"rm -fv {cifdir}/AF-{acc}-F*-model_v0.cif.gz {paedir}/AF-{acc}-F*-predicted_aligned_error_v0.json.gz {paramdir}/AF-{acc}-F*-alphafold_params.json")
        Log(f"obsolete CIF/PAE/params files deleted for acc", acc)

        # First, check if there are any best=1 mappings for this acc in table 'alphamap' (there shouldn't be)
        query = Query(f"SELECT * FROM {alphamap} WHERE type='{type}' AND version='{version}' AND afdb=0 AND map='{acc}' AND best=1")
        if Numrows(query) > 0:
            Die(f"Error: Found best=1 mappings for obsolete acc '{acc}' in table '{alphamap}' (shouldn't happen)")
        # Remove acc from table 'alphamap'
        Query(f"DELETE FROM {alphamap} WHERE type='{type}' AND version='{version}' AND afdb=0 AND map='{acc}'")
        Log(f"obsolete 'map' acc deleted from table '{alphamap}' for acc", acc)

        # Remove acc from table 'alphafrag'
        Query(f"DELETE FROM {alphafrag} WHERE afdb=0 AND acc='{acc}'")
        Log(f"obsolete acc deleted from table '{alphafrag}' for acc", acc)

        # Remove acc from table 'alphaseq'
        Query(f"DELETE FROM {alphaseq} WHERE afdb=0 AND acc='{acc}'")
        Log(f"obsolete acc deleted from table '{alphaseq}' for acc", acc)

        # Remove acc from table 'alphasa'
        Query(f"DELETE FROM {alphasa} WHERE afdb=0 AND acc='{acc}'")
        Log(f"obsolete acc deleted from table '{alphasa}' for acc", acc)

        # Remove acc from table 'alphacon'
        Query(f"DELETE FROM {alphacon} WHERE afdb=0 AND acc='{acc}'")
        Log(f"obsolete acc deleted from table '{alphacon}' for acc", acc)

    else:
        
        Log(f"obsolete CIF/PAE/params files & table rows would have been deleted (but -debug is active) for acc", acc)
        
        query = Query(f"SELECT * FROM {alphamap} WHERE type='{type}' AND version='{version}' AND afdb=0 AND map='{acc}' AND best=1")
        if Numrows(query) > 0:
            Die(f"Error: Found best=1 mappings for obsolete acc '{acc}' in table '{alphamap}' (shouldn't happen)")
        
        query = Query(f"SELECT * FROM {alphamap} WHERE type='{type}' AND version='{version}' AND afdb=0 AND map='{acc}'")
        if Numrows(query) == 0:
            Die(f"Error: No mappings found for obsolete acc '{acc}' in table '{alphamap}' (shouldn't happen)")
        
        query = Query(f"SELECT * FROM {alphafrag} WHERE afdb=0 AND acc='{acc}'")
        if Numrows(query) == 0:
            Die(f"Error: No rows found for obsolete acc '{acc}' in table '{alphafrag}' (shouldn't happen)")

        query = Query(f"SELECT * FROM {alphaseq} WHERE afdb=0 AND acc='{acc}'")
        if Numrows(query) == 0:
            Die(f"Error: No rows found for obsolete acc '{acc}' in table '{alphaseq}' (shouldn't happen)")

        query = Query(f"SELECT * FROM {alphasa} WHERE afdb=0 AND acc='{acc}'")
        if Numrows(query) == 0:
            Die(f"Error: No rows found for obsolete acc '{acc}' in table '{alphasa}' (shouldn't happen)")

        query = Query(f"SELECT * FROM {alphacon} WHERE afdb=0 AND acc='{acc}'")
        if Numrows(query) == 0:
            Die(f"Error: No rows found for obsolete acc '{acc}' in table '{alphacon}' (shouldn't happen)")

Show(lim=0)

print("\nDone!")
