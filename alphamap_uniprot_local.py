#!/usr/bin/env python3
"""
alphamap_uniprot.py:
Uses perfect sequence matching to map identifiers for a given version of UniProt to accessions in AlphaSync tables such as 'alphaseq'.
"""

# Initialize
# import pandas as pd
# import numpy as np
# import networkx as nx
# import seaborn as sns
# import matplotlib.pyplot as mp
import tarfile
import gzip
import io
# from Bio import SeqIO
from blang_mysql import *
from blang import *

type = "uniprot"

alphamap = "alphamap"
uniseq = "uniseq"

(version) = Args(1, "[Current UniProt version]\n -debug: Don't actually make any changes, just simulate", "2022_04")

if not Switch('debug'):
    # Clear(alphamap)
    print(f"Deleting type '{type}' version '{version}' rows from table '{alphamap}'...")
    query = Query(f"DELETE FROM {alphamap} WHERE type='{type}' AND version='{version}'")
    print(f"Rows affected: {Comma(Numrows(query))}")

# Start

State(f"Getting mapping from UniProt species mnemonics to NCBI taxon IDs from table 'unitax' and inserting into table '{alphamap}':")
tax = {}
for species, thistax in Fetch(Query(f"SELECT species, tax FROM unitax")):
    tax[species] = thistax

State(f"Filling table '{alphamap}' with type '{type}' mappings for UniProt version '{version}' based on perfect sequence matching in table '{uniseq}':")
for acc, species, seq in Fetch(Query(f"SELECT acc, species, seq FROM {uniseq} WHERE type IN ('UniProt', 'UniIso') ORDER BY species='human' DESC, species, acc")):
    Log("total accs", acc)
    Log("total species", species)
    Log("total sequences", seq)

    # Get all alphaseq accs that have this exact sequence (we don't care about the species here)
    query = Query(f"SELECT s.acc, s.afdb, AVG(a.plddt) AS avg_plddt FROM alphaseq s, alphasa a WHERE s.seq='{seq}' AND a.acc=s.acc AND a.afdb=s.afdb GROUP BY s.acc ORDER BY avg_plddt DESC, s.afdb ASC, s.acc")
    if Numrows(query) == 0:
        # Insert NULL mapping into alphamap
        Query(f"INSERT INTO {alphamap} SET type='{type}', version='{version}', species='{species}', tax='{tax[species]}', value='{acc}', map=NULL, afdb=NULL, avg_plddt=NULL, best=NULL")

        Log("sequence not found in alphaseq for acc (skipped)", acc)
        Log("sequence not found in alphaseq for seq (skipped)", seq)
        # Check if sequence contains non-AA characters
        if not Aa(seq):
            Log("sequence not found in alphaseq & contains non-AA characters for acc (skipped)", acc)
            Log("sequence not found in alphaseq & contains non-AA characters for seq (skipped)", seq)
        continue
    best = 1
    for alphacc, afdb, avg_plddt in query:
        # Insert mapping into alphamap
        Query(f"INSERT INTO {alphamap} SET type='{type}', version='{version}', species='{species}', tax='{tax[species]}', value='{acc}', map='{alphacc}', afdb='{afdb}', avg_plddt='{avg_plddt}', best='{best}'")
        best = 0

        Log("successfully mapped for acc|alphacc", f"{acc}|{alphacc}")
        Log("successfully mapped for acc", acc)
        Log("successfully mapped for alphacc", alphacc)
        Log("successfully mapped for seq", seq)
        Log("successfully mapped for species", species)
        Log("successfully mapped for tax", tax[species])
        if afdb == 0:
            Log("successfully mapped to non-AFDB AlphaSync acc for acc|alphacc", f"{acc}|{alphacc}")
            Log("successfully mapped to non-AFDB AlphaSync acc for acc", acc)
            Log("successfully mapped to non-AFDB AlphaSync acc for alphacc", alphacc)
            Log("successfully mapped to non-AFDB AlphaSync acc for seq", seq)
            Log("successfully mapped to non-AFDB AlphaSync acc for species", species)
            Log("successfully mapped to non-AFDB AlphaSync acc for tax", tax[species])

Show(lim=0, sort=True)

State("Average number of AlphaFold DB structures per unique sequence: " + str(round(len(Get("successfully mapped for acc|alphacc")) / len(Get("successfully mapped for acc")), 2)))
State("Sequences that still need structure predictions for complete coverage: " + Comma(len(Get("sequence not found in alphaseq for seq (skipped)"))))

# Optimize(alphamap)

print("\nDone!")
