#!/usr/bin/env python3
"""
alphamap_uniprot.py:
Fills table 'alphamap'. Uses perfect sequence matching to map identifiers for a given version of UniProt (in table 'alphauniprot', via UniProt's API) to AlphaFold DB accessions in AlphaSync table 'alphaseq'.
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

alphamap = "alphamap"
alphasa = "alphasa"
alphaseq = "alphaseq"
alphauniprot = "alphauniprot"

(version) = Args(1, "[Current UniProt version]\n -debug: Don't actually make any changes, just simulate", "2025_01")

# if not Switch('debug'):
#     if Switch('alphasync'):
#         # AlphaSync: Only delete mappings for accs that are not in AFDB (afdb=0)
#         print(f"Getting type '{type}' version '{version}' AlphaSync accs (afdb=0) from table '{alphaseq}'...")
#         # Get all alphaccs that are not from AFDB
#         query = Query(f"SELECT DISTINCT acc FROM {alphaseq} WHERE afdb=0")
#         alphasync_accs = FetchSet(query)
#         # Delete mappings for these accs (to re-map them)
#         print(f"Setting these to map=NULL in table '{alphamap}'...")
#         affected = 0
#         for acc in alphasync_accs:
#             query = Query(f"UPDATE {alphamap} SET map=NULL WHERE type='{type}' AND version='{version}' AND map='{acc}'")
#             affected += Numrows(query)
#         print(f"Rows affected: {Numrows(query):,}")
#     else:
# The only clean option is to re-run the mapping completely, even for an AlphaSync update:
# Clear(alphamap)
print(f"Deleting type '{type}' version '{version}' rows from table '{alphamap}'...")
query = Query(f"DELETE FROM {alphamap} WHERE type='{type}' AND version='{version}'")
print(f"Rows affected: {Numrows(query):,}")
print()

# # Pre-fetch AlphaFold DB structures for each unique sequence in AlphaSync table 'alphaseq' and put them into a dict
# Too slow
# State("Pre-fetching sequence-matched AlphaFold DB structures for each row in table 'alphaseq':")
# structures = {}
# query = Query(f"SELECT s.seq, s.acc, s.afdb, AVG(a.plddt) AS avg_plddt FROM {alphaseq} s, {alphasa} a WHERE a.acc=s.acc AND a.afdb=s.afdb GROUP BY s.acc LIMIT 100000")
# # query = Query(f"SELECT s.seq, s.acc, s.afdb, AVG(a.plddt) AS avg_plddt FROM {alphaseq} s, {alphasa} a WHERE s.acc='Q12753' AND a.acc=s.acc AND a.afdb=s.afdb GROUP BY s.acc ORDER BY avg_plddt DESC, s.afdb ASC, s.acc")
# for seq, acc, afdb, avg_plddt in Fetch(query):
#     if seq not in structures:
#         structures[seq] = []
#     structures[seq].append((acc, afdb, avg_plddt))
# # Sort each list by descending avg_plddt
# State("Sorting AlphaFold DB structures by descending avg_plddt:")
# for seq in tq(structures):
#     structures[seq] = sorted(structures[seq], key=lambda x: x[2], reverse=True)


# Start
State(f"Filling table '{alphamap}' with type '{type}' mappings for UniProt version '{version}' based on perfect sequence matching between tables '{alphauniprot}' and '{alphaseq}':")
# Note: This might seem pretty slow (~12 hours), but besides being an all-against-all comparison for 15 million rows it actually involves ranking structures by average pLDDT etc. and there isn't really a good way to speed it up.
inserted = 0
# for acc, species, tax, seq in Fetch(Query(f"SELECT acc, species, tax, seq FROM {alphauniprot} ORDER BY species='human' DESC, species, acc")):
# No longer need to order this since "alphauniprot.py -compara_species" now fetches the species in the compara_species order
# for acc, species, tax, seq in Fetch(Query(f"SELECT acc, species, tax, seq FROM {alphauniprot} WHERE acc='Q12753'")):
# Only map the 48 model & global health taxa:
# for acc, species, tax, seq in Fetch(Query(f"SELECT acc, species, tax, seq FROM {alphauniprot} WHERE tax IN (1352, 3702, 3847, 4577, 5671, 6183, 6239, 6248, 6279, 6282, 6293, 7227, 7955, 9606, 10090, 10116, 36087, 36329, 39947, 44689, 71421, 83332, 83333, 85962, 86049, 93061, 99287, 100816, 171101, 185431, 192222, 208964, 237561, 242231, 243232, 272631, 284812, 300267, 318479, 353153, 447093, 502779, 559292, 1125630, 1133849, 1299332, 1391915, 1442368)")):
# Map all:
for acc, species, tax, seq in Fetch(Query(f"SELECT acc, species, tax, seq FROM {alphauniprot}")):

    # Replace non-standard amino acids (B/Z/U/X) in sequence (for compatibility with AlphaFold)
    if rx("[BZUX]", seq):
        seq = ReplaceNonstandardAAs(seq)

    Log("total accs", acc)
    Log("total species", species)
    Log("total taxa", tax)
    Log("total sequences", seq)

    # Get all alphaseq accs that have this exact sequence (any species is fine - the only input to AlphaFold 2 is a sequence)
    query = Query(f"SELECT s.acc, s.afdb, AVG(a.plddt) AS avg_plddt FROM {alphaseq} s, {alphasa} a WHERE s.seq='{seq}' AND a.acc=s.acc AND a.afdb=s.afdb GROUP BY s.acc ORDER BY avg_plddt DESC, s.afdb ASC, s.acc")
    # query = structures.get(seq, [])
    if Numrows(query) == 0:
    # if len(query) == 0:
    
        # Insert NULL mapping into alphamap
        q = f"INSERT INTO {alphamap} SET type='{type}', version='{version}', species='{species}', tax='{tax}', value='{acc}', map=NULL, afdb=NULL, avg_plddt=NULL, best=NULL"
        if not Switch('debug'):
            Query(q)
        else:
            State(q)

        Log(f"sequence not found in {alphaseq} for acc (skipped)", acc)
        Log(f"sequence not found in {alphaseq} for seq (skipped)", seq)
        continue

    best = 1
    for alphacc, afdb, avg_plddt in query:
    
        # Insert mapping into alphamap
        q = f"INSERT INTO {alphamap} SET type='{type}', version='{version}', species='{species}', tax='{tax}', value='{acc}', map='{alphacc}', afdb='{afdb}', avg_plddt='{avg_plddt}', best='{best}'"
        if not Switch('debug'):
            Query(q)
        else:
            State(q)
        inserted += 1
        if best == 1:
            Log("successfully mapped for best=1 acc|alphacc", f"{acc}|{alphacc}")
        best = 0

        Log("successfully mapped for acc|alphacc", f"{acc}|{alphacc}")
        Log("successfully mapped for acc", acc)
        Log("successfully mapped for alphacc", alphacc)
        Log("successfully mapped for seq", seq)
        Log("successfully mapped for species", species)
        Log("successfully mapped for tax", tax)
        if afdb == 0:
            Log("successfully mapped to non-AFDB AlphaSync acc for acc|alphacc", f"{acc}|{alphacc}")
            Log("successfully mapped to non-AFDB AlphaSync acc for acc", acc)
            Log("successfully mapped to non-AFDB AlphaSync acc for alphacc", alphacc)
            Log("successfully mapped to non-AFDB AlphaSync acc for seq", seq)
            Log("successfully mapped to non-AFDB AlphaSync acc for species", species)
            Log("successfully mapped to non-AFDB AlphaSync acc for tax", tax)

# # Would show millions of log items
# Show(lim=0, sort=True)
Show(lim=50, sort=True)

State("Rows inserted: " + Comma(inserted))
State("Average number of AlphaFold DB structures per unique sequence: " + str(round(len(Get("successfully mapped for acc|alphacc")) / len(Get("successfully mapped for acc")), 2)))
State("Sequences that still need structure predictions for complete coverage: " + Comma(len(Get(f"sequence not found in {alphaseq} for seq (skipped)"))))

# Optimize(alphamap)

print("\nDone!")
