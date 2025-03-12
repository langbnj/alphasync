#!/usr/bin/env python3
"""
Get the list of UniProt accessions for which AlphaFold needs to be re-run in order to synchronise AlphaFold DB with the latest UniProt release.

Use table 'alphamap' to find sequences for which no matching AlphaFold DB sequence was found.
Uses sequence sources such as the current 'alphauniprot' table, obtained from UniProt via API.
"""

# Note: Now superseded by alphasync.py

# Initialize
# import pandas as pd
# import numpy as np
# import networkx as nx
# import seaborn as sns
# import matplotlib.pyplot as mp
from blang_mysql import *
from blang import *

# Set SQL table names to use
alphamap = "alphamap"
alphauniprot = "alphauniprot"
alphaensembl = "alphaensembl"
# alphaseq = "alphaseq_v2"
alphaseq = "alphaseq"
uniseq = "uniseq"
ensembl = "ensembl"
compara = "comparafasta_einsi_tree_1para"

# infile = Args(1, "[input sequence file]", "../alphafold_db/input/UP000005640_9606_HUMAN_v2.tar.seqs.txt")
(type, rel) = Args(2, "[sequence database to sync to: alphauniprot/alphaensembl/uniprot/ensembl/compara] [release: e.g. 2024_05/112]", "alphauniprot 2024_05\nalphaensembl 112")

# Parse database type
if (type == "alphauniprot"):
    table = alphauniprot
elif (type == "alphaensembl"):
    table = alphaensembl
elif (type == "uniprot"):
    table = uniseq
elif (type == "ensembl"):
    table = ensembl
elif (type == "compara"):
    table = compara
else:
    Die("Unhandled database type: expecting 'alphauniprot', 'alphaensembl', 'uniprot', 'ensembl', or 'compara'")

# Parse release number
if (type == "uniprot"):
    if (rx(r"^20\d{2}_\d{2}$", rel) == False):
        Die(f"Unexpected '{type}' release number: expecting e.g. '2024_05'")
elif (type == "alphaensembl"):
    if (rx(r"^\d{2,3}$", rel) == False):
        Die(f"Unexpected '{type}' release number: expecting e.g. '112'")
elif (type == "ensembl"):
    if (rx(r"^\d{2,3}$", rel) == False):
        Die(f"Unexpected '{type}' release number: expecting e.g. '108'")
elif (type == "compara"):
    if (rx(r"^\d{2,3}$", rel) == False):
        Die(f"Unexpected '{type}' release number: expecting e.g. '108'")

# Start
print(f"\nGetting unique '{type}' release '{rel}' sequences from table '{table}'...")
Time(1)
if (type == "alphauniprot"):
    seqs = FetchSet(Query(f"SELECT seq FROM {table}"))
if (type == "alphaensembl"):
    seqs = FetchSet(Query(f"SELECT seq FROM {table}"))
if (type == "uniprot"):
    # SELECT type, COUNT(DISTINCT seq) AS seqs FROM uniseq WHERE type IN ('UniProt', 'UniIso') GROUP BY type;
    # UniProt	480281
    # UniIso	66029
    # >> UniIso is just a small subset (not many species have a significant number of isoform sequences), including it
    # All UniProt canonical and isoform sequences
    seqs = FetchSet(Query(f"SELECT seq FROM {table} WHERE type IN ('UniProt', 'UniIso')"))
    # Human only:
    # seqs = FetchSet(Query(f"SELECT seq FROM {table} WHERE type IN ('UniProt', 'UniIso') AND species='human'"))
elif (type == "ensembl"):
    # All Ensembl sequences
    seqs = FetchSet(Query(f"SELECT seq FROM {table}"))
elif (type == "compara"):
    # Only those sequences actually used in Ensembl Compara
    seqs = FetchSet(Query(f"SELECT seq FROM {table}"))
    # Remove gaps (-) from sequences
    seqs = [seq.replace("-", "") for seq in seqs]
Time(2)

print(f"\nGetting unique AlphaFold DB sequences from table '{alphaseq}':")
alphaseqs = FetchSet(Query(f"SELECT seq FROM {alphaseq}"))
# Human only:
# alphaseqs = FetchSet(Query(f"SELECT seq FROM {alphaseq} WHERE species='human'"))

# Get all accs that are in seqs, but not in alphaseqs
seqs_both = set(seqs) & set(alphaseqs)
seqs_only = set(seqs) - set(alphaseqs)

print(f"\n{Comma(len(seqs_both))} sequences are in both '{type}' and AlphaFold DB")
print(f"{Comma(len(seqs_only))} sequences are in '{type}' but not in AlphaFold DB")

Show(lim=10)

print("\nDone!")
