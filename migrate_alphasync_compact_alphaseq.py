#!/usr/bin/env python3
"""
Main: Main script
"""

# Initialize
# import pandas as pd
# import numpy as np
# import networkx as nx
# import seaborn as sns
# import matplotlib.pyplot as mp
from blang_mysql import *
from blang import *

source = "alphasync.alphaseq"
dest = "alphasync_compact.alphaseq"
# alphamap = "alphasync_compact.alphamap"
alphamap = "alphasync.alphamap"

# # Args(0, "", "")
# var = Args(1, "[species]", "human")
# infile = f"input/input.txt"

# Run queries

Starttime()
Query(f"""TRUNCATE TABLE {dest};""", loud=1)
Stoptime();

# Indexes needed for AlphaSync website queries for alphaseq: tax, acc, name, species, seq

# Example with duplicated rows:
# Duplicate entry 'A0A6P3VU72-1-1'
# SELECT * FROM alphaseq WHERE acc='A0A6P3VU72' AND site='1' AND afdb='1';

Starttime()
# UniProt
type = "uniprot"
version = get_local_uniprot_release()
# Query("""INSERT INTO {dest} SELECT acc, name, species, tax, frags, afdb, seq FROM alphasync.alphaseq;""")
# Only get UniProt proteins that are mapped in table 'alphamap' (that have structures) and that are the best match (according to average pLDDT) for this sequence
Query(f"""INSERT INTO {dest} SELECT acc, name, species, tax, frags, afdb, seq FROM {source} WHERE acc IN (SELECT DISTINCT map FROM {alphamap} WHERE type='{type}' AND version='{version}' AND map IS NOT NULL AND best=1);""")
Stoptime();

Starttime()
Query(f"""ANALYZE TABLE {dest};""")
Stoptime();

Starttime()
Query(f"""OPTIMIZE TABLE {dest};""")
Stoptime();


Show(lim=0)

print("\nDone!")

