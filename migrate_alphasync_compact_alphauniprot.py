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

source = "alphasync.alphauniprot"
dest = "alphasync_compact.alphauniprot"
# alphamap = "alphasync_compact.alphamap"
alphamap = "alphasync.alphamap"

# # Args(0, "", "")
# var = Args(1, "[species]", "human")
# infile = f"input/input.txt"

# Run queries

# Simply copy this table, since it's already in InnoDB format in schema blang

Starttime()
Query(f"""DROP TABLE IF EXISTS {dest};""", loud=1)
Stoptime();

# Indexes needed for AlphaSync website queries for alphauniprot: tax, acc, name, species, seq
Starttime()
Query(f"""CREATE TABLE {dest} LIKE {source};""")
Stoptime();

Starttime()
# UniProt
type = "uniprot"
version = get_local_uniprot_release()
# Query(f"""INSERT INTO {dest} SELECT * FROM {source};""")
# Only get UniProt proteins that are mapped in table 'alphamap' (that have structures) and that are the best match (according to average pLDDT) for this sequence
Query(f"""INSERT INTO {dest} SELECT * FROM {source} WHERE acc IN (SELECT DISTINCT value FROM {alphamap} WHERE type='{type}' AND version='{version}' AND map IS NOT NULL AND best=1);""")
Stoptime();

Starttime()
Query(f"""ANALYZE TABLE {dest};""")
Stoptime();

Starttime()
Query(f"""OPTIMIZE TABLE {dest};""")
Stoptime();


Show(lim=0)

print("\nDone!")

