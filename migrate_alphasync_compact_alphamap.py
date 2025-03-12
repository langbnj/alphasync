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

# table = "table"

# # Args(0, "", "")
# var = Args(1, "[species]", "human")
# infile = f"input/input.txt"

# Run queries

Starttime()
Query("""TRUNCATE TABLE alphasync_compact.alphamap;""", loud=1)
Stoptime();

# Example with duplicated rows:
# Duplicate entry 'A0A6P3VU72-1-1'
# SELECT * FROM alphamap WHERE acc='A0A6P3VU72' AND site='1' AND afdb='1';


# Query("""INSERT INTO alphasync_compact.alphamap SELECT DISTINCT type, version, value, best, afdb, map, species, tax, avg_plddt FROM alphasync.alphamap;""")
# Skip unmapped sequences (map IS NULL)
Starttime()
# UniProt
type = "uniprot"
version = get_local_uniprot_release()
# Query(f"""INSERT INTO alphasync_compact.alphamap SELECT DISTINCT type, version, value, best, afdb, map, species, tax, avg_plddt FROM alphasync.alphamap WHERE map IS NOT NULL AND type='{type}' AND version='{version}';""")
Query(f"""INSERT INTO alphasync_compact.alphamap SELECT DISTINCT type, version, value, best, afdb, map, species, tax, avg_plddt FROM alphasync.alphamap WHERE map IS NOT NULL AND type='{type}' AND version='{version}' AND best=1;""")
# # Ensembl
# type = "ensembl"
# version = "112"
# Query(f"""INSERT INTO alphasync_compact.alphamap SELECT DISTINCT type, version, value, best, afdb, map, species, tax, avg_plddt FROM alphasync.alphamap WHERE map IS NOT NULL AND type='{type}' AND version='{version}';""")
Stoptime();

Starttime()
Query("""ANALYZE TABLE alphasync_compact.alphamap;""")
Stoptime();

Starttime()
Query("""OPTIMIZE TABLE alphasync_compact.alphamap;""")
Stoptime();


Show(lim=0)

print("\nDone!")

