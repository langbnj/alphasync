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
Query("""TRUNCATE TABLE alphasync_compact.alphasa;""", loud=1)
Stoptime();

# Example with duplicated rows (test_remove_duplicates_alphasa.py):
# Duplicate entry 'A0A6P3VU72-1-1'
# SELECT * FROM alphasa WHERE acc='A0A6P3VU72' AND site='1' AND afdb='1';

Starttime()
# Query("""INSERT INTO alphasync_compact.alphasa SELECT DISTINCT acc, site, afdb, aa, plddt, plddt10, asa, asa10, relasa, relasa10, dis, dis10, surf, surf10, sec, membrane, iso, phi, psi, omega, chi1, chi2, chi3, chi4, chi5, tau FROM alphasync.alphasa;""")
# # DISTINCT shouldn't be necessary anymore, since I removed all duplicates now (test_remove_duplicates_alphasa.py):
# Using INSERT IGNORE instead (skipping duplicated rows):
# Query("""INSERT IGNORE INTO alphasync_compact.alphasa SELECT acc, site, afdb, aa, plddt, plddt10, asa, asa10, relasa, relasa10, dis, dis10, surf, surf10, sec, membrane, iso, phi, psi, omega, chi1, chi2, chi3, chi4, chi5, tau FROM alphasync.alphasa;""")
Query("""INSERT INTO alphasync_compact.alphasa SELECT acc, site, afdb, aa, plddt, plddt10, asa, asa10, relasa, relasa10, dis, dis10, surf, surf10, sec, iso, phi, psi, omega, chi1, chi2, chi3, chi4, chi5, tau FROM alphasync.alphasa;""")
Stoptime();

Starttime()
Query("""ANALYZE TABLE alphasync_compact.alphasa;""")
Stoptime();

Starttime()
Query("""OPTIMIZE TABLE alphasync_compact.alphasa;""")
Stoptime();


Show(lim=0)

print("\nDone!")

