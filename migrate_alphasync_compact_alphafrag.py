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

# Note: alphafrag is not actually used by the AlphaSync website itself, but it will eventually be used by the API.

Starttime()
Query("""TRUNCATE TABLE alphasync_compact.alphafrag;""", loud=1)
Stoptime();


# Example with duplicated rows:
# Duplicate entry 'A0A6P3VU72-1-1'
# SELECT * FROM alphafrag WHERE acc='A0A6P3VU72' AND site='1' AND afdb='1';


Starttime()
Query("""INSERT INTO alphasync_compact.alphafrag SELECT acc, name, species, tax, frag, fragstart, fragstop, source, afdb, seq FROM alphasync.alphafrag;""")
Stoptime();

Starttime()
Query("""ANALYZE TABLE alphasync_compact.alphafrag;""")
Stoptime();

Starttime()
Query("""OPTIMIZE TABLE alphasync_compact.alphafrag;""")
Stoptime();


Show(lim=0)

print("\nDone!")

