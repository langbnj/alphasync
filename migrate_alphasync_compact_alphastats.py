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

source = "alphasync.alphastats"
dest = "alphasync_compact.alphastats"

# # Args(0, "", "")
# var = Args(1, "[species]", "human")
# infile = f"input/input.txt"

# Run queries

# Simply copy this table, since it's already in InnoDB format in schema blang

Starttime()
Query(f"""DROP TABLE IF EXISTS {dest};""", loud=1)
Stoptime();

Starttime()
Query(f"""CREATE TABLE {dest} LIKE {source};""")
Stoptime();

# Copy table
Starttime()
Query(f"""INSERT INTO {dest} SELECT * FROM {source};""")
Stoptime();

Starttime()
Query(f"""ANALYZE TABLE {dest};""")
Stoptime();

Starttime()
Query(f"""OPTIMIZE TABLE {dest};""")
Stoptime();


Show(lim=0)

print("\nDone!")

