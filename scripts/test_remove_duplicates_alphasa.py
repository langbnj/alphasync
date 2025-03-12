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

# Note: Should no longer occur

Starttime()
query = Query("""SELECT DISTINCT acc, COUNT(*) AS c FROM alphasa GROUP BY acc, site, afdb HAVING c>1""")
for acc, c in query:
    print(f"{acc}\t{c}")
Stoptime();


# Example with duplicated rows:
# Duplicate entry 'A0A6P3VU72-1-1'
# SELECT * FROM alphasa WHERE acc='A0A6P3VU72' AND site='1' AND afdb='1';


# Complete list (result of this script):

# Probably (from alphacon):
# A0A2D0PRQ2	2
# A0A6P3VU72	2



Show(lim=0)

print("\nDone!")

