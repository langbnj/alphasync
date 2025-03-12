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
query = Query("""SELECT DISTINCT acc, COUNT(*) AS c FROM alphacon GROUP BY acc, site1, site2, aa1, aa2, atom1, atom2, type, afdb HAVING c>1""")
for acc, c in query:
    print(f"{acc}\t{c}")
Stoptime();


# Example with duplicated rows:
# # Error Code: 1062. Duplicate entry 'A0A3B4ERA5-2-4-O-N-PolarHBondContacts-1' for key 'alphacon.PRIMARY'
# SELECT * FROM alphacon WHERE acc='A0A3B4ERA5' AND site1=2 AND site2=4 AND atom1='O' AND atom2='N' AND type='PolarHBondContacts';
# SELECT t.acc, t.c FROM (SELECT acc, COUNT(*) AS c FROM alphacon WHERE acc='A0A3B4ERA5' GROUP BY acc, site1, site2, aa1, aa2, atom1, atom2, type, afdb ORDER BY c DESC) t GROUP BY t.acc ORDER BY t.c DESC;
# SELECT t.acc, t.c FROM (SELECT acc, COUNT(*) AS c FROM alphacon GROUP BY acc, site1, site2, aa1, aa2, atom1, atom2, type, afdb ORDER BY c DESC) t GROUP BY t.acc ORDER BY t.c DESC;
# SELECT DISTINCT acc, COUNT(*) AS c FROM alphacon GROUP BY acc, site1, site2, aa1, aa2, atom1, atom2, type, afdb HAVING c>1 LIMIT 2;

# Complete list (result of this script):
# A0A2D0PRQ2	2
# A0A3B4ERA5	2
# A0A4W5K2R3	2
# A0A6P3VU72	2
# D3ZJX0	2
# 'A0A2D0PRQ2', 'A0A3B4ERA5', 'A0A4W5K2R3', 'A0A6P3VU72', 'D3ZJX0'

# SELECT DISTINCT acc, COUNT(*) AS c FROM alphacon WHERE acc IN ('A0A2D0PRQ2', 'A0A3B4ERA5', 'A0A4W5K2R3', 'A0A6P3VU72', 'D3ZJX0') GROUP BY acc, site1, site2, aa1, aa2, atom1, atom2, type, afdb;



Show(lim=0)

print("\nDone!")

