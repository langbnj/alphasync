#!/usr/bin/env python3
"""
Fills the 'membrane' column in table 'alphasa', flagging residues that, according to UniProt (table unifeat), are within 'transmembrane region' or 'intramembrane region' features.
"""

# Initialize
# import pandas as pd
# import numpy as np
# import networkx as nx
# import seaborn as sns
# import matplotlib.pyplot as mp
from blang_mysql import *
from blang import *

alphasa = "alphasa"

# Args(0, "", "")
# var = Args(1, "[species]", "human")
# infile = f"input/input.txt"

# Start
print(f"\nClearing 'membrane' column in table '{alphasa}'...")

# # First, initialise column by setting it to NULL
# query = Query(f"UPDATE {alphasa} SET membrane=NULL")
# print(f"{Comma(Numrows(query))} rows affected")

# Start
print(f"\nFilling the 'membrane' column in table '{alphasa}', flagging residues that, according to UniProt (table unifeat), are within 'transmembrane region' or 'intramembrane region' features:")

# Cycle through all unifeat features that are 'transmembrane region' or 'intramembrane region'
# for (acc, start, stop) in tq(Query(f"SELECT acc, start, stop FROM unifeat WHERE description IN ('transmembrane region', 'intramembrane region')")):
for (uniacc, seq, start, stop) in tq(Query(f"SELECT f.acc, s.seq, f.start, f.stop FROM unifeat f, uniseq s WHERE f.type IN ('transmembrane region', 'intramembrane region') AND f.acc=s.acc AND s.type IN ('UniProt', 'UniIso')")):
    Log("total membrane features for uniacc|start|stop", f"{uniacc}|{start}|{stop}")
    Log("total uniaccs", uniacc)
    Log("total sequences", seq)

    # Get all alphaseq accs that have this exact sequence (we don't care about the species here)
    for (acc,) in Query(f"SELECT DISTINCT acc FROM alphaseq WHERE seq='{seq}'"):
        Log("total accs", acc)
        
        # Set membrane column in table 'alphasa' to 1 for each residue of this feature
        Query(f"UPDATE {alphasa} SET membrane=1 WHERE acc='{acc}' AND site BETWEEN {start} AND {stop}")
        Log("successfully set membrane column 1 for residues of acc", acc)

        # TODO
        # Ideally, this script would also set orthologous residues' membrane values to 0.
        # For now, I can simply check what the status of the human PTM site is and assume it's the same for others.

Show(lim=0)

print("\nCounting residues that have membrane=1...")
affected = FetchOne(Query(f"SELECT COUNT(*) FROM {alphasa} WHERE membrane=1"))
print(f"{Comma(affected)} rows affected")

print("\nDone!")
