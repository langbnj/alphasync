#!/usr/bin/env python3
"""
Fills the "best" column in table "alphaseq", flagging the best median pLDDT structure for each unique sequence (since AlphaFold DB frequently contains multiple predictions for the same exact sequence)
"""

# Initialize
# import pandas as pd
# import numpy as np
# import networkx as nx
# import seaborn as sns
# import matplotlib.pyplot as mp
from blang_mysql import *
from blang import *

alphaseq = "alphaseq"

# Args(0, "", "")
# var = Args(1, "[species]", "human")
# infile = f"input/input.txt"

# Start
# TODO alphamap_….py currently uses average pLDDT to define best - use median instead? However, average is easiest to calculate in SQL.
print(f"\nFilling the 'best' column in table '{alphaseq}', flagging the best median pLDDT structure for each unique sequence (since AlphaFold DB frequently contains multiple predictions for the same exact sequence):")

# First, initialise column by setting it to NULL
Query(f"UPDATE {alphaseq} SET best=NULL")

# Then, set column to 1 for all sequences that only occur once
Query(f"UPDATE {alphaseq} s, (SELECT seq, COUNT(*) AS c FROM {alphaseq} GROUP BY seq HAVING c = 1) t SET s.best=1 WHERE s.seq=t.seq")

# Then, set column to 0 for all sequences that occur more than once
Query(f"UPDATE {alphaseq} s, (SELECT seq, COUNT(*) AS c FROM {alphaseq} GROUP BY seq HAVING c > 1) t SET s.best=0 WHERE s.seq=t.seq")

# Cycle through all sequences that occur more than once
# for (seq,) in tq(Query(f"SELECT DISTINCT seq FROM {alphaseq}")):
for (seq, count) in tq(Query(f"SELECT seq, COUNT(*) AS c FROM {alphaseq} GROUP BY seq HAVING c > 1")):
    # Get pLDDT values for each accession that has this sequence and choose best median
    # TODO finish here
    pass

Show(lim=0)

print("\nDone!")
