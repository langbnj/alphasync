#!/usr/bin/env python3
"""
alphauniprot_symbols.py:
Fill AlphaSync table 'alphauniprot_symbols', a normalised mapping of UniProt gene symbols and synonyms to UniProt accessions for quick lookup.
"""

# Initialize
# import pandas as pd
# import numpy as np
# import networkx as nx
# import seaborn as sns
# import matplotlib.pyplot as mp
from blang_mysql import *
from blang import *

alphauniprot = "alphauniprot"
alphauniprot_symbols = "alphauniprot_symbols"

# Args(0, "", "")
# var = Args(1, "[species]", "human")
# infile = f"input/input.txt"

# Clear table
Clear(alphauniprot_symbols)

# Start
print(f"\nFilling table '{alphauniprot_symbols}' with (not necessarily unambiguous) gene symbols and synonyms from '{alphauniprot}'...")
# Insert gene symbols and synonyms
for acc, species, tax, symbols, synonyms in Fetch(Query(f"SELECT acc, species, tax, symbols, synonyms FROM {alphauniprot} WHERE symbols IS NOT NULL OR synonyms IS NOT NULL")):
    if symbols:
        for symbol in symbols.split("|"):
            Query(f"INSERT IGNORE INTO {alphauniprot_symbols} SET acc='{acc}', species='{species}', tax='{tax}', type='symbol', alias='{Esc(symbol.strip())}'")
            Log("inserted total symbols", symbol.strip())
            Log("inserted symbol for acc", acc)
            Log("inserted total entries", symbol.strip())
            Log("inserted entry for acc", acc)
            Log("inserted entry for species", species)
            Log("inserted entry for tax", tax)
    if synonyms:
        for synonym in synonyms.split("|"):
            Query(f"INSERT IGNORE INTO {alphauniprot_symbols} SET acc='{acc}', species='{species}', tax='{tax}', type='synonym', alias='{Esc(synonym.strip())}'")
            Log("inserted total synonyms", synonym.strip())
            Log("inserted synonym for acc", acc)
            Log("inserted total entries", synonym.strip())
            Log("inserted entry for acc", acc)
            Log("inserted entry for species", species)
            Log("inserted entry for tax", tax)

Show(lim=0)

Optimize(alphauniprot_symbols)

print("\nDone!")
