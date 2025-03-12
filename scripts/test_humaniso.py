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

# # Args(0, "", "")
# var = Args(1, "[species]", "human")
infile = f"test_humaniso/chess_structure_CDS_aa_v1.2.faa"

# uniprot_version = "2024_06"
uniprot_version = "2025_01"
# uniprot_version = get_local_uniprot_release()

# Load human CHESS sequences from FASTA
print(f"Loading human CHESS sequences from FASTA file '{infile}'...")
# seqs = GetFasta(infile)
seqs = set()
with open(infile) as f:
    for title, seq in Fasta(f):
        # print(f">{title}\n{seq}")
        # seqs[title] = seq
        seqs.add(seq)
print(f" >> Loaded {len(seqs):,} sequences")

# Check how many of the human isoforms we want are in this FASTA file
# alphasync.py uniprot 2024_06 -humanonly -iso -no_u -no_c -no_z -no_x -nofrag # (human isoforms)
print(f"Checking how many of the human isoforms that we want are in this FASTA file...")
# for acc, seq in Fetch(Query(f"SELECT m.value, au.seq FROM alphamap m, alphauniprot au WHERE m.value=au.acc AND (au.reviewed=1 OR au.refprotcanon=1) AND m.type='uniprot' AND m.version='{uniprot_version}' AND m.map IS NULL AND m.tax=9606 HAVING LENGTH(au.seq)>=16 AND LENGTH(au.seq)<=2699 AND au.seq NOT REGEXP '[UBZX]' ORDER BY m.tax, RAND(42)")):
# Get all the isoforms (these are now mapped, of course, since I've completed the human isoforms):
for acc, seq in Fetch(Query(f"SELECT m.value, au.seq FROM alphamap m, alphauniprot au WHERE au.species='human' AND m.value=au.acc AND (au.reviewed=1 OR au.refprotcanon=1) AND acc!=canon AND m.type='uniprot' AND m.version='{uniprot_version}' AND m.map IS NOT NULL AND m.tax=9606")):
    if seq in seqs:
        Log(f"sequence in humaniso for acc", acc)
        Log(f"sequence in humaniso for seq", seq)
    else:
        Log(f"sequence not in humaniso for acc", acc)
        Log(f"sequence not in humaniso for seq", seq)

Show()

print("\nDone!")

# 2024_06:
# # Conclusion:
# # 12,576 human isoforms that we want are in the human CHESS FASTA file.
# # However, another 7,277 are missing.
# # This means only 12,576 / (12,576 + 7,277) = 12,576 / 19,853 = 63.3% of the human isoforms that we want are in the human CHESS FASTA file.
# # Since I'm already 14,000 CPU structures and 7,000 GPU structures into the human isoforms (out of 20,000), I'll skip this dataset and simply proceed with AlphaSync's pipeline for consistency.

# 2025_01:
# Conclusion:
# 14,379 human isoform accs that we want are in the human CHESS FASTA file.
# However, another 7,722 are missing.
# This means only 14,379 / (14,379 + 7,722) = 14,379 / 22,101 = 65.1% of the human isoforms that we want are in the human CHESS FASTA file.
# >> Sommer et al. covers 65.1% of sequences needed.
