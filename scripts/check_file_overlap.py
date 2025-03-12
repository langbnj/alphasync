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
# from blang_mysql import *
from blang import *

files = Return("ls -1 input/ftp/*.tar")
l = files.split("\n")

Run("Make directory for file lists", "mkdir -p input/diagnostic_lists/")
print("Getting file lists for all TAR files")
for f in tq(l):
    listfile = f"input/diagnostic_lists/{Basename(f)}.files.txt"
    if not Exists(listfile):
        Run("Get file list", fr"tar -tf {f} | perl -npe 's/\.pdb\.gz//' | sort | uniq > {listfile}")


print("Comparing contents of *.tar files in input/ftp/ (not showing pairs with zero overlap):")
files = Return("ls -1 input/diagnostic_lists/*.tar.files.txt")
l = files.split("\n")
for f1 in l:

    # if f == "input/swissprot_cif_v2.tar.files.txt":
    #     continue

    for f2 in l:

        if f1.lower() >= f2.lower():
            continue

        overlap = Return(f"comm -12 {f1} {f2} | wc -l")

        if overlap != "0":
            # print(f" >> {f1} vs. {f2} >> overlap {overlap}")
            print(f"{f1}\t{f2}\t{overlap}")
    
# d()

Show(lim=10)

print("\nDone!")
