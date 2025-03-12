#!/usr/bin/env python3
"""
Extract FASTA sequence from PDB file. Ultimately not used (too slow).
"""

# Initialize
# import pandas as pd
# import numpy as np
# import networkx as nx
# import seaborn as sns
# import matplotlib.pyplot as mp
# from Bio.PDB import *
from Bio import SeqIO
# from blang_mysql import *
from blang import *

# table = "table"

# infile = Args(1, "[file]", "../alphafold/test/input/UP000005640_9606_HUMAN_v2/AF-Q15008-F1-model_v2.cif.gz")
infile = Args(1, "[file]", "../alphafold/test/ranked_0.pdb")
# infile = f"input/input.txt"

# Start
print(f"\nReading '{infile}' and extracting FASTA sequence:")

# Define function to extract sequence from PDB file
# https://biopython.org/wiki/SeqIO:
# def pdbfile_to_seq(infile, format="pdb-seqres"):        # PDB: doesn't work with our local AlphaFold output, which doesn't have any header (although it does work with the AlphaFold DB PDB files, which do)
def pdbfile_to_seq(infile, format="pdb-atom"):        # PDB: works
# def pdbfile_to_seq(infile, format="cif-seqres"):      # CIF: works
# def pdbfile_to_seq(infile, format="cif-atom"):        # CIF: works

    # Filter out a warning that occurs in all AlphaFold PDB result files
    warnings.filterwarnings("ignore", message="'HEADER' line not found; can't determine PDB ID.")
    # Filter out a warning that occurs in all AlphaFold mmCIF result files
    warnings.filterwarnings("ignore", message="Could not determine the PDB ID.")

    with open(infile) as f:
        # for record in SeqIO.parse(f, 'pdb-atom'):
        #     print(record.seq)
    
        # Get "records" (chains)
        records = list(SeqIO.parse(f, format))
    
        # Check number of records
        if len(records) != 1:
            Die(f"Expected 1 record, but got {len(records)} records in '{infile}'")
    
        seq = records[0].seq
    
        return(seq)

# Extract sequence from input file
seq = pdbfile_to_seq(infile)
print(seq)




Show(lim=10)

print("\nDone!")
