#!/usr/bin/env python3
"""
alphamap_ensembl.py:
Uses perfect sequence matching to map identifiers for a given version of Ensembl to accessions in AlphaSync tables such as 'alphaseq'.
"""

# Initialize
# import pandas as pd
# import numpy as np
# import networkx as nx
# import seaborn as sns
# import matplotlib.pyplot as mp
import tarfile
import gzip
import io
# from Bio import SeqIO
from blang_mysql import *
from blang import *

type = "ensembl"

alphamap = "alphamap"
ensembl = "ensembl"

(version) = Args(1, "[Current Ensembl version] [-comparaonly]\n -comparaonly: Only use proteins from the ~200 species contained in Ensembl Compara (specifically, the comparafasta_… tables, as well as all human sequences)\n -debug: Don't actually make any changes, just simulate", "108 -comparaonly")

if not Switch('debug'):
    # Clear(alphamap)
    print(f"Deleting type '{type}' version '{version}' rows from table '{alphamap}'...")
    query = Query(f"DELETE FROM {alphamap} WHERE type='{type}' AND version='{version}'")
    print(f"Rows affected: {Comma(Numrows(query))}")

# Start

# Get set of wanted sequences and species
if Switch("comparaonly"):
    # Only sequences required for Ensembl Compara (as well as all human sequences)
    print(f"\nGetting set of wanted sequences from comparafasta_… tables (excluding comparafasta itself):")
    wanted_seqs = set()
    for comparafasta in FetchList(Query("SHOW TABLES LIKE 'comparafasta\\_%'")):
        print(f" >> {comparafasta}", end="")
        # print(f" >> {comparafasta}")
        tmpseqs = FetchSet(Query(f"SELECT REPLACE(seq, '-', '') FROM {comparafasta}"))
        print(f" >> {Comma(len(tmpseqs))} sequences", end="")
        # print(f"   >> {Comma(len(tmpseqs))} sequences")
        # Add tmpseqs to wanted_seqs
        wanted_seqs |= tmpseqs
        print(f" >> total {Comma(len(wanted_seqs))} sequences")
        # print(f"     >> total {Comma(len(wanted_seqs))} sequences")
    print(f" >> total {Comma(len(wanted_seqs))} sequences (plus all human sequences)")

    # # Add uniens sequences
    # print(f"\nGetting uniens sequences (successfully mapped between UniProtKB/Swiss-Prot and Ensembl)")
    # FetchSet(Query("SELECT e.seq FROM uniens ue, ensembl e WHERE ue.ensp=e.ensp"))
    # d()
    # Not necessary since many of these won't be in Compara anyway

    # Only species required for Ensembl Compara
    print(f"\nGetting set of wanted species from compara_species:")
    wanted_species = FetchSet(Query("SELECT DISTINCT species FROM compara_species"))
    print(f" >> total {Comma(len(wanted_species))} species")

State("Getting mapping from Ensembl fullspecies names to NCBI taxon IDs from table 'ensembl_species':")
tax = {}
for species, thistax in Fetch(Query(f"SELECT fullspecies, tax FROM ensembl_species")):
    tax[species] = thistax

State(f"Filling table '{alphamap}' with mappings for Ensembl version '{version}' based on perfect sequence matching in table '{ensembl}':")
for ensp, unispec, species, seq in Fetch(Query(f"SELECT ensp, species, LOWER(fullspecies), seq FROM {ensembl} ORDER BY species='human' DESC, fullspecies, ensp")):

    # Check if this sequence is in the wanted_seqs set or human (and skip otherwise)
    if Switch('comparaonly'):
        if species not in wanted_species or seq not in wanted_seqs:
        # # In addition to wanted_seqs, also retain all human sequences:
        # if tax[species] not in wanted_taxa or (seq not in wanted_seqs and unispec != 'HUMAN'):
            Log("skipped unwanted sequence for ensp", ensp)
            Log("skipped unwanted sequence for seq", seq)
            Log("skipped unwanted sequence for species", species)
            continue

    Log("total ensps", ensp)
    Log("total species", species)
    Log("total taxa", tax[species])
    Log("total sequences", seq)

    # Get all alphaseq accs that have this exact sequence (we don't care about the species here)
    query = Query(f"SELECT s.acc, s.afdb, AVG(a.plddt) AS avg_plddt FROM alphaseq s, alphasa a WHERE s.seq='{seq}' AND a.acc=s.acc AND a.afdb=s.afdb GROUP BY s.acc ORDER BY avg_plddt DESC, s.afdb ASC, s.acc")
    if Numrows(query) == 0:
        # Insert NULL mapping into alphamap
        Query(f"INSERT INTO {alphamap} SET type='{type}', version='{version}', species='{species}', tax='{tax[species]}', value='{ensp}', map=NULL, afdb=NULL, avg_plddt=NULL, best=NULL")

        Log("sequence not found in alphaseq for ensp (skipped)", ensp)
        Log("sequence not found in alphaseq for seq (skipped)", seq)
        # Check if sequence contains non-AA characters
        if not Aa(seq):
            Log("sequence not found in alphaseq & contains non-AA characters for ensp (skipped)", ensp)
            Log("sequence not found in alphaseq & contains non-AA characters for seq (skipped)", seq)
        continue
    best = 1
    for alphacc, afdb, avg_plddt in query:

        # Insert mapping into alphamap
        Query(f"INSERT INTO {alphamap} SET type='{type}', version='{version}', species='{species}', tax='{tax[species]}', value='{ensp}', map='{alphacc}', afdb='{afdb}', avg_plddt='{avg_plddt}', best='{best}'")
        best = 0

        Log("successfully mapped for ensp|alphacc", f"{ensp}|{alphacc}")
        Log("successfully mapped for ensp", ensp)
        Log("successfully mapped for alphacc", alphacc)
        Log("successfully mapped for seq", seq)
        Log("successfully mapped for species", species)
        Log("successfully mapped for tax", tax[species])
        if afdb == 0:
            Log("successfully mapped to non-AFDB AlphaSync acc for ensp|alphacc", f"{ensp}|{alphacc}")
            Log("successfully mapped to non-AFDB AlphaSync acc for ensp", ensp)
            Log("successfully mapped to non-AFDB AlphaSync acc for alphacc", alphacc)
            Log("successfully mapped to non-AFDB AlphaSync acc for seq", seq)
            Log("successfully mapped to non-AFDB AlphaSync acc for species", species)
            Log("successfully mapped to non-AFDB AlphaSync acc for tax", tax[species])

Show(lim=0, sort=True)

State("Average number of AlphaFold DB structures per unique sequence: " + str(round(len(Get("successfully mapped for ensp|alphacc")) / len(Get("successfully mapped for ensp")), 2)))
State("Sequences that still need structure predictions for complete coverage: " + Comma(len(Get("sequence not found in alphaseq for seq (skipped)"))))

# Optimize(alphamap)

print("\nDone!")
