#!/usr/bin/env python3
"""
alphaseq.py:
Combine AlphaFold fragment sequences from SQL table 'alphafrag' into SQL table 'alphaseq'.
This table serves as a list of sequences of interest for which DSSP (alphasa, accessible surface area) and Lahuta (alphacon, contacts) calculations will be made.
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

alphafrag = "alphafrag"     # SQL table with fragment protein sequences (>2700 aa proteins get split into 1400 aa fragments with a step size of 200 in AlphaFold DB, for human only - other species don't have results for >2700 aa proteins)
alphaseq = "alphaseq"       # SQL table with complete protein sequences

# Fragment length and step size used by DeepMind.
# Proteins longer than 2700 residues are split into windows of width 1400 with a step size of 200.
fraglen = 1400
fragstep = 200

Args(0, " \n -alphasync: Syncing: only refresh updated proteins from AlphaSync (those with afdb=0)\n -comparaonly: Only use proteins from the ~200 species contained in Ensembl Compara (specifically, the comparafasta_… tables, as well as all 'model organism' and 'global health proteomes' sequences from the AlphaFold Protein Structure Database)\n -debug: Don't actually make any changes, just simulate", " -comparaonly")



# Start

# # Get set of sequences already in alphaseq (could also simply clear the table first instead, but this makes the pipeline robust to re-running the script)
# Leads to issues where subsequent accs with the same sequence will be skipped, and therefore incomplete data. Disabling.
# print(f"\nGetting sequences already in table 'alphaseq' (to be skipped):")
# seqs_already_in_alphaseq = FetchSet(Query(f"SELECT seq FROM alphaseq"))
# print(f" >> {Comma(len(seqs_already_in_alphaseq))} sequences")

if not Switch('debug') and not Switch('debug2'):
    if Switch('alphasync'):
        # If syncing (switch -alphasync active): only clear AlphaSync fragments from table 'alphaseq'
        print(f"Clearing AlphaSync fragments from table '{alphaseq}'...")
        query = Query(f"DELETE FROM {alphaseq} WHERE afdb=0")
        print(f" >> {Comma(Numrows(query))} rows affected")
    else:
        # Clear entire table
        Clear(alphaseq)

# Get set of wanted sequences and taxa
# if Switch("comparaonly"):

# Only sequences required for Ensembl Compara (as well as all human sequences)
print(f"\nGetting set of wanted sequences from comparafasta_… tables (excluding comparafasta itself, which contains all paralogs and lower-quality orthologs):")
wanted_seqs = set()
for comparafasta in FetchList(Query("SHOW TABLES LIKE 'comparafasta\\_%'")):
    print(f" >> {comparafasta}", end="")
    # print(f" >> {comparafasta}")
    tmpseqs = FetchSet(Query(f"SELECT REPLACE(seq, '-', '') FROM {comparafasta}"))
    print(f" >> {Comma(len(tmpseqs))} sequences", end="")
    # print(f"   >> {Comma(len(tmpseqs))} sequences")
    wanted_seqs.update(tmpseqs)
    print(f" >> total {Comma(len(wanted_seqs))} sequences")
    # print(f"     >> total {Comma(len(wanted_seqs))} sequences")
# print(f" >> total {Comma(len(wanted_seqs))} sequences (plus all human sequences)")
print(f" >> total {Comma(len(wanted_seqs))} sequences (plus all human, model organism and global health proteome sequences)")

# # Add uniens sequences
# print(f"\nGetting uniens sequences (successfully mapped between UniProtKB/Swiss-Prot and Ensembl)")
# FetchSet(Query("SELECT e.seq FROM uniens ue, ensembl e WHERE ue.ensp=e.ensp"))
# d()
# Not necessary since many of these won't be in Compara anyway

# Only taxa required for Ensembl Compara
print(f"\nGetting set of wanted taxa from compara_species:")
wanted_taxa = FetchSet(Query("SELECT DISTINCT tax FROM compara_species"))
print(f" >> {Comma(len(wanted_taxa))} taxa")

# ...plus "Model organisms" and "Global health proteomes" from AlphaFold Protein Structure Database
print(f"\nGetting set of completely wanted taxa (model organisms and global health proteomes):")
wanted_complete_taxa = FetchSet(Query("SELECT DISTINCT tax FROM alphafrag WHERE source IN (SELECT DISTINCT source FROM alphafrag HAVING source LIKE 'UP%');"))
print(f" >> {Comma(len(wanted_complete_taxa))} taxa")
# Add wanted_complete_taxa to wanted_taxa
wanted_taxa.update(wanted_complete_taxa)
print(f" >> total {Comma(len(wanted_taxa))} taxa")



# print(f"\nWanted sequences that aren't already in alphaseq:")
# print(f" >> {Comma(len(wanted_seqs - seqs_already_in_alphaseq))} sequences")

# Main query: Loop over individual proteins
# This shows that grouping by acc alone is sufficient:
# SELECT COUNT(DISTINCT acc), COUNT(DISTINCT acc, name), COUNT(DISTINCT acc, name, species, tax) FROM alphafrag;
Starttime()
print(f"\nGetting fragment sequences from table '{alphafrag}' and combining them into complete protein sequences in table '{alphaseq}':")
# wanted_seqs = ("MVVSAGPWSSEKAEMNILEINETLRPQLAEKKQQFRNLKEKCFLTQLAGFLANRQKKYKYEECKDLIKFMLRNERQFKEEKLAEQLKQAEELRQYKVLFHSQERELTQLREKLREGRDASRSLNEHLQALLTPDEPDKSQGQDLQEQLAEGCRLAQHLVQKLSPENDNDDDEDVQVEVAEKVQKSSAPREMQKAEEKEVPEDSLEECAITYSNSHGPYDSNQPHRKTKITFEEDKVDSTLIGSSSHVEWEDAVHIIPENESDDEEEEEKGPVSPRNLQESEEEEVPQESWDEGYSTLSIPPEMLASYQSYSSTFHSLEEQQVCMAVDIGRHRWDQVKKEDQEATGPRLSRELLDEKGPEVLQDSLDRCYSTPSGCLELTDSCQPYRSAFYVLEQQRVGLAIDMDEIEKYQEVEEDQDPSCPRLSRELLDEKEPEVLQDSLDRCYSIPSGYLELPDLGQPYSSAVYSLEEQYLGLALDVDRIKKDQEEEEDQDPPCPRLSRELVEVVEPEVLQDSLDRCYSTPSSCLEQPDSCQPYGSSFYALEEKHVGFSLDVGEIEKKGKGKKRRGRRSKKERRRGRKEGEEDQNPPCPRLSRELLDEKGPEVLQDSLDRCYSTPSGCLELTDSCQPYRSAFYILEQQCVGLAVDMDEIEKYQEVEEDQDPSCPRLSRELLDEKEPEVLQDSLDRCYSIPSGYLELPDLGQPYSSAVYSLEEQYLGLALDVDRIKKDQEEEEDQDPPCPRLSRELVEVVEPEVLQDSLDRCYSTPSSCLEQPDSCQPYGSSFYALEEKHVGFSLDVGEIEKKGKGKKRRGRRSKKERRRGRKEGEEDQNPPCPRLSRELLDEKGPEVLQDSLDRCYSTPSGCLELTDSCQPYRSAFYILEQQCVGLAIDMDEIEKYQEVEEDQDPSCPRLSRELLDEKEPEVLQDSLDRCYSIPSGYLELPDLGQPYSSAVYSLEEQYLGLALDVDRIKKDQEEEEDQDPPCPRLSRELVEVVEPEVLQDSLDRCYSTPSSCLEQPDSCQPYGSSFYALEEKHVGFSLDVGEIEKKGKGKKRRGRRSKKERRRGRKEGEEDQNPPCPRLSRELLDEKGPEVLQDSLDRCYSTPSGCLELTDSCQPYRSAFYILEQQCVGLAVDMDEIEKYQEVEEDQDPSCPRLSRELLDEKEPEVLQDSLDRCYSIPSGYLELPDLGQPYSSAVYSLEEQYLGLALDVDRIKKDQEEEEDQDPPCPRLSRELVEVVEPEVLQDSLDRCYSTPSSCLEQPDSCQPYGSSFYALEEKHVGFSLDVGEIEKKGKGKKRRGRRSKKERRRGRKEGEEDQNPPCPRLSRELLDEKGPEVLQDSLDRCYSTPSGCLELTDSCQPYRSAFYILEQQCVGLAVDMDEIEKYQEVEEDQDPSCPRLSRELLDEKEPEVLQDSLDRCYSIPSGYLELPDLGQPYSSAVYSLEEQYLGLALDVDRIKKDQEEEEDQDPPCPRLSRELVEVVEPEVLQDSLDRCYSTPSSCLEQPDSCQPYGSSFYALEEKHVGFSLDVGEIEKKGKGKKRRGRRSKKERRRGRKEGEEDQNPPCPRLSRELLDEKGPEVLQDSLDRCYSTPSGCLELTDSCQPYRSAFYILEQQCVGLAVDMDEIEKYQEVEEDQDPSCPRLSRELLDEKEPEVLQDSLDRCYSIPSGYLELPDLGQPYSSAVYSLEEQYLGLALDVDRIKKDQEEEEDQDPPCPRLSRELLEVVEPEVLQDSLDRCYSTPSSCLEQPDSCQPYGSSFYALEEKHVGFSLDVGEIEKKGKGKKRRGRRSKKERRRGRKEGEEDQNPPCPRLSRELLDEKGPEVLQDSLDRCYSTPSGCLELTDSCQPYRSAFYILEQQCVGLAVDMDEIEKYQEVEEDQDPSCPRLSRELLDEKEPEVLQDSLDRCYSTPSGYLELPDLGQPYSSAVYSLEEQYLGLALDVDRIKKDQEEEEDQGPPCPRLSRELLEVVEPEVLQDSLDRCYSTPSSCLEQPDSCQPYGSSFYALEEKHVGFSLDVGEIEKKGKGKKRRGRRSKKERRRGRKEGEEDQNPPCPRLSRELLDEKGPEVLQDSLDRCYSTPSGCLELTDSCQPYRSAFYILEQQCVGLAVDMDEIEKYQEVEEDQDPSCPRLSRELLDEKEPEVLQDSLDRCYSTPSGYLELPDLGQPYSSAVYSLEEQYLGLALDVDRIKKDQEEEEDQGPPCPRLSRELLEVVEPEVLQDSLDRCYSTPSSCLEQPDSCQPYGSSFYALEEKHVGFSLDVGEIEKKGKGKKRRGRRSKKERRRGRKEGEEDQNPPCPRLSRELLDEKGPEVLQDSLDRCYSTPSGCLELTDSCQPYRSAFYILEQQCVGLAVDMDEIEKYQEVEEDQDPSCPRLSRELLDEKEPEVLQDSLDRCYSTPSGYLELPDLGQPYSSAVYSLEEQYLGLALDVDRIKKDQEEEEDQGPPCPRLSRELLEVVEPEVLQDSLDRCYSTPSSCLEQPDSCQPYGSSFYALEEKHVGFSLDVGEIEKKGKGKKRRGRRSKKERRRGRKEGEEDQNPPCPRLSRELLDEKGPEVLQDSLDRCYSTPSGCLELTDSCQPYRSAFYILEQQCVGLAVDMDEIEKYQEVEEDQDPSCPRLSRELLDEKEPEVLQDSLDRCYSTPSGYLELPDLGQPYSSAVYSLEEQYLGLALDVDRIKKDQEEEEDQGPPCPRLSRELLEVVEPEVLQDSLDRCYSTPSSCLEQPDSCQPYGSSFYALEEKHVGFSLDVGEIEKKGKGKKRRGRRSKKERRRGRKEGEEDQNPPCPRLSRELLDEKGPEVLQDSLDRCYSTPSGCLELTDSCQPYRSAFYILEQQCVGLAVDMDEIEKYQEVEEDQDPSCPRLSRELLDEKEPEVLQDSLDRCYSTPSGYLELPDLGQPYSSAVYSLEEQYLGLALDVDRIKKDQEEEEDQGPPCPRLSRELLEVVEPEVLQDSLDRCYSTPSSCLEQPDSCQPYGSSFYALEEKHVGFSLDVGEIEKKGKGKKRRGRRSKKERRRGRKEGEEDQNPPCPRLSRELLDEKGPEVLQDSLDRCYSTPSGCLELTDSCQPYRSAFYILEQQCVGLAVDMDEIEKYQEVEEDQDPSCPRLSRELLAEKEPEVLQDPLDRCYSTPSGYLELPDLGQPYSSAVYSLEEQYLGLALDVDRIKKDQEEEEDQGPPCPRLSRELLEVVEPEVLQDSLDRCYSTPSSCLEQPDSCQPYGSSFYALEEKHVGFSLDVGEIEKKGKGKKRRGRRSKKERRRGRKEGEEDQNPPCPRLSRELLDEKGPEVLQDSLDRCYSTPSGCLELTDSCQPYRSAFYILEQQCVGLAVDMDEIEKYQEVEEDQDPSCPRLSRELLAEKEPEVLQDPLDRCYSTPSGYLELPDLGQPYSSAVYSLEEQYLGLALDVDRIKKDQEEEEDQGPPCPRLSRELLEVVEPEVLQDSLDRCYSTPSSCLEQPDSCQPYGSSFYALEEKHVGFSLDVGEIEKKGKGKKRRGRRSKKERRRGRKEGEEDQNPPCPRLSRELLDEKGPEVLQDSLDRCYSTPSGCLELTDSCQPYRSAFYILEQQCVGLAVDMDEIEKYQEVEEDQDPSCPRLSRELLAEKEPEVLQDPLDRCYSTPSGYLELPDLGQPYSSAVYSLEEQYLGLALDVDRIKKDQEEEEDQGPPCPRLSRELLEVVEPEVLQDSLDRCYSTPSSCLEQPDSCQPYGSSFYALEEKHVGFSLDVGEIEKKGKGKKRRGRRSKKERRRGRKEGEEDQNPPCPRLNSVLMEVEEPEVLQDSLDGCYSTPSMYFELPDSFQHYRSVFYSFEEEHISFALYLDNRFFTLTVTSLHLVFQMLVIFPQ",)
# for (acc, name, species, tax, frags, afdb, seqs) in Query(f"SELECT acc, name, species, tax, MAX(frag), MIN(afdb), GROUP_CONCAT(DISTINCT seq ORDER BY frag SEPARATOR '|') FROM {alphafrag} WHERE species='human' GROUP BY acc LIMIT 100"):
# for (acc, name, species, tax, frags, afdb, seqs) in Query(f"SELECT acc, name, species, tax, MAX(frag), MIN(afdb), GROUP_CONCAT(DISTINCT seq ORDER BY frag SEPARATOR '|') FROM {alphafrag} WHERE acc='A0A087WUL8' GROUP BY acc ORDER BY species='human' DESC, species, tax, acc"):
if not Switch('alphasync'):
    query = Query(f"SELECT acc, name, species, tax, MAX(frag), MIN(afdb), GROUP_CONCAT(DISTINCT seq ORDER BY frag SEPARATOR '|') FROM {alphafrag} GROUP BY acc ORDER BY species='human' DESC, species, tax, acc")
else:
    # Syncing (switch -alphasync active): only get AlphaSync re-predicted proteins
    query = Query(f"SELECT acc, name, species, tax, MAX(frag), MIN(afdb), GROUP_CONCAT(DISTINCT seq ORDER BY frag SEPARATOR '|') FROM {alphafrag} WHERE afdb=0 GROUP BY acc ORDER BY species='human' DESC, species, tax, acc")
for (acc, name, species, tax, frags, afdb, seqs) in tq(query, total=Numrows(query)):
    # if Switch('debug'):
    #     if frags == 1:
    #         continue
    #     print(f" >> {acc} >> {name} >> {species} >> {tax}")

    # Verify sequences
    seqs = seqs.split('|')
    if len(seqs) != frags:
        Die(f"Error: Expected {frags} fragments, but got {len(seqs)} for acc '{acc}'")

    # Combine sequences
    if frags == 1:
        seq = seqs[0]
    else:
        # ...if there are multiple fragments
        seq = ''
        for i in range(frags):
            if Switch('debug'):
                print(f" >> {i + 1}")
                print(f"   ++ {seqs[i]}")
            # start = i * fragstep
            # stop = start + fraglen
            # seq[start:stop] = seqs[i]
            if i + 1 < frags:
                # Initial fragments (1–200)
                seq += seqs[i][0:fragstep]
            else:
                # Final fragment (1–End)
                seq += seqs[i]
            if Switch('debug'):
                print(f"   >> {seq}\n")

    # Check if this sequence is in the wanted_seqs set or human (and skip otherwise)
    if Switch('comparaonly'):
        # if tax not in wanted_taxa or seq not in wanted_seqs:
        # # In addition to wanted_seqs, also retain all human sequences:
        # if tax not in wanted_taxa or (seq not in wanted_seqs and species != 'HUMAN'):
        # In addition to wanted_seqs, also retain all model organism and global health proteome sequences:
        # if (tax not in wanted_complete_taxa and tax not in wanted_taxa) or (tax in wanted_taxa and seq not in wanted_seqs):
        if tax not in wanted_taxa or (tax not in wanted_complete_taxa and seq not in wanted_seqs):
            Log("skipped unwanted sequence for acc", acc)
            Log("skipped unwanted sequence for acc|frags", f"{acc}|{frags}")
            Log("skipped unwanted sequence for name", name)
            Log("skipped unwanted sequence for species", species)
            Log("skipped unwanted sequence for tax", tax)
            Log("skipped unwanted sequence for species|tax", f"{species}|{tax}")
            Log("skipped unwanted sequence for seq", seq)
            continue

    # # Check if this sequence is already in alphaseq (and skip if yes)
    # if seq in seqs_already_in_alphaseq:
    #     Log("skipped sequence already in alphaseq for acc", acc)
    #     Log("skipped sequence already in alphaseq for acc|frags", f"{acc}|{frags}")
    #     Log("skipped sequence already in alphaseq for name", name)
    #     Log("skipped sequence already in alphaseq for species", species)
    #     Log("skipped sequence already in alphaseq for tax", tax)
    #     Log("skipped sequence already in alphaseq for species|tax", f"{species}|{tax}")
    #     Log("skipped sequence already in alphaseq for seq", seq)
    #     continue

    # Insert complete sequences into alphaseq SQL table
    if Switch('debug'):
        print(f"\nINSERT INTO {alphaseq} SET acc='{acc}', name='{name}', species='{species}', tax='{tax}', frags='{frags}', afdb='{afdb}', seq='{seq}'")
    elif not Switch('debug2'):
        Query(f"INSERT INTO {alphaseq} SET acc='{acc}', name='{name}', species='{species}', tax='{tax}', frags='{frags}', afdb='{afdb}', seq='{seq}'")
        # # Verify existing sequence (for debugging)
        # query = Query(f"SELECT seq FROM {alphaseq} WHERE acc='{acc}'")
        # if Numrows(query) > 0:
        #     tmpseq = FetchOne(query)
        #     if tmpseq != seq:
        #         Die(f"Error: Sequence should be\n\n{tmpseq}\n\naccording to the existing alphaseq entry, but it is:\n\n{seq}\n\n")

    Log(f"successfully inserted into table '{alphaseq}' for acc", acc)
    Log(f"successfully inserted into table '{alphaseq}' for acc|frags", f"{acc}|{frags}")
    Log(f"successfully inserted into table '{alphaseq}' for name", name)
    Log(f"successfully inserted into table '{alphaseq}' for species", species)
    Log(f"successfully inserted into table '{alphaseq}' for tax", tax)
    Log(f"successfully inserted into table '{alphaseq}' for species|tax", f"{species}|{tax}")
    Log(f"successfully inserted into table '{alphaseq}' for seq", seq)



Show(lim=20)

if not Switch('debug'):
    Optimize(alphaseq)

Stoptime()
print("\nDone!")
