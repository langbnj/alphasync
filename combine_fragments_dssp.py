#!/usr/bin/env python3
"""
Combine DeepMind AlphaFold fragment runs into a single output file (pLDDT and DSSP relASA values, and predicted secondary structure)

Ignores scores within ±200 aa of "artificial termini" introduced by the fragmentation, and averages values across fragments (see Methods).

"""

# Initialize
import pandas as pd
# import numpy as np
from numpy import mean
# import networkx as nx
# import seaborn as sns
# import matplotlib.pyplot as mp
from collections import Counter
from blang_mysql import *
from blang import *

# Fragment length and step size used by DeepMind.
# Proteins longer than 2700 residues are split into windows of width 1400 with a step size of 200.
fraglen = 1400
fragstep = 200

# Solvent-accessible surface area maxima table for residues
asafile = "../../input/asa/max_asa_for_residues.tsv"

# . {windowsize} {disthresh} {acc} {maxfrag} {af2file}
(inpath, windowsize, disthresh, acc, maxfrag, outfile) = Args(6, "[directory] [window size] [disorder threshold] [UniProt accession] [maximum fragment number] [output filename]\n\n -alphasync: Updating AlphaSync proteins (non-AFDB, i.e. afdb=0)", ". 10 0.55 A0A087WUL8 14 AF-A0A087WUL8-F14-model_v4.af2")

# Set afdb to 0 for AlphaSync proteins (not in AFDB)
afdb = 1
if Switch('alphasync'):
    afdb = 0



# Start

# Read reference maximum ASA values (skipping comment lines)
refasas = read_tsv(asafile, header=5)
# Reduce to columns that are actually used
refasas = refasas[["1-letter", "Empirical"]]
# Rename columns to something simpler (aa and asa)
refasas = refasas.rename(columns = {"1-letter": "aa", "Empirical": "asa"})

# Get expected fragment sequences from table 'alphafrag'
seqs = FetchMap(Query(f"SELECT frag, seq FROM alphafrag WHERE acc='{acc}' AND afdb='{afdb}'"))



# Start parsing pLDDT values from mmCIF files
print(f"\nParsing pLDDT values from {maxfrag} .cif fragment files in directory '{inpath}':")
# Loop through fragments (1..maxfrag)
aas = {}
plddts = {}
for frag in tq(range(maxfrag)):
    frag += 1

    # Check if file exists (with tolerance for version number updates)
    infiles = nsort(Return(f"ls -1U {inpath}/AF-{acc}-F{frag}-model_v*.cif").split("\n"))
    if len(infiles) != 1:
        Die(f"Error: Found {len(infiles)} mmCIF fragment files with fragment number '{frag}' (expected 1)")

    for ciffile in infiles:

        # Initialize
        seq = ''

        if Switch('debug'):
            print(f" >> {ciffile}")

        # Parse pLDDT scores
        with open(ciffile) as f:
            for line in f:
                
                # "_ma_qa_metric_local" is only present in AFDB .cif files:
                # loop_
                # _ma_qa_metric_local.label_asym_id
                # _ma_qa_metric_local.label_comp_id
                # _ma_qa_metric_local.label_seq_id
                # _ma_qa_metric_local.metric_id
                # _ma_qa_metric_local.metric_value
                # _ma_qa_metric_local.model_id
                # _ma_qa_metric_local.ordinal_id
                # A MET 1   2 53.51 1 1   
                # A MET 2   2 62.17 1 2   
                # ...
                # A VAL 180 2 52.79 1 180 
                # A ARG 181 2 54.09 1 181 
                # #

                # "_ma_qa_metric_local" is not present in AlphaFold 2.3.2 (AlphaSync) .cif files, hence using "_atom_site":
                # #
                # loop_
                # _atom_site.group_PDB
                # _atom_site.id
                # _atom_site.type_symbol
                # _atom_site.label_atom_id
                # _atom_site.label_alt_id
                # _atom_site.label_comp_id
                # _atom_site.label_asym_id
                # _atom_site.label_entity_id
                # _atom_site.label_seq_id
                # _atom_site.pdbx_PDB_ins_code
                # _atom_site.Cartn_x
                # _atom_site.Cartn_y
                # _atom_site.Cartn_z
                # _atom_site.occupancy
                # _atom_site.B_iso_or_equiv
                # _atom_site.auth_seq_id
                # _atom_site.auth_asym_id
                # _atom_site.pdbx_PDB_model_num
                # ATOM 1     N N   . MET A 0 1    . 33.411  28.750   -62.996 1.00 50.54 1    A 1 
                # ATOM 2     C CA  . MET A 0 1    . 32.096  28.066   -62.959 1.00 50.54 1    A 1 
                # ATOM 3     C C   . MET A 0 1    . 32.002  27.022   -61.835 1.00 50.54 1    A 1 
                # ATOM 4     C CB  . MET A 0 1    . 31.781  27.471   -64.342 1.00 50.54 1    A 1 
                # ATOM 5     O O   . MET A 0 1    . 30.903  26.653   -61.459 1.00 50.54 1    A 1 
                # ATOM 6     C CG  . MET A 0 1    . 30.281  27.415   -64.647 1.00 50.54 1    A 1 
                # ATOM 7     S SD  . MET A 0 1    . 29.962  27.178   -66.415 1.00 50.54 1    A 1 
                # ATOM 8     C CE  . MET A 0 1    . 28.149  27.132   -66.430 1.00 50.54 1    A 1 
                # ATOM 9     N N   . GLU A 0 2    . 33.123  26.579   -61.256 1.00 53.84 2    A 1 
                # ATOM 10    C CA  . GLU A 0 2    . 33.210  25.372   -60.409 1.00 53.84 2    A 1 
                # ATOM 11    C C   . GLU A 0 2    . 32.645  25.499   -58.981 1.00 53.84 2    A 1 
                # ATOM 12    C CB  . GLU A 0 2    . 34.694  24.985   -60.337 1.00 53.84 2    A 1 
                # ATOM 13    O O   . GLU A 0 2    . 32.414  24.491   -58.321 1.00 53.84 2    A 1 
                # ATOM 14    C CG  . GLU A 0 2    . 35.281  24.676   -61.728 1.00 53.84 2    A 1 
                # ATOM 15    C CD  . GLU A 0 2    . 36.807  24.699   -61.713 1.00 53.84 2    A 1 
                # ATOM 16    O OE1 . GLU A 0 2    . 37.403  23.725   -62.211 1.00 53.84 2    A 1 
                # ATOM 17    O OE2 . GLU A 0 2    . 37.330  25.746   -61.266 1.00 53.84 2    A 1 
                # ATOM 18    N N   . ALA A 0 3    . 32.370  26.714   -58.494 1.00 56.43 3    A 1 
                # ATOM 19    C CA  . ALA A 0 3    . 31.952  26.964   -57.107 1.00 56.43 3    A 1 
                # ATOM 20    C C   . ALA A 0 3    . 30.546  26.434   -56.725 1.00 56.43 3    A 1 
                # ATOM 21    C CB  . ALA A 0 3    . 32.095  28.469   -56.835 1.00 56.43 3    A 1 
                # ATOM 22    O O   . ALA A 0 3    . 30.090  26.683   -55.612 1.00 56.43 3    A 1 
                # [...]
                # ATOM 10880 N N   . SER A 0 1399 . -28.410 22.537   18.712  1.00 40.64 1399 A 1 
                # ATOM 10881 C CA  . SER A 0 1399 . -27.615 23.766   18.961  1.00 40.64 1399 A 1 
                # ATOM 10882 C C   . SER A 0 1399 . -26.092 23.522   19.069  1.00 40.64 1399 A 1 
                # ATOM 10883 C CB  . SER A 0 1399 . -28.060 24.546   20.220  1.00 40.64 1399 A 1 
                # ATOM 10884 O O   . SER A 0 1399 . -25.579 23.303   20.169  1.00 40.64 1399 A 1 
                # ATOM 10885 O OG  . SER A 0 1399 . -27.926 23.790   21.406  1.00 40.64 1399 A 1 
                # ATOM 10886 N N   . VAL A 0 1400 . -25.343 23.690   17.971  1.00 36.64 1400 A 1 
                # ATOM 10887 C CA  . VAL A 0 1400 . -23.936 24.156   17.986  1.00 36.64 1400 A 1 
                # ATOM 10888 C C   . VAL A 0 1400 . -23.731 25.158   16.866  1.00 36.64 1400 A 1 
                # ATOM 10889 C CB  . VAL A 0 1400 . -22.903 23.011   17.927  1.00 36.64 1400 A 1 
                # ATOM 10890 O O   . VAL A 0 1400 . -24.286 24.919   15.775  1.00 36.64 1400 A 1 
                # ATOM 10891 C CG1 . VAL A 0 1400 . -21.518 23.410   17.384  1.00 36.64 1400 A 1 
                # ATOM 10892 C CG2 . VAL A 0 1400 . -22.681 22.436   19.333  1.00 36.64 1400 A 1 
                # ATOM 10893 O OXT . VAL A 0 1400 . -23.048 26.157   17.199  1.00 36.64 1400 A 1 
                # #

                # From AFDB:
                # #
                # loop_
                # _atom_site.group_PDB
                # _atom_site.id
                # _atom_site.type_symbol
                # _atom_site.label_atom_id
                # _atom_site.label_alt_id
                # _atom_site.label_comp_id
                # _atom_site.label_asym_id
                # _atom_site.label_entity_id
                # _atom_site.label_seq_id
                # _atom_site.pdbx_PDB_ins_code
                # _atom_site.Cartn_x
                # _atom_site.Cartn_y
                # _atom_site.Cartn_z
                # _atom_site.occupancy
                # _atom_site.B_iso_or_equiv
                # _atom_site.pdbx_formal_charge
                # _atom_site.auth_seq_id
                # _atom_site.auth_comp_id
                # _atom_site.auth_asym_id
                # _atom_site.auth_atom_id
                # _atom_site.pdbx_PDB_model_num
                # _atom_site.pdbx_sifts_xref_db_acc
                # _atom_site.pdbx_sifts_xref_db_name
                # _atom_site.pdbx_sifts_xref_db_num
                # _atom_site.pdbx_sifts_xref_db_res
                # ATOM 1    N N   . MET A 1 1   ? -8.753  -5.120  -6.078  1.0 30.38 ? 1   MET A N   1 A0A061ACK4 UNP 1   M 
                # ATOM 2    C CA  . MET A 1 1   ? -8.959  -4.233  -7.251  1.0 30.38 ? 1   MET A CA  1 A0A061ACK4 UNP 1   M 
                # ATOM 3    C C   . MET A 1 1   ? -9.229  -2.830  -6.739  1.0 30.38 ? 1   MET A C   1 A0A061ACK4 UNP 1   M 

                # # Skip ahead to pLDDT section
                # if line != "_ma_qa_metric_local.ordinal_id\n":
                # Skip ahead to PDB section and read pLDDT from B-factor columns
                if line != "_atom_site.pdbx_PDB_model_num\n":
                    continue
                # Read pLDDT from PDB section
                prevsite = 0
                for line in f:
                    # Skip extra _atom_site.… lines
                    if line.startswith("_atom_site."):
                        continue
                    # Stop parsing on '#'
                    if line == "#\n":
                        break

                    # Split line
                    a = re.split(r" +", line.rstrip())
                    
                    # # Parse pLDDT values
                    # aa3 = a[1]
                    # site = int(a[2])
                    # plddt = float(a[4])
                    # tmpsite = int(a[6])


                    # Parse pLDDT values
                    aa3 = a[5]
                    site = int(a[8])
                    plddt = float(a[14])

                    # Build sequence (for verification)
                    aa = ThreeToOne(aa3)

                    # if frag == 2:
                    #     d()

                    # Only process whenever site is increased compared to previous site (since each atom has its own row, leading to repeats)
                    if site > prevsite:
                        prevsite = site

                        # Check if seq[site - 1] is already defined and fail if aa is not as expected
                        if len(seq) >= site:
                            if seq[site - 1] != aa:
                                Die(f"Error: Amino acid mismatch at site '{site}': Expected '{seq[site - 1]}', but got '{aa}'")
                        else:
                            seq += aa

                        # Ignore values from dubious regions (artificial termini)
                        # Remove N-terminal 200 aa at an artificial N-terminus
                        if frag > 1:
                            if site <= fragstep:
                                if Switch('debug'):
                                    print(f" >> mmCIF >> Skipping residue {frag} {site} ({site + fragstep * (frag - 1)}) because it's too close to an artificial N-terminus")
                                continue
                        # Remove C-terminal 200 aa at an artificial C-terminus
                        if frag < maxfrag:
                            if site > len(seqs[frag]) - fragstep:
                                if Switch('debug'):
                                    print(f" >> mmCIF >> Skipping residue {frag} {site} ({site + fragstep * (frag - 1)}) because it's too close to an artificial C-terminus")
                                continue
                        
                        # Shift site according to fragment number
                        site += fragstep * (frag - 1)

                        # Store amino acid for this residue
                        if site in aas:
                            if aas[site] != aa:
                                Die(f"Error: Amino acid mismatch at site '{site}': Expected '{aas[site]}', but got '{aa}'")
                        else:
                            aas[site] = aa

                        # Store pLDDT value for this residue (in a list)
                        if site in plddts:
                            plddts[site].append(plddt)
                        else:
                            plddts[site] = [plddt]

                        if Switch('debug'):
                            print(f" >> mmCIF >> {acc} >> {frag} >> {site} >> {aa} >> {plddt}")

                # Finished parsing this fragment's mmCIF file:

                # Verify complete fragment sequence
                if seq != seqs[frag]:
                    Die(f"Error: Expected '{acc}' fragment '{frag}' sequence to be:\n\n{seqs[frag]}\n\n...but found this sequence in '{ciffile}':\n\n{seq}\n\n")
                else:
                    Log("sequence matches alphafrag sequence for acc|frag", f"{acc}|{frag}")



# Start parsing accessible surface areas and secondary structure calls from DSSP output files
print(f"\nParsing relASA values and secondary structure calls from {maxfrag} .dssp fragment files in directory '{inpath}':")
# Loop through fragments (1..maxfrag)
asas = {}
relasas = {}
secs = {}
for frag in tq(range(maxfrag)):
    frag += 1

    # Check if file exists (with tolerance for version number updates)
    infiles = nsort(Return(f"ls -1U {inpath}/AF-{acc}-F{frag}-model_v*.dssp").split("\n"))
    if len(infiles) != 1:
        Die(f"Error: Found {len(infiles)} DSSP fragment files with fragment number '{frag}' (expected 1)")

    for dsspfile in infiles:

        # Initialize
        seq = ''

        if Switch('debug'):
            print(f" >> {dsspfile}")

        # Parse DSSP file
        with open(dsspfile) as f:

            # Skip header
            for line in f:
                if line == "  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA\n":
                    break
            
            # Read content
            for line in f:
                #     1    1 A G              0   0   62      0, 0.0     2,-0.1     0, 0.0   994,-0.0   0.000 360.0 360.0 360.0 141.5  -93.6  -22.6   59.6
                #     2    2 A E        +     0   0  143    993,-0.1   991,-0.0   992,-0.1     0, 0.0   0.417 360.0 113.6  57.1 133.8  -91.8  -20.5   58.1
                # ...
                #  1399 1399 A L              0   0  161     -2,-0.5  -451,-0.1  -451,-0.0  -452,-0.0  -0.951 360.0 360.0-123.7 113.9  -26.8  -15.2   79.2
                #  1400 1400 A D              0   0  244     -2,-0.4    -2,-0.0     0, 0.0     0, 0.0  -0.948 360.0 360.0-151.9 360.0  -23.0  -14.7   79.8

                if line[13:14] == "!":
                    Log("warning ('!' in aa column) encountered in DSSP file for acc|frag (ignored)", f"{acc}|{frag}")
                    continue

                # Parse values (fixed positions)
                site = int(line[5:10])
                aa = line[13:14]
                sec = line[16:17]
                asa = float(line[35:38])

                # Calculate relative ASA
                relasa = asa / float(refasas[refasas["aa"] == aa]["asa"].iloc[0])

                # Build sequence (for verification)
                seq += aa

                # Ignore values from dubious regions (artificial termini)
                # Remove N-terminal 200 aa at an artificial N-terminus
                if frag > 1:
                    if site <= fragstep:
                        if Switch('debug'):
                            print(f" >> DSSP >> Skipping residue {frag} {site} ({site + fragstep * (frag - 1)}) because it's too close to an artificial N-terminus")
                        continue
                # Remove C-terminal 200 aa at an artificial C-terminus
                if frag < maxfrag:
                    if site > len(seqs[frag]) - fragstep:
                        if Switch('debug'):
                            print(f" >> DSSP >> Skipping residue {frag} {site} ({site + fragstep * (frag - 1)}) because it's too close to an artificial C-terminus")
                        continue
                
                # Shift site according to fragment number
                site += fragstep * (frag - 1)

                # Store ASA for this residue (in a list)
                if site in asas:
                    asas[site].append(asa)
                else:
                    asas[site] = [asa]

                # Store relative ASA for this residue (in a list)
                if site in relasas:
                    relasas[site].append(relasa)
                else:
                    relasas[site] = [relasa]

                # Store secondary structure predictions for this residue (in a list)
                if site in secs:
                    secs[site].append(sec)
                else:
                    secs[site] = [sec]

                if Switch('debug'):
                    print(f" >> DSSP >> {acc} >> {frag} >> {site} >> {aa} >> {asa} >> '{sec}'")
            
            # Finished parsing this fragment's DSSP file:

            # Verify complete fragment sequence
            if seq != seqs[frag]:
                Die(f"Error: Expected '{acc}' fragment '{frag}' sequence to be:\n\n{seqs[frag]}\n\n...but found this sequence in '{dsspfile}':\n\n{seq}\n\n")
            else:
                Log("sequence matches alphafrag sequence for acc|frag", f"{acc}|{frag}")



# Build complete protein sequence across fragments
seq = ''.join(aas.values())
# Verify complete protein sequence
if len(seq) != max(aas):
    Die(f"Error: Expected complete protein sequence of length {max(aas)}, but got {len(seq)}")
if len(seq) != max(plddts):
    Die(f"Error: Expected complete protein sequence of length {max(plddts)}, but got {len(seq)}")
if len(seq) != max(asas):
    Die(f"Error: Expected complete protein sequence of length {max(asas)}, but got {len(seq)}")
if len(seq) != max(secs):
    Die(f"Error: Expected complete protein sequence of length {max(secs)}, but got {len(seq)}")



# Average pLDDT scores and ASA values, and choose most common secondary structure, across fragments (ignoring 200 aa near artificial N- and C-termini introduced in fragments)
print(f"\nAveraging pLDDT scores and ASA values and choosing most common secondary structure across {maxfrag} fragments:")
with open(outfile, 'w') as out:

    # Print header
    print("#position	residue1	residue3	pLDDT	pLDDT_smoothed_10	DSSP_sec_struct	SASA_absolute	SASA_absolute_smoothed_10	SASA_relative	SASA_relative_smoothed_10	disordered", file=out)

    # Loop through residues
    for site in tq(range(len(seq))):
        site += 1

        if Switch('debug'):
            print(f" >> {site}")

        # Residue
        aa = aas[site]

        # Average pLDDT at this position
        plddt = mean(plddts[site])

        # Average ASA at this position
        asa = mean(asas[site])

        # Average relative ASA at this position
        # relasa = asa / float(refasas[refasas["aa"] == aa]["asa"])
        relasa = mean(relasas[site])



        # Smoothed values

        # Smoothed pLDDT (± windowsize aa)
        window_plddts = []
        # window_plddts_individual = []
        for i in range(site - windowsize, site + windowsize + 1):
            if i in plddts:
                # # Create a single list of all pLDDT values within the window
                # for p in plddts[i]:
                #     window_plddts_individual.append(p)
                # Average per residue, then average across residues (ensures that there aren't any artifacts across boundaries where one residue has fewer fragments covering it)
                window_plddts.append(mean(plddts[i]))
        plddt_smoothed = mean(window_plddts)
        # plddt_smoothed_individual = mean(window_plddts_individual)
        
        # Smoothed ASA (± windowsize aa)
        window_asas = []
        for i in range(site - windowsize, site + windowsize + 1):
            if i in asas:
                # Average per residue, then average across residues (ensures that there aren't any artifacts across boundaries where one residue has fewer fragments covering it)
                window_asas.append(mean(asas[i]))
        asa_smoothed = mean(window_asas)

        # Smoothed relASA (± windowsize aa)
        window_relasas = []
        for i in range(site - windowsize, site + windowsize + 1):
            if i in relasas:
                # Average per residue, then average across residues (ensures that there aren't any artifacts across boundaries where one residue has fewer fragments covering it)
                window_relasas.append(mean(relasas[i]))
        relasa_smoothed = mean(window_relasas)



        # Most common secondary structure type at this position (with pLDDT value as tiebreaker)
        
        # Get best secondary structures by count (most common)
        bestsecs_by_count = set()
        seccounts = Counter(secs[site])
        for sec in set(secs[site]):
            if seccounts[sec] == max(seccounts.values()):
                bestsecs_by_count.add(sec)
        
        # Get best secondary structures by pLDDT value
        # bestfrag = plddts[site].index(max(plddts[site]))
        # bestsecs_by_plddt = secs[site][bestfrag]
        bestsecs_by_plddt = set()
        for i in range(len(plddts[site])):
            if plddts[site][i] == max(plddts[site]):
                bestsecs_by_plddt.add(secs[site][i])

        # See if the two metrics agree
        if len(bestsecs_by_count) == 1:
            # If there is only one "best" secondary structure by count:
            # Convert set to iterator and get its only element
            bestsec = next(iter(bestsecs_by_count))
        else:
            # If there are multiple "best" secondary structures by count:
            # See if one of them is found at the residue with the best pLDDT value
            if len(bestsecs_by_count.intersection(bestsecs_by_plddt)) >= 1:
                # If there are secondary structure types that are best according to count as well as the maximum-pLDDT residues:
                # Convert set to iterator and get the first element (choosing a type at random if there are multiple, which should be extremely rare)
                bestsec = next(iter(bestsecs_by_count.intersection(bestsecs_by_plddt)))
            elif len(bestsecs_by_count.intersection(bestsecs_by_plddt)) == 0:
                # If there is no agreement:
                # Convert set to iterator and get the first element (choosing a type at random if there are multiple, which should be extremely rare)
                bestsec = next(iter(bestsecs_by_count))
        
        # Most common secondary structure type across fragments, with pLDDT as tiebreaker
        sec = bestsec
        
        # Disordered yes/no
        if relasa_smoothed >= disthresh:
            dis = 1
        else:
            dis = 0

        # Print to output file
        # #position	residue1	residue3	pLDDT	pLDDT_smoothed_10	DSSP_sec_struct	SASA_absolute	SASA_absolute_smoothed_10	SASA_relative	SASA_relative_smoothed_10	disordered
        # 1	M	MET	36.39	44.29	C	236	129.909	1	0.739	1
        # 2	V	VAL	39.44	44.884	C	129	132.333	0.782	0.735	1
        # 3	V	VAL	42.10	45.579	C	137	129.077	0.83	0.735	1
        s = f"{site}\t{aa}\t{OneToThree(aa)}\t{plddt}\t{plddt_smoothed}\t{sec}\t{asa}\t{asa_smoothed}\t{relasa}\t{relasa_smoothed}\t{dis}"
        if Switch('debug'):
            print(s)
        print(s, file=out)

        # if site == 2632:
        #     d()




Show(lim=0)

print(f"\nWrote to '{outfile}'")

print("\nDone!")
