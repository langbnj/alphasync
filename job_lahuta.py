#!/usr/bin/env python3
"""
Job script (runs on a given protein accession): Get contacts from mmCIF file using the lahuta module (https://github.com/bisejdiu/lahuta) and insert them into the 'alphacon' MySQL table
"""

# Initialize
import pandas as pd
import numpy as np
from blang_mysql import *
from blang import *
np.set_printoptions(suppress=True)

# Currently using an old Lahuta (a pre-release v0.6 version). Keeping this version for reproducibility.
# 
# To install Lahuta (current, old version):
# conda install -c conda-forge openbabel mdanalysis
# cd ./bin/lahuta
# pip install -e . 
# 
# Here are the original versions of these packages that I used in 2023-04-12:
# ~/miniconda3/lib/python3.10/site-packages/
# MDAnalysis-2.4.3
# openbabel 3.1.0
# Using the latest now (2024-03-18).
# 
# Installing these using pip (since conda proposes a crazy number of updates):
# conda install mdanalysis=2.4.3
# conda install openbabel=3.1.0
# pip install mdanalysis==2.4.3
# pip install openbabel
# 
# For updating Lahuta to the latest, current version (not done yet):
# conda install -c bisejdiu lahuta

# Import Lahuta
from lahuta import AtomGroup
from lahuta.config.defaults import CONTACTS
# from lahuta.contacts import (
#     AromaticContacts,
#     CarbonylContacts,
#     CovalentContacts,
#     HBondContacts,
#     HydrophobicContacts,
#     IonicContacts,
#     MetalContacts,
#     PolarHBondContacts,
#     WeakHBondContacts,
#     WeakPolarHBondContacts,
# )
# from lahuta.contacts.vdw import VanDerWaalsContacts
# from lahuta.contacts.atomplane import AtomPlaneContacts
# from lahuta.contacts.planeplane import PlanePlaneContacts
from lahuta.core.universe import Universe
import lahuta.contacts



# Define contact types to get
# Note that these are "hardcoded" as options for an ENUM field in the MySQL table to save space
types = (
    "AromaticContacts",
    "CarbonylContacts",
    "CovalentContacts",
    "HBondContacts",
    "HydrophobicContacts",
    "IonicContacts",
    "MetalContacts",
    "PolarHBondContacts",
    "VanDerWaalsContacts",
    "WeakHBondContacts",
    "WeakPolarHBondContacts",
)

alphafrag = "alphafrag"     # SQL table with fragment protein sequences (>2700 aa proteins get split into 1400 aa fragments with a step size of 200 in AlphaFold DB, for human only - other species don't have results for >2700 aa proteins)
alphaseq = "alphaseq"       # SQL table with complete protein sequences (used for retrieving additional protein information here)
alphacon = "alphacon"       # SQL table with residue-level contacts from Lahuta

# Fragment length and step size used by DeepMind.
# Proteins longer than 2700 residues are split into windows of width 1400 with a step size of 200.
fraglen = 1400
fragstep = 200

# Get arguments
(acc, maxfrag) = Args(2, "[UniProt accession] [Number of fragments]\n\n -alphasync: Updating AlphaSync proteins (non-AFDB, i.e. afdb=0)", "A0A087WUL8 14")

# inpath = f"../{source}/{acc}"
inpath = f"../{acc}"
# inpath = "."
logpath = "tmp/_logs"

# Set afdb to 0 for AlphaSync proteins (not in AFDB)
afdb = 1
if Switch('alphasync'):
    afdb = 0



# Functions

# Utility function that returns a data frame of contacts (with a column specifying the contact type)
def GetContacts(type, frag, universe, neighbors):

    # Get data frame of contacts
    df = getattr(lahuta.contacts, type)(universe, neighbors).contacts("dataframe", "expanded")

    # Add contact type as first column
    df.insert(0, "type", type)

    # Add current fragment number as second column
    df.insert(1, "frag", frag)

    # Return data frame
    return df

# Format Lahuta data frame for insertion into SQL table
def FormatContacts(df):
    
    # Target structure (SQL table 'alphacon'):
    # acc	species	tax	frags	afdb	site1	site2	aa1	aa2	type	atom1	atom2  dist

    # Add extra columns on left
    df.insert(0, "acc", acc)
    df.insert(1, "species", species)
    df.insert(2, "tax", tax)
    df.insert(3, "frags", maxfrag)
    df.insert(4, "afdb", afdb)

    # Replace AA3 (e.g. PHE) with AA (e.g. F)
    # df = df.assign(aa1 = ThreeToOne(df["residue1_resnames"]))
    # df = df.assign(aa2 = ThreeToOne(df["residue2_resnames"]))
    df["aa1"] = df["residue1_resnames"].apply(ThreeToOne)
    df["aa2"] = df["residue2_resnames"].apply(ThreeToOne)
    # df = df.drop("residue1_resnames", axis=1)
    # df = df.drop("residue2_resnames", axis=1)

    # Rename existing columns
    # residue1_resids  residue2_resids residue1_resnames residue2_resnames residue1_names residue2_names  residue1_indices  residue2_indices  distances
    #          26               19               PHE               PHE            CD2            CD1               204               149   3.650668
    #          19               26               PHE               PHE            CD1            CE2               149               206   3.996618
    #          19               26               PHE               PHE            CD1             CG               149               202   3.698249
    df = df.rename(columns={
        "residue1_resids": "site1",
        "residue2_resids": "site2",
        "residue1_names": "atom1",
        "residue2_names": "atom2",
        # "residue1_indices": "groupid1",
        # "residue2_indices": "groupid2",
        "distances": "dist",
    })

    # Select final columns in final order
    df = df[["acc", "species", "tax", "frag", "frags", "afdb", "site1", "site2", "aa1", "aa2", "type", "atom1", "atom2", "dist"]]

    # Flip cases where site1 > site2 (so that site1 will always be < site2)
    df.loc[df["site1"] > df["site2"], ["site1", "site2", "aa1", "aa2", "atom1", "atom2"]] = df.loc[df["site1"] > df["site2"]].rename(columns={"site1": "site2", "site2": "site1", "aa1": "aa2", "aa2": "aa1", "atom1": "atom2", "atom2": "atom1"})[["site2", "site1", "aa2", "aa1", "atom2", "atom1"]]

    # Return
    return df

# Combine fragments
# Ignore values from dubious regions (artificial termini)
def CombineFragmentContacts(df):
    
    # Ignore values from dubious regions (artificial termini)
    # if frag > 1:
    #     if site <= fragstep:
    #         if Switch('debug'):
    #             print(f" >> Skipping residue {frag} {site} ({site + fragstep * (frag - 1)}) because it's too close to an artificial N-terminus")
    #         continue
    # if frag < maxfrag:
    #     if site > len(seqs[frag]) - fragstep:
    #         if Switch('debug'):
    #             print(f" >> Skipping residue {frag} {site} ({site + fragstep * (frag - 1)}) because it's too close to an artificial C-terminus")
    #         continue
    # Remove N-terminal 200 aa at an artificial N-terminus
    # Remove C-terminal 200 aa at an artificial C-terminus
    # df.loc[(df["frag"] > 1 & df["site"] > fragstep) | (df["frag"] < maxfrag & df["site"] <= len(seqs[frag]) - fragstep)]
    # df = df.loc[(
    #     (((df["frag"] > 1) & (df["site1"] > fragstep)) & 
    #      ((df["frag"] < maxfrag) & (df["site1"] <= fraglen - fragstep))) & 
    #     (((df["frag"] > 1) & (df["site2"] > fragstep))) & 
    #      ((df["frag"] < maxfrag) & (df["site2"] <= fraglen - fragstep)))
    #     )]

    # Ignore values from dubious regions (artificial termini)

    # Remove N-terminal 200 aa at artificial N-termini (frag > 1 will have an artificial N-terminus)
    # df.loc[(df["frag"] > 1) & (df["site1"] <= fragstep)]
    # df.loc[(df["frag"] > 1) & (df["site2"] <= fragstep)]
    # Remove any residue in fragment 2 or above that is between 1-200
    df = df.loc[np.invert((df["frag"] > 1) & (df["site1"] <= fragstep))] # site1
    df = df.loc[np.invert((df["frag"] > 1) & (df["site2"] <= fragstep))] # site2

    # Remove C-terminal 200 aa at artificial C-termini (frag < maxfrag will have an artificial C-terminus)
    # df.loc[(df["frag"] < maxfrag) & (df["site1"] > fraglen - fragstep)]
    # df.loc[(df["frag"] < maxfrag) & (df["site2"] > fraglen - fragstep)]
    # Remove any residue except in the last fragment that is between 1201-1400
    df = df.loc[np.invert((df["frag"] < maxfrag) & (df["site1"] > fraglen - fragstep))] # site1
    df = df.loc[np.invert((df["frag"] < maxfrag) & (df["site2"] > fraglen - fragstep))] # site2

    # Shift site according to fragment number (i.e. residue 1 in fragment 2 will become 1401)
    df["site1"] += fragstep * (df["frag"] - 1)
    df["site2"] += fragstep * (df["frag"] - 1)

    # Verify that all amino acid positions listed are the correct residue
    for i, row in pd.concat([df[["site1", "aa1"]].rename(columns={"site1": "site", "aa1": "aa"}), df[["site2", "aa2"]].rename(columns={"site2": "site", "aa2": "aa"})]).iterrows():
        tmpsite = row["site"]
        tmpaa = row["aa"]
        # Expected amino acid
        expaa = seq[tmpsite-1:tmpsite]
        # Verify that the amino acid is correct (using the alphaseq sequence retrieved earlier)
        if tmpaa != expaa:
            Die(f"Error: Expected residue '{expaa}' at position '{tmpsite}' in acc '{acc}', but got '{tmpaa}'")

    # # Select final columns in final order (leaving out "frag" since fragments have been combined already here)
    # df = df[["acc", "species", "tax", "frags", "afdb", "site1", "site2", "aa1", "aa2", "type", "atom1", "atom2", "dist"]]

    # Remove "frag" column
    df = df.drop("frag", axis=1)

    # Get unique contacts (union across fragments, averaging the distance across fragments)
    # TODO: Could implement "majority vote" here, as for secondary structure type, but for presence/absence of contact? We do trust these parts of the fragments, though (>200 aa away from artificial termini). Currently using "union" of contacts across fragments.
    # df = df.drop_duplicates()
    # df = df.drop_duplicates(subset=["acc", "species", "tax", "frags", "afdb", "site1", "site2", "aa1", "aa2", "type", "atom1", "atom2"])
    # df = df.pivot_table(index='org', columns='cluster', values='time', aggfunc='mean').mean()
    # df = df.groupby(["acc", "species", "tax", "frags", "afdb", "site1", "site2", "aa1", "aa2", "type", "atom1", "atom2"], as_index=False).count()
    df = df.groupby(["acc", "species", "tax", "frags", "afdb", "site1", "site2", "aa1", "aa2", "type", "atom1", "atom2"], as_index=False).mean()

    # Sort rows by values (site1, site2, etc.)
    df = df.sort_values(["acc", "site1", "site2", "type", "dist"], ignore_index=True)
    # Return
    return df

# Insert formatted Lahuta data frame into SQL table
def InsertContacts(df):
    # Select final columns in final order
    # df = df[["acc", "species", "tax", "frag", "frags", "afdb", "site1", "site2", "aa1", "aa2", "type", "atom1", "atom2", "dist"]]
    # acc	site1	site2	aa1	aa2	atom1	atom2	type	afdb	dist	pae
    # A0A021WW64	3	5	V	V	CG1	CG2	HydrophobicContacts	1	4.40445	
    df = df[["acc", "site1", "site2", "aa1", "aa2", "atom1", "atom2", "type", "afdb", "dist"]]
    InsertPanda(df, alphacon)



# Start

# Check if data already exists for this acc in table 'alphacon'
query = Query(f"SELECT * FROM {alphacon} WHERE acc='{acc}' AND afdb='{afdb}' LIMIT 1")
# ...and exit if yes
if Numrows(query) == 1:
    print(f"Data already existed for acc '{acc}' (afdb={afdb}) in table '{alphacon}', exiting (skip)!")
    sys.exit()
else:
    print(f"No data yet for acc '{acc}' (afdb={afdb}) in table '{alphacon}', starting!")



# # Get additional information on this UniProt accession (species etc.) from table 'alphafrag'
# query = Query(f"SELECT DISTINCT species, tax, MAX(frag), afdb FROM {alphafrag} WHERE acc='{acc}'")
# (species, tax, tmpmaxfrag, afdb) = FetchOne(query)
# Get additional information on this acc from table 'alphaseq'
query = Query(f"SELECT DISTINCT name, species, tax, frags, seq FROM {alphaseq} WHERE acc='{acc}' AND afdb='{afdb}'")
(name, species, tax, tmpmaxfrag, seq) = FetchOne(query)

if maxfrag != tmpmaxfrag:
    # Die(f"Error: Expected {maxfrag} fragments for source '{source}' acc '{acc}', but got {tmpmaxfrag}")
    Die(f"Error: Expected {maxfrag} fragments for acc '{acc}', but got {tmpmaxfrag}")

# # Get additional information on this UniProt accession (sequence)
# query = Query(f"SELECT DISTINCT seq FROM {alphaseq} WHERE acc='{acc}'")
# (seq) = FetchOne(query)

# # Get expected fragment sequences from table 'alphafrag'
# seqs = FetchMap(Query(f"SELECT frag, seq FROM {alphafrag} WHERE acc='{acc}'"))



# Start

print(f"\nRunning Lahuta on '{inpath}' (acc '{acc}', {maxfrag} fragments):")

# # Change directory
# Not necessary now since the job moves to this directory before running this script (to avoid delays from Python's module import due to the huge number of log files in tmp/_logs/)
# os.chdir(inpath)

# infiles = nsort(Return(f"ls -1U {inpath}").split("\n"))
infiles = nsort(Return(f"ls -1U *.cif").split("\n"))
if len(infiles) != maxfrag:
    Die(f"Error: Expected {maxfrag} fragments in '{inpath}', but found {len(infiles)}")

contacts = None
# tmpfiles = set()
# Minimum length: no non-neighbor contacts possible until length ≥ 3
if len(seq) > 2:
    for ciffile in tq(infiles):
        # print(f" >> {ciffile}")

        # # Get only .cif files
        # if not rx(r"\.cif$", ciffile):
        #     # Non-mmCIF file
        #     Die("Error: Unexpected non-mmCIF file found: {inpath}/{ciffile}")
        #     continue

        Log("cif files processed", ciffile)
        # Don't delete CIF files here (they will ultimately be removed by job.py)
        # if not Switch('debug'):
        #     # In debug mode, keep the mmCIF file (don't remove it after running)
        #     tmpfiles.add(ciffile)

        # Get fragment number for this mmCIF file
        m = rx(r"^AF-"+acc+r"-F(\d+)-model_v\d+\.cif$", ciffile)
        if m:
            frag = m[0]
        else:
            Die(f"Couldn't parse filename '{ciffile}'")
            
        # Open mmCIF file

        # Filter out a specific warning

        # *** Open Babel Warning  in PerceiveBondOrders
        #   Failed to kekulize aromatic bonds in OBMol::PerceiveBondOrders (title is AF-A0A087WUL8-F12)
        # 
        # ~/miniconda3/lib/python3.10/site-packages/MDAnalysis/lib/util.py:664: RuntimeWarning: Constructed NamedStream from a NamedStream
        #   warnings.warn("Constructed NamedStream from a NamedStream",

        # print(f" >> Loading '{ciffile}'")
        # with warnings.catch_warnings():
        #     warnings.filterwarnings("ignore", message=".+Open Babel Warning +in PerceiveBondOrders.+")
        #     warnings.filterwarnings("ignore", message=".+Failed to kekulize aromatic bonds in OBMol::PerceiveBondOrders.+")
        #     warnings.filterwarnings("ignore", message=".+Constructed NamedStream from a NamedStream.+")
        # universe = Universe("test_alphafold_db_cif_Q7NBS4/AF-Q7NBS4-F1-model_v2.cif")
        universe = Universe(ciffile)

        # Compute neighboring residues
        # print(f" >> Computing neighbors")
        neighbors = universe.compute_neighbors()

        # Get contacts (of all the types specified above)
        for type in types:
            # print(f" >> {type}")

            # Get contacts using Lahuta
            type_contacts = GetContacts(type, frag, universe, neighbors)

            # Format data frame
            type_contacts = FormatContacts(type_contacts)

            # Append
            # contacts = pd.concat([contacts, type_contacts])
            contacts = pd.concat([contacts, type_contacts], ignore_index=True)

            # # Reset index
            # contacts = contacts.reset_index()

        # # Order rows across contact types
        # c = c.sort_values(["acc", "site1", "site2", "type", "groupid1", "groupid2", "dist"], ignore_index=True)

        Log(f"ran lahuta for acc|frag", f"{acc}|{frag}")
        


    # Verify that all fragments have been run
    if (frag != maxfrag):
        Die(f"Error: Expected to have processed {maxfrag} fragments, but only processed {frag}")

    # Accession is complete (Lahuta run for all fragments):
    # Combine output across fragments (by using the union of all contacts (ignoring any from dubious regions within 200 aa of artificial termini), and averaging distances)
    contacts = CombineFragmentContacts(contacts)

    # Insert into table
    # if not Switch('debug'):
    InsertContacts(contacts)

    Show(lim=20)

if contacts is not None and len(contacts) > 0:
    Query(f"UPDATE {alphaseq} SET nocon=0 WHERE acc='{acc}' AND afdb='{afdb}'")
    print(f"Successfully inserted {len(contacts):,} contacts into table '{alphacon}'")
    print(f" >> Setting column 'nocon'=0 for acc '{acc}' (afdb={afdb}) in table '{alphaseq}'")
else:
    Query(f"UPDATE {alphaseq} SET nocon=1 WHERE acc='{acc}' AND afdb='{afdb}'")
    print(f"No contacts found for acc '{acc}' (afdb={afdb})")
    print(f" >> Setting column 'nocon'=1 for acc '{acc}' (afdb={afdb}) in table '{alphaseq}'")


print("\nDone!")
