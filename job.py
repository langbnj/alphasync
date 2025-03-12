#!/usr/bin/env python3
"""
job.py: Job script that sequentially runs DSSP and Lahuta on individual fragment files.

- Calculates residue-level accessible surface area (dssp.py) and contacts (lahuta.py).
- Combines output across fragments using combine_fragments_dssp.py and combine_fragments_lahuta.py.
- Parses its results into the 'alphasa' and 'alphacon' MySQL tables.
- Removes temporary tmp/{acc} directory once complete.

"""

# Initialize
from blang_mysql import *
from blang import *

# alphafrag = "alphafrag"     # SQL table with fragment protein sequences (>2700 aa proteins get split into 1400 aa fragments with a step size of 200 in AlphaFold DB, for human only - other species don't have results for >2700 aa proteins)
alphaseq = "alphaseq"       # SQL table with complete protein sequences
alphasa = "alphasa"         # SQL table with residue-level accessible surface area values from DSSP
alphacon = "alphacon"       # SQL table with residue-level contacts from Lahuta

(acc, maxfrag) = Args(2, f"[UniProt accession] [Number of fragments]\n\n -alphasync: Treat this as an updated protein from AlphaSync (i.e. not in AFDB, afdb=0) and delete its data from tables '{alphasa}' and '{alphacon}' (and reset its 'nocon' (no contacts) flag in table '{alphaseq}') before re-running it\n -alphakeep: Keep existing AlphaSync data: don't delete from tables '{alphasa}' and '{alphacon}', and don't reset 'nocon' (no contacts) flag in table '{alphaseq}'. Use this for re-runs of main.py -alphasync jobs.", "A0A087WUL8 14")

inpath = f"../{acc}"

afdb = 1
tmp_alphasync = ""
if Switch('alphakeep'):
    SetSwitch('alphasync')
if Switch('alphasync'):
    afdb = 0
    tmp_alphasync = " -alphasync"

    if not Switch('debug'):
        if not Switch('alphakeep'):
            print(f"\nClearing AlphaSync data for acc '{acc}' before re-running it:")

            print(f"\n >> Clearing AlphaSync data for acc '{acc}' from table '{alphasa}'...")
            query = Query(f"DELETE FROM {alphasa} WHERE acc='{acc}' AND afdb={afdb}")
            print(f"   >> Rows affected: {Numrows(query):,}")

            print(f"\n >> Clearing AlphaSync data for acc '{acc}' from table '{alphacon}'...")
            query = Query(f"DELETE FROM {alphacon} WHERE acc='{acc}' AND afdb={afdb}")
            print(f"   >> Rows affected: {Numrows(query):,}")

            print(f"\n >> Clearing AlphaSync data for acc '{acc}' from table '{alphaseq} (resetting its 'nocon' (no contacts) flag)'...")
            query = Query(f"UPDATE {alphaseq} SET nocon=NULL WHERE acc='{acc}' AND afdb={afdb}")
            print(f"   >> Rows affected: {Numrows(query):,}")

            print()
        else:
            print(f"\nSwitch -alphakeep is active: Not clearing AlphaSync data for acc '{acc}'")

            print(f"\n >> AlphaSync data for acc '{acc}' in table '{alphasa}':")
            query = Query(f"SELECT * FROM {alphasa} WHERE acc='{acc}' AND afdb={afdb}")
            print(f"   >> Rows: {Numrows(query):,}")

            print(f"\n >> AlphaSync data for acc '{acc}' in table '{alphacon}'...")
            query = Query(f"SELECT * FROM {alphacon} WHERE acc='{acc}' AND afdb={afdb}")
            print(f"   >> Rows: {Numrows(query):,}")

            print(f"\n >> AlphaSync data for acc '{acc}' in table '{alphaseq} (where 'nocon' (no contacts) flag is not NULL)'...")
            query = Query(f"SELECT * FROM {alphaseq} WHERE acc='{acc}' AND afdb={afdb} AND nocon IS NOT NULL")
            print(f"   >> Rows: {Numrows(query):,}")

            print()
    
# Check if data already exists for this acc in table 'alphasa'
# query_alphasa = Query(f"SELECT id FROM {alphasa} WHERE acc='{acc}' LIMIT 1")
# alphasa also contains dihedral angles and proline isomerization states, so require these as well
query_alphasa = Query(f"SELECT * FROM {alphasa} WHERE acc='{acc}' AND afdb='{afdb}' AND iso IS NOT NULL LIMIT 1")

# Check if data already exists for this acc in table 'alphacon'
query_alphacon = Query(f"SELECT * FROM {alphacon} WHERE acc='{acc}' AND afdb='{afdb}' LIMIT 1")

# ...and exit if both exist already
if (Numrows(query_alphacon) == 1) and (Numrows(query_alphasa) == 1):
    print(f"Data already existed for acc '{acc}' (afdb={afdb}) in tables '{alphasa}' and '{alphacon}', exiting (skip)!")
    sys.exit()
elif (Numrows(query_alphacon) == 0) and (Numrows(query_alphasa) == 0):
    print(f"No data yet for acc '{acc}' (afdb={afdb}) from in tables '{alphasa}' and '{alphacon}', starting!")
else:
    # ...throw an error if only one of the two exists (^ is XOR)...
    if (Numrows(query_alphacon) == 1) ^ (Numrows(query_alphasa) == 1):
        Die(f"Error: acc '{acc}' (afdb={afdb}) is present in one, but not both of alphacon and alphasa: alphacon {Numrows(query_alphacon)}, alphasa {Numrows(query_alphasa)}")
        # # For dihedral angles, make this a print statement only:
        # print(f"Warning: acc '{acc}' is present in one, but not both of alphacon and alphasa: alphacon {Numrows(query_alphacon)}, alphasa {Numrows(query_alphasa)}")



# Start

# # Change directory
# Not necessary now since the job moves to this directory before running this script (to avoid delays from Python's module import due to the huge number of log files in tmp/_logs/)
# os.chdir(inpath)

# Preparations

# Get list of temporary CIF files that are present in this temporary directory (inpath)
infiles = nsort(Return(f"ls -1U *.cif").split("\n"))
if len(infiles) != maxfrag:
    Die(f"Error: Expected {maxfrag} fragments in '{inpath}', but found {len(infiles)}")

tmpfiles = set()
for ciffile in tq(infiles):
    Log("temporary cif files removed", ciffile)
    tmpfiles.add(ciffile)

    # Get fragment number for this mmCIF file
    m = rx(r"^AF-"+acc+r"-F(\d+)-model_v\d+\.cif$", ciffile)
    if m:
        frag = m[0]
    else:
        Die(f"Couldn't parse filename '{ciffile}'")

# Verify that all fragments are present
if (frag != maxfrag):
    Die(f"Error: Expected to find temporary CIF files for {maxfrag} fragments, but only found {frag}")



# Run individual tasks (DSSP, Lahuta, dihedral angles)

# 1. Run DSSP
Run("Getting relative accessible surface areas using DSSP (for SQL table 'alphasa')", f"../../job_dssp.py {acc} {maxfrag}{tmp_alphasync}", silent=False)

# 2. Run Lahuta
# Filter out specific warnings
# *** Open Babel Warning  in PerceiveBondOrders
#   Failed to kekulize aromatic bonds in OBMol::PerceiveBondOrders (title is AF-A0A087WUL8-F12)
# 
# ~/miniconda3/lib/python3.10/site-packages/MDAnalysis/lib/util.py:664: RuntimeWarning: Constructed NamedStream from a NamedStream
#   warnings.warn("Constructed NamedStream from a NamedStream",
# Note that this requires running bash (sh's redirect syntax doesn't seem to support filtering only STDERR)
Run("Getting contacts using Lahuta (for SQL table 'alphacon')", f"bash -c \"../../job_lahuta.py {acc} {maxfrag}{tmp_alphasync} 2> >(grep -viP '(Open Babel Warning +in PerceiveBondOrders|Failed to kekulize aromatic bonds in OBMol::PerceiveBondOrders|Constructed NamedStream from a NamedStream|^==============================$|^$)'>&2)\"", silent=False)

# 3. Run Bio.PDB dihedral angle calculation
# Run("Getting dihedral angles (for SQL table 'alphasa')", f"../../job_dihedral_angles.py {acc} {maxfrag}", silent=False)
Run("Getting dihedral angles (for SQL table 'alphasa')", f"bash -c \"../../job_dihedral_angles.py {acc} {maxfrag}{tmp_alphasync} 2> >(grep -viP '(Open Babel Warning +in PerceiveBondOrders|Failed to kekulize aromatic bonds in OBMol::PerceiveBondOrders|Constructed NamedStream from a NamedStream|^==============================$|^$)'>&2)\"", silent=False)



# Clean up

# Delete temporary CIF files
if not Switch('debug'):
    print()
    for tmpfile in nsort(tmpfiles):
        print(f"Removing temporary file '{tmpfile}'")
        # This will occasionally produce "NotADirectoryError: [Errno 20] Not a directory" on nfs, but the file will still be deleted correctly
        os.remove(tmpfile)

    # Remove temporary directory for this accession (will throw an error if not empty)
    os.rmdir(inpath)

print("\nDone!")
