#!/usr/bin/env python3
"""
alphacon_add_pae.py:
Parse AlphaFold PAE scores from JSON files into SQL table 'alphacon'.
Note: No fragment support since PAE scores aren't particularly useful for fragments (they're intended for full-length proteins and would be inaccurate).
"""

# TODO Should add PAE parsing to job_lahuta.py instead and get the PAE scores for a given protein there - would be much faster than adding them to a giant table later on, as done here

# Initialize
# import pandas as pd
# import numpy as np
# import networkx as nx
# import seaborn as sns
# import matplotlib.pyplot as mp
import tarfile
import gzip
import io
import json
import pickle
# from Bio import SeqIO
from blang_mysql import *
from blang import *

alphafrag = "alphafrag"
alphaseq = "alphaseq"
alphacon = "alphacon"

Args(0, " \n -refresh: Reset wanted list (do not use temporary file)\n -debug: Don't actually make any changes, just simulate", "")
paepath = "input/pae/besian"
tmpfile = "input/pae/besian_diagnostic/tmp-filelist-seqlen-accs.pkl"

# AlphaSync updated structures PAE scores
alphasync_paepath = "input/alphasync/pae"

afdb = 1
if Switch('alphasync'):
    # AlphaSync updated structures only (not in AlphaFold Protein Structure Database)
    afdb = 0
    paepath = alphasync_paepath

# Note - This script contains many experiments for speeding up import. Ultimately the fastest way was a join with a temporary JSON table (used here).

# fifofile = "tmp/json_fifo"
# fifofile = "../../mysql/backup/tmp-json-fifo.txt"
# fifofile_server = "/mnt/backup/tmp/json_fifo" # Server refuses due to SHOW VARIABLES LIKE 'secure_file_priv';
# fifofile_server = "/mnt/backup/tmp-json-fifo.txt"
# fifofile_server = "json_fifo"
# fifofile_server = "/var/lib/mysql-files/tmp/"

# # Takes very long to run (affects 3,423,869,501 rows)
# # (3,419,234,111 will have PAE scores after this script is run)
# if not Switch('debug'):
#     print(f"Clearing 'pae' column in table '{alphacon}'...")
#     query = Query(f"UPDATE {alphacon} SET pae=NULL")
#     print(f"Rows affected: {Comma(Numrows(query))}")

# If I wanted an "alphapae" table, MySQL has a JSON data type that sounds great:
# https://dev.mysql.com/doc/refman/8.3/en/json.html
# "Optimized storage format. JSON documents stored in JSON columns are converted to an internal format that permits quick read access to document elements. When the server later must read a JSON value stored in this binary format, the value need not be parsed from a text representation. The binary format is structured to enable the server to look up subobjects or nested values directly by key or array index without reading all values before or after them in the document."
# >> Can avoid text parsing
# >> However: Simply adding PAE scores to alphacon is sufficient.
# >> Actually using JSON in a temporary table now below to update alphacon faster through JOINs.


# Start
Starttime()

# Get list of .json.gz files to parse (should take ~5 minutes)
print(f"\nGetting list of .json.gz archives in '{paepath}'...")
Time(1)
if not Exists(tmpfile) or Switch('refresh') or Switch('alphasync'):
    
    # First 10 files only
    # infiles = nsort(Return(f"find {paepath} -name '*.json.gz' | head").split("\n"))
    # First 100,000 files only
    # infiles = nsort(Return(f"find {paepath} -name '*.json.gz' | head -n 100000").split("\n"))
    # All files
    infiles = nsort(Return(f"find {paepath} -name '*.json.gz'").split("\n"))

    # Get wanted accs from alphaseq
    Time(2)
    print(f"\nGetting list of 'wanted' accessions (that consist of a single fragment and have contacts) and their sequence lengths (for verifying that the PAE matrix dimensions equal protein length) from table '{alphaseq}':")
    seqlen = {}
    wanted_accs = set()
    # Require frags=1 since PAE scores aren't meaningful for fragmented proteins (>2700 aa)
    for acc, length in tq(Query(f"SELECT acc, LENGTH(seq) FROM {alphaseq} WHERE frags=1 AND afdb={afdb} AND nocon!=1")):
        # Only add to wanted list if alphacon already has PAE scores for every contact for this protein
        # if Numrows(Query(f"SELECT 1 FROM {alphacon} WHERE acc='{acc}' AND pae IS NULL LIMIT 1")) > 0:
        # https://dev.mysql.com/doc/refman/8.4/en/exists-and-not-exists-subqueries.html
        if Switch('alphasync'):
            # If updating: also update rows where pae is already set (not NULL)
            if FetchOne(Query(f"SELECT EXISTS (SELECT * FROM {alphacon} WHERE acc='{acc}' AND afdb={afdb})")) == 1:
                wanted_accs.add(acc)
                seqlen[acc] = length
        else:
            if FetchOne(Query(f"SELECT EXISTS (SELECT * FROM {alphacon} WHERE acc='{acc}' AND afdb={afdb} AND pae IS NULL)")) == 1:
                wanted_accs.add(acc)
                seqlen[acc] = length
    Time(2)

    # Filter infiles to only include wanted_accs
    print(f"\nFiltering input file list to only retain wanted accessions:")
    Time(3)
    tmpinfiles = infiles
    infiles = []
    accs = {}
    for infile in tq(tmpinfiles):

        # Get acc from source file name (e.g. input/pae/besian/store/afdb/data/x8772/AF-P53587-F1-predicted_aligned_error_v4.json.gz)
        m = rx(r"AF-(\w+(-\d+)?)-F(\d+)-predicted_aligned_error_v\d+\.json\.gz$", infile)
        if not m:
            Die(f"Error: Couldn't parse accession from '{infile}'")

        acc = m[0]
        frag = int(m[2])

        if acc not in wanted_accs:
            Log(f"unwanted infile skipped for infile (skipped)", infile)
            Log(f"unwanted infile skipped for acc (skipped)", acc)
            continue
        else:
            Log(f"wanted infile kept for infile (kept)", infile)
            Log(f"wanted infile kept for acc (kept)", acc)
            infiles.append(infile)
            accs[infile] = acc

        if frag != 1:
            Die(f"Error: Fragment is '{frag}' for acc '{acc}' (no support for fragment numbers higher than 1 in this script)")
    Time(3)

    print(f"Number of wanted accessions: {len(wanted_accs):,}\n")

    if not Switch('alphasync'):
        # Write to file
        with open(tmpfile, "wb") as f:
            pickle.dump([infiles, seqlen, accs], f)
        print(f" >> Wrote to '{tmpfile}'")

else:

    # Load from file
    with open(tmpfile, "rb") as f:
        (infiles, seqlen, accs) = pickle.load(f)
    print(f" >> Loaded from '{tmpfile}'")

Time(1)

if Switch('debug'):
    infiles = ["input/pae/besian/store/afdb/data/x0006/AF-A0A1S3ETA8-F1-predicted_aligned_error_v4.json.gz"] * 1000
    wanted_accs = ["A0A1S3ETA8"]
    # infiles = infiles[:100]

# Create temporary table
# print("\n >> Create temp table\t", end="") if Switch('debug') else None
# Starttime() if Switch('debug') else None
tmptable = f"alphacon_tmp_update"
# # PRIMARY and indexes
# q = f"""CREATE TEMPORARY TABLE {tmptable} (
# `id` int unsigned NOT NULL AUTO_INCREMENT,
# `acc` char(10) DEFAULT NULL,
# `site1` mediumint DEFAULT NULL,
# `site2` mediumint DEFAULT NULL,
# `pae` float DEFAULT NULL,
# PRIMARY KEY (`id`),
# KEY `Acc` (`acc`),
# KEY `Site1` (`site1`),
# KEY `Site2` (`site2`)
# );"""
# # Indexes
# q = f"""CREATE TEMPORARY TABLE {tmptable} (
# `acc` char(10) DEFAULT NULL,
# `site1` mediumint DEFAULT NULL,
# `site2` mediumint DEFAULT NULL,
# `pae` float DEFAULT NULL,
# KEY `Acc` (`acc`),
# KEY `Site1` (`site1`),
# KEY `Site2` (`site2`)
# );"""
# # Indexes, without acc
# q = f"""CREATE TEMPORARY TABLE {tmptable} (
# `site1` mediumint DEFAULT NULL,
# `site2` mediumint DEFAULT NULL,
# `pae` float DEFAULT NULL,
# KEY `Site1` (`site1`),
# KEY `Site2` (`site2`)
# );"""
# # Composite index
# q = f"""CREATE TEMPORARY TABLE {tmptable} (
# `acc` char(10) DEFAULT NULL,
# `site1` mediumint DEFAULT NULL,
# `site2` mediumint DEFAULT NULL,
# `pae` float DEFAULT NULL,
# KEY `Acc` (`acc`, `site1`, `site2`)
# );"""
# # Composite index, without acc
# q = f"""CREATE TEMPORARY TABLE {tmptable} (
# `site1` mediumint DEFAULT NULL,
# `site2` mediumint DEFAULT NULL,
# `pae` float DEFAULT NULL,
# KEY `Sites` (`site1`, `site2`)
# )"""
# # No indexes
# q = f"""CREATE TEMPORARY TABLE {tmptable} (
# `acc` char(10) DEFAULT NULL,
# `site1` mediumint DEFAULT NULL,
# `site2` mediumint DEFAULT NULL,
# `pae` float DEFAULT NULL
# );"""
# # No indexes, without acc
# q = f"""CREATE TEMPORARY TABLE {tmptable} (
# `site1` mediumint DEFAULT NULL,
# `site2` mediumint DEFAULT NULL,
# `pae` float DEFAULT NULL
# );"""
# JSON (can't be used with MEMORY engine, but performance hit should be small)
# InnoDB is slightly faster here than MyISAM, I think because it does more in-memory.
q = f"""CREATE TEMPORARY TABLE {tmptable} (
`pae` json DEFAULT NULL
) engine=InnoDB"""
# ) engine=MyISAM"""
Query(q)
# Stoptime() if Switch('debug') else None

# # Create FIFO pipe so MySQL can read compressed files through it
# if Exists(fifofile):
#     Run("Remove existing FIFO pipe", f"rm -fv {fifofile}")
# if not Exists(fifofile):
#     Run("Create FIFO pipe", f"mkfifo {fifofile}")
#     # Run("Chmod FIFO pipe", f"chmod 777 {fifofile}")



# Main loop: Parse .json.gz files
print(f"\nParsing PAE scores in .json.gz archives in '{paepath}' and updating table '{alphacon}':")
affected = 0
for infile in tq(infiles):
    # print(f" >> {infile}") if Switch('debug') else None
    
    # # Get acc from source file name (e.g. input/pae/besian/store/afdb/data/x8772/AF-P53587-F1-predicted_aligned_error_v4.json.gz)
    # m = rx(r"AF-(\w+(-\d+)?)-F(\d+)-predicted_aligned_error_v\d+\.json\.gz$", infile)
    # if not m:
    #     Die(f"Error: Couldn't parse accession from '{infile}'")
    # acc = m[0]

    acc = accs[infile]

    # Simply overwrite instead
    # # Check if acc has any contacts where PAE has not been set yet ('pae' IS NULL) in table 'alphacon'
    # if Numrows(Query(f"SELECT * FROM {alphacon} WHERE acc='{acc}' AND pae IS NULL LIMIT 1")) == 0:
    #     Log(f"acc already had PAE scores for every contact for acc (skipped)", acc)
    #     continue

    # # Skip unwanted accessions (i.e. those that are not in table 'alphaseq')
    # if acc not in wanted_accs:
    #     Log(f"unwanted acc skipped for acc", acc)
    #     continue
    # else:
    #     Log(f"found file for wanted acc for acc", acc)

    # print(f" >> {infile}") if Switch('debug') else None

    # # Parse JSON from .json.gz file
    # print("\n >> Load JSON\t", end="") if Switch('debug') else None
    # Starttime() if Switch('debug') else None
    # with gzip.open(infile, "rt") as f:
    #     data = json.load(f)
    # Stoptime() if Switch('debug') else None
    # # paes = data[0]["predicted_aligned_error"]
    # # maxpae = data[0]["max_predicted_aligned_error"]
    # # if maxpae != 31.75:
    # #     Die(f"Error: Expected max max_predicted_aligned_error to be '31.75', but got '{maxpae}' instead in '{infile}'")
    # # Not sure what maxpae is actually good for:
    # # It appears to always be 31.75.
    # # https://alphafold.ebi.ac.uk/faq
    # # "predicted_aligned_error: The PAE value of the residue pair, rounded to the closest integer. For the PAE value at position (i, j), i is the residue on which the structure is aligned, j is the residue on which the error is predicted."
    # # https://alphafold.ebi.ac.uk/faq
    # # "max_predicted_aligned_error: A number that denotes the maximum possible value of PAE. The smallest possible value of PAE is 0."
    # # https://alphafold.ebi.ac.uk/entry/G5EBL2#help
    # # Section: "How to interpret the Predicted Aligned Error"
    # # [PAE scores are] "the expected distance error in Ångströms (Å), ranging from 0 Å to an arbitrary cut-off of 31 Å."
    # # >> I've seen values of 32, so this isn't even correct (I suppose ≥31.5 gets rounded to 32). The actual maximum must be 31.75.
    # # Looking directly at AF 2.3.2 PAE results files, the diagonal has values like 0.3, so the function used must indeed be a normal round().

    # # Load JSON directly from .json.gz file through FIFO pipe
    # print("\n >> Load JSON\t", end="") if Switch('debug') else None
    # Starttime() if Switch('debug') else None
    # # with gzip.open(infile, "rt") as f:
    # #     data = json.load(f)
    # # gunzip in background (will terminate as soon as LOAD DATA LOCAL INFILE is done)
    # os.system(f"gunzip -c {infile} > {fifofile} &")
    # Query(f"TRUNCATE TABLE {tmptable}; LOAD DATA LOCAL INFILE '{fifofile}' INTO TABLE {tmptable} (pae)")
    # # Query(f"TRUNCATE TABLE {tmptable}; LOAD DATA INFILE '{fifofile_server}' INTO TABLE {tmptable} (pae)")
    # Stoptime() if Switch('debug') else None

    # Load JSON into memory, but don't parse it
    # print(" >> Load JSON\t", end="") if Switch('debug') else None
    # Starttime() if Switch('debug') else None
    with gzip.open(infile, "rt") as f:
        data = f.read()
    # Stoptime() if Switch('debug') else None

    # # Verify PAE matrix size
    # if seqlen[acc] != len(paes):
    #     Die(f"Error: Expected PAE matrix height '{seqlen[acc]}', but got '{len(paes)}' instead in '{infile}'")
    # if seqlen[acc] != len(paes[0]):
    #     Die(f"Error: Expected PAE matrix width '{seqlen[acc]}', but got '{len(paes[0])}' instead in '{infile}'")

# Previous method, with one UPDATE query per site1|site2 combination:
#     # Insert PAE scores for alphacon contacts
#     for site1, site2 in Query(f"SELECT DISTINCT site1, site2 FROM {alphacon} WHERE acc='{acc}'"):
# 
#         # Verify PAE matrix size
#         if seqlen[acc] != len(paes[site1-1]):
#             Die(f"Error: Expected PAE matrix width '{seqlen[acc]}', but got '{len(paes[site1-1])}' instead in '{infile}'")
# 
#         # # Verify that PAE matrix is symmetric
#         # It's not symmetric:
#         # https://alphafold.ebi.ac.uk/faq
#         # "predicted_aligned_error: The PAE value of the residue pair, rounded to the closest integer. For the PAE value at position (i, j), i is the residue on which the structure is aligned, j is the residue on which the error is predicted."
#         # https://alphafold.ebi.ac.uk/entry/G5EBL2#help
#         # Section: "How to interpret the Predicted Aligned Error"
#         # "Note that the PAE scores are asymmetrical, meaning there might be variations in PAE values between (x,y) and (y,x) positions. This is particularly relevant for loop regions with highly uncertain orientations, as seen on the DNA topoisomerase 3 (Q8T2T7)."
#         # >> Use maximum of the two values
#         # if paes[site1][site2] != paes[site2][site1]:
#         #     d()
#         #     Die(f"Error: PAE matrix is not symmetric for acc '{acc}' site1 '{site1}' site2 '{site2}' in '{infile}'")
#         mypae = max(paes[site1-1][site2-1], paes[site2-1][site1-1])
# 
#         # Update alphacon with PAE score
#         if not Switch('debug'):
#             query = Query(f"UPDATE {alphacon} SET pae='{mypae}' WHERE acc='{acc}' AND site1='{site1}' AND site2='{site2}'")
#             affected += Numrows(query)
#         else:
#             State(f"UPDATE {alphacon} SET pae='{mypae}' WHERE acc='{acc}' AND site1='{site1}' AND site2='{site2}'")
#             affected += 1
    


#     # Combined method using CASE, with one UPDATE query for all site1|site2 combinations
#     
#     # Insert PAE scores for alphacon contacts
#     # print("\n >> Build pae_values\t", end="")
#     # Starttime()
#     pae_values = {}
#     # Initialise pae_values dictionary
#     for i in range(0, 33):
#         pae_values[i] = []
#     # Fill it
#     for site1 in range(0, seqlen[acc]):
#         for site2 in range(0, seqlen[acc]):
#             if site1 == site2:
#                 continue
#             # mypae = max(paes[site1][site2], paes[site2][site1])
#             mypae = max(data[0]["predicted_aligned_error"][site1][site2], data[0]["predicted_aligned_error"][site2][site1])
#             pae_values[mypae].append((site1 + 1, site2 + 1))
#     # Stoptime()
#     
# 
# 
#     # Update alphacon with PAE score in a single query per acc
#     # Example:
#     # UPDATE alphacon
#     # SET pae = (
#     #   CASE 
#     #     WHEN (site1, site2) IN ((site1_value1, site2_value1), (site1_value2, site2_value2), ...) 
#     #     THEN 'value1' -- for (site1_value1, site2_value1)
#     #     WHEN (site1, site2) IN ((site1_value3, site2_value3), (site1_value4, site2_value4), ...) 
#     #     THEN 'value2' -- for (site1_value3, site2_value3)
#     #     -- Add more WHEN clauses for additional combinations with their corresponding pae values
#     #     ELSE NULL -- optional: handle unmatched combinations (if needed)
#     #   END
#     # )
#     # WHERE acc='A0A1S3EUJ9';
#     # Build query
#     # print("\n >> Build query\t", end="")
#     # Starttime()
#     cases = "\n".join(["WHEN (site1,site2) IN (" + ",".join([f"({site1},{site2})" for site1, site2 in pae_values[pae]]) + f") THEN {pae}" for pae in range(0, 33) if len(pae_values[pae]) > 0])
#     cases += " ELSE NULL END\n"
#     q = f"UPDATE {alphacon} SET pae=(CASE {cases}) WHERE acc='{acc}'"
#     # Stoptime()
#     
#     # print("\n >> Run query\t", end="")
#     # Starttime()
#     # if not Switch('debug'):
#     query = Query(q)
#     affected += Numrows(query)
#     # else:
#     #     State(q)
#     #     affected += 1
#     # Stoptime()
    


    # # Simpler method using individual queries for each site1|site2 combination
    # print("\n >> Run queries\t", end="")
    # Starttime()
    # for site1 in range(0, seqlen[acc]):
    #     for site2 in range(0, seqlen[acc]):
    #         if site1 == site2:
    #             continue
    #         # mypae = max(paes[site1][site2], paes[site2][site1])
    #         mypae = max(data[0]["predicted_aligned_error"][site1][site2], data[0]["predicted_aligned_error"][site2][site1])
    #         # Build query
    #         q = f"UPDATE {alphacon} SET pae='{mypae}' WHERE acc='{acc}' AND site1='{site1+1}' AND site2='{site2+1}'"
    #         # Update alphacon with PAE score
    #         if not Switch('debug'):
    #             query = Query(q)
    #             affected += Numrows(query)
    #         else:
    #             State(q)
    #             affected += 1
    # Stoptime()
    
    

    
    # # Simpler method using individual queries for each site1|site2 combination, but sent as a single multi-statement
    # print("\n >> Run queries\t", end="")
    # Starttime()
    # q = ""
    # for site1 in range(0, seqlen[acc]):
    #     for site2 in range(0, seqlen[acc]):
    #         if site1 == site2:
    #             continue
    #         # mypae = max(paes[site1][site2], paes[site2][site1])
    #         mypae = max(data[0]["predicted_aligned_error"][site1][site2], data[0]["predicted_aligned_error"][site2][site1])
    #         # Build query
    #         q += f"UPDATE {alphacon} SET pae='{mypae}' WHERE acc='{acc}' AND site1='{site1+1}' AND site2='{site2+1}';\n"
    # # Update alphacon with PAE scores
    # if not Switch('debug'):
    #     query = Query(q)
    #     # affected += Numrows(query)
    # else:
    #     State(q)
    #     # affected += 1
    # Stoptime()




#     # Combined method using a temporary table and joining it with alphacon
#     # Default engine for temporary tables is MEMORY (which is fast)
# #     print("\n >> Clear temp table\t", end="")
# #     Starttime()
# #     Query(f"TRUNCATE TABLE {tmptable}")
# #     Stoptime()
# 
#     print(" >> Build query\t", end="") if Switch('debug') else None
#     Starttime() if Switch('debug') else None
#     # q = f"INSERT INTO {tmptable} (id, acc, site1, site2, pae) VALUES "
#     # q = f"INSERT INTO {tmptable} (acc, site1, site2, pae) VALUES "
#     q = f"TRUNCATE TABLE {tmptable};INSERT INTO {tmptable} (site1,site2,pae) VALUES "
# 
#     # Build query
#     for site1 in range(0, seqlen[acc]):
#         # site1 is always < site2, hence we can start site2 at site1 + 1:
#         for site2 in range(site1 + 1, seqlen[acc]):
#             if site1 == site2:
#                 Die("Error: site1 == site2")
#             # Use the "worst possible" (highest) PAE value since site1 and site2 are equally important for a given contact
#             # (As mentioned elsewhere, the PAE matrix can be slightly asymmetric since it's the predicted aligned error at residue x for residue y, and vice versa)
#             mypae = max(data[0]["predicted_aligned_error"][site1][site2], data[0]["predicted_aligned_error"][site2][site1])
#             # with id column
#             # q += f"(NULL, '{acc}', {site1+1}, {site2+1}, {mypae}),"
#             # without id column
#             # q += f"('{acc}', {site1+1}, {site2+1}, {mypae}),"
#             # without id and acc columns
#             q += f"({site1+1},{site2+1},{mypae}),"
#     # Remove last character of q (the trailing comma)
#     q = q[:-1]
#     Stoptime() if Switch('debug') else None
#     
#     # Fill temporary table
#     print("\ >> Fill temp table\t", end="") if Switch('debug') else None
#     Starttime() if Switch('debug') else None
#     query = Query(q)
#     Stoptime() if Switch('debug') else None
#     
#     # Use temporary table to set PAE values in alphacon
#     print(" >> Update alphacon\t", end="") if Switch('debug') else None
#     Starttime() if Switch('debug') else None
#     # if not Switch('debug'):
#     # query = Query(f"UPDATE {alphacon} c, {tmptable} t SET c.pae=t.pae WHERE c.acc='{acc}' AND c.acc=t.acc AND c.site1=t.site1 AND c.site2=t.site2")
#     query = Query(f"UPDATE {alphacon} c, {tmptable} t SET c.pae=t.pae WHERE c.acc='{acc}' AND c.site1=t.site1 AND c.site2=t.site2")
#     affected += Numrows(query)
#     # else:
#     #     State(f"UPDATE {alphacon} c, {tmptable} t SET c.pae=t.pae WHERE c.acc='{acc}' AND c.acc=t.acc AND c.site1=t.site1 AND c.site2=t.site2")
#     #     affected += 1
#     Stoptime() if Switch('debug') else None



    # Use JSON directly
    # SELECT c.acc, c.site1, c.site2, c.pae, GREATEST(JSON_EXTRACT(t.pae, CONCAT('$[0].predicted_aligned_error[', site1-1, '][', site2-1, ']')), JSON_EXTRACT(t.pae, CONCAT('$[0].predicted_aligned_error[', site2-1, '][', site1-1, ']'))) AS t_pae FROM alphacon c, alphacon_tmp t WHERE c.acc='A0A1S3ETA8' HAVING c.pae!=t_pae
    # query = Query(f"SELECT c.acc, c.site1, c.site2, c.pae, GREATEST(JSON_EXTRACT(t.pae, CONCAT('$[0].predicted_aligned_error[', site1-1, '][', site2-1, ']')), JSON_EXTRACT(t.pae, CONCAT('$[0].predicted_aligned_error[', site2-1, '][', site1-1, ']'))) AS t_pae FROM alphacon c, alphacon_tmp t WHERE c.acc='{acc}'")
    # print(" >> Update alphacon\t", end="") if Switch('debug') else None
    # Starttime() if Switch('debug') else None
    # query = Query(f"UPDATE {alphacon} c, {tmptable} t SET c.pae=GREATEST(JSON_EXTRACT(t.pae, CONCAT('$[0].predicted_aligned_error[', site1-1, '][', site2-1, ']')), JSON_EXTRACT(t.pae, CONCAT('$[0].predicted_aligned_error[', site2-1, '][', site1-1, ']'))) WHERE c.acc='{acc}'")
    # Escape : in data
    # data = Esc(data)
    # Slightly faster
    data = data.replace(":", "\\:")
    Query(f"TRUNCATE {tmptable}")
    Query(f"INSERT INTO {tmptable} SET pae='{data}'")
    # query = Query(f"UPDATE {alphacon} c, {tmptable} t SET c.pae=GREATEST(CAST(JSON_UNQUOTE(JSON_EXTRACT(t.pae, CONCAT('$[0].predicted_aligned_error[', site1-1, '][', site2-1, ']'))) AS SIGNED), CAST(JSON_UNQUOTE(JSON_EXTRACT(t.pae, CONCAT('$[0].predicted_aligned_error[', site2-1, '][', site1-1, ']'))) AS SIGNED)) WHERE c.acc='{acc}' AND c.afdb={afdb}")
    # For AlphaSync structures (afdb=0): need to use ROUND() to round the values. Also removed unnecessary CAST to SIGNED that didn't work with FLOAT values.
    query = Query(f"UPDATE {alphacon} c, {tmptable} t SET c.pae=ROUND(GREATEST(JSON_UNQUOTE(JSON_EXTRACT(t.pae, CONCAT('$[0].predicted_aligned_error[', site1-1, '][', site2-1, ']'))), JSON_UNQUOTE(JSON_EXTRACT(t.pae, CONCAT('$[0].predicted_aligned_error[', site2-1, '][', site1-1, ']'))))) WHERE c.acc='{acc}' AND c.afdb={afdb}")
    # # Shorter JSON operator notation in MySQL 8.3 that includes unquoting (->>), see https://dev.mysql.com/doc/refman/8.3/en/json-search-functions.html
    # Still on MySQL 8.0 though
    # query = Query(f"UPDATE {alphacon} c, {tmptable} t SET c.pae=GREATEST(t.pae->>CONCAT('$[0].predicted_aligned_error[', site1-1, '][', site2-1, ']'), t.pae->>CONCAT('$[0].predicted_aligned_error[', site2-1, '][', site1-1, ']')) WHERE c.acc='{acc}'")
    affected += Numrows(query)
    # Stoptime() if Switch('debug') else None
    

    
    # print(" >> Log\t", end="") if Switch('debug') else None
    # Starttime() if Switch('debug') else None
    if Numrows(query) > 0:
        Log(f"successfully added PAE scores for acc", acc)
    # Stoptime() if Switch('debug') else None

# # Remove FIFO pipe
# Run("Remove FIFO pipe", f"rm -fv {fifofile}")

Show(lim=20)

print(f"\nRows affected: {Comma(affected)}")

# Optimize(alphacon)
Stoptime()

print("\nDone!")
