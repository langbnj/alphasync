#!/usr/bin/env python3
"""
Main: Main script
"""

# Initialize
import pandas as pd
# import numpy as np
# import networkx as nx
# import seaborn as sns
# import matplotlib.pyplot as mp
from blang_mysql import *
from blang import *

physsourcedb = "blang"
sourcedb = "alphasync"
destdb = "alphasync_compact"

# # Args(0, "", "")
# var = Args(1, "[species]", "human")
# infile = f"input/input.txt"

# Run individual table scripts (for the larger, more complex ones where I changed the indexes for InnoDB):

# In parallel:
Time(1)
Run("Migrate alphacon             to alphasync_compact", "~/scripts/qsub.sh migrate_alphasync_compact_alphacon.py")   # alphasync.alphacon is now the only version (already compact). As this is ~280 GB, I won't be copying it to alphasync_compact here. I'll only create a view.
Run("Migrate alphafrag            to alphasync_compact", "~/scripts/qsub.sh migrate_alphasync_compact_alphafrag.py")  # alphafrag is not used by the website, but it will eventually be used by the API.
Run("Migrate alphamap             to alphasync_compact", "~/scripts/qsub.sh migrate_alphasync_compact_alphamap.py")
Run("Migrate alphasa              to alphasync_compact", "~/scripts/qsub.sh migrate_alphasync_compact_alphasa.py")
Run("Migrate alphaseq             to alphasync_compact", "~/scripts/qsub.sh migrate_alphasync_compact_alphaseq.py")
Run("Migrate alphastats           to alphasync_compact", "~/scripts/qsub.sh migrate_alphasync_compact_alphastats.py")
Run("Migrate alphauniprot         to alphasync_compact", "~/scripts/qsub.sh migrate_alphasync_compact_alphauniprot.py")
Run("Migrate alphauniprot_species to alphasync_compact", "~/scripts/qsub.sh migrate_alphasync_compact_alphauniprot_species.py")
Run("Migrate alphauniprot_symbols to alphasync_compact", "~/scripts/qsub.sh migrate_alphasync_compact_alphauniprot_symbols.py")

Waitforjobs()
Time(1)

# # Other tables (compara_species etc.)
# 
# # Check which tables exist in schema alphasync (which actually only consists of views), but not yet in alphasync_compact
# source = FetchSet(Query(f"SELECT TABLE_NAME FROM information_schema.TABLES WHERE TABLE_SCHEMA='{sourcedb}'"))
# dest = FetchSet(Query(f"SELECT TABLE_NAME FROM information_schema.TABLES WHERE TABLE_SCHEMA='{destdb}'"))
# 
# for table in nsort(source-dest):
#   print(f" >> {table}")
# 
#   # Starttime()
#   # Query(f"""DROP TABLE IF EXISTS alphasync_compact.{table};""", loud=1)
#   # Stoptime();
# 
#   # Get CREATE statement
#   create = FetchOne(Query(f"""SHOW CREATE TABLE {physsourcedb}.{table};""", loud=1))[1]
#   print(f"   >> {create}")
#   # Add {destdb} to CREATE
#   create = re.sub(r"CREATE TABLE ", f"CREATE TABLE {destdb}.", create)
#   # Remove AUTO_INCREMENT
#   create = re.sub(r"AUTO_INCREMENT=\d+ ", "", create)
#   # Replace ENGINE with InnoDb
#   create = re.sub(r"ENGINE=MyISAM ", "ENGINE=InnoDB ", create)
#   # print(f"   >> {create}")
#   
#   # Run CREATE
#   Query(create)
# 
#   Starttime()
#   Query(f"""INSERT INTO {destdb}.{table} SELECT * FROM {physsourcedb}.{table};""", loud=1)
#   Stoptime();
# 
#   Starttime()
#   Query(f"""ANALYZE TABLE {destdb}.{table};""", loud=1)
#   Stoptime();
# 
#   Starttime()
#   Query(f"""OPTIMIZE TABLE {destdb}.{table};""", loud=1)
#   Stoptime();


Show(lim=0)

print("\nDone!")

