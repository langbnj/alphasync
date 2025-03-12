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

# table = "table"

# # Args(0, "", "")
# var = Args(1, "[species]", "human")
# infile = f"input/input.txt"

# Run queries

Starttime()
Query("""CREATE OR REPLACE VIEW alphasync_compact.alphacon AS SELECT * FROM alphasync.alphacon;""", loud=1)
Stoptime();

# Starttime()
# Query("""DROP TABLE IF EXISTS alphasync_compact.alphacon;""", loud=1)
# Stoptime();
# 
# Starttime()
# Query("""CREATE TABLE alphasync_compact.alphacon (
#   `acc` char(10),
#   `site1` mediumint,
#   `site2` mediumint,
#   `aa1` char(1) DEFAULT NULL,
#   `aa2` char(1) DEFAULT NULL,
#   `atom1` char(3),
#   `atom2` char(3),
#   `type` enum('AromaticContacts','CarbonylContacts','CovalentContacts','HBondContacts','HydrophobicContacts','IonicContacts','MetalContacts','PolarHBondContacts','VanDerWaalsContacts','WeakHBondContacts','WeakPolarHBondContacts'),
#   `afdb` char(1),
#   `dist` float DEFAULT NULL,
#   `pae` float DEFAULT NULL,
#   PRIMARY KEY (`acc`, `site1`, `site2`, `atom1`, `atom2`, `type`, `afdb`)
# ) ENGINE=InnoDB DEFAULT CHARSET=latin1 COMMENT='AlphaSync contacts';""")
# Stoptime();
# 
# Starttime()
# # Query("""INSERT INTO alphasync_compact.alphacon SELECT DISTINCT acc, site1, site2, aa1, aa2, atom1, atom2, type, afdb, dist, pae FROM alphasync.alphacon;""")
# # # DISTINCT shouldn't be necessary anymore, since I removed all duplicates now (test_remove_duplicates_alphacon.py):
# # Query("""INSERT INTO alphasync_compact.alphacon SELECT acc, site1, site2, aa1, aa2, atom1, atom2, type, afdb, dist, pae FROM alphasync.alphacon;""")
# # Using INSERT IGNORE instead (skipping duplicated rows):
# # Query("""INSERT IGNORE INTO alphasync_compact.alphacon SELECT acc, site1, site2, aa1, aa2, atom1, atom2, type, afdb, dist, pae FROM alphasync.alphacon;""")
# # The compact version is now the only one:
# Query("""INSERT INTO alphasync_compact.alphacon SELECT * FROM alphasync.alphacon;""")
# Stoptime();

Starttime()
# Query("""ANALYZE TABLE alphasync_compact.alphacon;""")
Query("""ANALYZE TABLE alphasync.alphacon;""")
Stoptime();

Starttime()
# Query("""OPTIMIZE TABLE alphasync_compact.alphacon;""")
Query("""OPTIMIZE TABLE alphasync.alphacon;""")
Stoptime();
    
Show(lim=0)

print("\nDone!")

