#!/usr/bin/env python3
"""
migrate_alphasync_tar_archives.py:
Compress all CIF/PAE/params output files into .tar archives
"""

# Initialize
import pandas as pd
# import numpy as np
# import networkx as nx
# import seaborn as sns
# import matplotlib.pyplot as mp
import tarfile
import gzip
import io
import requests
import json
# from Bio import SeqIO
from blang_mysql import *
from blang import *

# Paths for deleting outdated AlphaSync structure predictions (where an accession has since been made obsolete, or a sequence was changed, in UniProt)
cifdir = "input/alphasync/cif"
paedir = "input/alphasync/pae"
paramdir = "input/alphasync/params"

(version) = Args(1, "[Current UniProt version]\n -debug: Don't actually make any changes, just simulate", "2025_01")

# Output files
alphasync_cif_tar = f"alphasync_cif_{version}.tar"
alphasync_pae_tar = f"alphasync_pae_{version}.tar"
alphasync_params_tar = f"alphasync_params_{version}.tar"

# Compress all CIF/PAE/params output files into .tar archives
if not Switch('debug') and not Switch('debug2'):
    if glob(f"{cifdir}/AF-*.cif.gz"):
        # Run("Compress all CIF files into a single .tar file", f"cd {cifdir}; tar -cf {alphasync_tar} ../AF-*.cif.gz")
        # Run("Compress all CIF files into a single .tar file", f"cd {cifdir}; find . -name 'AF-*.cif.gz' -print0 | tar --null -T - -cf ../{alphasync_tar}")

        # Run sequentially
        Run("Compress all CIF files into a single .tar file (naturally sorted)",            f"cd {cifdir};   find . -name 'AF-*.cif.gz' -printf '%P\\n'  | sort -V | tar -T - -cf ../{alphasync_cif_tar}")
        Run("Compress all PAE JSON files into a single .tar file (naturally sorted)",       f"cd {paedir};   find . -name 'AF-*.json.gz' -printf '%P\\n' | sort -V | tar -T - -cf ../{alphasync_pae_tar}")
        Run("Compress all parameter JSON files into a single .tar file (naturally sorted)", f"cd {paramdir}; find . -name 'AF-*.json' -printf '%P\\n'    | sort -V | tar -T - -cf ../{alphasync_params_tar}")

        # # Run in parallel
        # Run("Compress all CIF files into a single .tar file (naturally sorted)",            f"""~/scripts/qsub.sh sh -c 'cd {cifdir};   find . -name "AF-*.cif.gz"  -printf '%P\\n'  | sort -V | tar -T - -cf ../{alphasync_cif_tar}'""")
        # Run("Compress all PAE JSON files into a single .tar file (naturally sorted)",       f"""~/scripts/qsub.sh sh -c 'cd {paedir};   find . -name "AF-*.json.gz" -printf '%P\\n' | sort -V | tar -T - -cf ../{alphasync_pae_tar}'""")
        # Run("Compress all parameter JSON files into a single .tar file (naturally sorted)", f"""~/scripts/qsub.sh sh -c 'cd {paramdir}; find . -name "AF-*.json"    -printf '%P\\n'    | sort -V | tar -T - -cf ../{alphasync_params_tar}'""")
        # Waitforjobs()

        # Run manually, sequentially (for 2025_01):
        # cd input/alphasync/cif;    find . -name 'AF-*.cif.gz' -printf '%P\n' | sort -V | tar -T - -cf ../alphasync_cif_2025_01.tar  ; echo " >> alphasync_cif_2025_01.tar complete";    cd input/alphasync; 
        # cd input/alphasync/pae;    find . -name 'AF-*.json.gz' -printf '%P\n' | sort -V | tar -T - -cf ../alphasync_pae_2025_01.tar ; echo " >> alphasync_pae_2025_01.tar complete";    cd input/alphasync; 
        # cd input/alphasync/params; find . -name 'AF-*.json' -printf '%P\n' | sort -V | tar -T - -cf ../alphasync_params_2025_01.tar ; echo " >> alphasync_params_2025_01.tar complete"; cd input/alphasync; 
        # In one line (for 2025_01):
        # cd input/alphasync/cif;    find . -name 'AF-*.cif.gz' -printf '%P\n' | sort -V | tar -T - -cf ../alphasync_cif_2025_01.tar  ; echo " >> alphasync_cif_2025_01.tar complete";    cd input/alphasync; cd input/alphasync/pae;    find . -name 'AF-*.json.gz' -printf '%P\n' | sort -V | tar -T - -cf ../alphasync_pae_2025_01.tar ; echo " >> alphasync_pae_2025_01.tar complete";    cd input/alphasync; cd input/alphasync/params; find . -name 'AF-*.json' -printf '%P\n' | sort -V | tar -T - -cf ../alphasync_params_2025_01.tar ; echo " >> alphasync_params_2025_01.tar complete"; cd input/alphasync; 
        # Run manually, in parallel (for 2025_01):
        # ~/scripts/qsub.sh sh -c 'cd input/alphasync/cif;    find . -name "AF-*.cif.gz"  -printf '%P\n' | sort -V | tar -T - -cf ../alphasync_cif_2025_01.tar'
        # ~/scripts/qsub.sh sh -c 'cd input/alphasync/pae;    find . -name "AF-*.json.gz" -printf '%P\n' | sort -V | tar -T - -cf ../alphasync_pae_2025_01.tar'
        # ~/scripts/qsub.sh sh -c 'cd input/alphasync/params; find . -name "AF-*.json"    -printf '%P\n' | sort -V | tar -T - -cf ../alphasync_params_2025_01.tar'
        # In one line (for 2025_01):
        # ~/scripts/qsub.sh sh -c 'cd input/alphasync/cif;    find . -name "AF-*.cif.gz"  -printf '%P\n' | sort -V | tar -T - -cf ../alphasync_cif_2025_01.tar'; ~/scripts/qsub.sh sh -c 'cd input/alphasync/pae;    find . -name "AF-*.json.gz" -printf '%P\n' | sort -V | tar -T - -cf ../alphasync_pae_2025_01.tar'; ~/scripts/qsub.sh sh -c 'cd input/alphasync/params; find . -name "AF-*.json"    -printf '%P\n' | sort -V | tar -T - -cf ../alphasync_params_2025_01.tar'


# # Compress all CIF files into per-species .tar files (as for all other AlphaFold Protein Structure Database files)
# I think it'd be too many...perhaps for human only?
# TODO: Add a separate file for human?
# if not Switch('debug') and not Switch('debug2'):
#     if glob(f"{cifdir}/AF-*.cif.gz"):
#         q1 = "SELECT DISTINCT m.tax"
#         q = q1 + q2
#         for tax in Fetch(Query(q)):
#             # Make directory for this taxon
#             taxdir = f"{cifdir}/{tax}"
#         
#             # Run("Compress all CIF files into a single .tar file", f"cd {cifdir}; tar -cf {alphasync_tar} ../AF-*.cif.gz")
#             # Run("Compress all CIF files into a single .tar file", f"cd {cifdir}; find . -name 'AF-*.cif.gz' -print0 | tar --null -T - -cf ../{alphasync_tar}")
#             Run("Compress all CIF files into a single .tar file (naturally sorted)", f"cd {cifdir}; find . -name 'AF-*.cif.gz' -printf '%P\\n' | sort -V | tar -T - -cf ../{alphasync_tar}")

Show(lim=0)

print("\nDone!")
