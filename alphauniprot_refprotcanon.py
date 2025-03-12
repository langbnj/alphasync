#!/usr/bin/env python3
"""
alphauniprot.py:
Get canonical reference proteome status for via UniProt's FTP server and update AlphaSync table 'alphauniprot' (column 'refprotcanon').

'refprotcanon' complements the 'refproteome' column: 'refprotcanon' are the "canonical" reference proteome sequences.
Mouse, for example, has ~55,000 reference proteome sequences, but only ~22,000 canonical reference proteome sequences (one per gene).
Details: https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/README

"These files, composed of canonical and additional sequences, are non-redundant
FASTA sets for the sequences of each reference proteome.
The additional set contains isoform/variant sequences for a given gene, and its
FASTA header indicates the corresponding canonical sequence ("Isoform of ...").
The FASTA format is the standard UniProtKB format."
"""

# Initialize
import pandas as pd
# import numpy as np
# import networkx as nx
# import seaborn as sns
# import matplotlib.pyplot as mp
# import tarfile
# import gzip
# import io
import os
import sys
# from ftplib import FTP
# import re
import requests
from requests.adapters import HTTPAdapter, Retry
import gzip
import io
# import tempfile
# from Bio import SeqIO
from blang_mysql import *
from blang import *

# Variables

alphauniprot = "alphauniprot"       # Table name
refprotcanon = "refprotcanon"       # Column name
alphafrag = "alphafrag"

# Start
Args(0, "[-debug]", "-debug: Don't actually make any changes, just simulate", "")

# Clear table
if not Switch('debug'):
    print(f"Clearing column '{refprotcanon}' in table '{alphauniprot}'...")
    # Query(f"UPDATE {alphauniprot} SET {refprotcanon}=NULL")
    query = Query(f"UPDATE {alphauniprot} SET {refprotcanon}=0")
    print(f"{Numrows(query):,} rows affected")

# Start

# Download and parse UniProt reference proteome FASTA files
State(f"Downloading and parsing UniProt reference proteome FASTA files and updating '{refprotcanon}' column in table '{alphauniprot}':")
# https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000000589/UP000000589_10090.fasta.gz
# proteomes = FetchList(Query("SELECT DISTINCT source FROM alphafrag HAVING source LIKE 'UP%'"))
# Stream FASTA file from FTP server over https
retries = Retry(total=5, backoff_factor=1, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount('https://', HTTPAdapter(max_retries=retries))
if not Switch('debug'):
    # Get the 48 "Model organisms" and "Global health proteomes" from the AlphaFold Protein Structure Database proteome identifiers (which start with UP...) only
    # This is actually pretty fast (~30 minutes). Might want to add the 200 Ensembl Compara species later (or all other AlphaSync species)
    query = Query(f"SELECT DISTINCT source FROM {alphafrag} HAVING source LIKE 'UP%'")
else:
    # For debugging: Get human only
    query = Query(f"SELECT DISTINCT source FROM {alphafrag} HAVING source='UP000005640_9606_HUMAN_v4'")
# query = Query(f"SELECT source, tax FROM {alphafrag} GROUP BY source HAVING source LIKE 'UP%'")
# query = Query(f"SELECT DISTINCT source, tax FROM {alphafrag} WHERE source IN (SELECT DISTINCT source FROM {alphafrag} HAVING source LIKE 'UP%')")
# Source: e.g. UP000005640_9606_HUMAN_v4, where proteome=UP000005640, tax=9606
i = 0
for source, in query:
    i += 1
    print(f" >> {i} / {Numrows(query)} >> {source}")

    Log("total sources", source)

    # Split
    proteome = source.split('_')[0]
    tax = source.split('_')[1]

    # Canonical reference proteome

    # Check if the file exists in /Eukaryota/, /Archaea/, /Bacteria/ or /Viruses/
    categories = ['Eukaryota', 'Archaea', 'Bacteria', 'Viruses']
    url = None
    for category in categories:
        # e.g. https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz
        test_url = f"https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/{category}/{proteome}/{proteome}_{tax}.fasta.gz"
        response = session.head(test_url)
        if response.status_code == 200:
            url = test_url
            break

    if url is None:
        # # Expected crashes for UP000020681_1299332_MYCUL_v4 and UP000325664_1352_ENTFC_v4
        # if source != 'UP000020681_1299332_MYCUL_v4' and source != 'UP000325664_1352_ENTFC_v4':
        #     # Otherwise, crash
        #     Die(f"Error: Canonical file not found for {source} in any category (Eukaryota, Archaea, Bacteria, Viruses)")
        print(f"SKIP (canonical file not found for {source} in any category (Eukaryota, Archaea, Bacteria, Viruses))")
        Log(f"canonical file not found in any subfolder (Eukaryota etc.) for source (skipped)", source)
        continue

    print(f"   >> canonical >> {url}")
    # print(f"   >> canonical >>  ", end="")
    response = session.get(url)
    if response.status_code == 200:
        # Decompress and parse FASTA file
        with gzip.open(io.BytesIO(response.content), 'rt') as f:
            lines = f.read().strip().split('\n')
            for line in tq(lines):
                if line.startswith('>'):
                    acc = line.split('|')[1]
                    if not Switch('debug'):
                        Query(f"UPDATE {alphauniprot} SET {refprotcanon}=1 WHERE tax='{tax}' AND acc='{acc}'")
                    else:
                        State(f"UPDATE {alphauniprot} SET {refprotcanon}=1 WHERE tax='{tax}' AND acc='{acc}'")
                    Log(f"successfully updated row for acc", acc)
                    Log(f"successfully updated row for canonical source", source)
    else:
        print(f" >> Error downloading {url}: {response.status_code}")
        Log(f"error: could not download url for canonical source (skipped)", source)
        Log(f"error: could not download url for canonical url (skipped)", url)

# This would set refprotcanon to 0 for isoforms and variants, but it only makes sense if we used NULL as the default (rather than 0).
# Using only 0 and 1 now for clarity: either an accession is in the canonical reference proteome for a given species (one per gene) or it is not.
#     # Additional sequences (isoforms/variants)
# 
#     # Check if the file exists in /Eukaryota/, /Archaea/, /Bacteria/ or /Viruses/
#     categories = ['Eukaryota', 'Archaea', 'Bacteria', 'Viruses']
#     url = None
#     for category in categories:
#         test_url = f"https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/{category}/{proteome}/{proteome}_{tax}_additional.fasta.gz"
#         response = session.head(test_url)
#         if response.status_code == 200:
#             url = test_url
#             break
# 
#     if url is None:
#         # print(f"SKIP (additional file not found for {source} in any category (Eukaryota, Archaea, Bacteria, Viruses))")
#         print(f"SKIP (additional file not found)")
#         Log(f"additional file not found in any subfolder (Eukaryota etc.) for source (skipped)", source)
#         continue
# 
#     # print(f"   >> additional >> {url}")
#     print(f"   >> additional >> ", end="")
#     response = session.get(url)
#     if response.status_code == 200:
#         # Decompress and parse FASTA file
#         with gzip.open(io.BytesIO(response.content), 'rt') as f:
#             lines = f.read().strip().split('\n')
#             for line in tq(lines):
#                 if line.startswith('>'):
#                     acc = line.split('|')[1]
#                     if not Switch('debug'):
#                         Query(f"UPDATE {alphauniprot} SET {refprotcanon}=0 WHERE tax='{tax}' AND acc='{acc}'")
#                     else:
#                         State(f"UPDATE {alphauniprot} SET {refprotcanon}=0 WHERE tax='{tax}' AND acc='{acc}'")
#                     Log(f"successfully updated row for acc", acc)
#                     Log(f"successfully updated row for additional source", source)
#     else:
#         print(f"Error downloading {url}: {response.status_code}")
#         Log(f"error: could not download url for additional source (skipped)", source)
#         Log(f"error: could not download url for additional url (skipped)", url)

Show(lim=20)

# if not Switch('debug'):
#     Optimize(alphauniprot)

print("\nDone!")
