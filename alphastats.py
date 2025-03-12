#!/usr/bin/env python3
"""
alphastats.py:
Fill AlphaSync table 'alphastats', a table containing precalculated statistics that are displayed on the AlphaSync website.
"""

# Initialize
# import pandas as pd
# import numpy as np
# import networkx as nx
# import seaborn as sns
# import matplotlib.pyplot as mp
from blang_mysql import *
from blang import *

alphastats = "alphastats"

# Args(0, "", "")
# var = Args(1, "[species]", "human")
# infile = f"input/input.txt"

# Clear table
Clear(alphastats)

# Start

uniprot_version = get_local_uniprot_release()
ensembl_version = "108"

# current_accs
Starttime()
print(f"\nGetting number of UniProt accessions successfully mapped to structures in table 'alphasync_compact.alphamap'...")
current_accs = FetchOne(Query(f"SELECT COUNT(DISTINCT value) FROM alphasync_compact.alphamap WHERE type='uniprot' AND version='{uniprot_version}' AND map IS NOT NULL AND best=1"))
print(f" >> current_accs {current_accs:,}")
Query(f"INSERT INTO {alphastats} SET stat='current_accs', value='{current_accs}'")
Stoptime()

# current_accs_tax_[tax]
Starttime()
print(f"\nGetting number of UniProt accessions successfully mapped to structures in table 'alphasync_compact.alphamap' for individual completed taxa...")
for tax, accs in Query(f"SELECT DISTINCT tax, COUNT(DISTINCT value) FROM alphasync_compact.alphamap WHERE type='uniprot' AND version='{uniprot_version}' AND map IS NOT NULL AND best=1 AND tax IN (SELECT DISTINCT tax FROM alphasync_compact.alphauniprot_species WHERE complete=1) GROUP BY tax"):
    print(f" >> current_accs_tax_{tax} >> {accs:,}")
    Query(f"INSERT INTO {alphastats} SET stat='current_accs_tax_{tax}', value='{accs}'")
# for tax, in Fetch(Query(f"SELECT DISTINCT tax FROM alphamap WHERE type='uniprot' AND version='{uniprot_version}' AND map IS NOT NULL AND best=1")):
#     acc_count = FetchOne(Query(f"SELECT COUNT(DISTINCT value) FROM alphamap WHERE type='uniprot' AND version='{uniprot_version}' AND map IS NOT NULL AND best=1 AND tax='{tax}'"))
#     print(f" >> current_accs_{tax} {acc_count}")
#     Query(f"INSERT INTO {alphastats} SET stat='current_accs_{tax}', value='{acc_count}'")
Stoptime()

# current_residues
Starttime()
print(f"\nGetting number of residues in table 'alphasync_compact.alphasa'...")
# # Faster but slightly inaccurate due to afdb=1 and afdb=0 being present in some cases, but only afdb=0 being used
# current_residues = FetchOne(Query(f"SELECT COUNT(*) FROM alphasync_compact.alphasa"))
current_residues = FetchOne(Query(f"SELECT COUNT(*) FROM alphasync_compact.alphamap m, alphasync_compact.alphasa a WHERE m.type='uniprot' AND m.version='{uniprot_version}' AND m.map=a.acc AND m.afdb=a.afdb"))
print(f" >> current_residues {current_residues:,}")
Query(f"INSERT INTO {alphastats} SET stat='current_residues', value='{current_residues}'")
Stoptime()

# current_contacts
Starttime()
print(f"\nGetting number of contacts in table 'alphasync.alphacon'...")
# # Faster but slightly inaccurate due to afdb=1 and afdb=0 being present in some cases, but only afdb=0 being used
# current_contacts = FetchOne(Query(f"SELECT COUNT(*) FROM alphasync_compact.alphacon"))
current_contacts = FetchOne(Query(f"SELECT COUNT(*) FROM alphasync_compact.alphamap m, alphasync.alphacon c WHERE m.type='uniprot' AND m.version='{uniprot_version}' AND m.map=c.acc AND m.afdb=c.afdb"))
print(f" >> current_contacts {current_contacts:,}")
Query(f"INSERT INTO {alphastats} SET stat='current_contacts', value='{current_contacts}'")
Stoptime()

# current_data
Starttime()
print(f"\nGetting total size of data in schema 'alphasync_compact' (also including 'alphasync.alphacon', which is only present as a view in 'alphasync_compact')...")
total_gb = FetchOne(Query(f"SELECT SUM(t.total_gb) AS `total_GB` FROM (SELECT table_schema AS `db`, table_name AS `table`, table_rows/1000000000 AS `billion_rows`, table_rows AS `rows`, round(((data_length) / 1024 / 1024), 2) AS `data_mb`, round(((index_length) / 1024 / 1024), 2) AS `index_mb`, round(((data_length + index_length) / 1024 / 1024), 2) AS `total_mb`, round(((data_length) / 1024 / 1024 / 1024), 2) AS `data_gb`, round(((index_length) / 1024 / 1024 / 1024), 2) AS `index_gb`, round(((data_length + index_length) / 1024 / 1024 / 1024), 2) AS `total_gb` FROM information_schema.TABLES WHERE (table_schema='alphasync_compact' AND table_name LIKE 'alpha%') OR (table_schema='blang' AND table_name='alphacon')) AS t"))
current_data = str(round(total_gb)) + " GB"
print(f" >> current_data {current_data}")
Query(f"INSERT INTO {alphastats} SET stat='current_data', value='{current_data}'")
Stoptime()

# ensembl_version
Starttime()
print(f"\nCurrent Ensembl version (Note: hardcoded here in alphastats.py):")
print(f" >> ensembl_version {ensembl_version}")
Query(f"INSERT INTO {alphastats} SET stat='ensembl_version', value='{ensembl_version}'")
Stoptime()

# current_predictions
Starttime()
print(f"\nGetting number of AlphaSync re-predicted structures (accessions) in table 'alphaseq'...")
# current_predictions = FetchOne(Query(f"SELECT COUNT(DISTINCT acc) FROM alphaseq WHERE afdb=0"))
current_predictions = FetchOne(Query(f"SELECT COUNT(DISTINCT map) FROM alphasync_compact.alphamap WHERE afdb=0"))
print(f" >> current_predictions {current_predictions:,}")
Query(f"INSERT INTO {alphastats} SET stat='current_predictions', value='{current_predictions}'")
Stoptime()

# current_predictions_frags
Starttime()
print(f"\nGetting number of AlphaSync re-predicted structures (accessions & fragments) in table 'alphaseq'...")
current_predictions_frags = FetchOne(Query(f"SELECT COUNT(DISTINCT acc, frag) FROM alphafrag WHERE afdb=0"))
print(f" >> current_predictions_frags {current_predictions_frags:,}")
Query(f"INSERT INTO {alphastats} SET stat='current_predictions_frags', value='{current_predictions_frags}'")
Stoptime()

# current_predictions_nofrag
Starttime()
print(f"\nGetting number of AlphaSync re-predicted structures (accessions & fragments) in table 'alphaseq'...")
current_predictions_nofrag = FetchOne(Query(f"SELECT COUNT(DISTINCT acc) FROM alphaseq WHERE afdb=0 AND frags=1"))
print(f" >> current_predictions_nofrag {current_predictions_nofrag:,}")
Query(f"INSERT INTO {alphastats} SET stat='current_predictions_nofrag', value='{current_predictions_nofrag}'")
Stoptime()

# current_isoforms
Starttime()
print(f"\nGetting number of UniProt non-canonical isoforms successfully mapped to structures in table 'alphasync_compact.alphamap'...")
current_isoforms = FetchOne(Query(f"SELECT COUNT(DISTINCT value) FROM alphasync_compact.alphamap WHERE type='uniprot' AND version='{uniprot_version}' AND map IS NOT NULL AND value REGEXP '-[0-9]+$'"))
print(f" >> current_isoforms {current_isoforms:,}")
Query(f"INSERT INTO {alphastats} SET stat='current_isoforms', value='{current_isoforms}'")
Stoptime()

# current_isoform_predictions
Starttime()
print(f"\nGetting number of AlphaSync re-predicted structures (accessions) for UniProt non-canonical isoforms in table 'alphaseq'...")
current_isoform_predictions = FetchOne(Query(f"SELECT COUNT(DISTINCT acc) FROM alphaseq WHERE afdb=0 AND acc REGEXP '-[0-9]+$'"))
print(f" >> current_isoform_predictions {current_isoform_predictions:,}")
Query(f"INSERT INTO {alphastats} SET stat='current_isoform_predictions', value='{current_isoform_predictions}'")
Stoptime()

# current_isoform_predictions_frags
Starttime()
print(f"\nGetting number of AlphaSync re-predicted structures (accessions) for UniProt non-canonical isoforms in table 'alphaseq'...")
current_isoform_predictions_frags = FetchOne(Query(f"SELECT COUNT(DISTINCT acc, frag) FROM alphafrag WHERE afdb=0 AND acc REGEXP '-[0-9]+$'"))
print(f" >> current_isoform_predictions_frags {current_isoform_predictions_frags:,}")
Query(f"INSERT INTO {alphastats} SET stat='current_isoform_predictions_frags', value='{current_isoform_predictions_frags}'")
Stoptime()

# current_species
Starttime()
print(f"\nGetting number of taxa (species) successfully mapped to structures in table 'alphasync_compact.alphamap'...")
current_species = FetchOne(Query(f"SELECT COUNT(DISTINCT tax) FROM alphasync_compact.alphamap WHERE type='uniprot' AND version='{uniprot_version}' AND map IS NOT NULL"))
print(f" >> current_species {current_species:,}")
Query(f"INSERT INTO {alphastats} SET stat='current_species', value='{current_species}'")
Stoptime()

# current_species_completed
Starttime()
# # Query to get the list of 39 fully finished, albeit ((au.reviewed=1 OR au.refprotcanon=1) AND acc=canon), taxa:
# print("\nGetting list of completed taxa using table 'alphamap'...")
# complete_taxa = FetchSet(Query(f"""SELECT DISTINCT t.mapped_tax FROM (SELECT *, mapped_tax IN (SELECT DISTINCT tax FROM alphafrag WHERE source IN (SELECT DISTINCT source FROM alphafrag HAVING source LIKE 'UP%')) AS is_model, mapped.mapped_seqs + unmapped.unmapped_seqs AS total_seqs, mapped.mapped_seqs / (mapped.mapped_seqs + unmapped.unmapped_seqs) AS mapped_fraction FROM
# 	(SELECT au.tax AS mapped_tax, au.species AS mapped_species, COUNT(DISTINCT au.seq) AS mapped_seqs FROM alphauniprot au, alphamap m WHERE au.acc=m.value AND m.type='uniprot' AND m.version='{uniprot_version}' AND m.map IS NOT NULL AND best=1 GROUP BY au.tax) mapped
#     LEFT OUTER JOIN
#     (SELECT au.tax AS unmapped_tax, au.species AS unmapped_species, COUNT(DISTINCT au.seq) AS unmapped_seqs FROM alphauniprot au, alphamap m WHERE (au.reviewed=1 OR au.refprotcanon=1) AND acc=canon AND au.acc=m.value AND m.type='uniprot' AND m.version='{uniprot_version}' AND m.map IS NULL AND au.seq NOT REGEXP '[BZX]' AND LENGTH(au.seq)>=16 GROUP BY au.tax) unmapped
#     ON mapped_tax=unmapped_tax
# LEFT OUTER JOIN compara_species s ON s.tax=mapped_tax GROUP BY mapped_tax HAVING unmapped_seqs IS NULL AND is_model=1 ORDER BY s.id, mapped_fraction DESC) t ORDER BY t.mapped_seqs DESC;"""))
# complete_taxa = FetchSet(Query(f"""SELECT DISTINCT t.mapped_tax FROM (SELECT *, mapped_tax IN (SELECT DISTINCT tax FROM alphafrag WHERE source IN (SELECT DISTINCT source FROM alphafrag HAVING source LIKE 'UP%')) AS is_model, mapped.mapped_seqs + unmapped.unmapped_seqs AS total_seqs, mapped.mapped_seqs / (mapped.mapped_seqs + unmapped.unmapped_seqs) AS mapped_fraction FROM
# 	(SELECT au.tax AS mapped_tax, au.species AS mapped_species, COUNT(DISTINCT au.seq) AS mapped_seqs FROM alphauniprot au, alphamap m WHERE au.acc=m.value AND m.type='uniprot' AND m.version='{uniprot_version}' AND m.map IS NOT NULL AND best=1 GROUP BY au.tax) mapped
#     LEFT OUTER JOIN
#     (SELECT au.tax AS unmapped_tax, au.species AS unmapped_species, COUNT(DISTINCT au.seq) AS unmapped_seqs FROM alphauniprot au, alphamap m WHERE (au.reviewed=1 OR au.refprotcanon=1) AND acc=canon AND au.acc=m.value AND m.type='uniprot' AND m.version='{uniprot_version}' AND m.map IS NULL GROUP BY au.tax) unmapped
#     ON mapped_tax=unmapped_tax
# LEFT OUTER JOIN alphauniprot_species s ON s.tax=mapped_tax GROUP BY mapped_tax HAVING unmapped_seqs IS NULL AND is_model=1 ORDER BY s.id, mapped_fraction DESC) t ORDER BY t.mapped_seqs DESC;"""))
print("\nGetting list of completed taxa from table 'alphauniprot_species'...")
complete_taxa = FetchSet(Query(f"""SELECT DISTINCT tax FROM alphauniprot_species WHERE complete=1"""))
print(f" >> current_species_completed {len(complete_taxa):,}")
current_species_completed = len(complete_taxa)
Query(f"INSERT INTO {alphastats} SET stat='current_species_completed', value='{current_species_completed}'")
Stoptime()

# current_species_completed_iso
Starttime()
print("\nGetting list of completed taxa (including isoforms) from table 'alphauniprot_species'...")
complete_iso_taxa = FetchSet(Query(f"""SELECT DISTINCT tax FROM alphauniprot_species WHERE complete_iso=1"""))
print(f" >> current_species_completed_iso {len(complete_iso_taxa):,}")
current_species_completed_iso = len(complete_iso_taxa)
Query(f"INSERT INTO {alphastats} SET stat='current_species_completed_iso', value='{current_species_completed_iso}'")
Stoptime()

# current_tablecount
Starttime()
print(f"\nGetting number of tables in schema 'alphasync_compact'...")
current_tablecount = FetchOne(Query(f"SELECT COUNT(*) FROM information_schema.TABLES WHERE TABLE_SCHEMA='alphasync_compact' AND TABLE_NAME LIKE 'alpha%'"))
print(f" >> current_tablecount {current_tablecount:,}")
Query(f"INSERT INTO {alphastats} SET stat='current_tablecount', value='{current_tablecount}'")
Stoptime()

# uniprot_version
Starttime()
print(f"\nGetting local UniProt version...")
# uniprot_version = get_local_uniprot_release()
print(f" >> uniprot_version {uniprot_version}")
Query(f"INSERT INTO {alphastats} SET stat='uniprot_version', value='{uniprot_version}'")
Stoptime()

Show(lim=0)

# Optimize(alphastats)

print("\nDone!")
