#!/usr/bin/env python3
"""
alphauniprot_species.py:
Fill AlphaSync table 'alphauniprot_species', a summary table of UniProt taxon IDs, latin names and common species names for efficient lookup.

The "completed" column also indicates whether the taxon is fully completed (i.e. all of its non-isoform reviewed or "canonical reference proteome" UniProt sequences are mapped to a structure in AlphaSync table 'alphamap').
The "completed_iso" column indicates whether the taxon is fully completed (including isoforms).
"""

# Initialize
# import pandas as pd
# import numpy as np
# import networkx as nx
# import seaborn as sns
# import matplotlib.pyplot as mp
from blang_mysql import *
from blang import *

alphauniprot = "alphauniprot"
alphamap = "alphamap"
ensembl_species = "ensembl_species"
compara_species = "compara_species"
alphauniprot_species = "alphauniprot_species"

# Args(0, "", "")
# var = Args(1, "[species]", "human")
# infile = f"input/input.txt"

uniprot_version = get_local_uniprot_release()

# Clear table
Clear(alphauniprot_species)

# Start



# Query to get the list of fully finished, albeit ((au.reviewed=1 OR au.refprotcanon=1) AND acc=canon) (canonical reference proteome or reviewed, and no isoforms) taxa:
print(f"Getting list of completed taxa in table '{alphauniprot}' using table '{alphamap}'...")
complete_taxa = FetchSet(Query(f"""SELECT DISTINCT t.mapped_tax FROM (SELECT *, mapped_tax IN (SELECT DISTINCT tax FROM alphafrag WHERE source IN (SELECT DISTINCT source FROM alphafrag HAVING source LIKE 'UP%')) AS is_model, mapped.mapped_seqs + unmapped.unmapped_seqs AS total_seqs, mapped.mapped_seqs / (mapped.mapped_seqs + unmapped.unmapped_seqs) AS mapped_fraction FROM
	(SELECT au.tax AS mapped_tax, au.species AS mapped_species, COUNT(DISTINCT au.seq) AS mapped_seqs FROM {alphauniprot} au, {alphamap} m WHERE au.acc=m.value AND m.type='uniprot' AND m.version='{uniprot_version}' AND m.map IS NOT NULL AND best=1 GROUP BY au.tax) mapped
    LEFT OUTER JOIN
    (SELECT au.tax AS unmapped_tax, au.species AS unmapped_species, COUNT(DISTINCT au.seq) AS unmapped_seqs FROM {alphauniprot} au, {alphamap} m WHERE (au.reviewed=1 OR au.refprotcanon=1) AND acc=canon AND au.acc=m.value AND m.type='uniprot' AND m.version='{uniprot_version}' AND m.map IS NULL GROUP BY au.tax) unmapped
    ON mapped_tax=unmapped_tax
LEFT OUTER JOIN alphauniprot_species s ON s.tax=mapped_tax GROUP BY mapped_tax HAVING unmapped_seqs IS NULL AND is_model=1 ORDER BY s.id, mapped_fraction DESC) t ORDER BY t.mapped_seqs DESC;"""))
# Concatenate for query
complete_taxa_string = "(" + ", ".join(str(tax) for tax in complete_taxa) + ")"

print(f"Completed taxa: {len(complete_taxa)}")
if len(complete_taxa) == 0:
    Die(f"Error: No completed taxa found in table '{alphamap}'")



# Query to get the list of fully finished, albeit (au.reviewed=1 OR au.refprotcanon=1) (canonical reference proteome or reviewed, including isoforms) taxa:
print(f"Getting list of completed taxa (including isoforms) in table '{alphauniprot}' using table '{alphamap}'...")
complete_iso_taxa = FetchSet(Query(f"""SELECT DISTINCT t.mapped_tax FROM (SELECT *, mapped_tax IN (SELECT DISTINCT tax FROM alphafrag WHERE source IN (SELECT DISTINCT source FROM alphafrag HAVING source LIKE 'UP%')) AS is_model, mapped.mapped_seqs + unmapped.unmapped_seqs AS total_seqs, mapped.mapped_seqs / (mapped.mapped_seqs + unmapped.unmapped_seqs) AS mapped_fraction FROM
	(SELECT au.tax AS mapped_tax, au.species AS mapped_species, COUNT(DISTINCT au.seq) AS mapped_seqs FROM {alphauniprot} au, {alphamap} m WHERE au.acc=m.value AND m.type='uniprot' AND m.version='{uniprot_version}' AND m.map IS NOT NULL AND best=1 GROUP BY au.tax) mapped
    LEFT OUTER JOIN
    (SELECT au.tax AS unmapped_tax, au.species AS unmapped_species, COUNT(DISTINCT au.seq) AS unmapped_seqs FROM {alphauniprot} au, {alphamap} m WHERE (au.reviewed=1 OR au.refprotcanon=1) AND au.acc=m.value AND m.type='uniprot' AND m.version='{uniprot_version}' AND m.map IS NULL GROUP BY au.tax) unmapped
    ON mapped_tax=unmapped_tax
LEFT OUTER JOIN alphauniprot_species s ON s.tax=mapped_tax GROUP BY mapped_tax HAVING unmapped_seqs IS NULL AND is_model=1 ORDER BY s.id, mapped_fraction DESC) t ORDER BY t.mapped_seqs DESC;"""))
# Concatenate for query
complete_iso_taxa_string = "(" + ", ".join(str(tax) for tax in complete_iso_taxa) + ")"

print(f"Completed taxa (including isoforms): {len(complete_iso_taxa)}")
if len(complete_iso_taxa) == 0:
    Die(f"Error: No completed taxa (including isoforms) found in table '{alphamap}'")



# Insert species
print(f"\nFilling table '{alphauniprot_species}' with UniProt taxon IDs, latin names and common species names from '{alphauniprot}'...")
# Order: human first, then completed taxa, then Ensembl Compara species, then sorted by UniProt species mnemonic
# query = Query(f"INSERT INTO {alphauniprot_species} SELECT DISTINCT NULL, au.tax, au.species, au.species_latin, REGEXP_REPLACE(REGEXP_REPLACE(au.species_latin, '\\\\s*\\\\([^)]+\\\\)$', ''), ' (subsp\\\\.|serotype|CBS) .+$', '') AS species_latin_short, au.species_common, e.fullspecies, c.id IS NOT NULL, au.tax IN {complete_taxa_string} FROM {alphauniprot} au LEFT OUTER JOIN {ensembl_species} e ON e.tax=au.tax LEFT OUTER JOIN {compara_species} c ON c.tax=au.tax GROUP BY au.tax ORDER BY au.tax=9606 DESC, au.tax IN {complete_taxa_string} DESC, c.id IS NOT NULL DESC, c.id, au.species")
# query = Query(f"INSERT INTO {alphauniprot_species} SELECT DISTINCT NULL, au.tax, au.species, au.species_latin, REGEXP_REPLACE(species_latin, '^(\\\\w+) (\\\\w+) .+', '$1 $2') AS species_latin_short, au.species_common, e.fullspecies, c.id IS NOT NULL, au.tax IN {complete_taxa_string} FROM {alphauniprot} au LEFT OUTER JOIN {ensembl_species} e ON e.tax=au.tax LEFT OUTER JOIN {compara_species} c ON c.tax=au.tax GROUP BY au.tax ORDER BY au.tax=9606 DESC, c.id IS NOT NULL DESC, c.id, au.tax IN {complete_taxa_string} DESC, species_latin_short")
query = Query(f"INSERT INTO {alphauniprot_species} SELECT DISTINCT NULL, au.tax, au.species, au.species_latin, REGEXP_REPLACE(species_latin, '^(\\\\w+) (\\\\w+) .+', '$1 $2') AS species_latin_short, au.species_common, e.fullspecies, c.id IS NOT NULL, au.tax IN {complete_taxa_string}, au.tax IN {complete_iso_taxa_string} FROM {alphauniprot} au LEFT OUTER JOIN {ensembl_species} e ON e.tax=au.tax LEFT OUTER JOIN {compara_species} c ON c.tax=au.tax GROUP BY au.tax ORDER BY au.tax=9606 DESC, c.id IS NOT NULL DESC, c.id, species_latin")
State(f"{Numrows(query):,} rows inserted")

Show(lim=0)

# Optimize(alphauniprot_species)

print("\nDone!")
