#!/usr/bin/env python3
"""
Run: Run entire AlphaSync pipeline
"""

# Initialize

import os
import sys
import requests

from blang_mysql import *
from blang import *

Args(0, "[-alphasync] [-refresh]", " -alphasync: Synchronize with latest UniProt only. Otherwise, load and process AlphaFold Protein Structure Database (AFDB) structures.\n -refresh: Run even if the local UniProt release is still current")

# Functions

# Check if remote UniProt has been updated, and exit if not (unless -refresh is active)
check_uniprot_release()
local_uniprot_release = get_local_uniprot_release()



# Start

if not Switch('alphasync'):

    # Download AlphaFold DB data (https://alphafold.ebi.ac.uk and https://console.cloud.google.com/marketplace/product/bigquery-public-data/deepmind-alphafold)

    Run("Download AlphaFold structure predictions from EBI FTP server", "download_ftp.py")
    # Run("Download AlphaFold structure predictions for all species from Google Cloud Services (~17 TB)", "download_gcs.py")
    Run("Download AlphaFold structure predictions for the ~200 Ensembl Compara species from Google Cloud Storage", "download_gcs.py -comparaonly")

    Clear("alphafrag")
    Clear("alphauniprot")
    Clear("alphaseq")
    Clear("alphamap")
    Clear("alphasa")
    Clear("alphacon")


    # alphafrag.py: Updates alphafrag (fragment-level sources and sequences)
    Run("Parse fragment sequences into table 'alphafrag'", "alphafrag.py")

    # alphaseq.py: Updates alphaseq (protein-level sequences)
    Run("Combine fragment sequences from table 'alphafrag' into table 'alphaseq'", "alphaseq.py")

    # main.py: Updates alphasa and alphacon (alphasa: accessible surface areas from DSSP) (alphacon: residue-residue contacts from Lahuta)
    Run("Main: Submit jobs that run DSSP and Lahuta on fragments and combine output across fragments, insert ASA values into table 'alphasa', and contacts into table 'alphacon'", "main.py")
    # # Use -checkseqs switch when updating existing alphaseq and alphasa tables (slow) (parses sequences)
    # Run("Main: Submit jobs that run DSSP and Lahuta on fragments and combine output across fragments, insert ASA values into table 'alphasa', and contacts into table 'alphacon'", "main.py -checkseqs")

    # Update alphamap (mapping from latest UniProt via API to AlphaFold accessions)
    Run("Update alphamap for UniProt (latest, via UniProt API)", f"alphamap_uniprot.py {local_uniprot_release}")

    # # Update alphamap (mapping from UniProt and Ensembl to AlphaFold accessions)
    # Run("Update alphamap for UniProt (local)", "alphamap_uniprot_local.py 2022_04")
    # Run("Update alphamap for Ensembl (local)", "alphamap_ensembl_local.py 108 -comparaonly")      # Note: Requires slightly over 16 GB of memory (32 GB works fine)

    # # Wait for jobs to finish
    # Waitforjobs()



    # Table optimizations (will be very slow, feel free to comment these out)

    # # Fill the "best" column in table "alphaseq", flagging the best median pLDDT structure for each unique sequence (since AlphaFold DB frequently contains multiple predictions for the same exact sequence)
    # Run("Add 'best' column to table 'alphaseq'", "alphaseq_add_best.py")

    # Optimize the smaller tables (which are already sorted, placing human first)
    Time(1)
    Optimize("alphafrag")
    Time(1)
    Optimize("alphaseq")
    Time(1)

    # Create a sorted ASA table (since inserts into it were concurrent from different jobs), placing human first
    Time(2)
    State("Creating sorted table 'alphasa_sorted'...")
    Query("DROP TABLE IF EXISTS alphasa_sorted")
    Query("CREATE TABLE alphasa_sorted LIKE alphasa")
    Query("INSERT INTO alphasa_sorted SELECT NULL, acc, species, tax, frags, afdb, site, aa, plddt, plddt10, asa, asa10, relasa, relasa10, dis, dis10, surf, surf10, sec FROM alphasa ORDER BY species='HUMAN' DESC, species, acc, afdb, site")
    Time(2)
    Optimize("alphasa_sorted")
    Time(2)
    State("Replacing table 'alphasa' with sorted table 'alphasa_sorted'...")
    Query("ALTER TABLE alphasa RENAME TO alphasa_raw")
    Query("ALTER TABLE alphasa_sorted RENAME TO alphasa")
    Query("DROP TABLE alphasa_raw")
    Time(2)
    # # Will need a large /var/tmp area (~500 GB) for filesort
    # Sort existing ASA table (since inserts into it were concurrent from different jobs) using ALTER TABLE
    # Time(2)
    # State("Sorting table 'alphasa'...")
    # Query("ALTER TABLE alphasa ORDER BY species, tax, acc, afdb, site")
    # Time(2)
    # Optimize("alphasa")
    # Time(2)

    # Create a sorted contacts table (since inserts into it were concurrent from different jobs), placing human first
    Time(3)
    State("Creating sorted table 'alphacon_sorted'...")
    Query("DROP TABLE IF EXISTS alphacon_sorted")
    Query("CREATE TABLE alphacon_sorted LIKE alphacon")
    Query("INSERT INTO alphacon_sorted SELECT NULL, acc, species, tax, frags, afdb, site1, site2, aa1, aa2, type, atom1, atom2, dist FROM alphacon ORDER BY species='HUMAN' DESC, species, acc, afdb, site1, site2, type, atom1, atom2")
    Time(3)
    Optimize("alphacon_sorted")
    Time(3)
    State("Replacing table 'alphacon' with sorted table 'alphacon_sorted'...")
    Query("ALTER TABLE alphacon RENAME TO alphacon_raw")
    Query("ALTER TABLE alphacon_sorted RENAME TO alphacon")
    Query("DROP TABLE alphacon_raw")
    Time(3)
    # # Sort existing contacts table (since inserts into it were concurrent from different jobs) using ALTER TABLE
    # # Will need a large /var/tmp area (~500 GB) for filesort
    # Time(3)
    # State("Sorting table 'alphacon'...")
    # Query("ALTER TABLE alphacon ORDER BY species, tax, acc, afdb, site1, site2, type, atom1, atom2")
    # Time(3)
    # Optimize("alphacon")
    # Time(3)

    # Fill PAE column in alphacon
    Time(4)
    Run("Parse AlphaFold PAE scores from JSON files into SQL table 'alphacon'", "alphacon_add_pae.py")
    Time(4)
    Optimize("alphacon")
    Time(4)



    # Verification queries for table completeness etc.:
    State("Running verification queries...")

    # Verify that alphasa proteins are complete (each residue is present)
    Time(11)
    q = "SELECT a.acc, COUNT(DISTINCT a.site) AS sites, LENGTH(s.seq) AS seqlen FROM alphasa a LEFT OUTER JOIN alphaseq s ON a.acc=s.acc WHERE s.acc IS NULL GROUP BY a.acc LIMIT 1;"
    print(q)
    query = Query(q)
    if Numrows(query) > 0:
        Warn("Warning: Some alphasa proteins are incomplete (each residue is not present)")
    else:
        print(" >> OK")
    Time(11)
    # >> None have s.acc NULL, i.e. the LEFT OUTER JOIN always works. Great!

    Time(12)
    q = "SELECT a.acc, COUNT(DISTINCT a.site) AS sites, LENGTH(s.seq) AS seqlen FROM alphasa a LEFT OUTER JOIN alphaseq s ON a.acc=s.acc GROUP BY a.acc HAVING sites!=seqlen LIMIT 1;"
    print(q)
    query = Query(q)
    if Numrows(query) > 0:
        Warn("Warning: Some alphasa proteins have the same site present twice (no COUNT DISTINCT), i.e. did multiple inserts occur for some accs?")
    else:
        print(" >> OK")
    Time(12)
    # >> None! alphasa is always complete. Awesome!

    # Do alphasa proteins ever have the same site present twice (no COUNT DISTINCT), i.e. did multiple inserts occur for some accs?
    Time(13)
    q = "SELECT a.acc, COUNT(a.site) AS sites, LENGTH(s.seq) AS seqlen FROM alphasa a LEFT OUTER JOIN alphaseq s ON a.acc=s.acc GROUP BY a.acc HAVING sites!=seqlen LIMIT 1;"
    print(q)
    query = Query(q)
    if Numrows(query) > 0:
        Warn("Warning: Some alphasa proteins have the same site present twice, i.e. multiple inserts occurred for some accs")
    else:
        print(" >> OK")
    Time(13)
    # >> None! Nothing was inserted twice in alphasa.

    # How about alphacon? Is the same site1|site2|type ever present twice?
    Time(14)
    q = "SELECT a.acc, COUNT(a.id) AS sites, COUNT(DISTINCT site1, site2, type) AS distinctsites FROM alphacon a LEFT OUTER JOIN alphaseq s ON a.acc=s.acc GROUP BY a.acc HAVING sites!=distinctsites LIMIT 1;"
    print(q)
    query = Query(q)
    if Numrows(query) > 0:
        Warn("Warning: Some alphacon proteins have the same site1|site2|type present twice")
    else:
        print(" >> OK")
    Time(14)
    # >> None! Nothing was inserted twice in alphacon.

    # How many accs still have pae=NULL after running alphacon_add_pae.py? (Expecting the 208 fragmented human proteins)
    query = Query(f"SELECT DISTINCT acc FROM alphacon WHERE pae IS NULL ORDER BY acc")
    print(f"\nNote: {Numrows(query)} accs in table 'alphacon' still have PAE score NULL (expected number: the 208 fragmented human proteins)\n")
    if (Numrows(query) != 208):
        Warn(f"{Numrows(query)} accs in table 'alphacon' still have PAE score NULL (expected number: the 208 fragmented human proteins)")
    for acc, in query:
        Log(f"acc still has PAE score NULL for acc", acc)
    print()
    Show("acc still has PAE score NULL for acc")
    print()



else:

    # AlphaSync: Synchronize with latest UniProt release

    # Update alphauniprot (latest UniProt release via API)
    # Run("Update alphauniprot (latest UniProt release, via UniProt API)", "alphauniprot.py -compara_species")
    Run("Update alphauniprot (latest UniProt release, via UniProt API)", "alphauniprot.py -compara_species -refresh") # Re-download protein annotation and sequences from API even if local UniProt release is still current
    Run("Get canonical reference proteome status (latest UniProt release, via UniProt FTP server) and update table 'alphauniprot' (column 'refprotcanon')", "alphauniprot_refprotcanon.py")
    Run("Update table 'alphauniprot_symbols', a normalised mapping of UniProt gene symbols and synonyms to UniProt accessions for efficient lookup", "alphauniprot_symbols.py")

    # alphasync.py: Run AlphaFold on proteins whose sequences have been updated, in species and proteins of interest

    # Synchronize (update) human proteins
    Run("AlphaSync: Run on human", f"alphasync.py uniprot {local_uniprot_release} all -humanonly -nofrag -cpufirst")
    Run("AlphaSync: Run on human (including >2700 aa, which will be fragmented)", f"alphasync.py uniprot {local_uniprot_release} all -humanonly -cpufirst")
    Run("AlphaSync: Run on human (including isoforms)", f"alphasync.py uniprot {local_uniprot_release} all -humanonly -iso -nofrag -cpufirst")
    Run("AlphaSync: Run on human (including >2700 aa, which will be fragmented) (including isoforms)", f"alphasync.py uniprot {local_uniprot_release} all -humanonly -iso -cpufirst")

    # Synchronize (update) model organisms & global health proteomes (includes human) (see AlphaFold DB), as long as there are ≤1000 proteins that require an updated structure for a given species
    Run("AlphaSync: Run on model organisms & global health proteomes", f"alphasync.py uniprot {local_uniprot_release} 1000 -modelhealthonly -nofrag -cpufirst")
    Run("AlphaSync: Run on model organisms & global health proteomes (including >2700 aa, which will be fragmented)", f"alphasync.py uniprot {local_uniprot_release} 1000 -modelhealthonly -cpufirst")
    Run("AlphaSync: Run on model organisms & global health proteomes (including isoforms)", f"alphasync.py uniprot {local_uniprot_release} 1000 -modelhealthonly -iso -nofrag -cpufirst")
    Run("AlphaSync: Run on model organisms & global health proteomes (including >2700 aa, which will be fragmented) (including isoforms)", f"alphasync.py uniprot {local_uniprot_release} 2500 -modelhealthonly -iso -cpufirst")

    # # Synchronize (update) Ensembl Compara species, as long as there are ≤100 proteins that require an updated structure for a given species
    # # Run("AlphaSync: Run on Ensembl Compara species", f"alphasync.py uniprot {local_uniprot_release} 100 -comparaonly -nofrag -cpufirst")
    # # Run("AlphaSync: Run on Ensembl Compara species (including >2700 aa, which will be fragmented)", f"alphasync.py uniprot {local_uniprot_release} 100 -comparaonly -cpufirst")
    # # Run("AlphaSync: Run on Ensembl Compara species (including isoforms)", f"alphasync.py uniprot {local_uniprot_release} 100 -comparaonly -iso -nofrag -cpufirst")
    # Run("AlphaSync: Run on Ensembl Compara species (including >2700 aa, which will be fragmented) (including isoforms)", f"alphasync.py uniprot {local_uniprot_release} 100 -comparaonly -iso -cpufirst")

    # Parse new AlphaSync structures

    Run("Compress all CIF/PAE/params output files into .tar archives", f"scripts/migrate_alphasync_tar_archives.py {local_uniprot_release}")
    Run("Parse updated fragment sequences into table 'alphafrag'", "alphafrag.py -alphasync")
    Run("Combine updated fragment sequences from table 'alphafrag' into table 'alphaseq'", "alphaseq.py -alphasync")
    Run("Main: Submit jobs that run DSSP and Lahuta on fragments and combine output across fragments, insert ASA values into table 'alphasa', and contacts into table 'alphacon'", "main.py -alphasync")
    Run("Parse AlphaFold PAE scores from JSON files into SQL table 'alphacon'", "alphacon_add_pae.py -alphasync")

    Run("Update table 'alphamap' that maps UniProt sequences to previously existing AFDB or new AlphaSync structures", f"alphamap_uniprot.py {local_uniprot_release}")
    Run("Remove AlphaSync prediction (CIF, PAE and params files) that are no longer necessary since there are better structures available for their sequences", f"alphasync_cleanup.py {local_uniprot_release}")
    Run("Compress all CIF/PAE/params output files into .tar archives (now minus the obsolete accessions)", f"scripts/migrate_alphasync_tar_archives.py {local_uniprot_release}")

    Run("Update table 'alphauniprot_species', a summary table of UniProt taxon IDs, latin names and common species names for efficient lookup (requires alphamap to be updated first)", "alphauniprot_species.py")
    Run("Migrate to optimized SQL tables (alphasync_compact) for tablespace export to web server", f"scripts/migrate_alphasync_compact.py")
    Run("Update precalculated statistics in SQL table 'alphastats'", f"alphastats.py")
    Run("Migrate table 'alphastats'", f"scripts/migrate_alphasync_compact_alphastats.py")

    # Successfully finished: update locally recorded UniProt release version
    # update_local_uniprot_release()
    check_uniprot_release(update = 1)

print("Done!")
