#!/usr/bin/env python3
"""
Extract structure files from AlphaFold Protein Structure Database TAR archives
"""

# Initialize
# import pandas as pd
# import numpy as np
# import networkx as nx
# import seaborn as sns
# import matplotlib.pyplot as mp
import tarfile
import gzip
import io
from blang_mysql import *
from blang import *

type = Args(1, "[structure type: cif/pdb]\n\n -nofrag: Exclude structures of proteins that needed to be fragmented due to protein size (>2700 aa)", "cif")

alphamap = "alphamap"
outpath = f"output_{type}"

# Verify that type is either "cif" or "pdb"
if type not in ("cif", "pdb"):
    Die("Error: Unhandled structure type (expecting 'cif' or 'pdb')")
if type == "pdb":
    Warn(f"Warning: Not all structures are available as 'pdb' files in the AlphaFold Protein Structure Database. Use 'cif' format to ensure structure availability.")

# Start

Run("Make temporary directory", f"mkdir -p {outpath}", silent=True)

# Get source .tar files and list of AlphaSync accessions (via table 'alphamap') for all human proteins in table 'uniprot'
wanted_accs = FetchSet(Query(f"SELECT DISTINCT m.map FROM uniprot p, {alphamap} m WHERE p.species='human' AND p.acc=m.value AND m.type='uniprot' AND m.version='2022_04' AND m.best=1"))
sources = FetchSet(Query(f"SELECT DISTINCT f.source FROM uniprot p, {alphamap} m, alphafrag f WHERE p.species='human' AND p.acc=m.value AND m.type='uniprot' AND m.version='2022_04' AND m.best=1 AND m.map=f.acc"))

if Switch('nofrag'):
    # Leave out fragmented files
    # (Fragmenting may introduce artefacts in analyses, due to artificially introduced termini)
    fragmented_accs = FetchSet(Query(f"SELECT DISTINCT acc FROM alphafrag WHERE species='human' AND frag!=1"))
    wanted_accs -= fragmented_accs



print(f"\nGetting .{type} files for {len(wanted_accs)} human proteins from {len(sources)} source archives:")

for source in tq(nsort(sources)):
    print(f" >> source {source}") if Switch('debug') else None
    # Reconstruct file path
    # If source starts with "proteome-tax_id-...":
    if rx(r"^proteome-tax_id-", source):
        # from GCS
        infile = f"input/gcs/{source}.tar"
    else:
        # from FTP
        if type == "pdb":
            # replace swissprot_cif (as it is called in table 'alphafrag') with swissprot_{type} to get the correct archive
            source = re.sub(r"swissprot_cif", f"swissprot_{type}", source)
            infile = f"input/ftp/{source}.tar"
    print(f"   >> infile {infile}") if Switch('debug') else None

    # Open TAR file
    # mode "r:": "Open for reading exclusively without compression" (https://docs.python.org/3/library/tarfile.html)
    tar = tarfile.open(infile, mode="r:")
    acc = None

    # Stream TAR archive
    for member in tar:
        print(f"     >> {member.name}") if Switch('debug') else None
        
        # frag = None

        if type == "pdb":
            # Skip .cif.gz (mmCIF) files
            if rx(r"\.cif\.gz$", member.name):
                print(f"       >> skipped cif file") if Switch('debug') else None
                Log("skipped cif file", member.name)
                continue
        elif type == "cif":
            # Skip .pdb.gz files
            if rx(r"\.pdb\.gz$", member.name):
                print(f"       >> skipped pdb file") if Switch('debug') else None
                Log("skipped pdb file", member.name)
                continue

        # Skip JSON files (metadata)
        if rx(r"\.json\.gz$", member.name):
            print(f"       >> skipped json file") if Switch('debug') else None
            Log("skipped json file", member.name)
            continue

        # Get only .{type}.gz files
        if not rx(rf"\.{type}\.gz$", member.name):
            # Non-PDB file
            # print(f"     >> non-cif, skipped")
            # Log("skipped non-cif file", member.name)
            Die(f"Error: Unexpected non-PDB, non-mmCIF file: {member.name}")
            # continue
            
        # Get accession (for handling fragments)
        # AF-A0A009IHW8-F1-model_v2.{type}.gz
        m = rx(rf"^AF-(\w+(-\d+)?)-F(\d+)-model_v\d+\.{type}\.gz$", member.name)
        if m:
            acc = m[0]
            frag = int(m[2])
        else:
            Die(f"Couldn't parse '{member.name}'")

        print(f"       >> {acc}") if Switch('debug') else None
            
        # Keep only wanted accessions
        if acc not in wanted_accs:
            print(f"         >> not in wanted list (skipped)") if Switch('debug') else None
            Log(f"skipped acc since it was not in the wanted_accs list for acc", acc)
            continue
        else:
            print(f"         >> in wanted list (kept)") if Switch('debug') else None
            Log(f"proceeded with acc since it was in the wanted_accs list for acc", acc)

        # Excluding fragment proteins:
        if Switch('nofrag'):
            # There should be no fragment number other than 1
            if frag > 1:
                Die(f"Fragment is '{frag}' for acc '{acc}'")

        # Extract .{type}.gz to .{type} using gzip (streaming, no temporary files)
        gz = tar.extractfile(member)
        # Open extracted .{type} in text mode (streaming, no temporary files)
        struct = gzip.open(gz, mode="rt")

        # Write .{type} to temporary file for DSSP
        structfile = f"{outpath}/{member.name}"

        # Remove .gz from file name
        structfile = re.sub(r"\.gz$", "", structfile)

        # Replace e.g. "_v4" with "-v4" for compatibility with the dMaSIF pipeline
        structfile = re.sub(r"_v(\d+)\.", "-v\\1.", structfile)

        # with open(structfile, "x") as structout:    # mode "x" would fail if file already exists (not desirable for re-runs for randomly failed jobs)
        with open(structfile, "w") as structout:
            print(struct.read(), file=structout)

        # Mark accession as retrieved
        wanted_accs -= {acc}
        Log("successfully retrieved structure file for acc", acc)
        Log("successfully retrieved structure file for archive", infile)
        print(f"           >> extracted (success!)") if Switch('debug') else None

Show(lim=0, sort=True)

# Check if all wanted accessions were retrieved
if len(wanted_accs) > 0:
    Warn(f"Warning: {len(wanted_accs)} wanted accessions were not retrieved:\n\n{wanted_accs}\n\n\n")

print("\nDone!")
