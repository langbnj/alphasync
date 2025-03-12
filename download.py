#!/usr/bin/env python3
"""
Download: Download source data
"""

# Initialize

# from blang_mysql import *
from blang import Run, Cd

# Download

# # From Bludau et al.:
# Run("Make IDR definition input file directory", "mkdir -p input/idrs")
# Run("Download", "wget -NP input/idrs 'https://raw.githubusercontent.com/MannLabs/structuremap_analysis/master/data/disorder_data/IDR_ROC_curve.csv'")
# Run("Show directory", "ls -lah input/idrs")
# # Note: Bludau et al. was ultimately not used (see discussion in job_dssp.py).

# Human CHESS sequence dataset from Sommer et al. (https://elifesciences.org/articles/82556)
# https://figshare.com/articles/dataset/Structure-guided_isoform_identification_for_the_human_transcriptome/21802476/1
# Notes:
# Possible human isoform structure source
# Used OpenFold with 2021 datasets (uniref30 etc.)
# Scanned for ORFs
# Only predicted proteins ≤1000 aa
# Could be useful, but not sure (since I've already got many human isoforms predicted by AlphaSync)
# Not sure how many actual human isoforms are included here
# Format is PDB, so I'd need to convert to mmCIF etc.
# They also filtered by pLDDT and aborted if scores were below 60 or so, so coverage is not complete. They also stopped early if there was a high-pLDDT structure, but it might not be the optimal one.
# All in all, I'd rather just proceed with AlphaSync.
# See also ~/scripts/test_humaniso.py.
Run("Make Sommer et al. input file directory", "mkdir -p test_humaniso")
# Download only the FASTA file:
Run("Download", "wget -NO test_humaniso/chess_structure_CDS_aa_v1.2.faa 'https://figshare.com/ndownloader/files/38695569'")

# ASpdb alternative splicing database
# https://pubmed.ncbi.nlm.nih.gov/39530217/
# https://biodataai.uth.edu/ASpdb/
# https://biodataai.uth.edu/ASpdb/tables/alternative_seq_tab.txt
# See also ~/scripts/test_aspdb.py.
Run("Make ASpdb input file directory", "mkdir -p test_aspdb")
# Download only the FASTA file:
Run("Download", "wget -NO test_aspdb/alternative_seq_tab.txt 'https://biodataai.uth.edu/ASpdb/tables/alternative_seq_tab.txt'")


# BFVD (viral structures)
# https://academic.oup.com/nar/article/53/D1/D340/7906834
# https://bfvd.steineggerlab.workers.dev
Run("Make BFVD viral structures input file directory", "mkdir -p input/bfvd")
Run("Download", "wget -O input/bfvd/bfvd.tar.gz 'https://bfvd.steineggerlab.workers.dev/bfvd.tar.gz'")
# Run("Unpack", "tar -xzf input/bfvd/bfvd.tar.gz -C input/bfvd")
# 351,243 viral structures
# One per acc
# .tar head: A0A2I6SBR6_unrelaxed_rank_001_alphafold2_ptm_model_2_seed_000.pdb
# .tar tail: A0A8S5U9K4_unrelaxed_rank_001_alphafold2_ptm_model_3_seed_000.pdb
# File names seem to keep track of OpenFold model used (1-5), rank always seems to be 001, and seed always seems to be 000

# Not used
# # Get mkdssp 4.4.0
# # https://github.com/PDB-REDO/dssp/releases/download/v4.4.0/mkdssp-4.4.0-linux-x64

# Get gemmi and command line program (for AF 2.3.2 mmCIF output fixing by conversion to PDB and back)
# pip install gemmi
# pip install gemmi-program

# Show directory
Run("Show directory", "ls -lah input/humaniso")

print("Done!")
