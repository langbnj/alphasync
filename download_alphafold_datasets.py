#!/usr/bin/env python3
"""
Download: Download latest versions of datasets used by AlphaFold (sequences, template PDB structures etc.)
"""

# Initialize

# from blang_mysql import *
from blang import *

# Create directories in input/alphasync
Run("Create 'alphafold_tools' directory", "mkdir -p input/alphasync/alphafold_tools")
Run("Create 'alphafold_datasets' directory", "mkdir -p input/alphasync/alphafold_datasets")

# Download

# Tools
# Distributed with repo

# hhblits v3.0-beta.3 (still used in AlphaFold 3)
# https://github.com/soedinglab/hh-suite/releases/download/v3.0-beta.3/hhsuite-3.0-beta.3-Linux.tar.gz
Run("Download hhblits v3.0-beta.3", "wget -nv -N -P input/alphasync/alphafold_tools 'https://github.com/soedinglab/hh-suite/releases/download/v3.0-beta.3/hhsuite-3.0-beta.3-Linux.tar.gz'")
Run("Decompress", "tar -xzf input/alphasync/alphafold_tools/hhsuite-3.0-beta.3-Linux.tar.gz -C input/alphasync/alphafold_tools")

# Datasets

# Small BFD
Run("Make directory", "mkdir -p input/alphasync/alphafold_datasets/small_bfd")
Run("Download MGnify clusters", "wget -nv -N -O input/alphasync/alphafold_datasets/small_bfd/bfd-first_non_consensus_sequences.fasta.gz 'https://storage.googleapis.com/alphafold-databases/reduced_dbs/bfd-first_non_consensus_sequences.fasta.gz'")
Run("Decompress", "gunzip -f input/alphasync/alphafold_datasets/small_bfd/bfd-first_non_consensus_sequences.fasta.gz")

# MGnify 2024_04 (latest)
Run("Make directory", "mkdir -p input/alphasync/alphafold_datasets/mgnify")
Run("Download MGnify clusters", "wget -nv -N -O input/alphasync/alphafold_datasets/mgnify/mgy_clusters_2024_04.fa.gz 'https://ftp.ebi.ac.uk/pub/databases/metagenomics/peptide_database/2024_04/mgy_clusters.fa.gz'")
Run("Decompress", "gunzip -f input/alphasync/alphafold_datasets/mgnify/mgy_clusters_2024_04.fa.gz")

# pdb70_from_mmcif_220313 (latest)
Run("Make directory", "mkdir -p input/alphasync/alphafold_datasets/pdb70_from_mmcif_220313")
Run("Download pdb70", "wget -nv -N -O input/alphasync/alphafold_datasets/pdb70_from_mmcif_220313/pdb70_from_mmcif_220313.tar.gz 'https://wwwuser.gwdguser.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/pdb70_from_mmcif_220313.tar.gz'")
Run("Decompress", "tar -xzf input/alphasync/alphafold_datasets/pdb70_from_mmcif_220313/pdb70_from_mmcif_220313.tar.gz -C input/alphasync/alphafold_datasets/pdb70")
Run("Clean up", "rm -fv input/alphasync/alphafold_datasets/pdb70_from_mmcif_220313/pdb70_from_mmcif_220313.tar.gz")

# pdb100 2023_05_17 (latest)
# pdb70 alternative from https://colabfold.mmseqs.com/ (ColabFold uses pdb100 instead of pdb70, and renames it to pdb70 for compatibility, as mentioned on www.colabfold.com):
Run("Make directory", "mkdir -p input/alphasync/alphafold_datasets/pdb100_foldseek_230517")
# https://wwwuser.gwdg.de/~compbiol/colabfold/pdb100_230517.fasta.gz (not needed I think)
Run("Download pdb100 fasta", "wget -nv -N -O input/alphasync/alphafold_datasets/pdb100_foldseek_230517/pdb100_230517.fasta.gz 'https://wwwuser.gwdg.de/~compbiol/colabfold/pdb100_230517.fasta.gz'")
Run("Decompress", "gunzip -f input/alphasync/alphafold_datasets/pdb100_foldseek_230517/pdb100_230517.fasta.gz")
# https://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/pdb100_foldseek_230517.tar.gz
Run("Download pdb100", "wget -nv -N -O input/alphasync/alphafold_datasets/pdb100_foldseek_230517/pdb100_from_mmcif_230517.tar.gz 'https://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/pdb100_foldseek_230517.tar.gz'")
Run("Decompress", "tar -xzf input/alphasync/alphafold_datasets/pdb100_foldseek_230517/pdb100_from_mmcif_230517.tar.gz -C input/alphasync/alphafold_datasets/pdb100")

# UniRef30 2023_02 (latest) and prior versions
for version in ['2023_02', '2022_02', '2021_03', '2020_06']:
    dir = "input/alphasync/alphafold_datasets/uniref30"
    Run("Make directory", f"mkdir -p {dir}")
    Run("Download UniRef90", f"wget -nv -N -P {dir} 'https://wwwuser.gwdguser.de/~compbiol/uniclust/{version}/UniRef30_{version}_hhsuite.tar.gz'")
    Run("Decompress", f"tar -xzf {dir}/UniRef30_{version}_hhsuite.tar.gz -C {dir}")
    Run("Clean up", f"rm -fv {dir}/UniRef30_{version}_hhsuite.tar.gz")

# uniclust30 (what UniRef30 used to be called) 2018_08
for version in ['2018_08']:
    dir = "input/alphasync/alphafold_datasets/uniclust30"
    Run("Make directory", f"mkdir -p {dir}")
    Run("Download UniRef90", f"wget -nv -N -P {dir} 'https://wwwuser.gwdguser.de/~compbiol/uniclust/{version}/uniclust30_{version}_hhsuite.tar.gz'")
    Run("Decompress", f"tar -xzf {dir}/uniclust30_{version}_hhsuite.tar.gz -C {dir}")
    Run("Clean up", f"rm -fv {dir}/uniclust30_{version}_hhsuite.tar.gz")

# UniRef90 2024_05 (latest)
Run("Make directory", "mkdir -p input/alphasync/alphafold_datasets/uniref90")
Run("Download UniRef90", "wget -nv -N -O input/alphasync/alphafold_datasets/uniref90_2024_05.fasta.gz 'https://ftp.ebi.ac.uk/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz'")
Run("Decompress", "gunzip -f input/alphasync/alphafold_datasets/uniref90_2024_05.fasta.gz")



# PDB mmCIF files (pdb_mmcif) (latest)

root_dir = f"input/alphasync/alphafold_datasets/pdb_mmcif"
raw_dir = f"{root_dir}/raw"
mmcif_dir = f"{root_dir}/mmcif_files"

Run("Make directory", f"mkdir -p {raw_dir}")
Run("Make directory", f"mkdir -p {mmcif_dir}")

# Run("Download pdb_mmcif", "download_pdb_mmcif.sh input/alphasync/alphafold_datasets")

# Run("Get latest mmCIF files from EBI PDBe via rsync", f"rsync --recursive --links --perms --times --compress --info=progress2 --delete rsync.ebi.ac.uk::pub/databases/pdb/data/structures/divided/mmCIF/ '{raw_dir}'")
# Run("Get latest mmCIF files from EBI PDBe via rsync", f"rsync -vi --recursive --links --perms --times --compress --progress --delete rsync://ftp.ebi.ac.uk/pub/databases/pdb/data/structures/divided/mmCIF/ '{raw_dir}'")
Run("Get latest mmCIF files from RCSB PDB via rsync", f"rsync -vi -rlpt --progress --delete --port=33444 rsync.rcsb.org::ftp_data/structures/divided/mmCIF/**/* '{raw_dir}'")
Run("Download obsolete.dat (list of superseded PDB IDs) from RCSB PDB", f"wget -nv -N -P {root_dir} 'https://files.wwpdb.org/pub/pdb/data/status/obsolete.dat'")

# # FTP instead (much slower - only use if port 33444 is blocked)
# # Run("Download only updated or new mmCIF files (timestamping, -N) from PDB", f"wget -N -r -nv -nd --retr-symlinks -P '{raw_dir}' 'ftp://ftp.ebi.ac.uk/pub/databases/pdb/data/structures/divided/mmCIF/*'")
# # Run("Download only updated or new mmCIF files (timestamping, -N) from PDB, retrying until complete", f"while ! wget -N -r -nv -nd --retr-symlinks -P '{raw_dir}' 'ftp://ftp.ebi.ac.uk/pub/databases/pdb/data/structures/divided/mmCIF/*'; do :; done")
# Run("Download only updated or new mmCIF files (timestamping, -N) from PDB, retrying until complete", f"while ! wget -N -r -nv -nd --retr-symlinks -P '{raw_dir}' 'ftp://ftp.ebi.ac.uk/pub/databases/pdb/data/structures/divided/mmCIF/*'; do sleep 5; done")
# Run("Clean up", f"rm -fv {raw_dir}/.listing")

Run("Unpack mmCIF files", f"find '{raw_dir}' -type f -iname '*.gz' -exec gunzip -kf " + "{}" + " +")
Run("Move all unpacked mmCIF files to final destination", f"find '{raw_dir}' -type f -iname '*.cif' -exec mv " + "{}" + f" {mmcif_dir}" + " \;")

# Run("Unzipping all mmCIF files", f"find '{raw_dir}/' -type f -iname '*.gz' -exec gunzip " + "{}" + " +")
# Run("Deleting empty directories", f"find '{raw_dir}' -type d -empty -delete")
# # Flattening all mmCIF files
# # for subdir in "${RAW_DIR}"/*; do
# #   mv "${subdir}/"*.cif "${MMCIF_DIR}"
# # done
# Run("Flattening all mmCIF files", f"for subdir in '{raw_dir}'/*; do mv '$subdir/'*.cif '{mmcif_dir}'; done")
# Run("Delete empty download directory structure", f"find '{raw_dir}' -type d -empty -delete")



Run("Show directory", "ls -ltrh input/alphasync/alphafold_datasets")

Done()
