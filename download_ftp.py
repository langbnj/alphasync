#!/usr/bin/env python3
"""
Download: Download source data from FTP server
"""

# Initialize

# from blang_mysql import *
from blang import *

Args()

# Create input directory
Run("Create FTP input directory", "mkdir -p input/ftp")

# Download

# Get change log
Run("Download", "wget -nv -N -P input/ftp 'http://ftp.ebi.ac.uk/pub/databases/alphafold/CHANGELOG.txt'")

# # Get accession list
# Run("Download", "wget -N -P input/ftp 'http://ftp.ebi.ac.uk/pub/databases/alphafold/accession_ids.csv'")
# wc -l input/gcs/accession_ids.csv 
# >> Contains 214,683,829 accessions

# # Get FASTA sequences
# Run("Download", "wget -nv -N -P input/ftp 'http://ftp.ebi.ac.uk/pub/databases/alphafold/sequences.fasta'")
# Better: gsutil -m cp 'gs://public-datasets-deepmind-alphafold/sequences.fasta' input
# gsutil -m cp 'gs://public-datasets-deepmind-alphafold/sequences.fasta' -|head
# >> This is a different file, though. FTP is from 2022-10-20, GCS is older, from 2022-07-21.
# GCS does not contain fragments (see download_gcs script). It's incomplete.
# ~/update/alphafold/input/ftp >> g -A1 '^>AFDB:AF-\w+\-F2 ' sequences.fasta | head
# >> FTP does not contain fragments either. It's incomplete.
# The two sequences.fasta files contain the same sequences. They are only sorted differently. See test_sequences_fasta: l commands.sh.
# How many sequences are in sequences.fasta?
# grep -iPc '^>' input/gcs/sequences.fasta 
# >> 214,683,829 sequences (headers, "^>")
# >> Same as in accession_ids.csv
# How many per species?
# >> For the Ensembl Compara species, the counts in sequences.fasta match those in alphafrag for v4, except for human, where I have 208 more (the fragmented proteins).
# >> The input sequences didn't change from v3 to v4, so it's okay that the date on sequences.fasta is older.
# Is it complete for most species?
# >> Apparently yes. >1,000,000 species.
# Does it have all the species I'd expect?
# >> It definitely has all the 200 Ensembl Compara species.
# Conclusion: >> The sequences.fasta files aren't useful (since they exclude the 208 human fragmented proteins).

# Get metadata (also includes sequences)
# gsutil -m cp 'gs://public-datasets-deepmind-alphafold/metadata/*' input/metadata

if Switch('humanonly'):
    # Download human only

    # HTTPS
    Run("Download", "wget -nv -N -r -e robots=off -l 1 -P input/ftp -A tar --accept-regex='(swissprot|_HUMAN_).+\.tar$' --reject-regex='(swissprot_pdb).+\.tar$' -R 'mane_overlap_v*.tar' -nd 'http://ftp.ebi.ac.uk/pub/databases/alphafold/latest/'")
    
else:
    # Download all species

    # FTP (sometimes works, sometimes doesn't (all requests time out))
    # Run("Download", "wget -N -P input/ftp 'ftp://ftp.ebi.ac.uk/pub/databases/alphafold/latest/*.tar'")

    # HTTPS
    Run("Download", "wget -nv -N -r -e robots=off -l 1 -P input/ftp -A tar --reject-regex='(swissprot_pdb).+\.tar$' -R 'mane_overlap_v*.tar' -nd 'http://ftp.ebi.ac.uk/pub/databases/alphafold/latest/'")


# Use TAR indexer
# Run("Download tarindexer", "wget -O bin/tarindexer.py 'https://github.com/devsnd/tarindexer/raw/master/tarindexer.py'")
# Run("Download tar2index", "wget -O bin/tar2index.py 'https://github.com/lpryszcz/tar_indexer/raw/master/tar2index.py'")
# Run("Refactor tar2index to Python 3", "2to3 -w bin/tar2index.py")
# Run("Index TAR files", "bin/tar2index.py -v -i input/ftp/*.tar -d input/tar_index.db3")
# >> Actually: Both TAR indexers are too slow when retrieving data. They both take around 0.3s. For the ~500,000 entries in the Swiss-Prot file alone, it would take >40 hours simply to retrieve them one by one, prior to any processing.
# >> Not using TAR indexers. I need to decompress the data.
# >> The file system is much faster. It only takes around 0.02 seconds to access a file, resulting in ~3 hours to get all ~500,000 entries from the Swiss-Prot file in a specific sequence.
# >> cat is 3-4x faster than zcat. However, decompressing increases the size 5-fold.

# Only human has fragments in v4, for 208 proteins (207 in Swiss-Prot):
# ~/update/alphafold/input/ftp >> 1 *_v4.tar | xa tar -tf {} | g '^AF-\w+\-F2-model'
# >> 208 accessions. UniProt's ID retrieval shows that all of these are human (207 reviewed, 1 unreviewed) (207 in table 'uniprot' as well, all human).


Run("Show change log", "cat input/ftp/CHANGELOG.txt")
Run("Show directory", "ls -ltrh input/ftp")
Run("Count files", "ls -1 input/ftp/*.tar | wc -l")

Done()
