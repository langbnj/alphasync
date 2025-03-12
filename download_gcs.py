#!/usr/bin/env python3
"""
Download: Download proteomes for all species or for only the 200 Ensembl Compara species (in Ensembl release 108) from Google Cloud Storage via gsutil
"""

# Initialize

from blang_mysql import *
from blang import *

Args(0, " -comparaonly: Only download proteomes for the ~200 species contained in Ensembl Compara", " -comparaonly")

# Create input directory
Run("Create GCS directory", "mkdir -p input/gcs")


# Documentation:
# https://github.com/deepmind/alphafold/blob/main/afdb/README.md

# Google Cloud Storage (GCS) bucket (files):
# https://console.cloud.google.com/storage/browser/public-datasets-deepmind-alphafold-v4;tab=objects?prefix=&forceOnObjectsSortingFiltering=false

# Google BigQuery metadata table:
# https://console.cloud.google.com/bigquery?p=bigquery-public-data&d=deepmind_alphafold&page=dataset


# Dataset exploration

# Object queries
# ... (but not A...)
# gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/B*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/C*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/D*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/E*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/F*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/G*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/H*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/I*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/J*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/K*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/L*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/M*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/N*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/O*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/P*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/Q*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/R*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/S*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/T*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/U*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/V*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/W*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/X*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/Y*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/Z*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/a*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/b*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/c*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/d*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/e*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/f*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/g*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/h*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/i*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/j*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/k*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/l*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/m*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/n*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/o*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/p*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/q*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/r*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/s*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/t*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/u*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/v*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/w*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/x*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/y*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/z*'; 
# A... (but not AF...)
# gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/AA*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/AB*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/AC*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/AD*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/AE*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/AG*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/AH*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/AI*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/AJ*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/AK*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/AL*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/AM*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/AN*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/AO*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/AP*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/AQ*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/AR*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/AS*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/AT*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/AU*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/AV*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/AW*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/AX*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/AY*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/AZ*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/Aa*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/Ab*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/Ac*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/Ad*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/Ae*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/Af*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/Ag*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/Ah*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/Ai*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/Aj*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/Ak*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/Al*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/Am*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/An*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/Ao*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/Ap*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/Aq*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/Ar*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/As*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/At*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/Au*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/Av*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/Aw*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/Ax*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/Ay*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/Az*'

# Objects found
# gs://public-datasets-deepmind-alphafold-v4/AF-A0A000-F1-confidence_v4.json
# gs://public-datasets-deepmind-alphafold-v4/AF-A0A000-F1-model_v4.cif
# gs://public-datasets-deepmind-alphafold-v4/AF-A0A000-F1-predicted_aligned_error_v4.json
# gs://public-datasets-deepmind-alphafold-v4/AF-A0A001-F1-confidence_v4.json
# gs://public-datasets-deepmind-alphafold-v4/AF-A0A001-F1-model_v4.cif
# gs://public-datasets-deepmind-alphafold-v4/AF-A0A001-F1-predicted_aligned_error_v4.json
# ...
# gs://public-datasets-deepmind-alphafold-v4/AF-...-F1-confidence_v4.json
# gs://public-datasets-deepmind-alphafold-v4/AF-...-F1-model_v4.cif
# gs://public-datasets-deepmind-alphafold-v4/AF-...-F1-predicted_aligned_error_v4.json
# gs://public-datasets-deepmind-alphafold-v4/accession_ids.csv
# gs://public-datasets-deepmind-alphafold-v4/manifests/
# gs://public-datasets-deepmind-alphafold-v4/metadata/
# gs://public-datasets-deepmind-alphafold-v4/proteomes/
# gs://public-datasets-deepmind-alphafold-v4/sequences.fasta

# Check subfolders as above, for surprising files:
# Subfolders:
# gs://public-datasets-deepmind-alphafold-v4/manifests/
# gs://public-datasets-deepmind-alphafold-v4/metadata/
# gs://public-datasets-deepmind-alphafold-v4/proteomes/
# Object queries:
# gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/manifests/*'
# gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/metadata/*' | head
# gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/metadata/A*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/metadata/B*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/metadata/C*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/metadata/D*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/metadata/E*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/metadata/F*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/metadata/G*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/metadata/H*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/metadata/I*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/metadata/J*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/metadata/K*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/metadata/L*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/metadata/M*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/metadata/N*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/metadata/O*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/metadata/P*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/metadata/Q*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/metadata/R*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/metadata/S*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/metadata/T*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/metadata/U*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/metadata/V*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/metadata/W*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/metadata/X*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/metadata/Y*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/metadata/Z*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/metadata/a*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/metadata/b*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/metadata/c*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/metadata/d*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/metadata/e*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/metadata/f*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/metadata/h*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/metadata/i*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/metadata/j*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/metadata/k*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/metadata/l*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/metadata/m*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/metadata/n*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/metadata/o*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/metadata/p*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/metadata/q*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/metadata/r*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/metadata/s*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/metadata/t*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/metadata/u*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/metadata/v*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/metadata/w*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/metadata/x*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/metadata/y*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/metadata/z*'
# gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/proteomes/*' | head
# gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/proteomes/A*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/proteomes/B*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/proteomes/C*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/proteomes/D*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/proteomes/E*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/proteomes/F*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/proteomes/G*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/proteomes/H*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/proteomes/I*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/proteomes/J*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/proteomes/K*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/proteomes/L*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/proteomes/M*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/proteomes/N*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/proteomes/O*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/proteomes/P*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/proteomes/Q*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/proteomes/R*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/proteomes/S*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/proteomes/T*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/proteomes/U*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/proteomes/V*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/proteomes/W*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/proteomes/X*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/proteomes/Y*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/proteomes/Z*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/proteomes/a*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/proteomes/b*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/proteomes/c*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/proteomes/d*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/proteomes/e*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/proteomes/f*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/proteomes/g*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/proteomes/h*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/proteomes/i*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/proteomes/j*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/proteomes/k*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/proteomes/l*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/proteomes/m*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/proteomes/n*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/proteomes/o*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/proteomes/q*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/proteomes/r*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/proteomes/s*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/proteomes/t*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/proteomes/u*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/proteomes/v*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/proteomes/w*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/proteomes/x*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/proteomes/y*'; gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/proteomes/z*'
# >> No surprising files beyond the ones listed below:

# There are 1,007,292 taxon IDs:
# cd test_gcs_taxa
# \gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/proteomes/*' > gcs_proteomes.txt
# cat gcs_proteomes.txt | g -o "tax_id-\d+-" | g -o "\d+" | uniq | wcl
# cat gcs_proteomes.txt | g -o "tax_id-\d+-" | g -o "\d+" | uniq | suq | wcl


# manifests:
# cat manifest-confidence_v4_json-part-001.csv | g -v "^AF-\w{6,11}-F1-confidence_v4\.json$"
# cat manifest-model_v4_cif-part-001.csv | g -v "^AF-\w{6,11}-F1-model_v4\.cif$"
# cat manifest-predicted_aligned_error_v4_json-part-001.csv | g -v "^AF-\w{6,11}-F1-predicted_aligned_error_v4\.json$"
# >> These are just lists of files

# metadata:
# e.g. gs://public-datasets-deepmind-alphafold-v4/metadata/gcd_metadata-00000-of-10000.json
# Should contain (https://github.com/deepmind/alphafold/blob/main/afdb/README.md):
# Column name	Data type	Description
# allVersions	ARRAY<INT64>	An array of AFDB versions this prediction has had
# entryId	STRING	The AFDB entry ID, e.g. "AF-Q1HGU3-F1"
# fractionPlddtConfident	FLOAT64	Fraction of the residues in the prediction with pLDDT between 70 and 90
# fractionPlddtLow	FLOAT64	Fraction of the residues in the prediction with pLDDT between 50 and 70
# fractionPlddtVeryHigh	FLOAT64	Fraction of the residues in the prediction with pLDDT greater than 90
# fractionPlddtVeryLow	FLOAT64	Fraction of the residues in the prediction with pLDDT less than 50
# gene	STRING	The name of the gene if known, e.g. "COII"
# geneSynonyms	ARRAY<STRING>	Additional synonyms for the gene
# globalMetricValue	FLOAT64	The mean pLDDT of this prediction
# isReferenceProteome	BOOL	Is this protein part of the reference proteome?
# isReviewed	BOOL	Has this protein been reviewed, i.e. is it part of SwissProt?
# latestVersion	INT64	The latest AFDB version for this prediction
# modelCreatedDate	DATE	The date of creation for this entry, e.g. "2022-06-01"
# organismCommonNames	ARRAY<STRING>	List of common organism names
# organismScientificName	STRING	The scientific name of the organism
# organismSynonyms	ARRAY<STRING>	List of synonyms for the organism
# proteinFullNames	ARRAY<STRING>	Full names of the protein
# proteinShortNames	ARRAY<STRING>	Short names of the protein
# sequenceChecksum	STRING	CRC64 hash of the sequence. Can be used for cheaper lookups.
# sequenceVersionDate	DATE	Date when the sequence data was last modified in UniProt
# taxId	INT64	NCBI taxonomy id of the originating species
# uniprotAccession	STRING	Uniprot accession ID
# uniprotDescription	STRING	The name recommended by the UniProt consortium
# uniprotEnd	INT64	Number of the last residue in the entry relative to the UniProt entry. This is equal to the length of the protein unless we are dealing with protein fragments.
# uniprotId	STRING	The Uniprot EntryName field
# uniprotSequence	STRING	Amino acid sequence for this prediction
# uniprotStart	INT64	Number of the first residue in the entry relative to the UniProt entry. This is 1 unless we are dealing with protein fragments.

# Interesting values:
# modelCreatedDate: In file gcd_metadata-00000-of-10000.json at least, this is always 2022-06-01 (g -o "modelCreatedDate[^,]+" gcd_metadata-00000-of-10000.json | suq)
# In fact, using BigQuery, I can confirm modelCreatedDate is always 2022-06-01:
# https://console.cloud.google.com/bigquery?p=bigquery-public-data&d=deepmind_alphafold&page=dataset
# SELECT DISTINCT modelCreatedDate FROM `bigquery-public-data.deepmind_alphafold.metadata`;
# >> 2022-06-01 is the only modelCreatedDate value in AlphaFold DB v4.
# sequenceVersionDate: In file gcd_metadata-00000-of-10000.json at least, this is never newer than 2021-09-29 (probably UniProt 2021_03) (g -o "sequenceVersionDate[^,]+" gcd_metadata-00000-of-10000.json | sort -h)
# Using BigQuery:
# SELECT MIN(sequenceVersionDate), MAX(sequenceVersionDate) FROM `bigquery-public-data.deepmind_alphafold.metadata`;
# >> 2021-09-29 is the latest sequenceVersionDate in AlphaFold DB v4.

# Are there any v3 entries left?
# l gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/AF-*_v3.*'
# >> None. All entries must now be v4.
# BigQuery:
# SELECT DISTINCT latestVersion FROM `bigquery-public-data.deepmind_alphafold.metadata`;
# >> All entries are now v4.


# >> It definitely makes sense to build AlphaSync since the structures in AlphaFold DB v4 are already outdated. They are no newer than UniProt 2021_03 or 2021_04 (https://www.uniprot.org/release-notes).


# proteomes: I've downloaded these.


# There are no fragments in sequences.fasta:
# input/gcs >> g -A1 '^>AFDB:AF-\w+\-F2 ' sequences.fasta | head
# >> None

# Number of sequences in sequences.fasta (which is from 2022-07-21, i.e. older than all structure files):
#  input/gcs >> g -c "^>" sequences.fasta 
# >> 214,683,829 sequences

# Are there fragments in GCS?
# l gsutil ls -d 'gs://public-datasets-deepmind-alphafold-v4/A*-F2-*'
# >> None. No fragments.
# Using BigQuery:
# SELECT taxId, COUNT(*) AS c FROM `bigquery-public-data.deepmind_alphafold.metadata` WHERE entryId LIKE '%-F1' GROUP BY taxId;
# >> Plenty of entries that end in -F1.
# SELECT taxId, COUNT(*) AS c FROM `bigquery-public-data.deepmind_alphafold.metadata` WHERE entryId LIKE '%-F2' GROUP BY taxId;
# >> None. No entries that end in -F2 (confirming that the GCS version of AlphaFold DB does not contain fragments, and that I still need the EBI FTP data).

# >> Therefore, I'll still need the FTP files as well (at least for human)




# Get summary files

# Accessions list (accession_ids.csv) (~7 GB)
# This file contains a list of all the UniProt accessions that have predictions in AlphaFold DB. The file is in CSV format and includes the following columns, separated by a comma:
# UniProt accession, e.g. A8H2R3
# First residue index (UniProt numbering), e.g. 1
# Last residue index (UniProt numbering), e.g. 199
# AlphaFold DB identifier, e.g. AF-A8H2R3-F1
# Latest version, e.g. 4

# acc,start,stop,afid,version
# A0A2I1PIX0,1,200,AF-A0A2I1PIX0-F1,4
# A0A349FBR4,1,269,AF-A0A349FBR4-F1,4
# A0A556C1E4,1,730,AF-A0A556C1E4-F1,4

# Does not contain any fragments:
# input/gcs >> cat accession_ids.csv | cut -d, -f4 | g "\-F2$"
# >> no results

Run("Get list of accessions", "gsutil -m cp -n 'gs://public-datasets-deepmind-alphafold-v4/accession_ids.csv' 'input/gcs/accession_ids.csv'")

# Current sequences (~92 GB)
# This file contains sequences for all proteins in the current database version in FASTA format. The identifier rows start with ">AFDB", followed by the AlphaFold DB identifier and the name of the protein. The sequence rows contain the corresponding amino acid sequences. Each sequence is on a single line, i.e. there is no wrapping.

# >AFDB:afid fullname tags
# >AFDB:AF-A0A1H0ADK1-F1 Ligand-binding SRPBCC domain-containing protein UA=A0A1H0ADK1 UI=A0A1H0ADK1_9SPHI OS=Pedobacter steynii OX=430522 GN=SAMN05421820_10754
# MPEIRLETYINADINLVFDLSRSIDLHQISTADTEEKVVGGRFSGLILLGEQVTWQARHFGLVQQLTAQITALHAPDFFVDEMQKGAFKSFKHEHIFRYNGQQTVMTDVFTYVSPCWIFGRIADYLFLKAYMKRFLFKRNQVIKTYAENGKGSELLLRCQAGFSEA
# >AFDB:AF-A0A3B9RYK9-F1 ABC transporter ATP-binding protein UA=A0A3B9RYK9 UI=A0A3B9RYK9_9DELT OS=Deltaproteobacteria bacterium OX=2026735 GN=DCG23_03755
# MNIEEKQDMCRFIIDINEQYGTTIVLIEHDMGVVMDLSERLIVLDYGRKIAEGTPDEIRGNQVVIDAYLGVAHTN

Run("Get current sequences", "gsutil -m cp -n 'gs://public-datasets-deepmind-alphafold-v4/sequences.fasta' 'input/gcs/sequences.fasta'")




# Get Ensembl Compara species by their taxon IDs
if Switch('comparaonly'):
    for fullspecies, tax in Query("SELECT DISTINCT species, tax FROM compara_species ORDER BY species"):
        # Each shard contains up to 10,000 structures (according to https://github.com/deepmind/alphafold/blob/main/afdb/README.md)
        Run(f"Download species proteome shards for Ensembl Compara species '{fullspecies}'", f"gsutil -m cp -n 'gs://public-datasets-deepmind-alphafold-v4/proteomes/proteome-tax_id-{tax}-*_v*.tar' input/gcs/")
        # At one point I thought of checking hashes before using gsutil cp (to avoid unnecessary copying if the file hasn't been updated), but this is difficult to do with all these different shards.
        # Also, downloading these 200 species's proteomes from GCS only takes around 3.5 hours, so that's tolerable.
else:
    # (never did this, it'd be 17 TB or so)
    Run(f"Download species proteome shards for all species", f"gsutil -m cp -n 'gs://public-datasets-deepmind-alphafold-v4/proteomes/proteome-tax_id-*-*_v*.tar' input/gcs/")

# >> This retrieves 200 species. Great!

# Number of shards per species (taxon ID):
# input/gcs >> 1 | g -o 'tax_id-\d+' | g -o '\d+' | sort | uniq | xa sh -c "ls -1U *tax_id-{}-* | wc -l"|sort -g
# >> 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,7,7,7,7,7,8,8,8,8,8,8,9,10,10,10,11,12,13,19
# >> Most species have more than 10,000 structures, which is good.
# Number of structures in each first shard file (usually 10,000, but sometimes fewer):
# input/gcs >> 1 *-0_v4.tar | xa sh -c "echo -n '{}: '; tar -tf {} | grep -iP '\-model_v4.cif.gz$' | wc -l"
# >> Some of the 200 Ensembl Compara species have very few structures:
# BigQuery:
# SELECT taxId, COUNT(DISTINCT entryId) AS c FROM `bigquery-public-data.deepmind_alphafold.metadata` WHERE taxId IN (80966, 9646, 61819, 80972, 161767, 64144, 8840, 28377, 132585, 37293, 223781, 8154, 7994, 9771, 158456, 43346, 30521, 30522, 72004, 9913, 6239, 9483, 7868, 9838, 286419, 9615, 9925, 7957, 1868482, 51154, 10141, 2715852, 9531, 84702, 106734, 34839, 60711, 9358, 8478, 7719, 51511, 7950, 56716, 93934, 10029, 8502, 8103, 244447, 28743, 630221, 7955, 9361, 9749, 299321, 13489, 10020, 7227, 9371, 8005, 7764, 83772, 9796, 9365, 27687, 8010, 9685, 59894, 8078, 8049, 9031, 69293, 48883, 1825980, 9595, 8153, 10181, 109280, 9606, 62062, 7998, 43179, 51337, 37003, 56723, 215358, 8187, 8630, 7897, 7918, 445787, 9785, 9541, 9544, 9545, 9568, 9994, 205130, 106582, 9103, 10036, 30608, 79684, 13616, 40151, 68415, 9669, 10089, 10090, 10093, 10103, 10096, 59463, 586833, 35670, 1026970, 32507, 452646, 61853, 9315, 8663, 105023, 9978, 10160, 8019, 8022, 74940, 8128, 9258, 9986, 123683, 8090, 30732, 183150, 30611, 9940, 9689, 9691, 74533, 9597, 9598, 9555, 1676925, 9157, 13735, 230844, 7757, 38626, 42100, 9755, 64176, 48698, 48699, 8081, 9601, 9813, 1328070, 379532, 8673, 132908, 303518, 42514, 10116, 59479, 61621, 61622, 559292, 39432, 8030, 8032, 96440, 283035, 9305, 55149, 113540, 52904, 9135, 41447, 1841481, 75366, 42254, 8175, 8508, 144197, 2489341, 441894, 9823, 59729, 31033, 2587831, 99883, 37347, 9739, 9999, 9643, 29073, 30538, 29139, 9627, 8364, 8083) GROUP BY taxId ORDER BY c DESC;
# taxId,COUNT(DISTINCT entryId)
# 9606,186018	9598,128104	9823,110422	8032,103657	9601,96866	8090,96407	75366,93434	62062,84554	10090,79300	9544,76936	8030,76473	9541,75291	7957,73824	9940,70004	8128,69023	8010,64861	8175,63049	30522,62620	7868,61045	7955,56996	59479,56314	64144,56046	8364,53884	205130,53592	10029,53548	30732,52770	8187,52579	586833,51872	9483,49545	9646,48862	31033,48337	8022,46649	99883,46098	9796,45584	52904,45450	8005,45121	9615,44473	29073,43956	9545,43590	9913,43333	9555,42494	9531,42255	9595,41892	61622,41281	9597,40347	61621,40069	9986,39659	9739,39580	9691,39286	37293,38925	7227,38754	9685,38237	9568,38080	7998,37812	61853,37706	2715852,37422	8154,37103	7994,36933	1676925,36818	9755,36574	39432,36536	9669,36388	56723,35849	72004,35598	106582,35168	9749,35027	158456,34817	10116,34524	42514,34371	13616,33620	215358,33581	230844,33108	9305,32984	64176,32819	8078,32761	10089,32725	38626,32644	8083,32557	7950,32513	9925,32407	105023,32375	8153,31921	30538,31735	9031,31728	161767,31725	132908,31170	10160,31093	1841481,30789	80966,30447	41447,30326	32507,30201	379532,30087	303518,30079	244447,29963	28743,29824	28377,29808	80972,29677	59729,29553	9643,29490	9258,29478	9838,29374	48699,29324	144197,29195	8081,28899	56716,28860	2587831,28638	10036,28613	29139,28420	48698,28332	37003,28013	9627,27960	69293,27628	1868482,27354	6239,26483	109280,26221	59894,25668	43346,25568	9365,25204	10181,25131	8673,24880	9103,24834	61819,24574	8840,24499	8663,24267	10020,24052	8502,23947	9785,23938	10141,23840	43179,23424	223781,22187	7897,22128	2489341,21614	113540,21517	7918,20509	123683,20144	59463,19898	40151,19494	51511,19130	13735,18929	30611,18765	7719,18451	60711,18248	452646,16063	441894,13254	7757,13113	48883,7913	559292,6719	13489,2343	8049,1575	9157,1084	30521,732	84702,726	30608,693	93934,629	9315,603	1026970,580	7764,555	42254,520	9361,485	74533,433	74940,428	10096,401	37347,281	8019,271	8508,234	9689,198	9978,178	79684,157	68415,156	9999,140	51337,133	9135,130	55149,114	10103,109	9813,109	9771,99	35670,94	283035,94	9358,90	10093,86	34839,84	9371,82	183150,78	8103,73	8630,68	27687,55	1328070,52	630221,40	106734,32	8478,29	42100,28	286419,28	299321,27	445787,17	96440,16	9994,16	1825980,8	83772,6	132585,6	51154,1
# alphafrag query: SELECT tax, COUNT(DISTINCT acc) AS c FROM alphafrag WHERE tax IN (9606, 9598, 9823, 8032, 9601, 8090, 75366, 62062, 10090, 9544, 8030, 9541, 7957, 9940, 8128, 8010, 8175, 30522, 7868, 7955, 59479, 64144, 8364, 205130, 10029, 30732, 8187, 586833, 9483, 9646, 31033, 8022, 99883, 9796, 52904, 8005, 9615, 29073, 9545, 9913, 9555, 9531, 9595, 61622, 9597, 61621, 9986, 9739, 9691, 37293, 7227, 9685, 9568, 7998, 61853, 2715852, 8154, 7994, 1676925, 9755, 39432, 9669, 56723, 72004, 106582, 9749, 158456, 10116, 42514, 13616, 215358, 230844, 9305, 64176, 8078, 10089, 38626, 8083, 7950, 9925, 105023, 8153, 30538, 9031, 161767, 132908, 10160, 1841481, 80966, 41447, 32507, 379532, 303518, 244447, 28743, 28377, 80972, 59729, 9643, 9258, 9838, 48699, 144197, 8081, 56716, 2587831, 10036, 29139, 48698, 37003, 9627, 69293, 1868482, 6239, 109280, 59894, 43346, 9365, 10181, 8673, 9103, 61819, 8840, 8663, 10020, 8502, 9785, 10141, 43179, 223781, 7897, 2489341, 113540, 7918, 123683, 59463, 40151, 51511, 13735, 30611, 7719, 60711, 452646, 441894, 7757, 48883, 559292, 13489, 8049, 9157, 30521, 84702, 30608, 93934, 9315, 1026970, 7764, 42254, 9361, 74533, 74940, 10096, 37347, 8019, 8508, 9689, 9978, 79684, 68415, 9999, 51337, 9135, 55149, 10103, 9813, 9771, 35670, 283035, 9358, 10093, 34839, 9371, 183150, 8103, 8630, 27687, 1328070, 630221, 106734, 8478, 42100, 286419, 299321, 445787, 96440, 9994, 1825980, 83772, 132585, 51154) GROUP BY tax ORDER BY c DESC;
# >> These numbers match what's in alphafrag with the exception of human, where I have 208 more (the fragmented proteins).
# Not sure why e.g. taxon 96440 has so few structures (16, also confirmed in sequences.fasta): it has 29,552 in UniProt: https://www.uniprot.org/uniprotkb?query=(taxonomy_id:96440)
# >> Looking at "Date Created", these are from 2022-01-19, i.e. too new for AlphaFold DB (even v4). AFDB v4's newest "sequenceVersionDate" is 2021-09-29 (see BigQuery above).
# >> 16 do indeed have an older "Date Created": https://www.uniprot.org/uniprotkb?query=(taxonomy_id:96440)%20AND%20(date_created:[*%20TO%202022-01-18])
# >> So it's not an issue with the taxon ID.
# AFDB has this taxon (96440) as "9SAUR" (BigQuery SELECT * FROM `bigquery-public-data.deepmind_alphafold.metadata` WHERE taxId=96440;), but the current UniProt has it as SALMN (see above), but that's okay.

# 55 species have fewer than 10,000 structures:
# Slow:
# input/gcs >> 1 *-0_v4.tar | xa sh -c "echo -n '{}: '; tar -tf {} | grep -iP '\-model_v4.cif.gz$' | wc -l"
# Fast (BigQuery, see above):
# 48883,7913	559292,6719	13489,2343	8049,1575	9157,1084	30521,732	84702,726	30608,693	93934,629	9315,603	1026970,580	7764,555	42254,520	9361,485	74533,433	74940,428	10096,401	37347,281	8019,271	8508,234	9689,198	9978,178	79684,157	68415,156	9999,140	51337,133	9135,130	55149,114	10103,109	9813,109	9771,99	35670,94	283035,94	9358,90	10093,86	34839,84	9371,82	183150,78	8103,73	8630,68	27687,55	1328070,52	630221,40	106734,32	8478,29	42100,28	286419,28	299321,27	445787,17	96440,16	9994,16	1825980,8	83772,6	132585,6	51154,1

# >> There is nothing I can do about this, and we don't have the resources to predict complete proteomes for these new species.
# >> Simply continue with what AlphaFold DB v4 currently has.



# Are there any newer files than 2022-11-01 (I think that's the release date of v4)?
# l gsutil -m ls -ldh 'gs://public-datasets-deepmind-alphafold-v4/AF*'
# g -v '2022-11-01' log-output-update_alphafold_gsutil__m_ls__ldh_gs___public_datasets_deepmind_alphafold_v4_AF_.txt
# >> None.
# >> No, no newer files, and neither are there any older files. All files are from 2022-11-01.

# Shortcuts used: 1: ls -1U, g: grep -iP, suq: sort | uniq -c | sort -g

# Extracting a specific UniProt accession from an archive: e.g.
# tar -xvf proteome-tax_id-7227-3_v4.tar --wildcards "*Q9W138*"

Run("Show directory", "ls -ltrh input/gcs")
Run("Count files", "ls -1 input/gcs/*.tar | wc -l")
Run("Count unique taxon IDs", "ls -1 input/gcs/*.tar | perl -npe 's/^input\/gcs\/proteomes\/proteome-tax_id-(\d+)-\d+_v\d+\.tar$/$1/' | sort | uniq | wc -l")

Done()
