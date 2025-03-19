# alphasync

AlphaSync (https://alphasync.stjude.org) is an updated AlphaFold structure database synchronised with UniProt.

AlphaSync predicts new structures to stay up-to-date with the latest UniProt release, and it additionally enhances all structures with residue-level data such as solvent accessibility and atom-level non-covalent contacts. A preprint describing AlphaSync should appear on bioRxiv soon.

## Please note

This repository provides the structure prediction and processing pipeline behind AlphaSync. It is not yet optimised for local deployment. It currently requires a local SQL server to be set up and is optimised for an LSF HPC environment with a Singularity container of AlphaFold 2. We have tentative plans for a more portable Docker version in future.

## Initial setup
- Requires Lahuta for residue-residue contacts, which is not yet publicly released (but should be soon) (https://bisejdiu.github.io/lahuta)
- Install Python packages (see individual scripts for imports)
- Install DSSP to make sure mkdssp is available (https://github.com/PDB-REDO/dssp)
- Update blang_mysql.py with SQL connection details
- Create tables in sql/sql_create_statements.sql and import .sql files

## To update
- Run run.py
    - Downloads structures from AlphaFold Protein Structure Database (AFDB) via FTP and GCS
    - Parses structures
    - Calculates RSA/dihedrals/contacts
    - Maps sequences to structures
- Run run.py -alphasync
    - Refreshes protein sequence and proteome data from UniProt REST API and FTP
    - Submits AlphaFold structure prediction jobs as needed
    - Parses structures
    - Calculates RSA/dihedrals/contacts
    - Maps sequences to structures
    - Can then migrate alphasync_compact SQL tables to web server (code available on request)
    - Repeat for new UniProt releases

## Acknowledgements
The code in input/alphasync/alphafold_tools is modified slightly from https://github.com/google-deepmind/alphafold, licensed under the Apache 2.0 license. The main change is a split into CPU- and GPU-based steps for more efficient parallelisation.
