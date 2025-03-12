#!/usr/bin/env python3
"""
alphafrag.py:
Parse AlphaFold fragment sequences into SQL table 'alphafrag'.
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
# from Bio import SeqIO
from blang_mysql import *
from blang import *

# SQL table with fragment protein sequences (>2700 aa proteins get split into 1400 aa fragments with a step size of 200 in AlphaFold DB, for human only - other species don't have results for >2700 aa proteins)
alphafrag = "alphafrag"
# SQL table with UniProt protein annotation and sequences
# alphauniprot = "alphauniprot_backup_2024_05"
alphauniprot = "alphauniprot"

# Paths
alphasyncpath = "input/alphasync"
ftppath = "input/ftp"
gcspath = "input/gcs"

# Paths for deleting outdated AlphaSync structure predictions (where an accession has since been made obsolete, or a sequence was changed, in UniProt)
cifdir = "input/alphasync/cif"
paedir = "input/alphasync/pae"
paramdir = "input/alphasync/params"

Args(0, f"\n -alphasync: Syncing: only parse updated fragments from AlphaSync TAR file in '{alphasyncpath}'\n -humanonly: Parse only human TAR file\n -debug: Don't actually make any changes, just simulate", " -humanonly")
# if not Switch('debug'):
#     Clear(alphafrag)




# Start

print(f"\nGetting TAR archives from '{alphasyncpath}', '{ftppath}' and '{gcspath}' and parsing fragment sequences into table '{alphafrag}':")

# Get list of TAR files to parse
if Switch('alphasync'):
    # AlphaSync updated structures only
    # infiles = nsort(Return(f"ls -1 {ftppath}/*_HUMAN_*.tar {gcspath}/*-9606-*.tar").split("\n"))
    # infiles = nsort(glob(f"{ftppath}/*_HUMAN_*.tar") + glob(f"{gcspath}/*-9606-*.tar"))
    # Get latest AlphaSync archive only (tail -n 1)
    infiles = nsort(Return(f"ls -1 {alphasyncpath}/alphasync_cif_*.tar | tail -n 1").split("\n"))
elif Switch('humanonly'):
    # Human only (NCBI taxon ID 9606)
    # infiles = nsort(Return(f"ls -1 {ftppath}/*_HUMAN_*.tar {gcspath}/*-9606-*.tar").split("\n"))
    # infiles = nsort(glob(f"{ftppath}/*_HUMAN_*.tar") + glob(f"{gcspath}/*-9606-*.tar"))
    infiles = nsort(Return(f"ls -1 {alphasyncpath}/*.tar {ftppath}/*_HUMAN_*.tar {gcspath}/*-9606-*.tar").split("\n"))
else:
    # All TAR files
    # infiles = nsort(Return(f"ls -1 {ftppath}/*.tar {gcspath}/*.tar").split("\n"))
    # infiles = nsort(glob(f"{ftppath}/*.tar") + glob(f"{gcspath}/*.tar"))
    infiles = nsort(Return(f"ls -1 {alphasyncpath}/*.tar {ftppath}/*.tar {gcspath}/*.tar").split("\n"))
    # State("Input files before moving human to the front:")
    # for file in infiles:
    #     print(file)
    # Move AlphaSync and human to the front of the list
    # for file in nsort(infiles):
    for file in reversed(nsort(infiles)):
        if rx("alphasync", file) or rx("_HUMAN_", file) or rx("-9606-", file):
            infiles.insert(0, infiles.pop(infiles.index(file)))
    # State("Input files after moving human to the front:")
    # for file in infiles:
    #     print(file)
    # print()

# If syncing (switch -alphasync active): first, clear AlphaSync fragments from table 'alphafrag'
if Switch('alphasync'):
    print(f"\nClearing AlphaSync fragments from table '{alphafrag}'...")
    if not Switch('debug'):
        query = Query(f"DELETE FROM {alphafrag} WHERE afdb=0")
        print(f" >> {Numrows(query):,} rows affected")
    else:
        query = Query(f"SELECT * FROM {alphafrag} WHERE afdb=0")
        print(f" >> {Numrows(query):,} rows affected")

# Main loop: Parse TAR files
Starttime()
for infile in infiles:

    afdb = 1
    if infile.startswith(alphasyncpath):
        # This structure isn't in the AlphaFold Protein Structure Database (it was predicted by AlphaSync instead!)
        afdb = 0

    if afdb == 1:
        print(f" >> {infile}")
    else:
        print(f" >> AlphaSync >> {infile}")

    # Format source file name (remove .tar, e.g. UP000000589_10090_MOUSE_v2.tar to UP000000589_10090_MOUSE_v2)
    source = re.sub(r"\.tar$", "", Basename(infile))
    
    # Open TAR file
    tar = tarfile.open(infile, mode="r:")
    acc = None
    prevacc = None
    maxfrag = 0
    for member in tq(tar, total=int(Return(f"tar -tf {infile} | wc -l"))):
        # print(f"   >> {member.name}")
        
        # Skip .pdb (PDB) files
        if rx(r"\.pdb\.gz$", member.name):
            Log("skipped pdb file", member.name)
            continue

        if rx(r"\.json\.gz$", member.name):
            Log("skipped json file", member.name)
            continue

        # Get only .cif.gz files
        if not rx(r"\.cif\.gz$", member.name):
            # Non-mmCIF file
            # print(f"     >> non-cif, skipped")
            # Log("skipped non-cif file", member.name)
            Die(f"Error: Unexpected non-PDB, non-mmCIF file: {member.name}")
            # continue
            
        # mmCIF file
        # print(f"     >> OK")
        
        # if not rx("Q9Y4Z8", member.name):   # random example, 1 fragment
        #     continue
        # if not rx("A0A024RBG1", member.name):   # random example, 1 fragment
        #     continue
        # if not rx("A0A087WUL8", member.name):   # random example, 14 fragments
        #     continue
        # if not rx("A0A075B6H5", member.name):     # DSSP error
        #     continue
        # if not rx("-(A0A024RBG1|A0A087WUL8|A0A075B6H5)-", member.name):     # Skip everything but this test set
        #     continue
        
        # Get accession (for handling fragments)
        # AF-A0A009IHW8-F1-model_v2.cif.gz
        m = rx(r"^AF-(\w+(-\d+)?)-F(\d+)-model_v\d+\.cif\.gz$", member.name)
        if m:
            acc = m[0]

            frag = int(m[2])
            if frag > maxfrag:
                maxfrag = frag
            # if frag > 1:
            #     Die(f"Fragment is '{frag}' for acc '{acc}'")
        else:
            Die(f"Couldn't parse '{member.name}'")
        
        # Extract .cif.gz to .cif using gzip (streaming, no temporary files)
        gz = tar.extractfile(member)
        # Open extracted .cif in text mode (streaming, no temporary files)
        cif = gzip.open(gz, mode="rt")

        # # Write .cif to temporary file for DSSP
        # os.mkdir(f"tmp/{acc}")
        # Run("Make temporary directory for DSSP for this acc", f"mkdir -p tmp/{acc}", silent=True)
        ciffile = f"tmp/{member.name}"
        # Remove .gz from file name
        ciffile = re.sub(r"\.gz$", "", ciffile)
        # Store list of mmCIF files and DSSP files, so they can be deleted once all fragments are present
        
        # # Write to output file (diagnostic)
        # with open(ciffile, "w") as cifout:
        #     print(cif.read(), file=cifout)
        # # Close
        # cif.close()
        # gz.close()
        # # Reopen
        # # Extract .cif.gz to .cif using gzip (streaming, no temporary files)
        # gz = tar.extractfile(member)
        # # Open extracted .cif in text mode (streaming, no temporary files)
        # cif = gzip.open(gz, mode="rt")

        # Parse mmCIF file (custom parsing, since Bio.SeqIO would require temporary files and be much slower)
        tmpacc = None
        name = None
        species = None
        tax = None
        fragstart = None
        fragstop = None
        tmplen = None
        seq = None
        for line in cif:
            # Parse line

            # Parse annotation
            if rx(r"^_ma_target_ref_db_details\.", line):
                # _ma_target_ref_db_details.db_accession                 A0A087WUL8
                # _ma_target_ref_db_details.db_code                      NBPFJ_HUMAN
                # _ma_target_ref_db_details.db_name                      UNP
                # _ma_target_ref_db_details.gene_name                    NBPF19
                # _ma_target_ref_db_details.ncbi_taxonomy_id             9606
                # _ma_target_ref_db_details.organism_scientific          "Homo sapiens"
                # _ma_target_ref_db_details.seq_db_align_begin           201
                # _ma_target_ref_db_details.seq_db_align_end             1600
                # _ma_target_ref_db_details.seq_db_isoform               ?
                # _ma_target_ref_db_details.seq_db_sequence_checksum     24A59AA23097CB90
                # _ma_target_ref_db_details.seq_db_sequence_version_date 2014-10-29
                # _ma_target_ref_db_details.target_entity_id             1
                
                # Accession
                m = rx(r"^_ma_target_ref_db_details\.db_accession +(\w+)", line)
                if m:
                    tmpacc = m[0]
                    # Verify that this is the correct UniProt accession (as expected from the file name)
                    if tmpacc != acc:
                        Die(f"Error: Expected UniProt accession '{acc}', but found '{tmpacc}' in '{ciffile}'")
                else:
                # Name (UniProt ID)
                    m = rx(r"^_ma_target_ref_db_details.db_code +(\w+)", line)
                    if m:
                        name = m[0]
                        # Parse species from e.g. NUD4B_HUMAN
                        m = rx(r"[A-Z0-9]+_([A-Z0-9]+)", name)
                        if m:
                            species = m[0]
                        else:
                            Die(f"Error: Couldn't parse species from UniProt ID '{name}'")
                    else:
                # Taxonomy (NCBI ID)
                        m = rx(r"^_ma_target_ref_db_details.ncbi_taxonomy_id +(\d+)", line)
                        if m:
                            tax = m[0]
                        else:
                # Start coordinate of this fragment
                            m = rx(r"^_ma_target_ref_db_details.seq_db_align_begin +(\d+)", line)
                            if m:
                                fragstart = m[0]
                # Stop coordinate of this fragment
                            m = rx(r"^_ma_target_ref_db_details.seq_db_align_end +(\d+)", line)
                            if m:
                                fragstop = m[0]
            
            # Parse sequence from _entity_poly.pdbx_seq_one_letter_code_can lines
            # Documentation:
            # https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v40.dic/Items/_entity_poly.pdbx_seq_one_letter_code_can.html (uses canonical residues as far as possible) (doesn't make any difference for AlphaFold, but should be the best choice for PDB structures as well)
            # https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_entity_poly.pdbx_seq_one_letter_code.html
            # https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v40.dic/Items/_struct_ref.pdbx_seq_one_letter_code.html
            # _entity_poly.pdbx_seq_one_letter_code_can 
            # ;EDSLEECAITYSNSHGPYDSNQPHRKTKITFEEDKVDSTLIGSSSHVEWEDAVHIIPENESDDEEEEEKGPVSPRNLQES
            # EEEEVPQESWDEGYSTLSIPPEMLASYQSYSSTFHSLEEQQVCMAVDIGRHRWDQVKKEDQEATGPRLSRELLDEKGPEV
            # LQDSLDRCYSTPSGCLELTDSCQPYRSAFYVLEQQRVGLAIDMDEIEKYQEVEEDQDPSCPRLSRELLDEKEPEVLQDSL
            # DRCYSIPSGYLELPDLGQPYSSAVYSLEEQYLGLALDVDRIKKDQEEEEDQDPPCPRLSRELVEVVEPEVLQDSLDRCYS
            # TPSSCLEQPDSCQPYGSSFYALEEKHVGFSLDVGEIEKKGKGKKRRGRRSKKERRRGRKEGEEDQNPPCPRLSRELLDEK
            # GPEVLQDSLDRCYSTPSGCLELTDSCQPYRSAFYILEQQCVGLAVDMDEIEKYQEVEEDQDPSCPRLSRELLDEKEPEVL
            # QDSLDRCYSIPSGYLELPDLGQPYSSAVYSLEEQYLGLALDVDRIKKDQEEEEDQDPPCPRLSRELVEVVEPEVLQDSLD
            # RCYSTPSSCLEQPDSCQPYGSSFYALEEKHVGFSLDVGEIEKKGKGKKRRGRRSKKERRRGRKEGEEDQNPPCPRLSREL
            # LDEKGPEVLQDSLDRCYSTPSGCLELTDSCQPYRSAFYILEQQCVGLAIDMDEIEKYQEVEEDQDPSCPRLSRELLDEKE
            # PEVLQDSLDRCYSIPSGYLELPDLGQPYSSAVYSLEEQYLGLALDVDRIKKDQEEEEDQDPPCPRLSRELVEVVEPEVLQ
            # DSLDRCYSTPSSCLEQPDSCQPYGSSFYALEEKHVGFSLDVGEIEKKGKGKKRRGRRSKKERRRGRKEGEEDQNPPCPRL
            # SRELLDEKGPEVLQDSLDRCYSTPSGCLELTDSCQPYRSAFYILEQQCVGLAVDMDEIEKYQEVEEDQDPSCPRLSRELL
            # DEKEPEVLQDSLDRCYSIPSGYLELPDLGQPYSSAVYSLEEQYLGLALDVDRIKKDQEEEEDQDPPCPRLSRELVEVVEP
            # EVLQDSLDRCYSTPSSCLEQPDSCQPYGSSFYALEEKHVGFSLDVGEIEKKGKGKKRRGRRSKKERRRGRKEGEEDQNPP
            # CPRLSRELLDEKGPEVLQDSLDRCYSTPSGCLELTDSCQPYRSAFYILEQQCVGLAVDMDEIEKYQEVEEDQDPSCPRLS
            # RELLDEKEPEVLQDSLDRCYSIPSGYLELPDLGQPYSSAVYSLEEQYLGLALDVDRIKKDQEEEEDQDPPCPRLSRELVE
            # VVEPEVLQDSLDRCYSTPSSCLEQPDSCQPYGSSFYALEEKHVGFSLDVGEIEKKGKGKKRRGRRSKKERRRGRKEGEED
            # QNPPCPRLSRELLDEKGPEVLQDSLDRCYSTPSGCLELTD
            # ;
            elif rx(r"^_entity_poly.pdbx_seq_one_letter_code_can +", line):

                m = rx(r"^_entity_poly.pdbx_seq_one_letter_code_can +(\w+)$", line)
                if m:
                    # Single-line sequence:
                    seq = m[0]
                else:
                    # Multi-line sequence:
                    seq = ""

                    # Parse entire sequence block
                    for line in cif:

                        # Break on closing semicolon
                        if rx(r"^;$", line):
                            break

                        # Ignore initial semicolon
                        m = rx(r"^;?([A-Z]+)$", line)
                        if m:
                            # Build complete sequence
                            seq += m[0]
                        else:
                            Die(f"Error: Couldn't parse '_entity_poly.pdbx_seq_one_letter_code_can' line:\n\n'{line}'\n\n")

            # AlphaSync structure (re-predicted using AlphaFold 2.3.2):
            # Parse sequence from _entity_poly_seq.mon_id lines instead (three-letter code, can convert using ThreeToOne())
            # Always one line per residue, terminated by #
            # _entity_poly_seq.mon_id
            # 0 1   MET 
            # 0 2   GLY 
            # 0 3   ARG 
            # 0 4   VAL 
            # 0 5   ARG 
            # ...
            # 0 270 GLN 
            # 0 271 GLU 
            # 0 272 GLN 
            # #
            elif afdb == 0 and rx(r"^_entity_poly_seq.mon_id$", line):
                seq = ""

                # Parse entire sequence block
                for line in cif:

                    # Break on #
                    if rx(r"^#$", line):
                        break

                    # Ignore initial semicolon
                    m = rx(r"^0 \d+ +([A-Z]{3}) +$", line)
                    if m:
                        # Build complete sequence
                        seq += ThreeToOne(m[0])
                    else:
                        Die(f"Error: Couldn't parse '_entity_poly_seq.mon_id' line:\n\n'{line}'\n\n")
                

        Log("parsed cif file", ciffile)

        if afdb == 0:
            # For AlphaSync, acc, name, species and tax will not be present in the .cif file. Need to retrieve these from alphauniprot.
            query = Query(f"SELECT acc, name, species, tax, seq FROM {alphauniprot} WHERE acc='{acc}'")
            if Numrows(query) == 0:
                # Log(f"alphasync acc from alphasync .cif file no longer found in table 'alphauniprot' for acc (acc must have been made obsolete) (skipped)", acc)
                Run(f"Delete all output files for obsolete acc '{acc}' (if they exist)", f"rm -fv {cifdir}/AF-{acc}-F*-model_v0.cif.gz {paedir}/AF-{acc}-F*-predicted_aligned_error_v0.json.gz {paramdir}/AF-{acc}-F*-alphafold_params.json")
                Log(f"alphasync acc from alphasync .cif file no longer found in table 'alphauniprot' for acc (acc must have been made obsolete) (loose /cif, /pae and /params files deleted) (skipped)", acc)
                continue
            tmpacc, name, species, tax, tmpseq = FetchOne(query)
            # Calculate fragment start and stop coordinates from fragment number with step size 200
            fragstart = 1 + 200 * (frag - 1)
            fragstop = 1400 + 200 * (frag - 1)
            if frag == 1 and len(seq) != 1400:
                fragstop = len(seq)
            elif frag > 1 and len(seq) != 1400:
                fragstop = fragstart + len(seq) - 1
            tmpseq = tmpseq[fragstart-1:fragstop]

            # Replace non-standard amino acids (B/Z/U/X) in sequence (for compatibility with AlphaFold)
            tmpseq = ReplaceNonstandardAAs(tmpseq)
            
            if seq != tmpseq:
                # Die(f"Error: Sequence mismatch for acc '{acc}' between table '{alphauniprot}' and AlphaSync .cif file '{ciffile}':\n\n{tmpseq}\n\n{seq}\n\n")
                # Log(f"sequence mismatch between alphasync .cif file and table 'alphauniprot' for acc (sequence must have been updated) (skipped)", acc)
                # Delete all files for this acc (if they exist): e.g. for acc='
                # {cifdir}/AF-Q96EY7-2-F*-model_v0.cif.gz
                # {paedir}/AF-Q96EY7-2-F*-predicted_aligned_error_v0.json.gz
                # {paramdir}/AF-Q96EY7-2-F*-alphafold_params.json
                Run(f"Delete all output files for acc with updated sequence '{acc}' (if they exist)", f"rm -fv {cifdir}/AF-{acc}-F*-model_v0.cif.gz {paedir}/AF-{acc}-F*-predicted_aligned_error_v0.json.gz {paramdir}/AF-{acc}-F*-alphafold_params.json")
                Log(f"sequence mismatch between alphasync .cif file and table 'alphauniprot' for acc (sequence must have been updated) (loose /cif, /pae and /params files deleted) (skipped)", acc)
                continue

        # Verify that all expected fields got parsed
        if tmpacc == None:
            Die(f"Error: Couldn't parse UniProt accession (acc) in mmCIF file '{ciffile}'")
        if name == None:
            Die(f"Error: Couldn't parse UniProt ID (name) in mmCIF file '{ciffile}'")
        if species == None:
            Die(f"Error: Couldn't parse UniProt species mnemonic (species) in mmCIF file '{ciffile}'")
        if tax == None:
            Die(f"Error: Couldn't parse NCBI taxon ID (tax) in mmCIF file '{ciffile}'")
        if fragstart == None:
            Die(f"Error: Couldn't parse fragment start coordinate (start) in mmCIF file '{ciffile}'")
        if fragstop == None:
            Die(f"Error: Couldn't parse fragment stop coordinate (stop) in mmCIF file '{ciffile}'")
        if seq == None:
            Die(f"Error: Couldn't parse _entity_poly.pdbx_seq_one_letter_code_can sequence (seq) in mmCIF file '{ciffile}'")

        # Verify that sequence is an AA sequence
        if not Aa(seq):
            Die(f"Error: Sequence '{seq}' contains non-AA characters in mmCIF file '{ciffile}'")
        
        # Verify sequence length
        tmplen = fragstop - fragstart + 1
        if len(seq) != tmplen:
            Die(f"Error: Expected sequence of length '{tmplen}', but got '{len(seq)}' aa in mmCIF file '{ciffile}'")
            
        # Insert fragment sequences into fragment SQL table
        q = f"INSERT INTO {alphafrag} SET acc='{acc}', name='{name}', species='{species}', tax='{tax}', frag='{frag}', fragstart='{fragstart}', fragstop='{fragstop}', source='{source}', afdb={afdb}, seq='{seq}'"
        if not Switch('debug'):
            Query(q)
        else:
            print(f"\n{q}")

        Log(f"successfully inserted into table '{alphafrag}' for acc", acc)
        Log(f"successfully inserted into table '{alphafrag}' for acc|frag", f"{acc}|{frag}")
        Log(f"successfully inserted into table '{alphafrag}' for acc|fragstart|fragstop", f"{acc}|{fragstart}|{fragstop}")
        Log(f"successfully inserted into table '{alphafrag}' for source|acc", f"{source}|{acc}")
        Log(f"successfully inserted into table '{alphafrag}' for name", name)
        Log(f"successfully inserted into table '{alphafrag}' for species", species)
        Log(f"successfully inserted into table '{alphafrag}' for tax", tax)
        Log(f"successfully inserted into table '{alphafrag}' for species|tax", f"{species}|{tax}")
        Log(f"successfully inserted into table '{alphafrag}' for source", source)
        Log(f"successfully inserted into table '{alphafrag}' for seq", seq)
            

        
Show(lim=20)

if not Switch('debug'):

    # Optimize(alphafrag)
    
    # Verify that fragment numbers make sense (always complete from 1..n)
    q = f"SELECT acc, source, MAX(frag) AS maxfrag, COUNT(DISTINCT frag) AS tmpmaxfrag FROM {alphafrag} GROUP BY acc, source HAVING maxfrag!=tmpmaxfrag"
    query = Query(q)
    if Numrows(query) > 0:
        print("Error:")
        for acc, source, maxfrag, tmpmaxfrag in Fetch(query):
            print(f" >> acc '{acc}'\tsource '{source}'\tmaxfrag '{maxfrag}'\ttmpmaxfrag '{tmpmaxfrag}'")
        Die(f"Error: MAX(frag) and COUNT(DISTINCT frag) don't match for '{Numrows(query)}' accs in table '{alphafrag}':\n\n{q}\n\n")

Stoptime()
print("\nDone!")
