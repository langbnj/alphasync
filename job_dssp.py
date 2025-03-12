#!/usr/bin/env python3
"""
Job script (runs on a given protein accession): Run DSSP on individual fragment files to calculate residue-level accessible surface area, combine them using combine_fragments_dssp.py, parse its results into the 'alphaseq' and 'alphasa' MySQL tables, and then remove temporary tmp/{acc} directory once complete
"""

# Initialize
import pandas as pd
import gemmi
from blang_mysql import *
from blang import *

alphaseq = "alphaseq"       # SQL table with complete protein sequences
alphasa = "alphasa"         # SQL table with residue-level accessible surface area values

# relasa ≥ 0.55: Disordered
windowsize = 10
disthresh = 0.55
# windowsize = 5
# disthresh = 0.58

# Disorder prediction calibration results from Akdel et al. 2022 (https://doi.org/10.1038/s41594-022-00849-w)
# 
# Smoothing window size:
# https://github.com/normandavey/ProcessedAlphafold
# Window	AUC
# IUPRED	0.870
# 0			0.843
# 5			0.922	# For small IDRs, 5 is already good
# 10		0.931	# No real gain >10, but a potential loss of sensitivity to small disordered regions
# 15		0.933	# relASA ≥0.55 is calibrated for window size 15, not for other sizes.
# 20		0.934
# 25		0.933
# 30		0.933
# 
# From GitHub Summary Statistics table:
# https://github.com/normandavey/ProcessedAlphafold/blob/main/Summary_statistics_IUPred2_vs_AlphaFold2.xlsx
# 
# relASA-based:
#                                   AUC      TPR at FPR=0.05    score cutoff for FPR=0.05
# IUPred2                           0.870    0.586    0.56
# AlphaFold2 accessibility, w=0     0.843    0.370    0.86
# AlphaFold2 accessibility, w=5     0.922    0.761    0.58
# AlphaFold2 accessibility, w=10    0.931    0.799    0.55
# AlphaFold2 accessibility, w=15    0.933    0.803    0.55
# AlphaFold2 accessibility, w=20    0.934    0.800    0.56
# AlphaFold2 accessibility, w=25    0.933    0.797    0.56
# AlphaFold2 accessibility, w=30    0.933    0.797    0.56
# AlphaFold2 accessibility, w=35    0.932    0.790    0.56
# AlphaFold2 accessibility, w=40    0.930    0.787    0.56
# AlphaFold2 accessibility, w=45    0.929    0.780    0.56
# AlphaFold2 accessibility, w=50    0.927    0.777    0.56
# 
# pLDDT-based (doesn't come anywhere close to relASA's TPR of 0.8):
#                           AUC      TPR at FPR=0.05    score cutoff for FPR=0.05
# IUPred2                   0.870    0.586    0.56
# AlphaFold2 pLDDT, w=0     0.891    0.682    0.39
# AlphaFold2 pLDDT, w=5     0.896    0.712    0.35
# AlphaFold2 pLDDT, w=10    0.900    0.722    0.33
# AlphaFold2 pLDDT, w=15    0.902    0.726    0.33
# AlphaFold2 pLDDT, w=20    0.902    0.724    0.33
# AlphaFold2 pLDDT, w=25    0.902    0.720    0.33
# AlphaFold2 pLDDT, w=30    0.902    0.716    0.33
# AlphaFold2 pLDDT, w=35    0.900    0.705    0.33
# AlphaFold2 pLDDT, w=40    0.899    0.697    0.33
# AlphaFold2 pLDDT, w=45    0.897    0.694    0.33
# AlphaFold2 pLDDT, w=50    0.895    0.680    0.34

# Bludau et al. (https://doi.org/10.1371/journal.pbio.3001636): https://github.com/MannLabs/structuremap_analysis/blob/master/data/disorder_data/IDR_ROC_curve.csv
# 
# This paper actually used the same training data as in https://doi.org/10.1038/s41594-022-00849-w, hence it should have the same results for SASA/RSA.
# Fig. 1d and Supp Fig. 4B from Akdel et al. show SASA/RSA results. Window size 10 performs well. Essentially zero improvement with larger windows.
# Thresholds aren't given in this paper, but the AUC values in Supp Fig. 4B perfectly match those above, which include threshold values (0.55 for window size 10).
# >> Use threshold 0.55 and RSA smoothing window size 10.
# 
# Bludau et al. used default RSA reference values from BioPython 1.79 (https://biopython.org/docs/1.79/api/Bio.PDB.DSSP.html):
# Sander and Rost 1994 https://doi.org/10.1002/prot.340200303
# I'm using Tien et al. (2013), which avoids values >100% more successfully (and hence might underestimate RSA compared to Sander and Rost).
# 
# Filtering this for FPR < 0.05 and sorting by descending TPR:
# >> input/idrs/IDR_ROC_curve_filtered_for_FPR_below_0.05_sorted_by_TPR_desc.tsv
# The threshold for RSA with a window size of 10 would be 0.42. With smaller and larger windows, the thresholds are similar.
# This is fairly similar to the 0.55 threshold calibrated by Akdel et al. 2022 (https://doi.org/10.1038/s41594-022-00849-w), and an average 42% exposure conceptually appears quite low for an IDR.
# It's also surprising since both used the same training data. Bludau et al. report AUC=0.943 for RSA/SASA20, while Akdel et al. report AUC=0.934 for RSA/SASA20.
# Bludau et al. also report AUC=0.931 for RSA/SASA10 (very good, nearly identical to RSA/SASA20).
# Akdel et al. doesn't specify RSA reference values, but may not have used Sander and Rost's 1994 values.
# 
# >> Keeping 0.55 threshold from Akdel et al. 2022 (https://doi.org/10.1038/s41594-022-00849-w).



# Relative residue accessible surface area <= 25% means buried (Emmanuel Levy's definition)
buried_threshold = 0.25;

# (source, acc, maxfrag) = Args(3, "[Source] [UniProt accession] [Number of fragments]", "UP000005640_9606_HUMAN_v2 A0A087WUL8 14")
(acc, maxfrag) = Args(2, "[UniProt accession] [Number of fragments]\n\n -alphasync: Updating AlphaSync proteins (non-AFDB, i.e. afdb=0)", "A0A087WUL8 14")

# inpath = f"../{source}/{acc}"
inpath = f"../{acc}"
# inpath = "."
logpath = "tmp/_logs"

# Set afdb to 0 for AlphaSync proteins (not in AFDB)
afdb = 1
tmpalphasync = ""
if Switch('alphasync'):
    afdb = 0
    tmpalphasync = " -alphasync"



# Start

# Check if data already exists for this acc in table 'alphasa'
query = Query(f"SELECT id FROM {alphasa} WHERE acc='{acc}' AND afdb='{afdb}' LIMIT 1")
# ...and exit if yes
if Numrows(query) == 1:
    print(f"Data already existed for acc '{acc}' (afdb={afdb}) in table '{alphasa}', exiting (skip)!")
    sys.exit()
else:
    print(f"No data yet for acc '{acc}' (afdb={afdb}) in table '{alphasa}', starting!")



# Get additional information on this acc from table 'alphaseq'
query = Query(f"SELECT DISTINCT name, species, tax, frags, seq FROM {alphaseq} WHERE acc='{acc}' AND afdb='{afdb}'")
(name, species, tax, tmpmaxfrag, seq) = FetchOne(query)

if maxfrag != tmpmaxfrag:
    # Die(f"Error: Expected {maxfrag} fragments for source '{source}' acc '{acc}', but got {tmpmaxfrag}")
    Die(f"Error: Expected {maxfrag} fragments for acc '{acc}', but got {tmpmaxfrag}")



# Start

print(f"\nRunning DSSP on '{inpath}' (acc '{acc}', {maxfrag} fragments):")

# # Change directory
# Not necessary now since the job moves to this directory before running this script (to avoid delays from Python's module import due to the huge number of log files in tmp/_logs/)
# os.chdir(inpath)

# infiles = nsort(Return(f"ls -1U {inpath}").split("\n"))
infiles = nsort(Return(f"ls -1U *.cif").split("\n"))
if len(infiles) != maxfrag:
    Die(f"Error: Expected {maxfrag} fragments in '{inpath}', but found {len(infiles)}")

tmpfiles = set()
for ciffile in tq(infiles):
    # print(f" >> {ciffile}")

    # # Get only .cif files
    # if not rx(r"\.cif$", ciffile):
    #     # Non-mmCIF file
    #     Die("Error: Unexpected non-mmCIF file found: {inpath}/{ciffile}")
    #     continue

    Log("cif files processed", ciffile)
    # Don't delete CIF files here (they will still be used by lahuta.py, and ultimately removed by job.py)
    # if not Switch('debug'):
    #     # In debug mode, keep the mmCIF file (don't remove it after running)
    #     tmpfiles.add(ciffile)

    # Get fragment number for this mmCIF file
    m = rx(r"^AF-"+acc+r"-F(\d+)-model_v\d+\.cif$", ciffile)
    if m:
        frag = m[0]
    else:
        Die(f"Couldn't parse filename '{ciffile}'")

    dsspfile = re.sub(r"\.cif$", ".dssp", ciffile)
    if not Exists(dsspfile):

        # If this is an AlphaSync structure (AlphaFold 2.3.2): Use gemmi to convert mmCIF file to PDB.
        # Unlike AFDB structures, AF 2.3.2's mmCIF output files are missing sections such as _pdbx_poly_seq_scheme, _entity_poly_seq, and _struct_asym, which mkdssp needs.
        # I haven't found a method to add these sections automatically, so converting to PDB appears to be the simplest solution.
        if afdb == 0:
            # Read structure
            tmp_structure = gemmi.read_structure(ciffile)

            # Write structure
            pdbfile = re.sub(r"\.cif$", ".pdb", ciffile)
            tmp_structure.write_pdb(pdbfile)
            tmpfiles.add(pdbfile)

            infile = pdbfile
        else:
            infile = ciffile

        # Run("Run DSSP", f"mkdssp '{infile}' '{dsspfile}'", silent=False)
        
        # This would suppress all warnings
        # Run("Run DSSP", f"mkdssp '{infile}' '{dsspfile}' --quiet", silent=False)

        # Suppress only one specific AlphaFold-related warning about citations
        # Note that this requires running bash (sh's redirect syntax doesn't seem to support filtering only STDERR)
        Run("Run DSSP", f"bash -c \"mkdssp '{infile}' '{dsspfile}' 2> >(grep -viP '^(Links for citation_author:citation:1 are incomplete|  There are 33 items in citation_author that don.t have matching parent items in citation|Warning, the input file is not valid. Run with --verbose to see why\.)'>&2)\"", silent=True)

        Log(f"ran DSSP for acc|frag", f"{acc}|{frag}")
    else:
        Log(f"skipped running DSSP since output already existed for acc|frag", f"{acc}|{frag}")
    tmpfiles.add(dsspfile)



# Verify that all fragments have been run
if (frag != maxfrag):
    Die(f"Error: Expected to have processed {maxfrag} fragments, but only processed {frag}")



# Accession is complete (DSSP run for all fragments):
# Process DSSP output
outfile = re.sub(r"\.cif$", ".combined", ciffile)
print()
Run("Combine fragments", f"../../combine_fragments_dssp.py . {windowsize} {disthresh} {acc} {maxfrag} {outfile}{tmpalphasync}", silent=False)
tmpfiles.add(outfile)

# Parse combined fragment data
# with open(outfile) as f:
#     for line in f:
#         if rx(r"^#", line):
#             continue
# Read TSV
df = read_tsv(outfile)

# Verify sequence by comparing to table 'alphaseq'
# Get sequence (by concatenating the "residue1" column)
tmpseq = df["residue1"].str.cat()

# Verify that sequence is an AA sequence
if not Aa(tmpseq):
    Die(f"Error: Sequence '{tmpseq}' contains non-AA characters")

# Verify by comparing to table 'alphaseq'
if seq != tmpseq:
    Die(f"\n\nError: Expected sequence:\n\n{seq}\n\n...but got:\n\n{tmpseq}\n\n")

# Previously, this script here filled alphaseq
# # Insert sequence (combined across fragments) into alphaseq SQL table
# q = f"INSERT INTO {alphaseq} SET acc='{acc}', name='{name}', species='{species}', tax='{tax}', frags='{maxfrag}', afdb='{afdb}', seq='{seq}'"
# if not Switch('debug'):
#     Query(q)
# else:
#     State(q)



# Insert into alphasa
# Parse row-wise
q = f"INSERT INTO {alphasa} (acc, species, tax, frags, afdb, site, aa, plddt, plddt10, asa, asa10, relasa, relasa10, dis, dis10, surf, surf10, sec) VALUES "
for i, a in df.iterrows():
    # print(a)
    site = a["#position"]
    aa = a["residue1"]
    plddt = a["pLDDT"]
    plddt10 = a["pLDDT_smoothed_10"]
    asa = a["SASA_absolute"]
    asa10 = a["SASA_absolute_smoothed_10"]
    relasa = a["SASA_relative"]
    relasa10 = a["SASA_relative_smoothed_10"]
    sec = a["DSSP_sec_struct"]
    dis10 = a["disordered"]
    
    # Secondary structure codes
    # https://github.com/PDB-REDO/dssp
    # https://github.com/PDB-REDO/dssp/blob/trunk/doc/mkdssp.pdf
    # The DSSP algorithm assigns secondary structure based on the energy calculated for H-bonds.
    # Table 1. Secondary Structures recognized:
    # 
    # DSSP Code    mmCIF Code        Description
    # H            HELX_RH_AL_P    Alphahelix
    # B            STRN            Betabridge
    # E            STRN            Strand
    # G            HELX_RH_3T_P    Helix_3
    # I            HELX_RH_PI_P    Helix_5
    # P            HELX_LH_PP_P    Helix_PPII
    # T            TURN_TY1_P      Turn
    # S            BEND            Bend
    # ' ' (space)  OTHER         Loop
    # 
    # Alphabetic:
    # DSSP Code    mmCIF Code        Description
    # ' ' (space)  OTHER           Loop
    # B            STRN            Betabridge
    # E            STRN            Strand
    # G            HELX_RH_3T_P    Helix_3
    # H            HELX_RH_AL_P    Alphahelix
    # I            HELX_RH_PI_P    Helix_5
    # P            HELX_LH_PP_P    Helix_PPII
    # S            BEND            Bend
    # T            TURN_TY1_P        Turn
    # 
    # Change "C" back to DSSP's " " for "LOOP"
    
    if sec == "C":
        sec = " "
    
    # Assign raw (non-smoothed) disorder based on relasa ≥disthresh
    if relasa >= disthresh:
        dis = "*"
    else:
        dis= "."
    
    # Format disorder classification (*: disordered, .: structured, as in CASP)
    # (Note that this column is based on SASA_relative_smoothed_{windowsize}, i.e. with smoothing in a ±10 aa window, i.e. across 21 aa, being ≥{disthresh})
    if dis10 == 0:
        dis10 = "."
    elif dis10 == 1:
        dis10 = "*"
    else:
        Die("Error: 'disordered' (smooth10) is not 0 or 1: {dis10}")
    
    # Format surface/core definition (S: surface, C: core, core being ≤25% relASA)
    if relasa <= buried_threshold:
        surf = "C"
    else:
        surf = "S"
    
    # Format surface/core definition (S: surface, C: core, core being ≤25% relASA)
    if relasa10 <= buried_threshold:
        surf10 = "C"
    else:
        surf10 = "S"
    
    # Insert residue-level data (processed across fragments) into alphasa SQL table
    # q = f"INSERT INTO {alphasa} SET acc='{acc}', name='{name}', species='{species}', frags='{maxfrag}', source='{source}', afdb=1, site='{site}', aa='{aa}', plddt='{plddt}', plddt10='{plddt10}', asa='{asa}', asa10='{asa10}', relasa='{relasa}', relasa10='{relasa10}', sec='{sec}', dis='{dis}', dis10='{dis10}', surf='{surf}', surf10='{surf10}'"
    # q = f"INSERT INTO {alphasa} SET acc='{acc}', name='{name}', species='{species}', frags='{maxfrag}', afdb=1, site='{site}', aa='{aa}', plddt='{plddt}', plddt10='{plddt10}', asa='{asa}', asa10='{asa10}', relasa='{relasa}', relasa10='{relasa10}', dis='{dis}', dis10='{dis10}', surf='{surf}', surf10='{surf10}', sec='{sec}'"
    # q = f"INSERT INTO {alphasa} SET acc='{acc}', species='{species}', tax='{tax}', frags='{maxfrag}', afdb=1, site='{site}', aa='{aa}', plddt='{plddt}', plddt10='{plddt10}', asa='{asa}', asa10='{asa10}', relasa='{relasa}', relasa10='{relasa10}', dis='{dis}', dis10='{dis10}', surf='{surf}', surf10='{surf10}', sec='{sec}'"
    # q += f"INSERT INTO {alphasa} SET acc='{acc}', species='{species}', tax='{tax}', frags='{maxfrag}', afdb=1, site='{site}', aa='{aa}', plddt='{plddt}', plddt10='{plddt10}', asa='{asa}', asa10='{asa10}', relasa='{relasa}', relasa10='{relasa10}', dis='{dis}', dis10='{dis10}', surf='{surf}', surf10='{surf10}', sec='{sec}'; "
    q += f"ROW ('{acc}', '{species}', '{tax}', '{maxfrag}', '{afdb}', '{site}', '{aa}', '{plddt}', '{plddt10}', '{asa}', '{asa10}', '{relasa}', '{relasa10}', '{dis}', '{dis10}', '{surf}', '{surf10}', '{sec}'), "

# Insert all rows in a single multi-row insert query
q = re.sub(r", $", ";", q)
if not Switch('debug'):
    Query(q)
else:
    State(q)
            

# Delete temporary files
if not Switch('debug'):
    print()
    for tmpfile in nsort(tmpfiles):
        print(f"Removing temporary file '{tmpfile}'")
        os.remove(tmpfile)

#     # Move back to original location (tmp/_logs)
#     os.chdir(logpath)
# 
#     # Remove temporary directory for this accession (will throw an error if not empty)
#     os.rmdir(inpath)

# # Also try to remove the source's directory (e.g. tmp/UP000005640_9606_HUMAN_v2), if empty (it'll get recreated by main.py for the next DSSP job)
# try:
#     os.rmdir(f"../{source}")
# except:
#     pass

Show(lim=20)

print(f"Successfully inserted ASA values into table '{alphasa}'")

print("\nDone!")
