#!/usr/bin/env python3
"""
Job script (runs on a given protein accession):
- Get dihedral angles and proline isomerization states (cis/trans) from mmCIF file using BioPDB
- Update 'alphasa' MySQL table columns: 'iso', 'phi', 'psi', 'omega', 'chi1', 'chi2', 'chi3', 'chi4', 'chi5', 'tau'

"""

# Initialize
import pandas as pd
import numpy as np
import math
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.ic_rebuild import structure_rebuild_test
from blang_mysql import *
from blang import *
np.set_printoptions(suppress=True)  # Disable scientific format

alphafrag = "alphafrag"     # SQL table with fragment protein sequences (>2700 aa proteins get split into 1400 aa fragments with a step size of 200 in AlphaFold DB, for human only - other species don't have results for >2700 aa proteins)
alphaseq = "alphaseq"       # SQL table with complete protein sequences (used for retrieving additional protein information here)
alphasa = "alphasa"         # SQL table with residue-level accessible surface area values
tmptable = f"alphasa_tmp_update"    # Temporary SQL table for updating the 'alphasa' table efficiently

# Fragment length and step size used by DeepMind.
# Proteins longer than 2700 residues are split into windows of width 1400 with a step size of 200.
fraglen = 1400
fragstep = 200

# Proline isomerization states: omega angle thresholds
omega_cis_max = 50.0
omega_trans_min = 130.0

# Get arguments
(acc, maxfrag) = Args(2, "[UniProt accession] [Number of fragments]\n\n -alphasync: Updating AlphaSync proteins (non-AFDB, i.e. afdb=0)", "A0A087WUL8 14")

# inpath = f"../{source}/{acc}"
inpath = f"../{acc}"
# inpath = "."
logpath = "tmp/_logs"

# Set afdb to 0 for AlphaSync proteins (not in AFDB)
afdb = 1
if Switch('alphasync'):
    afdb = 0



# Create temporary table for updating the 'alphasa' table efficiently
# q = f"""DROP TABLE IF EXISTS {tmptable}"""
# Query(q)
# q = f"""CREATE TABLE {tmptable} (
q = f"""CREATE TEMPORARY TABLE {tmptable} (
`dihedrals` json DEFAULT NULL
) engine=InnoDB"""
# q = f"""CREATE TABLE {tmptable} (
#   `acc` char(10) DEFAULT NULL,
#   `site` mediumint DEFAULT NULL,
#   `iso` char(1) DEFAULT NULL,
#   `phi` float DEFAULT NULL,
#   `psi` float DEFAULT NULL,
#   `omega` float DEFAULT NULL,
#   `chi1` float DEFAULT NULL,
#   `chi2` float DEFAULT NULL,
#   `chi3` float DEFAULT NULL,
#   `chi4` float DEFAULT NULL,
#   `chi5` float DEFAULT NULL,
#   `tau` float DEFAULT NULL,
#   PRIMARY KEY (`id`),
#   KEY `Acc` (`acc`),
#   KEY `Species` (`species`),
#   KEY `Tax` (`tax`),
#   KEY `Frags` (`frags`),
#   KEY `Afdb` (`afdb`),
#   KEY `Site` (`site`),
#   KEY `Aa` (`aa`),
#   KEY `Dis` (`dis`),
#   KEY `Dis10` (`dis10`),
#   KEY `Surf` (`surf`),
#   KEY `Surf10` (`surf10`),
#   KEY `Sec` (`sec`),
#   KEY `Membrane` (`membrane`),
#   KEY `Iso` (`iso`)
# ) ENGINE=InnoDB DEFAULT CHARSET=latin1"""
Query(q)





# Functions

# Average angles
def average_angles(angles):

    # print(f" >> ANGLES >> {angles}")
    # # Filter out None and np.nan values
    # valid_angles = [angle for angle in angles if angle is not None and not np.isnan(angle)]
    # 
    # # If no valid angles remain, return None or handle as appropriate
    # if not valid_angles:
    #     return None  # or return np.nan, depending on how you want to handle this case
    # 
    # # Convert degrees to radians
    # angles_rad = [math.radians(angle) for angle in valid_angles]
    # d()
    angles_rad = [math.radians(angle) for angle in angles]

    # Sum Cartesian coordinates
    X = sum(math.cos(angle) for angle in angles_rad)
    Y = sum(math.sin(angle) for angle in angles_rad)

    # Compute the average angle in radians
    avg_angle_rad = math.atan2(Y, X)

    # Convert back to degrees
    avg_angle_deg = math.degrees(avg_angle_rad)

    return avg_angle_deg


# Combine fragments
# Ignore values from dubious regions (artificial termini)
def CombineFragments(df):
    
    # Ignore values from dubious regions (artificial termini)

    # Remove N-terminal 200 aa at artificial N-termini (frag > 1 will have an artificial N-terminus)
    # Remove any residue in fragment 2 or above that is between 1-200
    df = df.loc[np.invert((df["frag"] > 1) & (df["site"] <= fragstep))]

    # Remove C-terminal 200 aa at artificial C-termini (frag < maxfrag will have an artificial C-terminus)
    # Remove any residue except in the last fragment that is between 1201-1400
    df = df.loc[np.invert((df["frag"] < maxfrag) & (df["site"] > fraglen - fragstep))]

    # Shift site according to fragment number (i.e. residue 1 in fragment 2 will become 1401)
    df["site"] += fragstep * (df["frag"] - 1)

    # Verify that all amino acid positions listed are the correct residue
    for i, row in df.iterrows():
        tmpsite = row["site"]
        tmpaa = row["aa"]
        # Expected amino acid
        expaa = seq[tmpsite-1:tmpsite]
        # Verify that the amino acid is correct (using the alphaseq sequence retrieved earlier)
        if tmpaa != expaa:
            Die(f"Error: Expected residue '{expaa}' at position '{tmpsite}' in acc '{acc}', but got '{tmpaa}'")
        # if tmpaa != 'P':
        #     Die(f"Error: Expected residue 'P' at position '{tmpsite}' in acc '{acc}', but got '{tmpaa}'")

    # # Select final columns in final order (leaving out "frag" since fragments have been combined already here)
    # df = df[["acc", "species", "tax", "frags", "afdb", "site1", "site2", "aa1", "aa2", "type", "atom1", "atom2", "dist"]]

    # Remove "frag" column
    # d()
    df = df.drop("frag", axis=1)

    # Get average dihedral angles (averaging across fragments)
    # df = df.groupby(["site", "aa"], as_index=False).mean(numeric_only=True)
    # Using average_angles function:
    # d()
    df = df.groupby(["site", "aa"], as_index=False).agg({
        "iso":   lambda x: None,
        "phi":   lambda x: average_angles(x),
        "psi":   lambda x: average_angles(x),
        "omega": lambda x: average_angles(x),
        "chi1":  lambda x: average_angles(x),
        "chi2":  lambda x: average_angles(x),
        "chi3":  lambda x: average_angles(x),
        "chi4":  lambda x: average_angles(x),
        "chi5":  lambda x: average_angles(x),
        "tau":   lambda x: average_angles(x)
    })

    # Update isomerization state based on average omega angle
    # if np.abs(angle) < 50.0:
    #     states.append('cis')
    # else:
    #     states.append('trans')
    # Update data frame
    # df.loc[df["omega"] < 50.0, "iso"] = 'cis'
    # df.loc[df["omega"] >= 50.0, "iso] = 'trans'
    # df.loc[(df["omega"] < 50.0) & (df["iso"] == "trans")]
    # df.loc[(df["omega"] >= 50.0) & (df["iso"] == "cis")]
    # Update data frame with a lambda function
    # df["iso"] = df["omega"].apply(lambda x: "cis" if x < 50.0 else "trans")
    # df["iso"] = df["omega"].apply(lambda x: "cis" if abs(x) <= 50.0 else "trans" if abs(x) >= 130.0 else None)
    # df["iso"] = df["omega"].apply(lambda x: "c" if abs(x) <= omega_cis_max else "t" if abs(x) >= omega_trans_min else None)
    df["iso"] = df["omega"].apply(lambda x: "c" if abs(x) <= omega_cis_max else "t" if abs(x) >= omega_trans_min else " ")
    # # With handling for omega == None, returning None (which should not happen, but this way it'll be traceable):
    # df["iso"] = df["omega"].apply(lambda x: None if np.isnan(x) else "c" if abs(x) <= omega_cis_max else "t" if abs(x) >= omega_trans_min else " ")
    # d()
    # df.insert(loc=2, column="iso", value=df["omega"].apply(lambda x: "cis" if abs(x) <= 50.0 else "trans" if abs(x) >= 130.0 else None))

    # Sort rows by site
    df = df.sort_values(["site"], ignore_index=True)

    # Return
    # d()
    return df
        



# Start

# Verify that residue data already exists for this acc in table 'alphasa'
query = Query(f"SELECT id FROM {alphasa} WHERE acc='{acc}' AND afdb='{afdb}' LIMIT 1")
# ...and exit if not
if Numrows(query) == 0:
    Die(f"Error: No data yet for acc '{acc}' (afdb={afdb}) in table '{alphasa}'")
    sys.exit()
else:
    print(f"Found data for acc '{acc}' (afdb={afdb}) in table '{alphasa}', starting!")

# Check if dihedral angles and proline isomerization states already exist for this acc in table 'alphasa'
query = Query(f"SELECT id FROM {alphasa} WHERE acc='{acc}' AND afdb='{afdb}' AND iso IS NOT NULL LIMIT 1")
# ...and exit if yes
if Numrows(query) == 1:
    print(f"Dihedral angles and proline isomerization data already existed in column 'iso' for acc '{acc}' (afdb={afdb}) in table '{alphasa}', exiting (skip)!")
    sys.exit()
else:
    print(f"No dihedral angle and proline isomerization data yet for acc '{acc}' (afdb={afdb}) in table '{alphasa}', starting!")



# Get additional information on this acc from table 'alphaseq'
query = Query(f"SELECT DISTINCT name, species, tax, frags, seq FROM {alphaseq} WHERE acc='{acc}' AND afdb='{afdb}'")
(name, species, tax, tmpmaxfrag, seq) = FetchOne(query)

if maxfrag != tmpmaxfrag:
    # Die("Error: Expected {maxfrag} fragments for source '{source}' acc '{acc}', but got {tmpmaxfrag}")
    Die(f"Error: Expected {maxfrag} fragments for acc '{acc}' (afdb={afdb}), but got {tmpmaxfrag}")




if maxfrag != tmpmaxfrag:
    # Die("Error: Expected {maxfrag} fragments for source '{source}' acc '{acc}', but got {tmpmaxfrag}")
    Die(f"Error: Expected {maxfrag} fragments for acc '{acc}', but got {tmpmaxfrag}")

# Get additional information on this UniProt accession (sequence)
query = Query(f"SELECT DISTINCT seq FROM {alphaseq} WHERE acc='{acc}' AND afdb='{afdb}'")
(seq) = FetchOne(query)

# Get expected fragment sequences from table 'alphafrag'
seqs = FetchMap(Query(f"SELECT frag, seq FROM {alphafrag} WHERE acc='{acc}' AND afdb='{afdb}'"))



# Start

print(f"\nRunning dihedral angle and proline isomerization state detection on '{inpath}' (acc '{acc}', {maxfrag} fragments, afdb={afdb}):")

# # Change directory
# Not necessary now since the job moves to this directory before running this script (to avoid delays from Python's module import due to the huge number of log files in tmp/_logs/)
# os.chdir(inpath)

# infiles = nsort(Return(f"ls -1U {inpath}").split("\n"))
infiles = nsort(Return(f"ls -1U *.cif").split("\n"))
if len(infiles) != maxfrag:
    Die(f"Error: Expected {maxfrag} fragments in '{inpath}', but found {len(infiles)}")

dihedrals = []
affected = 0
# tmpfiles = set()
for ciffile in tq(infiles):
    # print(f" >> {ciffile}")

    # # Get only .cif files
    # if not rx(r"\.cif$", ciffile):
    #     # Non-mmCIF file
    #     Die("Error: Unexpected non-mmCIF file found: {inpath}/{ciffile}")
    #     continue

    Log("cif files processed", ciffile)
    # Don't delete CIF files here (they will ultimately be removed by job.py)
    # if not Switch('debug'):
    #     # In debug mode, keep the mmCIF file (don't remove it after running)
    #     tmpfiles.add(ciffile)

    # Get fragment number for this mmCIF file
    m = rx(r"^AF-" + acc + r"-F(\d+)-model_v\d+\.cif$", ciffile)
    if m:
        frag = m[0]
    else:
        Die(f"Couldn't parse filename '{ciffile}'")
        
    # Parse structure from mmCIF file
    # parser = MMCIFParser()
    # myProtein = parser.get_structure(acc, ciffile)
    # myChain = myProtein[0]["A"]
    structure = MMCIFParser().get_structure(acc, ciffile)

    for model in structure:
        for chain in model:    

            rows = []

            # Calculate dihedral angles and cis/trans isomerization states
            # chain.atom_to_internal_coordinates()
            chain.atom_to_internal_coordinates(verbose=True)
            # verbose=True prints harmless warnings about backbone continuity, highlighting bonds > 1.4 A:
            # "chain break at GLU  1228  due to MaxPeptideBond (1.4 angstroms) exceeded"
            # These are fairly rare, usually no more than one bond per 1400 aa protein.
            # >> Even if the peptide bond distance is slightly higher, the angles shouldn't be fully invalidated. Ignoring these.

            # Test whether the structure can be rebuilt from the internal coordinates using a BioPDB function
            # Roughly doubles the running time (takes around 5 seconds)
            rebuild_test = structure_rebuild_test(chain)
            assert rebuild_test["pass"] == True

            # Get angles
            for res in chain:
                site = res.id[1]
                aa3 = res.resname
                aa = ThreeToOne(aa3)
                c = res.internal_coord
                # print(f" >> {site} >> {aa}")

                # According to the https://biopython.org/docs/latest/api/Bio.PDB.internal_coords.html pick_angle() documentation, these are the supported dihedral angles:
                phi = c.get_angle("phi")
                psi = c.get_angle("psi")
                chi1 = c.get_angle("chi1")
                chi2 = c.get_angle("chi2")
                chi3 = c.get_angle("chi3")
                chi4 = c.get_angle("chi4")
                chi5 = c.get_angle("chi5")
                tau = c.get_angle("tau")
                omega = c.get_angle("omega")

                # Replace None with np.nan to avoid average_angle errors later on
                phi = np.nan if phi is None else phi
                psi = np.nan if psi is None else psi
                chi1 = np.nan if chi1 is None else chi1
                chi2 = np.nan if chi2 is None else chi2
                chi3 = np.nan if chi3 is None else chi3
                chi4 = np.nan if chi4 is None else chi4
                chi5 = np.nan if chi5 is None else chi5
                tau = np.nan if tau is None else tau
                omega = np.nan if omega is None else omega

                # Isomerization state based on omega angle (for proline, mainly)
                # d()
                # iso = None
                iso = " "
                if omega is not None: 
                    if np.abs(omega) <= omega_cis_max:
                        # Cis: rare state (~ 0 degrees)
                        # iso = "cis"
                        iso = "c"
                    # else:
                    #     # Trans: common state (~ 180/-180)
                    #     iso = "trans"
                    elif np.abs(omega) >= omega_trans_min:
                        # Trans: common state (~ 180/-180)
                        # iso = "trans"
                        iso = "t"

                # print(f"   >> iso   = {iso}")
                # print(f"   >> phi   = {phi}")
                # print(f"   >> psi   = {psi}")
                # print(f"   >> omega = {omega}")
                # print(f"   >> chi1  = {chi1}")
                # print(f"   >> chi2  = {chi2}")
                # print(f"   >> chi3  = {chi3}")
                # print(f"   >> chi4  = {chi4}")
                # print(f"   >> chi5  = {chi5}")
                # print(f"   >> tau   = {tau}")

                # # Replace None with np.nan
                # phi = np.nan if phi is None else phi
                # psi = np.nan if psi is None else psi
                # chi1 = np.nan if chi1 is None else chi1
                # chi2 = np.nan if chi2 is None else chi2
                # chi3 = np.nan if chi3 is None else chi3
                # chi4 = np.nan if chi4 is None else chi4
                # chi5 = np.nan if chi5 is None else chi5
                # tau = np.nan if tau is None else tau
                # omega = np.nan if omega is None else omega
                # state = np.nan if state is None else state
                # 
                # # Replace None with np.nan
                # phi = np.nan if phi is None else phi
                # psi = np.nan if psi is None else psi
                # chi1 = np.nan if chi1 is None else chi1
                # omega = np.nan if omega is None else omega
                # state = np.nan if state is None else state
                # 
                # # Round all to 3 decimals for comparison
                # phi = np.round(phi, 3)
                # psi = np.round(psi, 3)
                # chi1 = np.round(chi1, 3)
                # omega = np.round(omega, 3)
                # rows.append([site, aa, phi, psi, chi1, omega, state])

                rows.append([frag, site, aa, iso, phi, psi, omega, chi1, chi2, chi3, chi4, chi5, tau])


    # Create data frame to save
    dihedrals.append(pd.DataFrame(rows, columns=["frag", "site", "aa", "iso", "phi", "psi", "omega", "chi1", "chi2", "chi3", "chi4", "chi5", "tau"]))

    # # Verify that all AAs are proline
    # if not isos['aa'].eq('P').all():
    #     Die(f"Error: Expected all residues to be proline in acc '{acc}', but found other residues")

    # Save data frame
    # df.to_csv(f"{outfile}.csv", sep="\t", index=False)

    Log(f"successfully calculated dihedral angles for acc|frag", f"{acc}|{frag}")

# Verify that all fragments have been run
if (frag != maxfrag):
    Die(f"Error: Expected to have processed {maxfrag} fragments, but only processed {frag}")

# Accession is complete (proline isomerization state detection complete for all fragments):
# Combine list of data frames (for performance) into a single data frame
dihedrals = pd.concat(dihedrals)
# Combine output across fragments (by using the union of all isos (ignoring any from dubious regions within 200 aa of artificial termini), and averaging distances)
dihedrals = CombineFragments(dihedrals)

# Insert into table
data = dihedrals.to_json(orient='records').replace(":", "\\:")
# Query(f"TRUNCATE {tmptable}")
query = Query(f"INSERT INTO {tmptable} SET dihedrals='{data}'")
# d()
if not Switch('debug'):
    # Query("UPDATE alphasa_tmp SET iso=NULL, phi=NULL, psi=NULL, omega=NULL, chi1=NULL, chi2=NULL, chi3=NULL, chi4=NULL, chi5=NULL, tau=NULL")
    query = Query(f"""UPDATE {alphasa} a, {tmptable} t SET 
    a.iso=NULLIF(JSON_UNQUOTE(JSON_EXTRACT(dihedrals, CONCAT('$[', a.site-1, '].iso'))), 'null'),
    a.phi=NULLIF(JSON_UNQUOTE(JSON_EXTRACT(dihedrals, CONCAT('$[', a.site-1, '].phi'))), 'null'),
    a.psi=NULLIF(JSON_UNQUOTE(JSON_EXTRACT(dihedrals, CONCAT('$[', a.site-1, '].psi'))), 'null'),
    a.omega=NULLIF(JSON_UNQUOTE(JSON_EXTRACT(dihedrals, CONCAT('$[', a.site-1, '].omega'))), 'null'),
    a.chi1=NULLIF(JSON_UNQUOTE(JSON_EXTRACT(dihedrals, CONCAT('$[', a.site-1, '].chi1'))), 'null'),
    a.chi2=NULLIF(JSON_UNQUOTE(JSON_EXTRACT(dihedrals, CONCAT('$[', a.site-1, '].chi2'))), 'null'),
    a.chi3=NULLIF(JSON_UNQUOTE(JSON_EXTRACT(dihedrals, CONCAT('$[', a.site-1, '].chi3'))), 'null'),
    a.chi4=NULLIF(JSON_UNQUOTE(JSON_EXTRACT(dihedrals, CONCAT('$[', a.site-1, '].chi4'))), 'null'),
    a.chi5=NULLIF(JSON_UNQUOTE(JSON_EXTRACT(dihedrals, CONCAT('$[', a.site-1, '].chi5'))), 'null'),
    a.tau=NULLIF(JSON_UNQUOTE(JSON_EXTRACT(dihedrals, CONCAT('$[', a.site-1, '].tau'))), 'null')
    WHERE a.acc='{acc}' AND afdb='{afdb}'""")

    # # Shorter JSON operator notation in MySQL 8.3 that includes unquoting (->>), see https://dev.mysql.com/doc/refman/8.3/en/json-search-functions.html
    # Still on MySQL 8.0 though
    affected += Numrows(query)

Show(lim=20)

print(f"Successfully inserted dihedral angles and proline isomerization states (cis/trans) into table '{alphasa}'")

print(f"\nRows affected: {Comma(affected)}")

print("\nDone!")
