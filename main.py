#!/usr/bin/env python3
"""
main.py: Prepare and submit jobs that insert DSSP accessible surface area results into SQL table 'alphasa', and insert Lahuta contact results into SQL table 'alphacon'.
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

alphafrag = "alphafrag"     # SQL table with fragment protein sequences (>2700 aa proteins get split into 1400 aa fragments with a step size of 200 in AlphaFold DB, for human only - other species don't have results for >2700 aa proteins)
alphaseq = "alphaseq"       # SQL table with complete protein sequences
alphasa = "alphasa"         # SQL table with residue-level accessible surface area values from DSSP
alphacon = "alphacon"       # SQL table with residue-level contacts from Lahuta

# Maximum number of jobs to have running simultaneously on the LSF cluster (if more jobs than this are running, the script will sleep before submitting more jobs)
# mysql.connector seems fairly smart about thread usage, so there's actually no risk of hitting the MySQL connection limit of 1000.
# 5000 is the maximum number of overall jobs per user (according to "busers"/"busers all"), so it doesn't make sense to submit a lot more:
# "bqueues -l standard" shows current scheduling priority per user
# maxjobs = 4000
# # I did run into "OS errno: 24 - Too many open files", so I do need to set some kind of limit here
# Should be fixed now
maxjobs = 2000
# maxjobs = 500
# maxjobs = 100
maxpending = 2000
# maxpending = 1000
# maxpending = 500
# maxpending = 100
# maxpending = maxjobs
# Seconds to pause between LSF requests
sleeptime = 1

# Spread jobs between multiple queues (evenly)
queues = ["standard"]
# queues = ["short"]
# queues = ["priority"]
# queues = ["short", "standard"]
# queues = ["short", "priority"]



Args(0,
"""-alphasync: Syncing: Only re-run updated AlphaSync proteins
 -alphakeep: Syncing: Only re-run updated AlphaSync proteins, and keep existing AlphaSync data. Use this for re-runs of 'main.py -alphasync' jobs.
 -debug: Don't submit cluster jobs (only print the submit commands that would have been used)
 -humanonly: Parse only human TAR file
 -keepincompletes: Keep incomplete proteins in tables alphasa and alphacon (rather than deleting and re-running them).""",
" -debug -humanonly")

alphasyncpath = "input/alphasync"
ftppath = "input/ftp"
gcspath = "input/gcs"
diagpath = "input/diagnostic_lists"
tmppath = "tmp"
logdir = "_logs"
# Replace non-word characters with underscores
tmplogdir = re.sub(r"[^\w]", "_", logdir)
logpath = f"{tmppath}/{logdir}"
mainpath = os.getcwd()

if Switch("debug"):
    tmpdebug = " -debug"
    tmpdebug2 = "__debug"
else:
    tmpdebug = ""
    tmpdebug2 = ""

if Switch('alphakeep'):
    SetSwitch('alphasync')



# Start



# Get annotation from table 'alphafrag'

# Get fragment counts per acc|source
Starttime()
print(f"Initialize: Getting fragment counts per UniProt accession and source file from table '{alphafrag}'...")
frags = FetchMap(Query(f"SELECT CONCAT(acc, '|', source) AS accsource, MAX(frag) AS maxfrag FROM {alphafrag} GROUP BY acc, source"))

# # Get annotation per acc|frag|source (just for logging)
# print(f"Getting annotation (species, sequence etc.) per UniProt accession and source file from table '{alphafrag}'...")
# annotation = FetchMap(Query(f"SELECT CONCAT_WS('|', acc, frag, source), CONCAT_WS('|', name, species, tax, fragstart, fragstop, seq) FROM {alphafrag}"))



# Get list of running jobs from LSF
running_accs = set()
for line in ReturnList("bjobs -w 2> /dev/null | tail -n+2"): # Skip the bjobs header & redirect "No unfinished job found"

    # Parse job name from LSF output
    job = re.split(r" +", line)[6]
    # d()

    # UniProt accession regular expression from the UniProt help section:
    # https://www.uniprot.org/help/accession_numbers
    # ^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$

    # # e.g. update_alphasync_tmp__logs_job_py_swissprot_cif_v4_Q9Y6Y0_1
    # m = rx(r"job_py_(\w+)_(([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}))_\d+$", job)
    # e.g. update_alphasync_tmp__logs_cd____A0A6P6LMY9_________job_py_A0A6P6LMY9_1
    m = rx(r"^update_alphasync_.*job_py_(([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})_?(\d+)?)_\d+(__alphasync)?$", job)
    if m:
        # source = m[0]
        acc = m[1]
        if m[3] is not None:
            acc += f"-{m[3]}"
        running_accs.add(acc)
        # State(acc)

if (len(running_accs) > 0):
    print(f"\nInitialize: Currently running on LSF:")
    for acc in nsort(running_accs):
        print(f" >> {acc}")
    print(f"\nCurrently running on LSF:")
    print(f" >> {Comma(len(running_accs))} UniProt accessions (including pending jobs), considering these completed\n")
    # print(f"\nCurrently running on LSF:\n{Comma(len(running_accs))} UniProt accessions (including pending jobs)")
    # # Not waiting for jobs to finish can lead to problems e.g. with the alphaseq_accs/alphasa_accs/alphacon_accs queries below (they might briefly be out of sync)
    # Should be fine now since these accs get ignored below
    # Waitforjobs()




# Get list of existing directories (these are okay to re-run so long as they aren't in the running_accs set - if they aren't, these directories must be from failed jobs)
# dir_exists = ReturnSet(f"ls -1U {tmppath}")
dir_exists = ReturnSet(f"find {tmppath} -mindepth 1 -type d -printf '%P\n'")    # Directories only
# Remove _logs from dir_exists set
dir_exists -= {logdir}
# Subtract running_accs from this set (since we want failed accs only)
dir_exists -= running_accs



# Get sets of accs that already exist in alphasa and alphacon (so that I can skip submitting these)

alphasync_accs = set()
afdb = 1
if Switch('alphasync'):

    # If -alphasync active: Only re-run updated AlphaSync proteins
    afdb = 0

    print(f"Initialize: Getting AlphaSync-calculated (afdb = 0) accession list from table '{alphaseq}' (-alphasync active)... ", end='')
    alphasync_accs = FetchSet(Query(f"SELECT DISTINCT acc FROM {alphaseq} WHERE afdb='{afdb}'"))
    print(Comma(len(alphasync_accs)))
    # Only run AlphaSync accs
    alphaseq_accs = alphasync_accs

    # Consider all (AlphaSync) accs unfinished (to re-run them)
    alphasa_accs = set()
    alphacon_accs = set()

#     # alphasa contains accessible surface area results from DSSP
#     print(f"Initialize: Getting accession list from table '{alphasa}'... ", end='')
#     # alphasa_accs = FetchSet(Query(f"SELECT DISTINCT acc FROM {alphasa}"))
#     # alphasa now also contains dihedral angles and proline isomerisation states (cis/trans)
#     alphasa_accs = FetchSet(Query(f"SELECT DISTINCT acc FROM {alphasa} WHERE iso IS NOT NULL AND afdb='{afdb}'"))
#     print(Comma(len(alphasa_accs)))
# 
#     # alphacon contains contacts from Lahuta
#     print(f"Initialize: Getting accession list from table '{alphacon}'... ", end='')
#     alphacon_accs = FetchSet(Query(f"SELECT DISTINCT acc FROM {alphacon} WHERE afdb='{afdb}'"))
#     print(Comma(len(alphacon_accs)))
# 
#     # alphaseq also contains information on accs that do not have any contacts
#     print(f"Initialize: Getting accession list known to be without contacts from table '{alphaseq}'... ", end='')
#     alphanocon_accs = FetchSet(Query(f"SELECT DISTINCT acc FROM {alphaseq} WHERE nocon=1 AND afdb='{afdb}'"))
#     print(Comma(len(alphanocon_accs)))
#     # Add these as "completed" alphacon accs (they were previously run, but no contacts exist in these structures)
#     alphacon_accs |= alphanocon_accs

else:

    # Run all proteins

    # alphaseq contains the "desired" set of accs (alphaseq.py can be run using -comparaonly to restrict the set of UniProt accessions DSSP and Lahuta will be run on)
    print(f"Initialize: Getting accession list from table '{alphaseq}' with afdb={afdb}... ", end='')
    alphaseq_accs = FetchSet(Query(f"SELECT DISTINCT acc FROM {alphaseq} WHERE afdb='{afdb}'"))
    print(Comma(len(alphaseq_accs)))

    # To recalculate everything, use empty sets:
    # alphasa_accs = set()
    # alphacon_accs = set()

    # alphasa contains accessible surface area results from DSSP
    print(f"Initialize: Getting accession list from table '{alphasa}' with afdb={afdb}... ", end='')
    # alphasa_accs = FetchSet(Query(f"SELECT DISTINCT acc FROM {alphasa}"))
    # alphasa now also contains dihedral angles and proline isomerisation states (cis/trans)
    # alphasa_accs = FetchSet(Query(f"SELECT DISTINCT acc FROM {alphasa} WHERE iso IS NOT NULL AND afdb='{afdb}'"))
    # Much faster:
    alphasa_accs = FetchSet(Query(f"SELECT DISTINCT s.acc FROM {alphaseq} s, {alphasa} a WHERE a.iso IS NOT NULL AND s.afdb='{afdb}' AND s.acc=a.acc"))
    print(Comma(len(alphasa_accs)))

    # alphacon contains contacts from Lahuta
    print(f"Initialize: Getting accession list from table '{alphacon}' with afdb={afdb}... ", end='')
    # alphacon_accs = FetchSet(Query(f"SELECT DISTINCT acc FROM {alphacon} WHERE afdb='{afdb}'"))
    # Much faster:
    alphacon_accs = FetchSet(Query(f"SELECT DISTINCT s.acc FROM {alphaseq} s, {alphacon} c WHERE s.afdb={afdb} AND s.acc=c.acc"))
    print(Comma(len(alphacon_accs)))

    # alphaseq also contains information on accs that do not have any contacts
    print(f"Initialize: Getting accession list known to be without contacts from table '{alphaseq}' with afdb={afdb}... ", end='')
    alphanocon_accs = FetchSet(Query(f"SELECT DISTINCT acc FROM {alphaseq} WHERE nocon=1 AND afdb='{afdb}'"))
    print(Comma(len(alphanocon_accs)))
    # Add these as "completed" alphacon accs (they were previously run, but no contacts exist in these structures)
    alphacon_accs |= alphanocon_accs

Stoptime()

# Ignore running accs (meaning they also won't end up being flagged as incomplete_accs)
alphaseq_accs -= running_accs
alphasa_accs -= running_accs
alphacon_accs -= running_accs

# Remove any accs that are only partially finished so they can be re-run completely
total_accs = alphasa_accs.union(alphacon_accs)
finished_accs = alphasa_accs.intersection(alphacon_accs)
incomplete_accs = total_accs - finished_accs
# Note: This test for 'incompletes' will also catch a handful of very small proteins that do not have any non-neighbour contacts in alphacon.
# These are now flagged "nocon=1" in table 'alphaseq'.


if len(incomplete_accs) > 0:
    alphasa_accs_deleted = 0
    alphacon_accs_deleted = 0
    
    if not Switch('keepincompletes'):

        Starttime()
        print(f"\nInitialize: Deleting {len(incomplete_accs):,} incomplete accs from tables alphasa and alphacon:")

        if not Switch('debug'):
            query = Query(f"DELETE FROM {alphasa} WHERE acc IN ('" + "', '".join(incomplete_accs) + f"') AND afdb='{afdb}'")
        else:
            query = Query(f"SELECT * FROM {alphasa} WHERE acc IN ('" + "', '".join(incomplete_accs) + f"') AND afdb='{afdb}'")
        alphasa_accs_deleted = Numrows(query)
        print(f" >> Deleted {Comma(alphasa_accs_deleted)} rows from table '{alphasa}'")

        if not Switch('debug'):
            query = Query(f"DELETE FROM {alphacon} WHERE acc IN ('" + "', '".join(incomplete_accs) + f"') AND afdb='{afdb}'")
        else:
            query = Query(f"SELECT * FROM {alphacon} WHERE acc IN ('" + "', '".join(incomplete_accs) + f"') AND afdb='{afdb}'")
        alphacon_accs_deleted = Numrows(query)
        print(f" >> Deleted {Comma(alphacon_accs_deleted)} rows from table '{alphacon}'")
        Stoptime()
        print()

        # Get sets of accs again
        # alphaseq contains the "desired" set of accs (alphaseq.py can be run using -comparaonly to restrict the set of accessions DSSP and Lahuta will be run on)
        print(f"Initialize: Getting accession list from table '{alphaseq}'... ", end='')
        # alphaseq_accs = FetchSet(Query(f"SELECT DISTINCT acc FROM {alphaseq}"))
        print(Comma(len(alphaseq_accs)))
        # alphasa and alphacon contain the accessible surface area results from DSSP, and the contacts from Lahuta, respectively
        print(f"Initialize: Getting accession list from table '{alphasa}'... ", end='')
        # alphasa_accs = FetchSet(Query(f"SELECT DISTINCT acc FROM {alphasa}"))
        alphasa_accs -= incomplete_accs
        print(Comma(len(alphasa_accs)))
        print(f"Initialize: Getting accession list from table '{alphacon}'... ", end='')
        # alphacon_accs = FetchSet(Query(f"SELECT DISTINCT acc FROM {alphacon}"))
        alphacon_accs -= incomplete_accs
        print(Comma(len(alphacon_accs)))

# Ignore running accs
alphaseq_accs -= running_accs
alphasa_accs -= running_accs
alphacon_accs -= running_accs

# Verify
if alphasa_accs != alphacon_accs:
    if not Switch('keepincompletes'):
        Die(f"Initialize: Mismatch between acc sets in alphasa ({len(alphasa_accs):,} accs) and alphacon ({len(alphacon_accs):,} accs)")
    else:
        # If keeping incomplete proteins: warning only
        Warn(f"Initialize: Mismatch between acc sets in alphasa ({len(alphasa_accs):,} accs) and alphacon ({len(alphacon_accs):,} accs)")
finished_accs = alphasa_accs

# Get list of wanted accs
print(f"Initialize: Getting list of 'wanted' accessions... ", end='')
wanted_accs = alphaseq_accs - finished_accs
print(Comma(len(wanted_accs)))



# Parse AlphaFold DB TAR archives
Starttime()
if not Switch('alphasync'):
    # AlphaFold DB and AlphaSync archives
    print(f"\nGetting TAR archives from '{alphasyncpath}', '{ftppath}' and '{gcspath}' and submitting jobs (which run DSSP on fragments, combine its output using combine_fragments_dssp.py, insert ASA values into table 'alphasa', and insert residue-residue contacts from Lahuta into table 'alphacon', before cleaning up their temporary files:")
else:
    # AlphaSync archives only
    print(f"\nGetting TAR archives from '{alphasyncpath}' and submitting jobs (which run DSSP on fragments, combine its output using combine_fragments_dssp.py, insert ASA values into table 'alphasa', and insert residue-residue contacts from Lahuta into table 'alphacon', before cleaning up their temporary files:")

# Get list of TAR files to parse
if Switch('alphasync'):
    # AlphaSync updated structures only (for latest UniProt release, hence tail -n 1)
    # infiles = nsort(Return(f"ls -1 {alphasyncpath}/*.tar 2> /dev/null").split("\n"))
    infiles = nsort(Return(f"ls -1 {alphasyncpath}/alphasync_cif_*.tar | tail -n 1").split("\n"))
elif Switch('humanonly'):
    # Human only (NCBI taxon ID 9606)
    infiles = nsort(Return(f"ls -1 {alphasyncpath}/*.tar {ftppath}/*_HUMAN_*.tar {gcspath}/*-9606-*.tar 2> /dev/null").split("\n"))
else:
    # All TAR files
    infiles = nsort(Return(f"ls -1 {alphasyncpath}/*.tar {ftppath}/*.tar {gcspath}/*.tar").split("\n"))
    # Move AlphaSync and human to the front of the list
    for file in reversed(nsort(infiles)):
        if rx("alphasync", file) or rx("_HUMAN_", file) or rx("-9606-", file):
            infiles.insert(0, infiles.pop(infiles.index(file)))

# Make temporary directory for job log files
Run("Make temporary directory for job logs", f"mkdir -p {logpath}", silent=True)



# Main loop: Parse TAR files
i = 0
tmpfiles = set()
for infile in infiles:
    i += 1

    afdb = 1
    tmpalphasync = ""
    tmpalphasync2 = ""
    if infile.startswith(alphasyncpath):
        # This structure isn't in the AlphaFold Protein Structure Database (it was re-predicted by AlphaSync instead!)
        afdb = 0
        tmpalphasync = " -alphasync"
        tmpalphasync2 = "__alphasync"

        # Pass on switch -alphakeep to job.py
        if Switch('alphakeep'):
            tmpalphasync += " -alphakeep"
            tmpalphasync2 = "__alphakeep"

    if afdb == 1:
        print(f" >> {i} / {len(infiles)} >> {infile}")
    else:
        print(f" >> {i} / {len(infiles)} >> AlphaSync >> {infile}")

    # Format source file name (remove .tar, e.g. UP000000589_10090_MOUSE_v2.tar to UP000000589_10090_MOUSE_v2)
    source = re.sub(r"\.tar$", "", Basename(infile))
    
    # for member in tq(tar, total=int(Return(f"cat '{diagpath}/{source}.tar.files.txt' | wc -l"))):      # Get number of files in archive from diagnostic list (not really faster and would require those lists)
    # for member in tq(tar, total=int(Return(f"tar -tf {infile} | wc -l"))):                             # Get number of files in archive directly from tar -tf (quite slow)

    # Get number of accs contained in archive from table 'alphafrag' (should be efficient - COUNT(*) is much faster than even COUNT(id))
    tmptotal = FetchOne(Query(f"SELECT COUNT(*) FROM alphafrag WHERE source='{source}'"))
    # Skip archive if it contains no .cif.gz files (e.g. swissprot_pdb_v4)
    if tmptotal == 0:
        continue
    if rx(r'^swissprot', source):
        # swissprot_cif_v4 only has one file per acc (the mmCIF file):
        tmptotal = tmptotal * 1
    elif rx(r'^UP', source):
        # UP... proteomes have two files per acc:
        tmptotal = tmptotal * 2
    elif rx(r'^proteome', source):
        # GCS proteomes have three files per acc:
        tmptotal = tmptotal * 3
    elif rx(r'^alphasync', source):
        # AlphaSync updated acc archives have one file per acc:
        tmptotal = tmptotal * 1
    else:
        Warn(f"Unhandled source type: '{source}'")
        tmptotal = tmptotal * 1

    # Open TAR file
    # mode "r:": "Open for reading exclusively without compression" (https://docs.python.org/3/library/tarfile.html)
    tar = tarfile.open(infile, mode="r:")
    acc = None
    # prevacc = None
    myjobs = Myjobs()
    mypending = Pendingjobs()
    pending_per_queue = [Pendingjobs(queue) for queue in queues]
    queue = queues[pending_per_queue.index(min(pending_per_queue))]
    submitted = 0

    # Stream TAR archive
    for member in tq(tar, total=tmptotal):
        # print(f"   >> {member.name}")
        
        frag = None

        # Skip .pdb (PDB) files
        if rx(r"\.pdb\.gz$", member.name):
            Log("skipped pdb file", member.name)
            continue

        # Skip JSON files (metadata)
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
            
        # Get accession (for handling fragments)
        # AF-A0A009IHW8-F1-model_v2.cif.gz
        m = rx(r"^AF-(\w+(-\d+)?)-F(\d+)-model_v\d+\.cif\.gz$", member.name)
        if m:
            acc = m[0]
            frag = int(m[2])

            # if frag > maxfrag:
            #     maxfrag = frag
            # if frag > 1:
            #     Die(f"Fragment is '{frag}' for acc '{acc}'")
        else:
            Die(f"Couldn't parse '{member.name}'")
            
        # Skip acc if a job for it is already running (or pending)
        # if f"{source}/{acc}" in running_accs:
        if acc in running_accs:
            Log(f"skipped acc since a job for it is already running (or pending) for acc", acc)
            continue
        
        # Skip this acc if it already exists in tables 'alphaseq', 'alphasa' and 'alphacon
        if acc not in wanted_accs:
            Log(f"skipped acc since it was not in the wanted_accs list ({alphaseq}) for acc", acc)
            continue
        else:
            Log(f"proceeded with acc since it was in the wanted_accs list ({alphaseq}) for acc", acc)
            
        # Skip this acc if it already exists in tables 'alphasa' and 'alphacon
        if acc in finished_accs:
            Log(f"skipped acc since it was already in tables '{alphasa}' and '{alphacon}' for acc", acc)
            continue
        else:
            Log(f"proceeded with acc since it was not yet in tables '{alphasa}' and '{alphacon}' for acc", acc)
        
        # # Skip acc if temporary directory already exists for it
        # Not skipping, since these are failed jobs and should be rerun
        # if DirExists(f"{tmppath}/{acc}"):
        if acc in dir_exists:
            # Log(f"skipped acc since temporary directory already existed for acc", acc)
            # Log(f"skipped acc since temporary directory already existed for source|acc", f"{source}|{acc}")
            # continue
            Log(f"reran job where temporary directory already existed for acc (any files in it will be overwritten)", acc)

        # Replace - with _ for e.g. log filenames
        tmpacc = acc.replace("-", "_")

        # Get number of fragments for this accession
        # query = Query(f"SELECT MAX(frag), COUNT(DISTINCT frag) FROM {alphafrag} WHERE acc='{acc}' AND source='{source}'")
        # (maxfrag, tmpmaxfrag) = FetchOne(query)
        # if maxfrag != tmpmaxfrag:
        #     Die(f"Error: MAX(frag) and COUNT(DISTINCT frag) don't match for acc '{acc}' in table 'alphafrag'")
        maxfrag = frags[f"{acc}|{source}"]
        
        # Extract .cif.gz to .cif using gzip (streaming, no temporary files)
        gz = tar.extractfile(member)
        # Open extracted .cif in text mode (streaming, no temporary files)
        cif = gzip.open(gz, mode="rt")
        
        # Write .cif to temporary file for DSSP
        # os.mkdir(f"{tmppath}/{acc}")
        # Run("Make temporary directory for DSSP for this acc", f"mkdir -p {tmppath}/{source}/{acc}", silent=True)
        Run("Make temporary directory for DSSP for this acc", f"mkdir -p {tmppath}/{acc}", silent=True)
        # ciffile = f"{tmppath}/{source}/{acc}/{member.name}"
        ciffile = f"{tmppath}/{acc}/{member.name}"
        # Remove .gz from file name
        ciffile = re.sub(r"\.gz$", "", ciffile)
        
        # with open(ciffile, "x") as cifout:    # mode "x" would fail if file already exists (not desirable for re-runs for randomly failed jobs)
        with open(ciffile, "w") as cifout:
            print(cif.read(), file=cifout)

        tmpfiles.add(ciffile)

        # # Get name/species/tax/sequence for this fragment from table 'alphafrag'
        # # query = Query(f"SELECT name, species, tax, fragstart, fragstop, seq FROM {alphafrag} WHERE acc='{acc}' AND frag='{frag}' AND source='{source}'")
        # # (name, species, tax, fragstart, fragstop, seq) = FetchOne(query)
        # (name, species, tax, fragstart, fragstop, seq) = annotation[f"{acc}|{frag}|{source}"].split('|')

        # # Logging with extra annotation
        # Log(f"wrote temporary CIF file(s) for acc", acc)
        # Log(f"wrote temporary CIF file(s) for acc|frag", f"{acc}|{frag}")
        # Log(f"wrote temporary CIF file(s) for acc|fragstart|fragstop", f"{acc}|{fragstart}|{fragstop}")
        # Log(f"wrote temporary CIF file(s) for source|acc", f"{source}|{acc}")
        # Log(f"wrote temporary CIF file(s) for name", name)
        # Log(f"wrote temporary CIF file(s) for species", species)
        # Log(f"wrote temporary CIF file(s) for tax", tax)
        # Log(f"wrote temporary CIF file(s) for species|tax", f"{species}|{tax}")
        # Log(f"wrote temporary CIF file(s) for source", source)
        # Log(f"wrote temporary CIF file(s) for seq", seq)
        # 
        # # # Logging without extra annotation
        # # Log(f"wrote temporary CIF file(s) for acc", acc)
        # # Log(f"wrote temporary CIF file(s) for acc|frag", f"{acc}|{frag}")
        # # Log(f"wrote temporary CIF file(s) for source|acc", f"{source}|{acc}")
        # # Log(f"wrote temporary CIF file(s) for source", source)
        
        # Accession is complete according to alphafrag table (all fragments written to temporary CIF files):
        # Submit job, which runs DSSP and then combines its output across fragments using combine_fragments_dssp.py
        # if (acc != prevacc):
        # print(f"\n >> {len(tmpfiles)} / {maxfrag}\n\n")
        if (len(tmpfiles) == maxfrag):
            # print(f"\n >> {acc} >> {len(tmpfiles)} of {maxfrag} fragments")

            # Submit job for this completely-extracted acc
            if not Switch('debug'):

                # Change directory so log files end up in logdir
                os.chdir(logpath)

                # # Choose queue randomly
                # queue = random.choice(queues)
                # Choose queue with the fewest pending jobs
                # for queue in queues:
                #     mypending = Pendingjobs(queue)
                
                # Sleep if total job number is higher than maxjobs
                if myjobs >= maxjobs:
                    while myjobs >= maxjobs:
                        myjobs = Myjobs()
                        pending_per_queue = [Pendingjobs(queue) for queue in queues]
                        queue = queues[pending_per_queue.index(min(pending_per_queue))]
                        if myjobs >= maxjobs:
                            time.sleep(sleeptime)

                # Sleep if total pending job number is higher than maxpending
                if mypending >= maxpending:
                    while mypending >= maxpending:
                        mypending = Pendingjobs()
                        pending_per_queue = [Pendingjobs(queue) for queue in queues]
                        queue = queues[pending_per_queue.index(min(pending_per_queue))]
                        if mypending >= maxpending:
                            time.sleep(sleeptime)

                Run(f"Submit job (which runs DSSP and Lahuta, parses their results into the alphasa and alphacon MySQL tables, and cleans up its {tmppath}/acc directory once complete)", f"""bsub -P idr -J update_alphasync_tmp_{tmplogdir}_cd____{tmpacc}_________job_py_{tmpacc}_{maxfrag}{tmpalphasync2}{tmpdebug2} -L /bin/bash -env 'LSB_JOB_REPORT_MAIL=N' -q {queue} -n 1 -R "rusage[mem=4G]" "bash -c 'cd ../{acc}; ../../job.py {acc} {maxfrag}{tmpalphasync}{tmpdebug} > ../{logdir}/log-output-update_alphasync_tmp_{tmplogdir}_cd____{tmpacc}_________job_py_{tmpacc}_{maxfrag}{tmpalphasync2}{tmpdebug2}.txt 2> ../{logdir}/log-errors-update_alphasync_tmp_{tmplogdir}_cd____{tmpacc}_________job_py_{tmpacc}_{maxfrag}{tmpalphasync2}{tmpdebug2}.txt; if [[ ! -s ../{logdir}/log-errors-update_alphasync_tmp_{tmplogdir}_cd____{tmpacc}_________job_py_{tmpacc}_{maxfrag}{tmpalphasync2}{tmpdebug2}.txt ]]; then rm -f ../{logdir}/log-errors-update_alphasync_tmp_{tmplogdir}_cd____{tmpacc}_________job_py_{tmpacc}_{maxfrag}{tmpalphasync2}{tmpdebug2}.txt ../{logdir}/log-output-update_alphasync_tmp_{tmplogdir}_cd____{tmpacc}_________job_py_{tmpacc}_{maxfrag}{tmpalphasync2}{tmpdebug2}.txt; fi'" -e /dev/null -o /dev/null > /dev/null &""", silent=True)

                # print("+", end="")
                myjobs += 1
                mypending += 1
                submitted += 1
                pending_per_queue[queues.index(queue)] += 1
                queue = queues[pending_per_queue.index(min(pending_per_queue))]

                # Add submitted acc to running_accs
                # Note: This will prevent the same acc to be run in any subsequent archives.
                # This means that each acc gets run from the first archive encountered according to the order of infiles.
                # This means that the human Swiss-Prot proteome will be run first (taxon 9606), followed by human UniProt, followed by all of Swiss-Prot, followed by other species' individual UniProt proteomes.
                # ...none of which is actually important since internally, DeepMind stores all accessions' results in a single GCS folder/bucket.
                # ...hence the prediction should be exactly the same for a given acc no matter what archive it is encountered in here.
                running_accs.add(acc)
                tmpfiles = set()

                # Move back to script directory (top level)
                os.chdir(mainpath)
            else:
                # Print only
                # print(f"   >> Submit job (which runs DSSP, parses its results into the alphasa and alphacon MySQL tables, and cleans up its tmp/acc directory once complete): ~/scripts/qsub.sh ../../dssp.py {source} {acc} {maxfrag}")
                print(f"   >> Submit job (which runs DSSP and Lahuta, parses their results into the alphasa and alphacon MySQL tables, and cleans up its {tmppath}/acc directory once complete): " + f"""bsub -P idr -J update_alphasync_tmp_{tmplogdir}_cd____{tmpacc}_________job_py_{tmpacc}_{maxfrag}{tmpalphasync2}{tmpdebug2} -L /bin/bash -env 'LSB_JOB_REPORT_MAIL=N' -q {queue} -n 1 -R "rusage[mem=4G]" "bash -c 'cd ../{acc}; ../../job.py {acc} {maxfrag}{tmpalphasync}{tmpdebug} > ../{logdir}/log-output-update_alphasync_tmp_{tmplogdir}_cd____{tmpacc}_________job_py_{tmpacc}_{maxfrag}{tmpalphasync2}{tmpdebug2}.txt 2> ../{logdir}/log-errors-update_alphasync_tmp_{tmplogdir}_cd____{tmpacc}_________job_py_{tmpacc}_{maxfrag}{tmpalphasync2}{tmpdebug2}.txt; if [[ ! -s ../{logdir}/log-errors-update_alphasync_tmp_{tmplogdir}_cd____{tmpacc}_________job_py_{tmpacc}_{maxfrag}{tmpalphasync2}{tmpdebug2}.txt ]]; then rm -f ../{logdir}/log-errors-update_alphasync_tmp_{tmplogdir}_cd____{tmpacc}_________job_py_{tmpacc}_{maxfrag}{tmpalphasync2}{tmpdebug2}.txt ../{logdir}/log-output-update_alphasync_tmp_{tmplogdir}_cd____{tmpacc}_________job_py_{tmpacc}_{maxfrag}{tmpalphasync2}{tmpdebug2}.txt; fi'" -e /dev/null -o /dev/null > /dev/null""")
                submitted += 1
                tmpfiles = set()
        
            # # Logging with extra annotation
            # Log(f"submitted job for acc", acc)
            # Log(f"submitted job for source|acc", f"{source}|{acc}")
            # Log(f"submitted job for name", name)
            # Log(f"submitted job for species", species)
            # Log(f"submitted job for tax", tax)
            # Log(f"submitted job for species|tax", f"{species}|{tax}")
            # Log(f"submitted job for source", source)

            # Logging without extra annotation
            Log(f"submitted job for acc", acc)
            Log(f"submitted job for acc|frags", f"{acc}|{maxfrag}")
            Log(f"submitted job for source|acc", f"{source}|{acc}")
            Log(f"submitted job for source", source)
    
    # TAR file completely processed
    print(f"     >> Submitted {Comma(submitted)} jobs")
        
Show("submitted job for acc") if Switch('debug') else None
# Show()
Show(lim=20)

if not Switch('debug'):
    Waitforjobs()

Stoptime()

print("\nDone!")
