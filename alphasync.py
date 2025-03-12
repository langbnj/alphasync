#!/usr/bin/env python3
"""
Run AlphaFold on sequences that are needed to synchronise AlphaFold DB with the latest UniProt release.

Use table 'alphamap' to find sequences for which no matching AlphaFold DB sequence was found.
Uses sequence sources such as the current 'alphauniprot' table, obtained from UniProt via API.
"""

# Initialize
# import pandas as pd
# import numpy as np
# import networkx as nx
# import seaborn as sns
# import matplotlib.pyplot as mp
from blang_mysql import *
from blang import *

# Set SQL table names to use
alphamap = "alphamap"
alphaseq = "alphaseq"
alphauniprot = "alphauniprot"
ensembl = "ensembl"

# Paths
alphafolddir = "input/alphasync/alphafold"
cifdir = "input/alphasync/cif"
paedir = "input/alphasync/pae"
paramdir = "input/alphasync/params"
maindir = os.getcwd()
cpufirst_attempted_path = "cpufirst-attempted"

# Maximum number of jobs to have running simultaneously on the LSF cluster (if more jobs than this are running, the script will sleep before submitting more jobs)
# mysql.connector seems fairly smart about thread usage, so there's actually no risk of hitting the MySQL connection limit of 1000.
# 5000 is the maximum number of overall jobs per user (according to "busers"/"busers all"), so it doesn't make sense to submit a lot more:
# "bqueues -l standard" shows current scheduling priority per user
# maxjobs = 7500
# # I did run into "OS errno: 24 - Too many open files", so I do need to set some kind of limit here
# maxjobs = 25000
maxjobs = 10000
# maxpending_gpu = 25000
# maxpending_gpu = 2500
maxpending_gpu = 1000
# maxpending_gpu = 500
# maxpending_gpu = 100
# maxpending_gpu = 10
# maxpending_cpu = 25000
maxpending_cpu = 2000
# maxpending_cpu = 1000
# maxpending_cpu = 500
# Seconds to pause between LSF requests
sleeptime = 1

# Spread GPU jobs between multiple queues (submitting to queue with fewer pending)
queues = ["gpu"]

# Spread CPU jobs between multiple queues (submitting to queue with fewer pending)
cpu_queues = ["standard"]


# Note: Useful commands:

# Isoform structure count:
# input/alphasync/cif
# 1|g "\-\d+\-F"|wcl
# 1|wcl

# Verifying acc counts in input/alphasync output folders
# 1 cif|wcl # CIF count
# 1 params|wcl  # Params count, should match CIF count
# 1 cif|g "\-F1-"|g -o "^AF-\w+(-\d+)?-F"|suq|wcl    # Total accs
# 1 cif|g "\-F2-"|g -o "^AF-\w+(-\d+)?-F"|suq|wcl    # Minus fragmented accs
# 1 pae|wcl # PAE count, should match total CIF accs minus fragmented accs
# 1 cif|g "\-\d+\-F"|g -o "^AF-\w+(-\d+)?-F"|suq|wcl    # Isoform accs
# 1 cif|g "\-\d+\-F"|wcl    # Isoform acc|frags

# See reasons for pending jobs
# bjobs -p|suq

# Switch all standard jobs to priority:
# q|g standard|g -o "^\d{9,}"|perl -npe 's/\n/ /g'|xa -d' ' bswitch priority {}
# Switch all gpu jobs to gpu_priority:
# q|g gpu|g -o "^\d{9,}"|perl -npe 's/\n/ /g'|xa -d' ' bswitch gpu_priority {}

# Switch first 200 jobs from standard to priority:
# bjobs -p -q standard -o jobid|t1|head -n200|xa bswitch priority {}

# Switch the last 20 standard CPU jobs to priority:
# q|g standard|g -o "^\d{9,}"|tail -n20|perl -npe 's/\n/ /g'|xa -d' ' bswitch priority {}

# Switch first 40 jobs (that are pending because of gpu availability) from gpu to gpu_priority:
# bjobs -p -q gpu|g -B1 ngpus_phys|g -o "^\d+"|head -n40|xa bswitch gpu_priority {}

# Change parameters from 12 cores & 15 GB to 4 cores & 45 GB, for all jobs in 'priority':
# q|g " priority"|g " pend "|g -o "^\d+"|xa sh -c 'bmod -n 4 -R "span[hosts=1]" -R "rusage[mem=45G]" {} &'

# All structures per day:
# 1|xa lsla {}|perl -npe 's/Domain Users [\d\.]+K //'|nsort|g -o "(Sep|Oct|Nov|Dec|Jan|Feb|Mar) +\d+"|suq|perl -npe 's/^( *\d+) (\w+ +\d+)/$2$1/; s/^Sep  /2024-09-0/; s/^Sep /2024-09-/; s/^Oct  /2024-10-0/; s/^Oct /2024-10-/; s/^Nov  /2024-11-0/; s/^Nov /2024-11-/; s/^Dec  /2024-12-0/; s/^Dec /2024-12-/; s/^Jan  /2025-01-0/; s/^Jan /2025-01-/; s/^Feb  /2025-02-0/; s/^Feb /2025-02-/; s/^Mar  /2025-03-0/; s/^Mar /2025-03-/'|sort -g

# All structures per day in the last 10 days:
# find . -type f -mtime -10|xa lsla {}|perl -npe 's/Domain Users [\d\.]+K //'|nsort|g -o "(Sep|Oct|Nov|Dec|Jan|Feb|Mar) +\d+"|suq|perl -npe 's/^( *\d+) (\w+ +\d+)/$2$1/; s/^Sep  /2024-09-0/; s/^Sep /2024-09-/; s/^Oct  /2024-10-0/; s/^Oct /2024-10-/; s/^Nov  /2024-11-0/; s/^Nov /2024-11-/; s/^Dec  /2024-12-0/; s/^Dec /2024-12-/; s/^Jan  /2025-01-0/; s/^Jan /2025-01-/; s/^Feb  /2025-02-0/; s/^Feb /2025-02-/; s/^Mar  /2025-03-0/; s/^Mar /2025-03-/'|sort -g

# Isoform structures per day:
# cd input/alphasync/cif
# Sorted chronologically
# 1|g "\-\d+\-F"|xa lsla {}|perl -npe 's/Domain Users [\d\.]+K //'|nsort|g -o "(Sep|Oct|Nov|Dec|Jan|Feb|Mar) +\d+"|suq|perl -npe 's/^( *\d+) (\w+ +\d+)/$2$1/; s/^Sep  /2024-09-0/; s/^Sep /2024-09-/; s/^Oct  /2024-10-0/; s/^Oct /2024-10-/; s/^Nov  /2024-11-0/; s/^Nov /2024-11-/; s/^Dec  /2024-12-0/; s/^Dec /2024-12-/; s/^Jan  /2025-01-0/; s/^Jan /2025-01-/; s/^Feb  /2025-02-0/; s/^Feb /2025-02-/; s/^Mar  /2025-03-0/; s/^Mar /2025-03-/'|sort -g
# 1|g "\-\d+\-F"|xa lsla {}|perl -npe 's/Domain Users [\d\.]+K //'|nsort|g -o "(Jan|Feb|Mar) +\d+"|suq|perl -npe 's/^( *\d+) (\w+ +\d+)/$2$1/; s/^Jan /2025-01-/; s/^Feb  /2025-02-0/; s/^Feb /2025-02-/; s/^Mar  /2025-03-0/; s/^Mar /2025-03-/'|sort -g
# 1|g "\-\d+\-F"|xa lsla {}|perl -npe 's/Domain Users [\d\.]+K //'|nsort|g -o "(Jan|Feb|Mar) +\d+"|suq|perl -npe 's/^ *(\d+) (\w+ +\d+)/$2\t$1/; s/^Jan /2025-01-/; s/^Feb  /2025-02-0/; s/^Feb /2025-02-/; s/^Mar  /2025-03-0/; s/^Mar /2025-03-/'|sort -g
# Sorted by days with most structures
# 1|g "\-\d+\-F"|xa lsla {}|perl -npe 's/Domain Users [\d\.]+K //'|nsort|g -o "(Jan|Feb|Mar) +\d+"|suq
# 1|g "\-\d+\-F"|xa lsla {}|perl -npe 's/Domain Users [\d\.]+K //'|nsort
# All structures per day:
# lsla|perl -npe 's/Domain Users [\d\.]+K //'|nsort|g -o "(Oct|Nov|Dec|Jan) \d+"|suq

# To see why jobs are queued (pending_reason):
# bjobs -q gpu -pl|g 'PENDING REASONS' -A1|suq
# bjobs -q gpu -pl|g 'PENDING REASONS' -A6|suq

# Show disk usage:
# cd input/alphasync/alphafold
# du -hd0 .

# Show MSA file counts across all running & incomplete accs:
# cd input/alphasync/alphafold
# 1 */*/alphafold_params.json */*/msas/*|g -o "[^/]+$"|suq

# Print errors:
# cd input/alphasync/alphafold
# find . -wholename "*/log-err*"|xa cat {}|suq
# Show non-empty error logs
# find . -wholename "*/log-err*" -size +0c
# Print errors from previous 24 hours:
# find . -wholename "*/log-err*" -mtime -1|xa cat {}|suq
# Show non-empty error logs from previous 24 hours:
# find . -wholename "*/log-err*" -mtime -1 -size +0c

# To find directories with ≥1000 files:
# find . -mindepth 1 -type d -print0|xargs -0 -I{} echo '{}'|natsort|xargs -I{} sh -c 'echo -n "{}"; echo -ne "\t"; find "{}" -mindepth 1 | wc -l'|perl -npe 's/^(.+)\t(.+)/$2 $1/'|sort -n|awk '$1 >= 1000'

# Note: -cpufirst speedup:
# e.g. Q14C86-6_F1:
# CPU 3 h  45 mins  37.493 sec
# GPU 1 h  17 mins  17.216 sec
# >> Huge reduction in GPU time by using -cpufirst (should be ~25%).




# infile = Args(1, "[input sequence file]", "../alphafold_db/input/UP000005640_9606_HUMAN_v2.tar.seqs.txt")
(type, version, max_unmapped_seqs) = Args(3, "[sequence database to sync to: uniprot/ensembl] [version: e.g. 2025_01/113] [maximum number of unmapped sequences for a given species (otherwise it'll be skipped): e.g. 1000, or 'all'] [-humanonly] [-modelonly] [-modelhealthonly] [-iso] [-nofrag] [-nopeptides] [-no_b] [-no_z] [-no_u] [-no_x] [-submitnow] [-cpufirst] [-debug] [-debug2]",
f""" -iso: Include isoform accessions (e.g. P04637-4)
 -humanonly: Only run on human sequences (taxon 9606)
 -modelonly: Only run on model organism sequences (human, mouse, Drosophila, C. elegans, yeast) (taxa: 9606, 10090, 7227, 6239, 559292)
 -modelhealthonly: Only run on model organism sequences and AlphaFold DB Global Health Proteomes (48 species in total)
 -iso: Include isoform accessions (e.g. P04637-4)
 -nofrag: Skip sequences that would need to be fragmented (>=2700 aa)
 -nopeptides: Skip small peptides (<16 aa). AlphaFold DB excluded these as AlphaFold 2 can be unreliable for them.
 -no_b: Skip B-containing sequences (asparagine/aspartate), rather than keeping them by replacing B with asparagine (N)
 -no_z: Skip Z-containing sequences (glutamine/glutamate), rather than keeping them by replacing Z with glutamine (Q)
 -no_u: Skip U-containing sequences (selenocysteine), rather than keeping them by replacing U with cysteine (C)
 -no_x: Skip X-containing sequences (unknown amino acids), rather than keeping them by replacing terminal X with glycine and internal X with alanine (A)
 -submitnow: Submit all jobs immediately, instead of spacing them out by having a maximum number of jobs pending at a time (and {maxjobs:,} total jobs active)
 -cpufirst: To optimize cluster usage, first run AlphaFold in a CPU-only mode to generate multiple sequence alignments (MSAs). Requires re-running alphasync.py a second time later on to predict the final structures (GPU-based).
 -debug: Print debug output, and don't actually submit jobs
 -debug2: Do not print debug output, and don't actually submit jobs""", 
 "uniprot 2025_01 1000 -humanonly -iso")

# If switch '-submitnow' is active, submit all jobs immediately, instead of spacing them out by having a maximum number of jobs pending at a time
if Switch('submitnow'):
    # (Still have 100,000 as a reasonable limit - above this I'm not sure LSF would be happy)
    maxjobs = 100000
    maxpending_gpu = 100000
    maxpending_cpu = 100000

if not Switch('iso'):
    # By default, only include canonical accessions (e.g. P04637) (no isoforms such as P04637-4)
    criteria = "(au.reviewed=1 OR au.refprotcanon=1) AND au.acc=au.canon"
else:
    # '-iso' is active: Include isoforms.
    criteria = "(au.reviewed=1 OR au.refprotcanon=1)"


if type == "uniprot":
    # q = f"SELECT m.value, au.seq FROM {alphamap} m, {alphauniprot} au WHERE m.value=au.acc AND m.type='{type}' AND m.version='{version}' AND m.map IS NULL"
    # refproteome=1: Still far too many
    # q = f"SELECT m.value, au.seq FROM {alphamap} m, {alphauniprot} au WHERE m.value=au.acc AND au.refproteome=1 AND m.type='{type}' AND m.version='{version}' AND m.map IS NULL"
    # # acc=canon (no isoforms) and reviewed=1 (Swiss-Prot) only: manageable number for human (376 proteins)
    # q = f"SELECT m.value, au.seq FROM {alphamap} m, {alphauniprot} au WHERE m.value=au.acc AND au.acc=au.canon AND au.reviewed=1 AND m.type='{type}' AND m.version='{version}' AND m.map IS NULL"
    # acc=canon (no isoforms) and (reviewed=1 (Swiss-Prot) only OR refprotcanon=1): manageable number for human (376 proteins)
    q1 = "SELECT m.value, au.seq"
    q2 = f" FROM {alphamap} m, {alphauniprot} au WHERE m.value=au.acc AND {criteria} AND m.type='{type}' AND m.version='{version}' AND m.map IS NULL"
elif type == "ensembl":
    q1 = "SELECT m.value, e.seq"
    q2 = f" FROM {alphamap} m, {ensembl} e WHERE m.value=e.ensp AND m.type='{type}' AND m.version='{version}' AND m.map IS NULL"
else:
    Die(f"Unhandled sequence database '{type}'")


if Switch('humanonly'):
    q2 += " AND m.tax=9606"
    # -iso (with frags):
    # >> 4,502 jobs

if Switch('modelonly'):
    # HUMAN MOUSE DROME CAEEL YEAST
    q2 += " AND m.tax IN (9606, 10090, 7227, 6239, 559292)"

if Switch('modelhealthonly'):
    # Model organisms & Global Health Proteomes

    # All accs (without {criteria} = reviewed/refprotcanon/acc=canon)
    # q2 += " AND m.tax IN (SELECT DISTINCT tax FROM alphafrag WHERE source IN (SELECT DISTINCT source FROM alphafrag HAVING source LIKE 'UP%'))"
    # >> ~57,349 jobs (+ frags)

    # With -iso (rest of isoforms, other than human):
    # alphasync.py uniprot 2025_01 -modelhealthonly -iso -debug2
    # >> 425 jobs (only!)

    # Skip species that have more than max_unmapped_seqs proteins for which we would need to predict structures (some are very unreasonably high since most of the proteome has changed since AFDB's version)
    tmp_max_unmapped_seqs = f" HAVING unmapped_seqs<={max_unmapped_seqs}"
    q2 += f" AND m.tax IN (SELECT t.tax FROM (SELECT au.tax, au.species, COUNT(DISTINCT au.seq) AS unmapped_seqs FROM alphauniprot au, alphamap m WHERE {criteria} AND au.acc=m.value AND m.type='{type}' AND m.version='{version}' AND m.map IS NULL AND au.tax IN (SELECT DISTINCT tax FROM alphafrag WHERE source IN (SELECT DISTINCT source FROM alphafrag HAVING source LIKE 'UP%')) GROUP BY au.tax{tmp_max_unmapped_seqs}) t)"

if Switch('comparaonly'):
    # Ensembl Comparative Genomics pipeline species (200 representative species, mostly vertebrates) (www.ensembl.org)
    # q2 += " AND m.tax IN (SELECT DISTINCT tax FROM compara_species)"
    # >> ~35,580 jobs (+ frags)

    # Skip species that have more than max_unmapped_seqs proteins for which we would need to predict structures (some are very unreasonably high since most of the proteome has changed since AFDB's version)
    # max_unmapped_seqs = 1000
    # >> 6,385 -nofrag jobs
    # max_unmapped_seqs = 100
    # >> 1,800 -nofrag jobs
    tmp_max_unmapped_seqs = f" HAVING unmapped_seqs<={max_unmapped_seqs}"
    q2 += f" AND m.tax IN (SELECT t.tax FROM (SELECT au.tax, au.species, COUNT(DISTINCT au.seq) AS unmapped_seqs FROM alphauniprot au, alphamap m WHERE {criteria} AND au.acc=m.value AND m.type='{type}' AND m.version='{version}' AND m.map IS NULL AND au.tax IN (SELECT DISTINCT tax FROM compara_species) GROUP BY au.tax{tmp_max_unmapped_seqs}) t)"


# Also ensure that sequences are unique (to avoid running the same sequence multiple times)
# # This would be nice and efficient (alphamap is clever enough to do it, see e.g. https://alphasync.stjude.org/display/A0A8C8L8Z5 and https://alphasync.stjude.org/display/A0PJE2-4), but sadly, I'll need complete CIF, PAE and Params folders for archives etc.
# # alphamap and even the Display download buttons will still use the best avg_plddt structure for these cases, but the Downloads TAR archives will contain all the accessions.
# Actually: the archives are already not complete. I only predict new structures where there isn't already a match to an existing sequence, even if that sequence is from a different acc or even species. So the archives are never fully complete.
# They only provide structures for sequences that are otherwise missing. They don't provide structures for all missing/updated accs.
# So: ensuring I only predict one structure per sequence makes good sense here.
if type == "uniprot":
    q2 += " GROUP BY au.seq"
elif type == "ensembl":
    q2 += " GROUP BY e.seq"
else:
    Die(f"Unhandled sequence database '{type}'")

# # Also filter out sequences that will be fragmented, are too short, or contain X
# q2 += f" HAVING au.seq NOT REGEXP 'X' AND LENGTH(au.seq)>=16 AND LENGTH(au.seq)<=2699"

# Also filter out small peptides (<16 aa) as AlphaFold 2 reportedly is unreliable on them (and AFDB does not contain any peptides shorter than 16 aa)
# As there are only very few (32 in human, 184 in all of modelhealth), I am now allowing them in order to reach full proteome coverage. I assume the pLDDT scores will reflect that these are not good predictions.
if Switch('nofrag'):
    # -nofrag: Filter out sequences that would need to be fragmented (>=2700 aa)
    if Switch('nopeptides'):
        # Exclude small peptides <16 aa
        q2 += f" HAVING LENGTH(au.seq)>=16 AND LENGTH(au.seq)<=2699"
    else:
        # Include small peptides <16 aa
        q2 += f" HAVING LENGTH(au.seq)<=2699"
else:
    if Switch('nopeptides'):
        # Exclude small peptides <16 aa
        q2 += f" HAVING LENGTH(au.seq)>=16"


# if Switch('debug'):
#     # q = f"SELECT m.value, au.seq FROM {alphamap} m, {alphauniprot} au WHERE m.value=au.acc AND au.reviewed=1 AND m.type='{type}' AND m.version='{version}'"
#     # q += " AND m.value='P20929-1'"
#     # q += " AND m.value='P20929-1'"
#     q += " ORDER BY m.value='Q92616' DESC" # length 2671

# Order by taxon (lower IDs first, which is likely to prioritise more standard model organisms), then order randomly within it (with fixed seed), rather than running accessions alphabetically
q2 += " ORDER BY m.tax, RAND(42)"

# Assemble query
q = q1 + q2

if Switch('debug'):
    print(f"Query: {q}")
    # sys.exit()

# Create main directories in case they didn't exist yet
if not Switch('debug') and not Switch('debug2'):
    os.makedirs(alphafolddir, exist_ok=True)
    os.makedirs(cifdir, exist_ok=True)
    os.makedirs(paedir, exist_ok=True)
    os.makedirs(paramdir, exist_ok=True)



# # Utility functions
# Now in blang.py
# 
# # Replace X stretches with alanines (A), used for short internal (≤3 aa) stretches of X
# def replace_with_alanines(match):
#     return "A" * len(match.group())
# 
# # Replace X stretches with glycine linkers (GGGGS repeated), used for terminal or long internal (≥4 aa) stretches of X
# def replace_with_ggggs(match):
#     n = len(match.group())
#     return ("GGGGS" * (n//5) + "GGGGS"[:n%5])



# Start

# Get list of running jobs from LSF
running_frags = set()
for line in ReturnList("bjobs -w 2> /dev/null | tail -n+2"): # Skip the bjobs header & ignore e.g. "No unfinished job found" error

    # UniProt accession regular expression from the UniProt help section:
    # https://www.uniprot.org/help/accession_numbers
    # ^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$

    m = rx(r"update_alphasync_input_alphasync_alphafold_(([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(_\d+)?)_(F\d+)", line)
    if m:
        acc = m[0]
        frag = m[4]
        running_frags.add(f"{acc}_{frag}")

if (len(running_frags) > 0):
    print(f"\nInitialize: Currently running on LSF:")
    # for frag in running_frags:
    for frag in nsort(running_frags):
        print(f" >> {frag}")
    print(f"\nCurrently running on LSF:")
    print(f" >> {Comma(len(running_frags))} fragments (including pending jobs), considering these completed\n")



print(f"\nGetting unmapped '{type}' version '{version}' sequences from table '{alphamap}' and submitting AlphaFold structure prediction jobs for them:")

# Get total jobs
myjobs = Myjobs()

# Get pending jobs for (shortest) CPU queue
if Switch('cpufirst'):
    pending_per_cpu_queue = [Pendingjobs(cpu_queue) for cpu_queue in cpu_queues]
    cpu_queue = cpu_queues[pending_per_cpu_queue.index(min(pending_per_cpu_queue))]
    mypending_cpu = Pendingjobs(cpu_queue)
    submitted_cpu = 0

# Get pending jobs for (shortest) GPU queue
pending_per_queue = [Pendingjobs(queue) for queue in queues]
queue = queues[pending_per_queue.index(min(pending_per_queue))]
mypending = Pendingjobs(queue)
submitted = 0

for value, seq in Fetch(Query(q)):

    if Switch('debug'):
        print()
        print(f" >> {value}")
        print(f"   >> Full sequence ({len(seq):,} aa):")
        print(f"     >> {seq}")

    Log("total accs", value)
    Log("total full sequences", seq)

    # Replace - with _ for e.g. log filenames
    tmpvalue = value.replace("-", "_")
    
    # Replace non-standard amino acids:
    # Could also use ReplaceNonstandardAAs() now, but here, we want reporting of what's replaced

    # Get non-standard amino acids
    seq_before = seq
    nonstd = re.findall("[^ACDEFGHIKLMNPQRSTVWY]", seq)
    for aa in nonstd:
        Log("non-standard amino acids found", aa)

    # if len(seq) < 2700:
    #     continue
    # if len(seq) < 100:
    #     continue
    # if len(seq) > 200:
    #     continue
    # # Skip sequences below 16 aa (AlphaFold DB minimum is 16: SELECT MIN(LENGTH(seq)), MAX(LENGTH(seq)) FROM alphafrag WHERE afdb=1;)
    # No need to skip them here as the above query already skips these
    # if len(seq) < 16:
    #     if Switch('debug'):
    #         print(f"       >> SKIP (sequence length <16 aa)")
    #     Log("skipped because sequence length <16 aa for acc (skipped)", value)
    #     # Log("skipped because sequence length <16 aa for sequence (skipped)", seq)
    #     continue

    # Selenocysteine-containing proteins: by default, keep them by replacing selenocysteine (U) with cysteine (C)
    if rx("U", seq):
        if not Switch('no_u'):
            # Replace selenocysteine (U) with cysteine (C)
            seq = seq.replace("U", "C")
            Log("selenocysteine (U) replaced with cysteine (C) for acc (kept)", value)
            # Log("selenocysteine (U) replaced with cysteine (C) for sequence (kept)", seq)
            # Log("selenocysteine (U) replaced with cysteine (C) for sequence (kept)", seq)

    # B/Z: by default, keep them by replacing B (D/N) with asparagine (N) and Z (E/Q) with glutamine (Q)
    if rx("B", seq):
        if not Switch('no_b'):
            # Replace B with asparagine (N)
            seq = seq.replace("B", "N")
            Log("ambiguous asparagine/aspartate (B) replaced with asparagine (N) for acc (kept)", value)
    if rx("Z", seq):
        if not Switch('no_z'):
            # Replace Z with glutamine (Q)
            seq = seq.replace("Z", "Q")
            Log("ambiguous glutamine/glutamate (Z) replaced with glutamine (Q) for acc (kept)", value)
            # Log("glutamine/glutamate (Z) replaced with glutamine (Q) for sequence (kept)", seq)

    # X: by default, keep them by replacing terminal or long internal stretches of X with glycine linkers (GGGGS repeated), and short internal stretches of X with alanines (A)
    if rx("X", seq):
        if not Switch('no_x'):
            # # Replace N- or C-terminal X stretches (any length) with glycine linkers (GGGGS repeated)
            # seq = re.sub(r'^X+|X+$', replace_with_ggggs, seq)
            # # Replace short internal X stretches ≤3 aa with alanines (A) to maintain secondary structure
            # seq = re.sub(r'(?<!^)X{1,3}(?!$)', replace_with_alanines, seq)
            # # Replace long internal X stretches ≥4 aa with glycine linkers (GGGGS repeated)
            # seq = re.sub(r'(?<!^)X{4,}(?!$)', replace_with_ggggs, seq)
            # # (?<!^) - Negative Lookbehind:
            # # ?<! means "make sure what's before this position is not ^" (start of the string)
            # # (?!$) - Negative Lookahead:
            # # ?! means "make sure what's after this position is not $" (end of the string)

            # Replace N-terminal X stretches (any length) with glycine linkers (GGGGS repeated) to maintain flexibility
            if rx(r"^X+", seq):
                seq = re.sub(r'^X+', replace_with_ggggs, seq)
                Log("unknown N-terminal amino acid stretch (X...) replaced with glycine linker (GGGGS... repeated) for acc (kept)", value)

            # Replace C-terminal X stretches (any length) with glycine linkers (GGGGS repeated) to maintain flexibility
            if rx(r"X+$", seq):
                seq = re.sub(r'X+$', replace_with_ggggs, seq)
                Log("unknown C-terminal amino acid stretch (...X) replaced with glycine linker (...GGGGS repeated) for acc (kept)", value)

            # Replace long internal X stretches ≥4 aa with glycine linkers (GGGGS repeated)
            if rx(r"(?<!X)X{4,}(?!X)", seq):
                seq = re.sub(r'(?<!X)X{4,}(?!X)', replace_with_ggggs, seq)
                Log("unknown long internal amino acid stretch (...XXXX... or longer) replaced with glycine linker (...GGGGS... repeated) for acc (kept)", value)

            # Replace short internal X stretches ≤3 aa with alanines (A) to maintain secondary structure
            if rx(r"(?<!X)X{1,3}(?!X)", seq):
                seq = re.sub(r'(?<!X)X{1,3}(?!X)', replace_with_alanines, seq)
                Log("unknown short internal amino acid stretch (...X..., ...XX..., ...XXX...) replaced with alanines (AAA) for acc (kept)", value)

            Log("unknown amino acid(s) (X) replaced with alanine(s) (A) or glycine linker(s) (GGGGS repeated) for acc (kept)", value)

    # If non-standard amino acids (U/B/Z/X) were found, show them and the processed sequence
    if nonstd != []:
        for aa in nonstd:
            Log("non-standard amino acids kept (kept)", aa)
        if Switch('debug'):
            print(f"   >> Processed sequence after replacing U/B/Z/X ({len(seq):,} aa):")
            print(f"     >> {seq}")

    # Verify that sequence length remained the same
    if len(seq) != len(seq_before):
        Die(f"Error: Sequence length changed from {len(seq_before):,} to {len(seq):,} for acc '{value}':\n\n{seq_before}\n\n{seq}\n\n")

    # Skip sequences still containing characters other than the standard 20 amino acids (after replacements above, some or all of which might be deactivated due to switches)
    if not rx("^[ACDEFGHIKLMNPQRSTVWY]+$", seq):
        if Switch('debug'):
            print(f"       >> SKIP (non-standard amino acids)")
        Log("skipped because non-standard amino acids found for acc (skipped)", value)
        # Log("skipped because non-standard amino acids found for sequence (skipped)", seq)

        # Get non-standard amino acids
        nonstd = re.findall("[^ACDEFGHIKLMNPQRSTVWY]", seq)
        for aa in nonstd:
            Log("non-standard amino acids skipped (skipped)", aa)

        continue

    
    
    # Get number of fragments (maxfrag)
    if len(seq) <= 2699:
        frag = 1
    else:
        frag = 1
        for start in range(0, len(seq), 200):
            end = min(start + 1400, len(seq))
            if end == len(seq):
                break
            frag += 1
    maxfrag = frag



    # Check if complete output CIF and PAE files (for all fragments) already exist for this acc in cifdir and paedir
    # e.g. AF-{acc}-F{frag}-model_v0.cif.gz
    # if (len(glob(f"{cifdir}/AF-{value}-F*-model_v0.cif.gz")) == maxfrag) and (len(glob(f"{paedir}/AF-{value}-F*-predicted_aligned_error_v0.json.gz")) == maxfrag):
    if len(glob(f"{cifdir}/AF-{value}-F*-model_v0.cif.gz")) == maxfrag:
        if Switch('debug'):
            print(f"   >> SKIP (CIF output files already complete for acc '{value}')")
        Log("CIF output files were already complete for all fragments for acc (skipped)", value)
        continue
    else:
        if Switch('debug'):
            print(f"   >> No output yet for acc '{value}' ({maxfrag} fragments), running AlphaFold!")
        # d()
    if Switch('debug'):
        print()



    # Make directory for this accession (AlphaFold will be run here)
    # Write FASTA files (fragmented if >= 2700 aa, into fragments of 1400 aa with a step size of 200 aa)
    # (AlphaFold DB maximum is 2699: SELECT MIN(LENGTH(seq)), MAX(LENGTH(seq)) FROM alphafrag WHERE afdb=1;)
    frag = 1
    if len(seq) <= 2699:

        Log("total acc|frags", f"{value}|{frag}")

        if f"{tmpvalue}_F{frag}" in running_frags:
            if Switch('debug'):
                print(f"       >> SKIP: {value}_F{frag} is already running on LSF")
            Log("already running on LSF for acc|frag (skipped)", f"{value}|{frag}")
            continue

        # Check if output CIF and PAE files for this fragment already exist for this acc in cifdir and paedir
        # if (Exists(f"{cifdir}/AF-{value}-F{frag}-model_v0.cif.gz")) and (Exists(f"{paedir}/AF-{value}-F{frag}-predicted_aligned_error_v0.json.gz")):
        if Exists(f"{cifdir}/AF-{value}-F{frag}-model_v0.cif.gz"):
            if Switch('debug'):
                print(f"     >> SKIP (CIF and PAE output files already existed for acc '{value}' fragment {frag})")
            Log("CIF and PAE output files already existed for fragment for acc|frag (skipped)", f"{value}|{frag}")
            continue
        else:
            if Switch('debug'):
                print(f"     >> No output yet for acc '{value}' fragment {frag}, running AlphaFold!")
            # 
            # # Check if output CIF file for this fragment already exists for this acc in cifdir and paedir
            # if Exists(f"{cifdir}/AF-{value}-F{frag}-model_v0.cif.gz"):
            #     Log("output CIF file, but not PAE file already existed for fragment for acc|frag (running anyway)", f"{value}|{frag}")
            # 
            # # Check if output PAE file for this fragment already exists for this acc in cifdir and paedir
            # if Exists(f"{paedir}/AF-{value}-F{frag}-predicted_aligned_error_v0.json.gz"):
            #     if Switch('debug'):
            #         print(f"     >> SKIP (output PAE file already existed for acc '{value}' fragment {frag})")
            #     Log("output PAE file, but not CIF file already existed for fragment for acc|frag (running anyway)", f"{value}|{frag}")

        fragdir = f"{alphafolddir}/{value}_F{frag}"

        if not Switch('debug') and not Switch('debug2'):
            os.makedirs(fragdir, exist_ok=True)
            os.makedirs(f"{fragdir}/fasta", exist_ok=True)

        acc = f"{value}_F{frag}"
        fastafile = f"{fragdir}/fasta/{value}_F{frag}.fasta"
        fragseq = seq
        start = 0
        end = len(fragseq)

        if not Switch('debug') and not Switch('debug2'):
            with open(fastafile, "w") as f:
                f.write(f">{acc}\n{seq}\n") 

        if Switch('debug'):
            print(f"     >> F{frag} >> {acc} >> {start+1}-{end} >> {len(fragseq)} aa")
            print(f"       >> '{fastafile}'")
            print(f"         >> {fragseq}")
        # Log("total fragment sequences", fragseq)

    else:
        for start in range(0, len(seq), 200):

            end = min(start + 1400, len(seq))
            # if f"{value}_F{frag}" == "Q23551_F7" or f"{value}_F{frag}" == "Q9N533_F1":
            #     d()

            Log("total acc|frags", f"{value}|{frag}")

            if f"{tmpvalue}_F{frag}" in running_frags:
                if Switch('debug'):
                    print(f"       >> SKIP: {value}_F{frag} is already running on LSF")
                Log("already running on LSF for acc|frag (skipped)", f"{value}|{frag}")
            else:
                # Check if output CIF and PAE files for this fragment already exist for this acc in cifdir and paedir
                # if (Exists(f"{cifdir}/AF-{value}-F{frag}-model_v0.cif.gz")) and (Exists(f"{paedir}/AF-{value}-F{frag}-predicted_aligned_error_v0.json.gz")):
                if Exists(f"{cifdir}/AF-{value}-F{frag}-model_v0.cif.gz"):
                    if Switch('debug'):
                        print(f"     >> SKIP (CIF output files already existed for acc '{value}' fragment {frag})")
                    Log("CIF output files already existed for fragment for acc|frag (skipped)", f"{value}|{frag}")
                else:
                    if Switch('debug'):
                        print(f"     >> No output yet for acc '{value}' fragment {frag}, running AlphaFold!")
                    #                 
                    # # Check if output CIF file for this fragment already exists for this acc in cifdir and paedir
                    # if Exists(f"{cifdir}/AF-{value}-F{frag}-model_v0.cif.gz"):
                    #     Log("output CIF file, but not PAE file already existed for fragment for acc|frag (running anyway)", f"{value}|{frag}")
                    # 
                    # # Check if output PAE file for this fragment already exists for this acc in cifdir and paedir
                    # if Exists(f"{paedir}/AF-{value}-F{frag}-predicted_aligned_error_v0.json.gz"):
                    #     if Switch('debug'):
                    #         print(f"     >> SKIP (output PAE file already existed for acc '{value}' fragment {frag})")
                    #     Log("output PAE file, but not CIF file already existed for fragment for acc|frag (running anyway)", f"{value}|{frag}")

                    fragdir = f"{alphafolddir}/{value}_F{frag}"
    
                    if not Switch('debug') and not Switch('debug2'):
                        os.makedirs(fragdir, exist_ok=True)
                        os.makedirs(f"{fragdir}/fasta", exist_ok=True)

                    acc = f"{value}_F{frag}"
                    fastafile = f"{fragdir}/fasta/{value}_F{frag}.fasta"

                    fragseq = seq[start:end]

                    # Skip fragments smaller than step size unless it's the last piece (shouldn't happen)
                    if len(fragseq) < 200 and start + 200 < len(seq):
                        Die(f"Error: Fragment smaller than step size for acc '{value}' (fragseq {start}:{end}):\n\n{seq[start:end]}\n\nseq\n\n")
                        # continue

                    if not Switch('debug') and not Switch('debug2'):
                        with open(fastafile, "w") as f:
                            # f.write(f">{acc}\n{fragseq}\n")
                            f.write(f">{acc}_{start+1}_{end}\n{fragseq}\n")

                    if Switch('debug'):
                        print(f"     >> F{frag} >> {acc} >> {start+1}-{end} >> {len(fragseq)} aa")
                        print(f"       >> '{fastafile}'")
                        print(f"         >> {fragseq}")
                    # Log("total fragment sequences", fragseq)

            if end == len(seq):
                break

            frag += 1

    for frag in range(1, maxfrag+1):
    
        if f"{tmpvalue}_F{frag}" in running_frags:
            if Switch('debug'):
                print(f"     >> SKIP: {value}_F{frag} is already running on LSF")
            Log("already running on LSF for acc|frag (skipped)", f"{value}|{frag}")
            continue

        fragdir = f"{alphafolddir}/{value}_F{frag}"

        # Check if output CIF and PAE files for this fragment already exist for this acc in cifdir and paedir
        # if (Exists(f"{cifdir}/AF-{value}-F{frag}-model_v0.cif.gz")) and (Exists(f"{paedir}/AF-{value}-F{frag}-predicted_aligned_error_v0.json.gz")):
        if Exists(f"{cifdir}/AF-{value}-F{frag}-model_v0.cif.gz"):
            if Switch('debug'):
                print(f"     >> SKIP (CIF output files already existed for acc '{value}' fragment {frag})")
            Log("CIF output files already existed for fragment for acc|frag (skipped)", f"{value}|{frag}")
            continue
        else:
            if Switch('debug'):
                print(f"     >> No output yet for acc '{value}' fragment {frag}, running AlphaFold!")
            # 
            # # Check if output CIF file for this fragment already exists for this acc in cifdir and paedir
            # if Exists(f"{cifdir}/AF-{value}-F{frag}-model_v0.cif.gz"):
            #     Log("output CIF file, but not PAE file already existed for fragment for acc|frag (running anyway)", f"{value}|{frag}")
            # 
            # # Check if output PAE file for this fragment already exists for this acc in cifdir and paedir
            # if Exists(f"{paedir}/AF-{value}-F{frag}-predicted_aligned_error_v0.json.gz"):
            #     if Switch('debug'):
            #         print(f"     >> SKIP (output PAE file already existed for acc '{value}' fragment {frag})")
            #     Log("output PAE file, but not CIF file already existed for fragment for acc|frag (running anyway)", f"{value}|{frag}")

        # Run AlphaFold

        # # Choose queue randomly
        # queue = random.choice(queues)
        # Choose queue with the fewest pending jobs
        # for queue in queues:
        #     mypending = Pendingjobs(queue)
        # pending_per_queue = [Pendingjobs(queue) for queue in queues]
        # queue = queues[pending_per_queue.index(min(pending_per_queue))]
        
        # Run locally
        # cmd = f"../../../../alphafold.py {value} {maxfrag}"

        # Check if cpufirst-attempted, "msas" and alphafold_params.json don't exist yet for this fragment (and run CPU-only MSA generation job if none of them do)
        cpu_job_submitted = 0
        # if Switch('cpufirst') and not Exists(f"{fragdir}/{cpufirst_attempted_path}") and not Exists(f"{fragdir}/{value}_F{frag}/msas") and not Exists(f"{fragdir}/{value}_F{frag}/alphafold_params.json"):
        # Instead, check if all expected MSA files are already present, and run CPU job if not:
        if Switch('cpufirst') and (not Exists(f"{fragdir}/{value}_F{frag}/msas/bfd_uniref_hits.a3m") or not Exists(f"{fragdir}/{value}_F{frag}/msas/mgnify_hits.sto") or not Exists(f"{fragdir}/{value}_F{frag}/msas/pdb_hits.hhr") or not Exists(f"{fragdir}/{value}_F{frag}/msas/uniref90_hits.sto")):
            # Run CPU job (CPU-based MSA generation only, switch -cpu active)
            if Switch('debug'):
                print(f"     >> Running CPU-based MSA generation job (-cpufirst is active)")

            # Sleep if total job number is higher than maxjobs
            if myjobs >= maxjobs:
                while myjobs >= maxjobs:
                    myjobs = Myjobs()
                    pending_per_cpu_queue = [Pendingjobs(cpu_queue) for cpu_queue in cpu_queues]
                    cpu_queue = cpu_queues[pending_per_cpu_queue.index(min(pending_per_cpu_queue))]
                    if myjobs >= maxjobs:
                        time.sleep(sleeptime)

            # Sleep if total pending job number is higher than maxpending_cpu
            if mypending_cpu >= maxpending_cpu:
                while mypending_cpu >= maxpending_cpu:
                    # mypending_cpu = Pendingjobs()
                    pending_per_cpu_queue = [Pendingjobs(cpu_queue) for cpu_queue in cpu_queues]
                    cpu_queue = cpu_queues[pending_per_cpu_queue.index(min(pending_per_cpu_queue))]
                    mypending_cpu = Pendingjobs(cpu_queue)
                    if mypending_cpu >= maxpending_cpu:
                        time.sleep(sleeptime)

            # Submit job (4 cores, 45G RAM each = 180G, standard queue - just to ensure that we have enough memory)
            cmd = f"""bsub -P idr -J "update_alphasync_input_alphasync_alphafold_{tmpvalue}_F{frag}_____________alphafold_py_{tmpvalue}_{frag}_{maxfrag}__cpu" -L /bin/bash -env 'LSB_JOB_REPORT_MAIL=N' -n 4 -q {cpu_queue} -R "span[hosts=1]" -R "rusage[mem=45G]" "bash -c '../../../../alphafold.py {value} {frag} {maxfrag} -cpu > log-output-update_alphasync_input_alphasync_alphafold_{tmpvalue}_F{frag}_____________alphafold_py_{tmpvalue}_{frag}_{maxfrag}__cpu.txt 2> log-errors-update_alphasync_input_alphasync_alphafold_{tmpvalue}_F{frag}_____________alphafold_py_{tmpvalue}_{frag}_{maxfrag}__cpu.txt'" -e /dev/null -o /dev/null > /dev/null"""

            # To see progress of -cpu jobs:
            # input/alphasync/alphafold >> 1 */*/alphafold_params.json */*/msas/*|g -o "[^/]+$"|suq            

            if not Switch('debug') and not Switch('debug2'):
                os.chdir(fragdir)
                Run(f"Submit job (which runs AlphaFold and processes its output) (CPU-based MSA generation only, switch -cpu active)", cmd, silent=True)
                os.chdir(maindir)
            else:
                if not Switch('debug2'):
                    State(f"""Submit job (which runs AlphaFold and processes its output) (CPU-based MSA generation only, switch -cpu active): {cmd}""")

            myjobs += 1
            mypending_cpu += 1
            submitted_cpu += 1
            pending_per_cpu_queue[cpu_queues.index(cpu_queue)] += 1
            cpu_queue = cpu_queues[pending_per_cpu_queue.index(min(pending_per_cpu_queue))]

            cpu_job_submitted = 1

            Log(f"submitted CPU-only job for acc|frag", f"{value}|{frag}")
            Log(f"submitted CPU-only job for acc", value)



        # Run GPU job (using existing MSAs)

        if not Switch('debug') and not Switch('debug2'):
            # Sleep if total job number is higher than maxjobs
            if myjobs >= maxjobs:
                while myjobs >= maxjobs:
                    myjobs = Myjobs()
                    pending_per_queue = [Pendingjobs(queue) for queue in queues]
                    queue = queues[pending_per_queue.index(min(pending_per_queue))]
                    if myjobs >= maxjobs:
                        time.sleep(sleeptime)

            # Sleep if total pending job number is higher than maxpending_gpu
            if mypending >= maxpending_gpu:
                while mypending >= maxpending_gpu:
                    # mypending = Pendingjobs()
                    pending_per_queue = [Pendingjobs(queue) for queue in queues]
                    queue = queues[pending_per_queue.index(min(pending_per_queue))]
                    mypending = Pendingjobs(queue)
                    if mypending >= maxpending_gpu:
                        time.sleep(sleeptime)

        # Submit job (A100 specs: 128 cores (16 per GPU), 1500 GB RAM (~187 GB per GPU), 8 GPUs)
        if cpu_job_submitted == 0:
            # GPU only: simply submit GPU job
            if Switch('debug'):
                print(f"     >> Running final GPU-based job (no dependency)")
            # 1 GPU, 4 CPUs, 4*30 = 120 GB RAM
            cmd = f"""bsub -P idr -J "update_alphasync_input_alphasync_alphafold_{tmpvalue}_F{frag}_____________alphafold_py_{tmpvalue}_{frag}_{maxfrag}" -L /bin/bash -env 'LSB_JOB_REPORT_MAIL=N' -n 4 -q {queue} -gpu "num=1/host" -R "span[hosts=1]" -R "rusage[mem=30G]" "bash -c '../../../../alphafold.py {value} {frag} {maxfrag} > log-output-update_alphasync_input_alphasync_alphafold_{tmpvalue}_F{frag}_____________alphafold_py_{tmpvalue}_{frag}_{maxfrag}.txt 2> log-errors-update_alphasync_input_alphasync_alphafold_{tmpvalue}_F{frag}_____________alphafold_py_{tmpvalue}_{frag}_{maxfrag}.txt'" -e /dev/null -o /dev/null > /dev/null"""
        else:
            # -cpufirst: wait for CPU job to finish before running GPU job (using bsub -w)
            if Switch('debug'):
                print(f"     >> Running final GPU-based job (dependent on CPU job finishing first)")
            # 1 GPU, 4 CPUs, 4*30 = 120 GB RAM
            cmd = f"""bsub -P idr -w "ended(update_alphasync_input_alphasync_alphafold_{tmpvalue}_F{frag}_____________alphafold_py_{tmpvalue}_{frag}_{maxfrag}__cpu)" -J "update_alphasync_input_alphasync_alphafold_{tmpvalue}_F{frag}_____________alphafold_py_{tmpvalue}_{frag}_{maxfrag}" -L /bin/bash -env 'LSB_JOB_REPORT_MAIL=N' -n 4 -q {queue} -gpu "num=1/host" -R "span[hosts=1]" -R "rusage[mem=30G]" "bash -c '../../../../alphafold.py {value} {frag} {maxfrag} > log-output-update_alphasync_input_alphasync_alphafold_{tmpvalue}_F{frag}_____________alphafold_py_{tmpvalue}_{frag}_{maxfrag}.txt 2> log-errors-update_alphasync_input_alphasync_alphafold_{tmpvalue}_F{frag}_____________alphafold_py_{tmpvalue}_{frag}_{maxfrag}.txt'" -e /dev/null -o /dev/null > /dev/null"""
    
        if not Switch('debug') and not Switch('debug2'):
            os.chdir(fragdir)
            Run(f"Submit job (which runs AlphaFold and processes its output)", cmd, silent=True)
            os.chdir(maindir)
        else:
            if not Switch('debug2'):
                State(f""" >> Submit job (which runs AlphaFold and processes its output) >> {cmd}""")

        myjobs += 1
        mypending += 1
        submitted += 1
        pending_per_queue[queues.index(queue)] += 1
        queue = queues[pending_per_queue.index(min(pending_per_queue))]

        Log(f"submitted GPU job for acc|frag", f"{value}|{frag}")
        Log(f"submitted GPU job for acc", value)

        # Logging without extra annotation
        Log(f"submitted job for acc|frag", f"{value}|{frag}")
        Log(f"submitted job for acc", value)
        # Log(f"submitted job for full sequence", seq)

    #     break
    # break
    
print()
Show("submitted job for acc|frag")
Show("submitted job for acc")
Show(lim=200, sort=True)

if not Switch('cpufirst'):
    State(f" >> Submitted {submitted:,} jobs")
else:
    State(f" >> Submitted {submitted_cpu:,} CPU-only jobs (will run first)")
    State(f" >> Submitted {submitted:,} GPU jobs (will run once the associated CPU job completes)")

# Wait for jobs to finish
if not Switch('debug') and not Switch('debug2'):
    if Switch('cpufirst'):
        for cpu_queue in cpu_queues:
            Waitforjobs(cpu_queue)
    for queue in queues:
        Waitforjobs(queue)

print("\nDone!")
