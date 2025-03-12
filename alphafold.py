#!/usr/bin/env python3
"""
alphafold.py: AlphaFold job script that runs AlphaFold on individual fragment files.

Copies CIF files to their final destination:    input/alphasync/cif
Copies PAE files to their final destination:    input/alphasync/pae
Copies params files to their final destination: input/alphasync/params

Note: Using "v0" as the "AlphaFold DB release version" for AlphaSync structures.
"""

# Initialize
# from blang_mysql import *
from blang import *

# Notes
# Check recent parameters used (-1 day):
# input/alphasync/params >> find . -type f -mtime -1 | xa cat {} | suq

# Versions

alphafold_version = "2.3.2"

# Versions of hhblits to use (sequentially, until successful)
# hhblits_versions = ['3.3.0_high_limits', '3.3.0', '3.3.0_low_limits', '3.0-beta.3_high_limits', '3.0-beta.3', '3.0-beta.3_low_limits']
# hhblits_versions = ['3.3.0_high_limits', '3.3.0_low_limits', '3.0-beta.3_high_limits', '3.0-beta.3_low_limits']
# high_limits leads to high RAM usage (hitting a 192 GB limit for Q9UJX3_F1), hence using high_cpus instead (which only increases CPUs from 4 to 12 from the default):
# high_cpu (12 CPUs max) and low_limits (4 CPUs, lower max seq numbers)
# hhblits_versions = ['3.3.0_high_cpus', '3.3.0_low_limits', '3.0-beta.3_high_cpus', '3.0-beta.3_low_limits']
# Default (4 CPUs max) and low_limits (4 CPUs, lower max seq numbers)
hhblits_versions = ['3.3.0', '3.3.0_low_limits', '3.0-beta.3', '3.0-beta.3_low_limits']
# hhblits_versions = ['3.0-beta.3']

# Versions of uniref30 to use (sequentially, until successful)
uniref30_versions = ['2023_02', '2022_02', '2021_03', '2020_06', '2018_08']
# uniref30_versions = ['2023_02']

# Versions of other datasets
# PDB template cutoff
# max_template_date = "$(date +%F)"
max_template_date = "2024-11-25"
# BFD
bfd_version = "bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt"
# UniRef90
uniref90_version = "2024_05"
# MGnify clusters
mgnify_version = "2024_04"
# # pdb70/pdb100
pdb_types = ['pdb100', 'pdb70']
# # pdb70 2022_03_13
# pdb_type = "pdb70"
# pdb_version = "pdb70_from_mmcif_220313"
# Use pdb100 in place of pdb70, as ColabFold does:
# # pdb100 2023_05_17
# pdb_type = "pdb100"
# pdb_version = "pdb100_foldseek_230517"

# Version of jackhmmer.py to use
# # Default: built-in, 8 CPUs
# jackhmmer_version = []
# 4 CPUs instead
jackhmmer_version = ["-B ../../alphafold_tools/jackhmmer_4_cpus.py:/app/alphafold/alphafold/data/tools/jackhmmer.py",]

# Version of run_alphafold.py to use
# Default: built-in
run_alphafold_py_version = []

cpufirst_attempted = "cpufirst-attempted"



(acc, frag, maxfrag) = Args(3, "[UniProt accession] [Fragment number] [Total fragments] [-cpu]", " -cpu: Use 4 CPUs for jackhmmer and hhblits, and stop after MSA generation (before the GPU part)", "A0A087WUL8 1 14")

cifdir = f"../../cif"
paedir = f"../../pae"
paramdir = f"../../params"

if Switch('cpu'):

    # Use CPU only and stop after MSA generation (before the GPU part)

    print("Switch active: -cpu")
    Run(f"Create '{cpufirst_attempted}' file", f"touch {cpufirst_attempted}")

    # hhblits_versions = ['3.3.0', '3.0-beta.3']
    # # high_cpu (12 CPUs max) and low_limits (4 CPUs, lower max seq numbers)
    # hhblits_versions = ['3.3.0_high_cpus', '3.3.0_low_limits', '3.0-beta.3_high_cpus', '3.0-beta.3_low_limits']
    # Default (4 CPUs max) and low_limits (4 CPUs, lower max seq numbers)
    hhblits_versions = ['3.3.0', '3.3.0_low_limits', '3.0-beta.3', '3.0-beta.3_low_limits']

    # 4 CPUs for jackhmmer
    jackhmmer_version = ["-B ../../alphafold_tools/jackhmmer_4_cpus.py:/app/alphafold/alphafold/data/tools/jackhmmer.py",]
    # # 8 CPUs for jackhmmer (default)
    # jackhmmer_version = []
    # 12 CPUs for jackhmmer
    # jackhmmer_version = ["-B ../../alphafold_tools/jackhmmer_12_cpus.py:/app/alphafold/alphafold/data/tools/jackhmmer.py",]

    run_alphafold_py_version = ["-B ../../alphafold_tools/run_alphafold_msas_only.py:/app/alphafold/run_alphafold.py",]

# Set SINGULARITY_TMPDIR
# This may require as much as 71 GB of free space (e.g. for Q5SSE9_F19's jackhmmer run against uniref90_2024_05)
# tmpdir = "tmp"
user = os.environ['USER']
singularity_tmpdir = f"/scratch_local/{user}"

# Set AlphaFold output directory
outdir = f"{acc}_F{frag}"

# To remove giant .pkl files:
# find . \( -name "result_model_*.pkl" -o -name "features*.pkl" \) -exec rm -fv {} +

# # Check if output CIF and PAE files (for at least one fragment) already exist for this acc in cifdir and paedir
# # e.g. AF-{acc}-F{frag}-model_v0.cif.gz
# if glob(f"{cifdir}/AF-{acc}-F*-model_v0.cif.gz") or glob(f"{paedir}/*/AF-{acc}-F*-predicted_aligned_error_v0.json.gz"):
#     Die(f"Error: Data already exists for acc '{acc}' in directory '{cifdir}' or '{paedir}', exiting (skip)!")
# else:
#     print(f"No data yet for acc '{acc}' in directories '{cifdir}' and '{paedir}', starting!")

# Check if complete output CIF and PAE files for this fragment already exist for this acc in cifdir and paedir
# e.g. AF-{acc}-F{frag}-model_v0.cif.gz
# if (not Exists(f"{cifdir}/AF-{acc}-F{frag}-model_v0.cif.gz") and not Exists(f"{paedir}/AF-{acc}-F{frag}-predicted_aligned_error_v0.json.gz")):
if not Exists(f"{cifdir}/AF-{acc}-F{frag}-model_v0.cif.gz"):
    # No data yet
    print(f"No data yet for acc '{acc}' fragment '{frag}' in directories '{cifdir}' and '{paedir}', starting!")
# elif (Exists(f"{cifdir}/AF-{acc}-F{frag}-model_v0.cif.gz") and Exists(f"{paedir}/AF-{acc}-F{frag}-predicted_aligned_error_v0.json.gz")):
elif Exists(f"{cifdir}/AF-{acc}-F{frag}-model_v0.cif.gz"):
    # Data already exists
    print(f"Data already exists for acc '{acc}' fragment '{frag}' in directory '{cifdir}' and '{paedir}', exiting (skip)!")
    sys.exit()
# else:
#     # Incomplete data (shouldn't happen)
#     Warn(f"Incomplete data for acc '{acc}' fragment '{frag}' in directories '{cifdir}' and '{paedir}', starting!")
# 
#     # Check if output CIF file for this fragment already exists for this acc in cifdir and paedir
#     if Exists(f"{cifdir}/AF-{acc}-F{frag}-model_v0.cif.gz"):
#         Log("output CIF file, but no PAE file existed for fragment for acc|frag (running)", f"{acc}|{frag}")
# 
#     # Check if output PAE file for this fragment already exists for this acc in cifdir and paedir
#     if Exists(f"{paedir}/AF-{acc}-F{frag}-predicted_aligned_error_v0.json.gz"):
#         if Switch('debug'):
#             print(f"     >> SKIP (output PAE file already existed for acc '{acc}' fragment {frag})")
#         Log("output PAE file, but no CIF file existed for fragment for acc|frag (running)", f"{acc}|{frag}")



# Start

Starttime()

# Lists of temporary files and directories to clean up later on
tmpfiles = []
tmpdirs = []

# # Create temporary directory for Singularity container
# os.makedirs(singularity_tmpdir, exist_ok=True)
# tmpdirs.append(singularity_tmpdir)

# Get list of (fragment) FASTA files that are present in this temporary directory
fastafile = f"fasta/{outdir}.fasta"
if not Exists(fastafile):
    Die(f"Error: Expected to find FASTA file '{fastafile}' for acc '{acc}' fragment {frag}, but it doesn't exist")

# First, check if all fragments are present, and if we can parse the file names correctly
Log("temporary FASTA files removed", fastafile)
tmpfiles.append(fastafile)
tmpdirs.append("fasta")

# Get fragment number for this mmCIF file
m = rx(r"^fasta/"+acc+r"_F(\d+)\.fasta$", fastafile)
if m:
    frag = m[0]
else:
    Die(f"Couldn't parse filename '{fastafile}'")

# Run AlphaFold
if not Switch('cpu'):
    # GPU
    State("\nRunning AlphaFold:")
else:
    # CPU-only
    State("\nRunning AlphaFold (CPU only: switch -cpu active):")

print(f" >> {fastafile}")

# Run AlphaFold
# print(f"   >> Running AlphaFold...")
# Run("AlphaFold", f"../../../../alphafold.sh {fastafile} {outdir}", silent=False)



if maxfrag > 1:
    # Multi-fragment proteins: Only run monomer (PAE scores would not be meaningful)
    models = ['monomer']
else:
    # Single-fragment proteins: Run both monomer_ptm (which produces PAE scores) and monomer
    models = ['monomer_ptm', 'monomer']


for model in models:

    if not Exists(f"{outdir}/ranked_0.cif"):

#         # Downgrade hhblits first
# 
#         # Try hhblits 3.3.0 (latest), then try higher limits for maxseq, maxfilt and realign_max, then lower limits, then downgrade to 3.0-beta.3 (the version used in all published AlphaFold papers, including AlphaFold 3)
#         # for hhblits_version in ['3.3.0_high_limits', '3.3.0', '3.3.0_low_limits', '3.0-beta.3_high_limits', '3.0-beta.3', '3.0-beta.3_low_limits']:
#         for hhblits_version in ['3.3.0_high_limits', '3.3.0_low_limits', '3.0-beta.3_high_limits', '3.0-beta.3_low_limits']:
#         # for hhblits_version in ['3.0-beta.3']:
# 
#             # Try most recent version of uniref30/uniclust30 first, then keep downgrading until one works
#             for uniref30_version in ['2023_02', '2022_02', '2021_03', '2020_06', '2018_08']:
#             # for uniref30_version in ['2023_02']:

        # Try pdb100, then pdb70 (pdb100 only seems to fail for P9WMK7, due to a missing 5cvx.cif file, which seems to be an odd issue in obsolete.dat where the ID got recycled)
        for pdb_type in pdb_types:

            if pdb_type == 'pdb100':
                pdb_version = "pdb100_foldseek_230517"
            elif pdb_type == 'pdb70':
                pdb_version = "pdb70_from_mmcif_220313"
            else:
                Die(f"Unhandled pdb_type '{pdb_type}'")
        
            # Downgrade UniRef30 first

            # Try most recent version of uniref30/uniclust30 first, then keep downgrading until one works
            for uniref30_version in uniref30_versions:

                # Try hhblits 3.3.0 (latest), then try higher limits for maxseq, maxfilt and realign_max, then lower limits, then downgrade to 3.0-beta.3 (the version used in all published AlphaFold papers, including AlphaFold 3)
                for hhblits_version in hhblits_versions:
                        
                    # if model == 'monomer_ptm':
                    if model == models[0]:
                        # If this is the first model:
                        print(f"\n\n\n   >> Running AlphaFold with model '{model}' using hhblits '{hhblits_version}' and uniref30 '{uniref30_version}' and {pdb_type} '{pdb_version}'...")
                    else:
                        print(f"\n\n\n   >> Running AlphaFold with model '{model}' and existing alignments...")

                    if uniref30_version in ['2023_02', '2022_02']:
                        uniref30 = f"/data/alphafold_datasets_latest/uniref30/UniRef30_{uniref30_version}"
                    elif uniref30_version in ['2021_03', '2020_06']:
                        uniref30 = f"/data/alphafold_data/uniref30/UniRef30_{uniref30_version}"
                    elif uniref30_version in ['2018_08']:
                        uniref30 = f"/data/alphafold_data/uniclust30/uniclust30_{uniref30_version}/uniclust30_{uniref30_version}"
                    else:
                        Die(f"Unhandled uniref30_version '{uniref30_version}'")

                    singularity_args = [
                        f"module load alphafold/{alphafold_version} 2>&1;",
                        f"echo;",
                        f"hostname;",
                        f"echo LD_LIBRARY_PATH=$LD_LIBRARY_PATH;",
                        f"echo CUDA_LIBRARY_PATH=$CUDA_LIBRARY_PATH;",
                        # Exit if the path specified in CUDA_LIBRARY_PATH doesn't exist, and print a warning
                        f"if [ ! -d $CUDA_LIBRARY_PATH ]; then echo 'Warning: CUDA_LIBRARY_PATH does not exist'; exit 1; fi;",
                        f"echo;",
                        # f"SINGULARITY_TMPDIR={singularity_tmpdir}",
                        # f"TMPDIR={singularity_tmpdir}",
                        "singularity exec --no-home",
                        "--env LD_LIBRARY_PATH='$LD_LIBRARY_PATH:$CUDA_LIBRARY_PATH:/opt/conda/lib'",
                        "-B .:/etc",
                        "-B $CUDA_LIBRARY_PATH:$CUDA_LIBRARY_PATH",
                        # Note: $DATA_DIR is already set to /lustre_scratch/reference/public/alphafold_data, no need to set this manually
                        "-B $DATA_DIR:/data/alphafold_data",
                        f"-B {singularity_tmpdir}:/tmp",
                        "-B ../../../../input/alphasync/alphafold_datasets:/data/alphafold_datasets_latest",
                        ]

                    # hhblits version
                    if hhblits_version.startswith('3.3.0'):
                        singularity_args += [
                            # "-B ~/bin/hhblits:/usr/bin/hhblits",
                        ]
                    elif hhblits_version.startswith('3.0-beta.3'):
                        singularity_args += [
                            "-B ../../alphafold_tools/hhsuite-3.0-beta.3-Linux/bin/hhblits:/usr/bin/hhblits",
                            "-B ../../alphafold_tools/hhsuite-3.0-beta.3-Linux:/app/hhlib",
                            # # Use high cpus (cpu=12)
                            # "-B ../../alphafold_tools/hhblits_high_cpus.py:/app/alphafold/alphafold/data/tools/hhblits.py",
                            "--env HHLIB='/app/hhlib'",
                        ]
                    else:
                        Die(f"Unhandled hhblits_version '{hhblits_version}'")

                    # hhblits parameters (high_limits, high_cpus, low_limits)
                    if hhblits_version.endswith('_high_limits'):
                        hhblits_limits = "high_limits"
                        singularity_args += [
                            # Use high limits (cpu=12, maxseq=1_000_000, maxfilt=1_000_000, realign_max=1_000_000)
                            "-B ../../alphafold_tools/hhblits_high_limits.py:/app/alphafold/alphafold/data/tools/hhblits.py",
                        ]
                    elif hhblits_version.endswith('_high_cpus'):
                        hhblits_limits = "high_cpus"
                        singularity_args += [
                            # Use high limits (cpu=12, maxseq=1_000_000, maxfilt=100_000, realign_max=100_000)
                            "-B ../../alphafold_tools/hhblits_high_cpus.py:/app/alphafold/alphafold/data/tools/hhblits.py",
                        ]
                    elif hhblits_version.endswith('_low_limits'):
                        hhblits_limits = "low_limits"
                        singularity_args += [
                            # Use low limits (cpu=4 (default), maxseq=65535 (default, and as was still hardcoded in 3.0-beta.3), maxfilt=100_000 (default), realign_max=100_000 (default))
                            "-B ../../alphafold_tools/hhblits_low_limits.py:/app/alphafold/alphafold/data/tools/hhblits.py",
                        ]
                    else:
                        # Use default limits (cpu=4, maxseq=1_000_000, maxfilt=100_000, realign_max=100_000)
                        hhblits_limits = "default"

                    # Set jackhmmer_version (empty for default, special 4-CPU version if switch '-cpu' is active)
                    singularity_args += jackhmmer_version
                    
                    # Set run_alphafold.py version (empty for default, special MSA-only (CPU-only) version if switch '-cpu' is active)
                    singularity_args += run_alphafold_py_version
                    
                    singularity_args += [
                        "--nv",
                        "$SIF",
                        "python /app/alphafold/run_alphafold.py",
                        ]

                    alphafold_args = [
                        "--model_preset=" + model,
                        "--fasta_paths='/etc/" + fastafile + "'",
                        "--output_dir='/etc'",
                        "--max_template_date='" + max_template_date + "'",
                        "--use_precomputed_msas=True",
                        # "--db_preset=reduced_dbs",
                        # "--small_bfd_database_path=/data/alphafold_datasets_latest/small_bfd/bfd-first_non_consensus_sequences.fasta",
                        "--uniref90_database_path=/data/alphafold_datasets_latest/uniref90/uniref90_" + uniref90_version + ".fasta",
                        "--mgnify_database_path=/data/alphafold_datasets_latest/mgnify/mgy_clusters_" + mgnify_version + ".fa",
                        "--bfd_database_path=/data/alphafold_data/bfd/" + bfd_version,
                        # "--uniref30_database_path=/data/alphafold_datasets_latest/uniref30/UniRef30_2023_02",
                        # "--uniref30_database_path=/data/alphafold_data/uniclust30/uniclust30_2018_08/uniclust30_2018_08",
                        "--uniref30_database_path=" + uniref30,
                        "--pdb70_database_path=/data/alphafold_datasets_latest/" + pdb_version + "/" + pdb_type,
                        "--template_mmcif_dir=/data/alphafold_datasets_latest/pdb_mmcif/mmcif_files",
                        "--obsolete_pdbs_path=/data/alphafold_datasets_latest/pdb_mmcif/obsolete.dat",
                        "2>&1",
                        ]

                    Starttime(2)
                    Run("AlphaFold", " ".join(singularity_args + alphafold_args))
                    print()
                    print("     >> ", end="")
                    Stoptime(2)

                    # Check if output CIF file exists (if the run worked with this version of uniref30/uniclust30), or if all four MSA output files are present for a CPU-only run (bfd_uniref_hits.a3m, mgnify_hits.sto, pdb_hits.hhr, uniref90_hits.sto)
                    if Exists(f"{outdir}/ranked_0.cif") or (Switch('cpu') and Exists(f"{outdir}/msas/bfd_uniref_hits.a3m") and Exists(f"{outdir}/msas/mgnify_hits.sto") and Exists(f"{outdir}/msas/pdb_hits.hhr") and Exists(f"{outdir}/msas/uniref90_hits.sto")):
                        break
                # Check if output CIF file exists (if the run worked with this version of uniref30/uniclust30), or if all four MSA output files are present for a CPU-only run (bfd_uniref_hits.a3m, mgnify_hits.sto, pdb_hits.hhr, uniref90_hits.sto)
                if Exists(f"{outdir}/ranked_0.cif") or (Switch('cpu') and Exists(f"{outdir}/msas/bfd_uniref_hits.a3m") and Exists(f"{outdir}/msas/mgnify_hits.sto") and Exists(f"{outdir}/msas/pdb_hits.hhr") and Exists(f"{outdir}/msas/uniref90_hits.sto")):
                    break
            # Check if output CIF file exists (if the run worked with this version of uniref30/uniclust30), or if all four MSA output files are present for a CPU-only run (bfd_uniref_hits.a3m, mgnify_hits.sto, pdb_hits.hhr, uniref90_hits.sto)
            if Exists(f"{outdir}/ranked_0.cif") or (Switch('cpu') and Exists(f"{outdir}/msas/bfd_uniref_hits.a3m") and Exists(f"{outdir}/msas/mgnify_hits.sto") and Exists(f"{outdir}/msas/pdb_hits.hhr") and Exists(f"{outdir}/msas/uniref90_hits.sto")):
                break
    else:
        # print(f"   >> Skipping AlphaFold with model '{model}' since '{outdir}/ranked_0.cif' already exists")
        print(f"   >> Skipping AlphaFold run since '{outdir}/ranked_0.cif' already exists")
        sys.exit()

    if not Switch('cpu'):
        # GPU
        if not Exists(f"{outdir}/ranked_0.cif"):
            Die(f"Error: Expected AlphaFold output CIF file '{outdir}/ranked_0.cif' doesn't exist (all parameter combinations failed)")
    else:
        # CPU only
        if not Exists(f"{outdir}/msas/bfd_uniref_hits.a3m") or not Exists(f"{outdir}/msas/mgnify_hits.sto") or not Exists(f"{outdir}/msas/pdb_hits.hhr") or not Exists(f"{outdir}/msas/uniref90_hits.sto"):
            Die(f"Error: Expected AlphaFold output MSA files are incomplete in '{outdir}/msas' (all parameter combinations failed)")



    # Store the successful parameter combination in a JSON file in the output directory
    # Only for monomer_ptm, since monomer will simply load the existing alignment files and give incorrect parameters here
    # if model == 'monomer_ptm':
    # # If this is the first model:
    # if model == models[0]:
    paramfile = f"{outdir}/alphafold_params.json"
    # If this is the first time we're writing parameters:
    if not Exists(paramfile):
        with open(paramfile, "w") as f:
            json.dump({
                "alphafold_version": alphafold_version,
                "hhblits_version": hhblits_version.replace('_high_limits', '').replace('_high_cpus', '').replace('_low_limits', ''),
                "hhblits_limits": hhblits_limits,
                "uniref30_version": uniref30_version,
                "bfd_version": bfd_version,
                "uniref90_version": uniref90_version,
                "mgnify_version": mgnify_version,
                "pdb_type": pdb_type,
                "pdb_version": pdb_version,
                "max_template_date": max_template_date,
            }, f, indent=4)
            # Final newline
            f.write("\n")

    if not paramfile in tmpfiles:
        tmpfiles.append(paramfile)
    
    if Switch('cpu'):
        print("\n\n >> Exiting before GPU-only part (switch -cpu active)\n\nDone!\n\n")
        sys.exit()


    # Mark temporary files and directories for cleanup
    if not f"{outdir}/msas" in tmpdirs:
        tmpdirs.append(f"{outdir}/msas")
    if not outdir in tmpdirs:
        tmpdirs.append(outdir)
    for file in ["features.pkl", "ranked_0.cif", "ranked_0.pdb", "ranked_1.cif", "ranked_1.pdb", "ranked_2.cif", "ranked_2.pdb", "ranked_3.cif", "ranked_3.pdb", "ranked_4.cif", "ranked_4.pdb", "ranking_debug.json", "relax_metrics.json", "timings.json"]:
        file = f"{outdir}/{file}"
        if not file in tmpfiles:
            tmpfiles.append(file)

    # Rename some 'monomer_ptm' output files so they don't get overwritten by 'monomer'
    if model == 'monomer_ptm':
        suffix = "_ptm"
        for file in ["features.pkl", "ranked_0.cif", "ranked_0.pdb", "ranked_1.cif", "ranked_1.pdb", "ranked_2.cif", "ranked_2.pdb", "ranked_3.cif", "ranked_3.pdb", "ranked_4.cif", "ranked_4.pdb", "ranking_debug.json", "relax_metrics.json", "timings.json"]:
            file = f"{outdir}/{file}"
            base, ext = os.path.splitext(file)
            print(f"     >> Renaming '{file}' to '{base}{suffix}{ext}'")
            os.rename(file, f"{base}{suffix}{ext}")
            if not f"{base}{suffix}{ext}" in tmpfiles:
                tmpfiles.append(f"{base}{suffix}{ext}")




# Process AlphaFold output
State(f"Processing AlphaFold output (to find best model etc.) and copying final CIF, PAE and parameter files to '{cifdir}', '{paedir} and {paramdir}':")

# Final file names
finalcif = f"{cifdir}/AF-{acc}-F{frag}-model_v0.cif.gz"
finalpae = f"{paedir}/AF-{acc}-F{frag}-predicted_aligned_error_v0.json.gz"
finalparam = f"{paramdir}/AF-{acc}-F{frag}-alphafold_params.json"


# monomer_ptm
# Just for removing temporary files:
if 'monomer_ptm' in models:
    # Get top model number (1-5) from "ranking_debug_ptm.json"
    ranking = json.load(open(f"{outdir}/ranking_debug_ptm.json"))
    model = ranking["order"][0].split("_")[1]

    # Look for output CIF file
    ciffile = f"{outdir}/relaxed_model_{model}_ptm_pred_0.cif"
    if not Exists(ciffile):
        Die(f"Expected AlphaFold output CIF file '{ciffile}' doesn't exist")
    tmpfiles.append(ciffile)

    # Look for output PDB file
    pdbfile = f"{outdir}/relaxed_model_{model}_ptm_pred_0.pdb"
    if not Exists(pdbfile):
        Die(f"Expected AlphaFold output PDB file '{pdbfile}' doesn't exist")
    tmpfiles.append(pdbfile)



# monomer

# Get top model number (1-5) from "ranking_debug.json"
ranking = json.load(open(f"{outdir}/ranking_debug.json"))
model = ranking["order"][0].split("_")[1]

# Look for output CIF file
# ciffile = f"{outdir}/ranked_0.cif" # Gets overwritten by monomer_ptm
ciffile = f"{outdir}/relaxed_model_{model}_pred_0.cif"
if not Exists(ciffile):
    Die(f"Expected AlphaFold output CIF file '{ciffile}' doesn't exist")
tmpfiles.append(ciffile)

# Look for output PDB file
pdbfile = f"{outdir}/relaxed_model_{model}_pred_0.pdb"
if not Exists(pdbfile):
    Die(f"Expected AlphaFold output PDB file '{pdbfile}' doesn't exist")
tmpfiles.append(pdbfile)

# Compress with gzip (keeping original file as well)
Run("gzip", f"gzip -kf {ciffile}", silent=True)

# Copy CIF file to destination
# e.g. AF-A0A0A0MRZ7-F1-model_v0.cif.gz
print(f" >> Copy '{ciffile}.gz' >> '{finalcif}'")
shutil.copy(f"{ciffile}.gz", finalcif)
tmpfiles.append(f"{ciffile}.gz")



# Copy PAE files to destination (only if maxfrag == 1)
# 
# To delete files from multi-fragment proteins:
# input/alphasync/pae >> 1|g '\-F2-'|g -o "^AF-(\w+(-\d+)?)-F"|perl -npe 's/^AF-//; s/-F$//'|sort|uniq|xa sh -c "ls -1U AF-{}-F*"|xa rm -fv {}
# 
# - Get top model number (1-5) from "ranking_debug.json"
# {
#     "plddts": {
#         "model_1_pred_0": 60.15109852333826,
#         "model_2_pred_0": 56.950351625898364,
#         "model_3_pred_0": 55.12731358164582,
#         "model_4_pred_0": 58.30320976316117,
#         "model_5_pred_0": 57.839050144204144
#     },
#     "order": [
#         "model_1_pred_0",
#         "model_4_pred_0",
#         "model_5_pred_0",
#         "model_2_pred_0",
#         "model_3_pred_0"
#     ]
# }
# - Compress e.g. "pae_model_1_ptm_pred_0.json" and copy it to "{cifdir}/AF-O05625-F1-predicted_aligned_error_v0.json.gz"

if (maxfrag == 1):

    # Look for output PAE file
    paefile = f"{outdir}/pae_model_{model}_ptm_pred_0.json"
    if not Exists(paefile):
        Die(f"Expected AlphaFold output PAE file '{paefile}' doesn't exist")

    # Compress with gzip (keeping original file as well)
    Run("gzip", f"gzip -kf {paefile}", silent=True)

    # Copy PAE file to destination
    # e.g. AF-A0A0A0MRZ7-F1-predicted_aligned_error_v0.json.gz

    # # To avoid paedir containing too many files:
    # # Use last 3 characters of accession (ignoring isoform number, e.g. only use P20929 from acc='P20929-1')
    # # This should work well: SELECT SUBSTRING(acc, -3), COUNT(DISTINCT acc) AS accs FROM alphaseq GROUP BY SUBSTRING(acc, -3) ORDER BY accs DESC;
    # # This leads to around 300 accs per subdir, which is great.
    # accsubdir = acc.split('-')[0][-3:]
    # os.makedirs(f"{paedir}/{accsubdir}", exist_ok=True)
    # print(f" >> Copy '{paefile}.gz' >> '{paedir}/{accsubdir}/AF-{acc}-F{frag}-predicted_aligned_error_v0.json.gz'")
    # shutil.copy(f"{paefile}.gz", f"{paedir}/{accsubdir}/AF-{acc}-F{frag}-predicted_aligned_error_v0.json.gz")

    print(f" >> Copy '{paefile}.gz' >> '{finalpae}'")
    shutil.copy(f"{paefile}.gz", finalpae)
    tmpfiles.append(f"{paefile}.gz")
else:
    print(f" >> Not copying PAE file to output destination since this protein has multiple fragments (PAE scores would be unreliable)")




# Look for output param file
paramfile = f"{outdir}/alphafold_params.json"
if not Exists(paramfile):
    Die(f"Expected AlphaFold output param file '{paramfile}' doesn't exist")

# # Compress with gzip (keeping original file as well)
# Run("gzip", f"gzip -kf {paramfile}", silent=True)

# Copy param file to destination
# e.g. AF-A0A0A0MRZ7-F1-predicted_aligned_error_v0.json.gz

# Not used since it's not actually that many files:
# # To avoid paramdir containing too many files:
# # Use last 3 characters of accession (ignoring isoform number, e.g. only use P20929 from acc='P20929-1')
# # This should work well: SELECT SUBSTRING(acc, -3), COUNT(DISTINCT acc) AS accs FROM alphaseq GROUP BY SUBSTRING(acc, -3) ORDER BY accs DESC;
# # This leads to around 300 accs per subdir, which is great.
# accsubdir = acc.split('-')[0][-3:]
# os.makedirs(f"{paramdir}/{accsubdir}", exist_ok=True)
# print(f" >> Copy '{paramfile}.gz' >> '{paramdir}/{accsubdir}/AF-{acc}-F{frag}-predicted_aligned_error_v0.json.gz'")
# shutil.copy(f"{paramfile}.gz", f"{paramdir}/{accsubdir}/AF-{acc}-F{frag}-predicted_aligned_error_v0.json.gz")

print(f" >> Copy '{paramfile}' >> '{finalparam}'")
shutil.copy(paramfile, finalparam)



# Add additional output files to be cleaned up

# MSAs
tmpfiles.append(f"{outdir}/msas/bfd_uniref_hits.a3m")
tmpfiles.append(f"{outdir}/msas/mgnify_hits.sto")
tmpfiles.append(f"{outdir}/msas/pdb_hits.hhr")
tmpfiles.append(f"{outdir}/msas/uniref90_hits.sto")

# Others
for model in range(1, 6):
    tmpfiles.append(f"{outdir}/confidence_model_{model}_pred_0.json")
    tmpfiles.append(f"{outdir}/result_model_{model}_pred_0.pkl")
    tmpfiles.append(f"{outdir}/unrelaxed_model_{model}_pred_0.cif")
    tmpfiles.append(f"{outdir}/unrelaxed_model_{model}_pred_0.pdb")

    # Only if monomer_ptm was run:
    if maxfrag == 1:
        tmpfiles.append(f"{outdir}/confidence_model_{model}_ptm_pred_0.json")
        tmpfiles.append(f"{outdir}/pae_model_{model}_ptm_pred_0.json")
        tmpfiles.append(f"{outdir}/result_model_{model}_ptm_pred_0.pkl")
        tmpfiles.append(f"{outdir}/unrelaxed_model_{model}_ptm_pred_0.cif")
        tmpfiles.append(f"{outdir}/unrelaxed_model_{model}_ptm_pred_0.pdb")



# Clean up

print("\n\nTotal time:")
Stoptime()

# Delete temporary files and (empty) directories
tmpacc = acc.replace("-", "_")
tmpoutdir = outdir.replace("-", "_")
if not Switch('debug'):
    print()
    for tmpfile in tmpfiles:
        print(f"Removing temporary file '{tmpfile}'")
        os.remove(tmpfile)

    for tmpdir in tmpdirs:
        print(f"Removing (empty) temporary directory '{tmpdir}'")
        os.rmdir(tmpdir)
    
    # Remove log files if error log is empty, otherwise keep
    errlog = f"log-errors-update_alphasync_input_alphasync_alphafold_{tmpoutdir}_____________alphafold_py_{tmpacc}_{frag}_{maxfrag}.txt"
    outlog = f"log-output-update_alphasync_input_alphasync_alphafold_{tmpoutdir}_____________alphafold_py_{tmpacc}_{frag}_{maxfrag}.txt"
    errlog_cpu = f"log-errors-update_alphasync_input_alphasync_alphafold_{tmpoutdir}_____________alphafold_py_{tmpacc}_{frag}_{maxfrag}__cpu.txt"
    outlog_cpu = f"log-output-update_alphasync_input_alphasync_alphafold_{tmpoutdir}_____________alphafold_py_{tmpacc}_{frag}_{maxfrag}__cpu.txt"
    if os.path.getsize(errlog) == 0:
        # Also remove cpufirst_attempted file and cpu log files if successfully finished
        if os.path.exists(cpufirst_attempted) and os.path.exists(errlog_cpu) and os.path.getsize(errlog_cpu) == 0:
            print(f"Removing temporary file '{cpufirst_attempted}'")
            os.remove(cpufirst_attempted)
            print(f"Removing log file '{errlog_cpu}'")
            os.remove(errlog_cpu)
            print(f"Removing log file '{outlog_cpu}'")
            os.remove(outlog_cpu)
        print(f"Removing log file '{errlog}'")
        os.remove(errlog)
        print(f"Removing log file '{outlog}'")
        os.remove(outlog)

    # Remove main directory for this accession (will throw an error if not empty)
    os.rmdir(f"../{outdir}")

print("\n\nTotal time:")
Stoptime()

print("\nDone!")
