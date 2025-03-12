"""Utility functions & package imports"""

# System
import copy
from itertools import islice
from math import nan, isnan, floor, ceil
from glob import glob
import json
import os
# from pandas import concat, set_option
import random
import re
import requests
import shutil
from statistics import mean, median, stdev
import sys
import tabulate
import time
import warnings
from Bio.SeqIO.FastaIO import SimpleFastaParser as Fasta


# Set terminal width for some printing functions
# termwidth = 80
# termwidth = os.get_terminal_size()[0]
# termwidth = shutil.get_terminal_size(fallback=(318, 91))[0]   # 5K screen
# termwidth = shutil.get_terminal_size(fallback=(222, 64))[0]     # 16" MacBook Pro
termwidth = shutil.get_terminal_size(fallback=(220, 64))[0]     # 16" MacBook Pro

# # Make pandas print all rows of a data frame
# set_option('display.max_rows', None)
# # set_option('display.max_columns', None)
# # set_option('display.width', None)
# # set_option('display.max_colwidth', -1)

def rx(pattern, string):
    """Simpler regular expression match function"""
    m = re.search(pattern, string)
    if m:
        # g = m.groups()
        g = list(m.groups())
        if len(g) > 0:
            # Regexp matched and groups defined:
            
            # Convert all-numeric matches to int
            for i, s in enumerate(g):
                if s != None:
                    if re.fullmatch(r"^[0-9]+$", s):
                        g[i] = int(s)

            # Return groups
            return g
        else:
            # Regexp matched, but no groups defined:
            return True
    else:
        # Regexp not matched:
        # Return False
        return False

# Debugger (ipdb, IPython)
import ipdb
d = ipdb.set_trace


# Warn
from warnings import warn

# Natural sort
# from natsort import natsorted as nsort
from natsort import natsorted, ns
import functools
# https://pypi.org/project/natsort/
# https://natsort.readthedocs.io/en/stable/api.html#natsort.ns
# "GROUPLETTERS, G:     Tell natsort to group lowercase and uppercase letters together when sorting. For example, ['Banana', 'apple', 'banana', 'Apple'] would be sorted as ['Apple', 'apple', 'Banana', 'banana']. Useless when used with IGNORECASE; use with LOWERCASEFIRST to reverse the order of upper and lower case. Generally not needed with LOCALE."
# "PRESORT, PS:     Sort the input as strings before sorting with the nasort algorithm. This can help eliminate inconsistent sorting in cases where two different strings represent the same number. For example, “a1” and “a01” both are internally represented as (“a”, “1), so without PRESORT the order of these two values would depend on the order they appeared in the input (because Python’s sorted is a stable sorting algorithm)."
nsort = functools.partial(natsorted, alg=ns.GROUPLETTERS | ns.PRESORT)
natsort = nsort

# Unique
def Unique(a):
    return nsort(set(a))

# Progress
from tqdm import tqdm
tqd = functools.partial(tqdm, file=sys.stdout, unit='', smoothing=0, bar_format="{desc:<5.5}{percentage:3.0f}%[{bar:32}] {n:,} / {total:,} >> {rate_fmt}{postfix}, {elapsed} elapsed, {remaining} remaining")

# Function that automatically uses the rowcount of an sqlalchemy CursorResult as the total for tqdm, and the line count of a file as the total for tqdm
def tq(it, *args, **kwargs):
    # Unless "total" has been given explicitly:
    if "total" not in kwargs:
        if it.__class__.__name__ == "CursorResult":
            # For SQL CursorResult objects, use query.rowcount as the total
            kwargs["total"] = it.rowcount
        elif it.__class__.__name__ == "TextIOWrapper":
            # For files, use line count as the total
            kwargs["total"] = Lines(it.name)
    # Any others will simply use len() (the default)
    return tqd(it, *args, **kwargs)

# read_tsv and write_tsv (only if pandas has been loaded)
if "pandas" in sys.modules:
    import pandas as pd
    # Set column width to auto (pd.options.display.width defaults to 80)
    pd.set_option('display.width', None)
    import csv
    
    # Read TSV table
    read_tsv = functools.partial(pd.read_csv, sep="\t")

    # Write dataframe to TSV
    def write_tsv(data, outfile):
        data.to_csv(outfile, sep="\t", quoting=csv.QUOTE_NONE, index=False)

# Run
# import subprocess as sp
# import shlex

# Run function:
# Run("List directory contents", "ls -1")       # With description
# Run("ls -1")                                  # Without descriptions: arguments are swapped
def Run(command, description="", silent=False):

    # Swap variables for flexible syntax (one or two arguments)
    if description != "":
        command, description = description, command
    
    # Print command
    if silent == 0:

        if description != "":

            # With description
            # print("\n--------------------------------------------------------------------------------")
            # print(f"> {description} ({command})")
            # print("--------------------------------------------------------------------------------\n", flush=True)
            print()
            print("-" * termwidth)
            print(f"> {description} ({command})")
            print("-" * termwidth)
            print(flush=True)

        else:

            # Without description
            # print("\n--------------------------------------------------------------------------------")
            # print(f"> {command}")
            # print("--------------------------------------------------------------------------------\n", flush=True)
            print()
            print("-" * termwidth)
            print(f"> ({command})")
            print("-" * termwidth)
            print(flush=True)
    
    # Run command
    # sp.run(shlex.split(command), shell=True, capture_output=False)    # Doesn't work with pipes
    os.system(command)

def Return(command):
    """Return: Capture STDOUT from a shell command as a string (can be multiline) (and print STDERR as a warning)"""
    # r = sp.run(shlex.split(command), capture_output=True)    # Doesn't work with pipes
    # if r.stderr != b"":
    #     Warn(r.stderr.rstrip())
    # return(str(r.stdout.decode().rstrip()))
    r = os.popen(command).read().rstrip()

    # Convert all-numeric output to int
    if re.fullmatch(r"^[0-9]+$", r):
       r = int(r)

    return r

def ReturnList(command):
    """Return: Capture multi-line STDOUT from a shell command as a list (and print STDERR as a warning)"""
    lines = os.popen(command).readlines()
    
    # Process individual lines
    r = []
    for line in lines:

        # rstrip lines (removing linebreaks)
        line = line.rstrip()
        
        r.append(line)
    
    # # Convert single-line output from list to string
    # if len(r) == 1:
    #     r = str(r[0])
    #     # Convert all-numeric output to int
    #     if re.fullmatch(r"^[0-9]+$", r):
    #        r = int(r)
    
    return r

def ReturnSet(command):
    """Return: Capture multi-line STDOUT from a shell command as a set (and print STDERR as a warning)"""
    return set(ReturnList(command))

def Lines(f):
    """Get number of lines in a file using wc -l"""
    r = Return(f"cat '{f}' | wc -l")
    r = int(r)
    return r

def Files(f):
    """Get number of files in a TAR archive using wc -l"""
    r = Return(f"tar -tf '{f}' | wc -l")
    r = int(r)
    return r

# Use warnings.warn instead    
def Warn(warning):
    # print(warning, file=sys.stderr)
    # print(f"Warning: {warning}", file=sys.stderr)
    # print(f"{warning}", file=sys.stderr)
    print(f"\nWarning: {warning}\n\n", file=sys.stderr)
    # print(f"\n\nWarning: {warning}\n\n\n", file=sys.stderr)
    # print(f"\n\n\n{warning}\n", file=sys.stderr)

def Die(message=""):
    # exit(message)
    # exit(f"Error: {message}")
    # print("\n\n")
    # raise Exception(message)
    raise Exception(f"\n\n\nError: {message}\n")
    # raise Exception(f"\n\n\n{message}\n")

def q(message = ""):
    sys.exit(message)

# Set-based logging
def Log(myset, item):
    """Log an item to a set"""
    global blang_log

    # Create dictionary if it doesn't exist yet
    if "blang_log" not in globals():
        blang_log = {}

    # Add item to set
    if myset not in blang_log:
        # Create new set if it doesn't exist yet
        blang_log[myset] = set()

    # Add to set
    blang_log[myset].add(str(item))
# Add "Add" as an alias for "Log"
Add = Log

# Get a Log set
def Get(myset):
    """Return a log set"""
    global blang_log
    # # With nsort (slow)
    # return nsort(blang_log[myset])
    # # With sort
    # return sorted(blang_log[myset])
    # Unsorted
    if myset in blang_log.keys():
        return blang_log[myset]
    else:
        return set()

# Delete a Log set
def Delete(myset):
    """Delete a log set"""
    global blang_log
    del blang_log[myset]

# Show Log (individual set, or all sets)
def Show(myset=None, lim=-1, sort=False):
    """
    lim:  Sets below this limit will be shown in full
    sort: False shows log sets in order of logging, True shows them alphabetically, and 'count' shows them sorted by the number of items in each set
    """
    
    if myset != None:
        if "blang_log" in globals():
            # Show an individual set
            if myset not in blang_log.keys():
                # Warn(f"Show(): Set '{myset}' is empty (0)")
                print(myset + " (" + format(0, ',d') + "):\n" + "\n")
            else:
                # With nsort (slow)
                # print(myset + " (" + format(len(blang_log[myset]), ',d') + "):\n" + "\n".join(nsort(blang_log[myset])) + "\n")
                # With sort (slow)
                print(myset + " (" + format(len(blang_log[myset]), ',d') + "):\n" + "\n".join(sorted(blang_log[myset])) + "\n")
                # Unsorted
                # print(myset + " (" + format(len(blang_log[myset]), ',d') + "):\n" + "\n".join(blang_log[myset]) + "\n")
    else:
        # Show all sets
        if "blang_log" in globals():
            
            # Get longest set name for formatting (padding)
            maxlen = len(max(blang_log.keys(), key=len))

            # Get width of longest set item count for formatting (padding)
            maxcount = 0
            for myset in blang_log.keys():
                if maxcount < len(format(len(blang_log[myset]), ',d')):
                    maxcount = len(format(len(blang_log[myset]), ',d'))

            print()
            
            logkeys = blang_log.keys()
            if sort == True:
                logkeys = nsort(logkeys)
            elif sort == "count":
                logkeys = sorted(nsort(logkeys), key=lambda x: len(blang_log[x]), reverse=True)
            for myset in logkeys:
                if len(blang_log[myset]) <= lim or lim == 0:
                    # Up to the limit: show full set contents
                    # With sort
                    print(myset + " (" + format(len(blang_log[myset]), ',d') + "):\n" + "\n".join(sorted(blang_log[myset])) + "\n")
                else:
                    # Above limit: show abbreviated
                    # Unsorted (need to convert to list to retrieve the example, though)
                    print(myset + (' · ' * floor((maxlen - len(myset)) / 3) + ' ' * floor((maxlen - len(myset)) % 3)) + ": " + format(len(blang_log[myset]), ',d') + (' ' * (maxcount - len(format(len(blang_log[myset]), ',d')))) + " (e.g. " + list(blang_log[myset])[0] + ")")
                
def printn(s):
    """Print without newline"""
    print(s, end="")

def State(s):
    """Make a statement. (i.e. print something with nice newlines around it for dramatic effect)"""
    print("\n"+s+"\n")

def Comma(s):
    """Format number with thousands separator (comma)"""
    return format(s, ',d')
    
def Percent(s, digits=0):
    """Format number as percentage with 'digits' fractional digits"""
    return format(s, f'.{digits}%')
    
def Acc(s):
    """Check if a string is a valid UniProt accession"""
    # The regular expression is from the UniProt help section:
    # https://www.uniprot.org/help/accession_numbers
    if re.fullmatch("^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$", s.upper()):
        return True
    else:
        return False

def Args(needed=0, usage="", example="", silent=False):
    """Parse command line arguments and switches"""
    global blang_switches
    
    # Create "switches" dictionary if it doesn't exist yet
    if "blang_switches" not in globals():
        blang_switches = set()

    # Parse command line switches
    # for switch in filter(re.compile(r"^-").search, sys.argv[1:]):
    # for switch in filter(lambda v: re.search(r"^-(.+)", v), sys.argv[1:]):
    for arg in sys.argv[1+needed:]:
        # for m in re.finditer(r"^-(.+)", arg):
        arg = re.sub(r'^-', '', arg)
        # blang_switches.add(m[1])
        blang_switches.add(arg)
        if silent == False:
            print(f"Switch active: -{arg}")
    
    # Return args
    # return list(filter(re.compile(r"^[^\-]").search, sys.argv[1:]))
    # args = list(filter(lambda v: not re.search(r"^-", v), sys.argv[1:1+needed]))
    args = sys.argv[1:1+needed]
    
    # If there are fewer arguments than required: show usage instructions
    if len(args) < needed:
        # sys.exit(usage)
        s = f"\nError: {needed} arguments needed. ({len(args)} given)\n"
        i = 1
        for arg in args:
            s = s + f" >> argument {i}: '{arg}'\n"
            i = i + 1
        # Show usage
        s = s + "\nUsage:   " + re.sub(r"^\.\/", "", sys.argv[0]) + " " + usage
        # Show example
        s = s + "\nExample: " + re.sub(r"^\.\/", "", sys.argv[0]) + " " + example
        sys.exit(s)

    # Convert all-numeric values to int or float
    for i, v in enumerate(args):
        if re.fullmatch(r"^[0-9]+$", v):
            args[i] = int(v)
        elif re.fullmatch(r"^-?[0-9]+\.[0-9]+$", v):
            args[i] = float(v)

    # Convert tuple to str if there's only one argument being returned (to avoid having to do e.g. "(arg,) = Args(1)")
    if (len(args) == 1):
        (args,) = args

    # Return arguments
    return args

def Switch(s):
    """Return True if a given command line switch (parsed by Args()) is active"""
    global blang_switches
    if s in blang_switches:
        return True
    else:
        return False

def SetSwitch(s):
    """Set given command line switch to active"""
    global blang_switches
    blang_switches.add(s)

def DirExists(f):
    """Check if a directory exists"""
    if os.path.isdir(f):
        return True
    else:
        return False

def Exists(f):
    """Check if a file or directory exists"""
    if os.path.exists(f):
        return True
    else:
        return False

def Nonzero(f):
    """Check if a file exists and is non-zero in size"""
    if os.path.exists(f) and os.path.getsize(f) > 0:
        return True
    else:
        return False

def Basename(f):
    """Return base name of file"""
    return os.path.basename(f)

def Cd(dir, silent=False):
    """Change directory (and announce it)"""

    # Print
    if silent != True:
        print()
        print("( " + "-" * (termwidth - 4) + " )")
        print(f"( Changing directory to {dir} )")
        print("( " + "-" * (termwidth - 4) + " )")
        print(flush=True)

    # Change directory
    os.chdir(dir)

# FASTA functions
def Split60(seq):
    """Add linebreaks every 60 characters to a sequence (FASTA-style), remove trailing newline and return"""
    # return re.sub("(.{60})", "\\1\n", seq)
    return re.sub(r'\n$', '', re.sub("(.{60})", "\\1\n", seq))

# Alignment functions

# Utility function to get nth amino acid character in a sequence that contains gaps (i.e. to get the alignment position corresponding to the nth amino acid in the reference species)
def site_to_alnsite(seq, site, pattern=r"[ACDEFGHIKLMNPQRSTVWY]"):
    """Find the index of the nth occurrence of a regular expression pattern in a string.

    Args:
        seq (str): The input string to search within.
        site (int): The desired occurrence number (1-based).
        pattern (str, optional): The regular expression pattern to match. Defaults to amino acid characters [ACDEFGHIKLMNPQRSTVWY].

    Returns:
        int: The index of the nth occurrence of the pattern in the string (1-based).
             If the pattern is not found or site is invalid (0 or negative), -1 is returned.
    """
    occ = [m.start() + 1 for m in re.finditer(pattern, seq, re.IGNORECASE)]
    return occ[site-1] if len(occ) >= site and site > 0 else Die("Error: site_to_alnsite(): site out of range")

# Utility function to get the occurrence number of a regular expression pattern in a string at a specific index (i.e. to get the species position corresponding to the nth amino acid in the alignment)
def alnsite_to_site(seq, alnsite, pattern=r"[ACDEFGHIKLMNPQRSTVWY]"):
    """Find the occurrence number of a regular expression pattern in a string at a specific index.

    Args:
        seq (str): The input string to search within.
        alnsite (int): The index of the pattern in the string (1-based).
        pattern (str, optional): The regular expression pattern to match. Defaults to amino acid characters [ACDEFGHIKLMNPQRSTVWY].

    Returns:
        int: The occurrence number of the pattern in the string (1-based).
             If the pattern is not found or site is invalid (0 or negative), -1 is returned.
    """
    occ = [m.start() + 1 for m in re.finditer(pattern, seq, re.IGNORECASE)]
    # return occ.index(alnsite) + 1 if alnsite in occ else Die("Error: alnsite_to_site(): alnsite out of range")
    return occ.index(alnsite) + 1 if alnsite in occ else None

# Amino acid functions (Biopython probably has these, but I prefer having full control over these)
def Aa(s):
    """Check if string is made up exclusively of the standard 20 amino acids (FASTA format)"""
    if re.fullmatch(r"^[ACDEFGHIKLMNPQRSTVWY]+$", s):
        return True
    else:
        return False

def Aax(s):
    """Check if string is made up exclusively of the standard 20 amino acids, plus X (FASTA format)"""
    if re.fullmatch(r"^[ACDEFGHIKLMNPQRSTVWXY]+$", s):
        return True
    else:
        return False

def Aau(s):
    """Check if string is made up exclusively of the standard 20 amino acids, plus selenocysteine (U) (FASTA format)"""
    if re.fullmatch(r"^[ACDEFGHIKLMNPQRSTUVWY]+$", s):
        return True
    else:
        return False

def Aaux(s):
    """Check if string is made up exclusively of the standard 20 amino acids, plus selenocysteine (U) and X (FASTA format)"""
    if re.fullmatch(r"^[ACDEFGHIKLMNPQRSTUVWXY]+$", s):
        return True
    else:
        return False

def ThreeToOne(s):
    # Check sequence length
    if len(s) % 3 != 0:
        Die(f"Error: Sequence length is not a multiple of 3 ({len(s)}):\n\n{s}\n\n")
    
    # Convert to upper case
    s = s.upper()

    # Go through triplets
    m = re.findall(r"([A-Z]{3})", s)
    if not m:
        Die("Error: Couldn't parse triplets:\n\n{s}\n\n")
    res = ""
    for three in m:
        if three == "ALA": res += "A"
        if three == "ARG": res += "R"
        if three == "ASN": res += "N"
        if three == "ASP": res += "D"
        if three == "CYS": res += "C"
        if three == "GLN": res += "Q"
        if three == "GLU": res += "E"
        if three == "GLY": res += "G"
        if three == "HIS": res += "H"
        if three == "ILE": res += "I"
        if three == "LEU": res += "L"
        if three == "LYS": res += "K"
        if three == "MET": res += "M"
        if three == "PHE": res += "F"
        if three == "PRO": res += "P"
        if three == "SEC": res += "U"   # Selenocysteine
        if three == "SER": res += "S"
        if three == "THR": res += "T"
        if three == "TRP": res += "W"
        if three == "TYR": res += "Y"
        if three == "VAL": res += "V"
        
    return res

def OneToThree(s):
    # Convert to upper case
    s = s.upper()

    # Go through residues
    m = re.findall(r"([A-Z])", s)
    if not m:
        Die("Error: Couldn't parse residues:\n\n{s}\n\n")
    res = ""
    for one in m:
        if one == "A": res += "ALA"
        if one == "R": res += "ARG"
        if one == "N": res += "ASN"
        if one == "D": res += "ASP"
        if one == "C": res += "CYS"
        if one == "Q": res += "GLN"
        if one == "E": res += "GLU"
        if one == "G": res += "GLY"
        if one == "H": res += "HIS"
        if one == "I": res += "ILE"
        if one == "L": res += "LEU"
        if one == "K": res += "LYS"
        if one == "M": res += "MET"
        if one == "F": res += "PHE"
        if one == "P": res += "PRO"
        if one == "U": res += "SEC"   # Selenocysteine
        if one == "S": res += "SER"
        if one == "T": res += "THR"
        if one == "W": res += "TRP"
        if one == "Y": res += "TYR"
        if one == "V": res += "VAL"
        
    return res

def Strandflip(s):
    """Flips DNA strand (reverse complement)"""
    
    m = rx(r"^[ACGT]+$", s)
    if not m:
        Die(f"Error while strand-flipping: Sequence contains characters other than A, C, G, T: '{s}'") ;
    
    # Reverse string
    res = ''.join(reversed(s));
    # Convert to lower case (for transliteration)
    # res =~ tr/ACGT/TGCA/;
    res = res.lower()
    # Transliterate (complementary strand)
    res = re.sub("a", "T", res)
    res = re.sub("c", "G", res)
    res = re.sub("g", "C", res)
    res = re.sub("t", "A", res)
    
    return res

# AlphaSync functions

# Replace X stretches with alanines (A), used for short internal (≤3 aa) stretches of X
def replace_with_alanines(match):
    return "A" * len(match.group())

# Replace X stretches with glycine linkers (GGGGS repeated), used for terminal or long internal (≥4 aa) stretches of X
def replace_with_ggggs(match):
    n = len(match.group())
    return ("GGGGS" * (n//5) + "GGGGS"[:n%5])

# Replace non-standard amino acids (B/Z/U/X) with standard AAs in sequence (for compatibility with AlphaFold)
def ReplaceNonstandardAAs(seq):
    # Replace U with C (selenocysteine)
    seq = seq.replace("U", "C")
    # Replace B with asparagine (N)
    seq = seq.replace("B", "N")
    # Replace Z with glutamine (Q)
    seq = seq.replace("Z", "Q")
    # Replace N-terminal X stretches (any length) with glycine linkers (GGGGS repeated) to maintain flexibility
    seq = re.sub(r'^X+', replace_with_ggggs, seq)
    # Replace C-terminal X stretches (any length) with glycine linkers (GGGGS repeated) to maintain flexibility
    seq = re.sub(r'X+$', replace_with_ggggs, seq)
    # Replace long internal X stretches ≥4 aa with glycine linkers (GGGGS repeated)
    seq = re.sub(r'(?<!X)X{4,}(?!X)', replace_with_ggggs, seq)
    # Replace short internal X stretches ≤3 aa with alanines (A) to maintain secondary structure
    seq = re.sub(r'(?<!X)X{1,3}(?!X)', replace_with_alanines, seq)
    # Return standard-AA sequence
    return seq

# UniProt functions

# Get local UniProt release (check release_file)
def get_local_uniprot_release(release_file = "input/uniprot/local_release.txt"):
    if os.path.exists(release_file):
        with open(release_file, 'r') as f:
            return f.read().strip()
    return None

# Get current UniProt release (check UniProt FTP)
def get_current_uniprot_release(url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/reldate.txt"):
    response = requests.get(url)
    # Parse release number 2024_05 from this response:
    # UniProt Knowledgebase Release 2024_05 consists of:
    # UniProtKB/Swiss-Prot Release 2024_05 of 02-Oct-2024
    # UniProtKB/TrEMBL Release 2024_05 of 02-Oct-2024
    if response.status_code == 200:
        lines = response.text.strip().split('\n')
        for line in lines:
            if "UniProt Knowledgebase Release" in line:
                return line.split()[3]
    return None
    
# Update local UniProt release, by default to current FTP release (update release_file)
def update_local_uniprot_release(release = "", release_file = "input/uniprot/local_release.txt"):
    if release == "":
        release = get_current_uniprot_release()
    os.makedirs(os.path.dirname(release_file), exist_ok=True)
    with open(release_file, 'w') as f:
        f.write(release)
    return None

# Check if UniProt has been updated, and exit if not (unless -refresh is active)
# Specify update = 1 to update the locally recorded release in release_file
def check_uniprot_release(update = 0):
    local_release = get_local_uniprot_release()
    current_release = get_current_uniprot_release()
    local_release_display = local_release
    if local_release is None:
        local_release_display = "unknown"
    print(f" >> UniProt release:")
    print(f"   >> Local: \t{local_release_display}")
    print(f"   >> Current: \t{current_release} (checked via FTP)")
    if local_release == current_release:
        if Switch("refresh"):
            print(" >> Local UniProt version is still current, but running anyway (-refresh is active)...")
        else:
            print(" >> Local UniProt version is still current, exiting (use -refresh to run anyway)!")
            sys.exit(0)
    if update == 0:
        print(f" >> Updating to UniProt release: {current_release}...")
        update_local_uniprot_release(release = current_release)
    else:
        # Update release file
        print(f" >> Updated to UniProt release: {current_release}")
        update_local_uniprot_release(release = current_release)
    return None

# Cluster functions (LSF)
def Myjobs():
    """Get number of LSF jobs currently running for my user (interactive and non-interactive, running and pending)"""
    n = None
    while not re.fullmatch(r"\d+", str(n)):
        # Keep running this until it returns a number (it returns 0 even if there are no jobs, and even if the user doesn't )
        # # Number of cores
        # n = Return(r"busers $USER | awk '{print $4}' | grep -v 'NJOBS'")
        # Number of jobs
        n = Return(r"bjobs -u $USER 2>/dev/null | grep -v '^JOBID\t' | wc -l")
    return n
    
def Pendingjobs(queue = ""):
    """Get number of LSF jobs currently running for my user (interactive and non-interactive, pending only)"""
    pending = None
    while not re.fullmatch(r"\d+", str(pending)):
        # Keep running this until it returns a number (it returns 0 even if there are no jobs, and even if the user doesn't )
        if queue == "":
            # All queues
            # # Number of cores:
            # pending = Return(r"busers $USER | awk '{print $5}' | grep -v 'PEND'")
            # Number of jobs
            pending = Return(r"bjobs -u $USER 2>/dev/null | grep -v '^JOBID\t' | awk '{print $3}' | grep '^PEND$' | wc -l")
        else:
            # Specific queue
            # Number of jobs
            pending = Return(r"bjobs -u $USER -q "+queue+r" 2>/dev/null | grep -v '^JOBID\t' | awk '{print $3}' | grep '^PEND$' | wc -l")
    return pending
    
def Runningjobs():
    """Get number of LSF jobs currently running for my user (interactive and non-interactive, running only)"""
    n = None
    while not re.fullmatch(r"\d+", str(n)):
        # Keep running this until it returns a number (it returns 0 even if there are no jobs, and even if the user doesn't )
        # # Number of cores
        # n = Return(r"busers $USER | awk '{print $6}' | grep -v 'RUN'")
        # Number of jobs
        n = Return(r"bjobs -u $USER 2>/dev/null | grep -v '^JOBID\t' | awk '{print $3}' | grep '^RUN$' | wc -l")
    return n
    
def Interactivejobs():
    """Get number of LSF jobs currently running for my user (interactive only)"""
    n = None
    while not re.fullmatch(r"\d+", str(n)):
        # Keep running this until it returns a number (it returns 0 even if there are no jobs, and even if the user doesn't )
        n = Return(r"bjobs -w -u $USER -q interactive 2>/dev/null | awk '{print $7}' | grep -v 'JOB_NAME' | wc -l")
    return n

def Thesejobs(queue = ""):
    """Get number of LSF jobs currently running for my user, in the current directory"""
    locale = Return("~/scripts/locale.sh")
    n = None
    while not re.fullmatch(r"\d+", str(n)):
        # Keep running this until it returns a number (it returns 0 even if there are no jobs, and even if the user doesn't )
        if queue == "":
            n = Return(r"bjobs -w -u $USER 2>/dev/null | awk '{print $7}' | grep -v 'JOB_NAME' | grep -P '^"+locale+r"_' | wc -l")
        else:
            n = Return(r"bjobs -w -u $USER -q "+queue+r" 2>/dev/null | awk '{print $7}' | grep -v 'JOB_NAME' | grep -P '^"+locale+r"_' | wc -l")
            # print(r"bjobs -w -u $USER -q "+queue+r" 2>/dev/null | awk '{print $7}' | grep -v 'JOB_NAME' | grep -P '^"+locale+r"_' | wc -l")
    return n

def Nodetype():
    """Get node type (submit host or compute node)"""
    hostname = Return("hostname")
    if rx("^splprhpc", hostname):
        return "lsf"
    else:
        return "node"

def Waitforjobs(queue = ""):
    """Wait for jobs from a particular script pipeline directory to finish"""
    if Nodetype() != "node" or queue != "":
        # If script is running on submit host, or if we're looking at a specific queue:
        # Expect 0+ jobs
        minjobs = 0 + Interactivejobs()
    else:
        # Expect 1+ jobs
        minjobs = 1 + Interactivejobs()
    
    # Give LSF a few seconds to show the job
    time.sleep(5)

    if queue == "":
        State(f"Waiting for these jobs to finish (currently {Thesejobs(queue)}, expecting {minjobs})...")
    else:
        State(f"Waiting for these jobs to finish (queue '{queue}') (currently {Thesejobs(queue)}, expecting {minjobs})...")
    
    while (Thesejobs(queue) > minjobs):
        time.sleep(3)

def Waitforalljobs():
    """Wait for all of my jobs to finish"""
    if Nodetype() != "node":
        # If script is running on submit host:
        minjobs = 0 + Interactivejobs()
    else:
        minjobs = 1 + Interactivejobs()
    
    # Give LSF a few seconds to show the job
    time.sleep(5)

    State(f"Waiting for all jobs to finish (currently {Myjobs()}, expecting {minjobs})...")
    
    while (Myjobs() > minjobs):
        time.sleep(3)

def Starttime(timer=0):
    """Start a timer"""
    global blang_timer

    # Create timer dictionary if it doesn't exist yet
    if "blang_timer" not in globals():
        blang_timer = {}

    # Set timer
    blang_timer[timer] = time.perf_counter()
# Alias
Timestart = Starttime

def Stoptime(timer=0):
    """Stop a timer and show elapsed time"""
    global blang_timer

    if timer in blang_timer:
        # Get seconds elapsed
        s = time.perf_counter() - blang_timer[timer]
        
        # Process
        m = 0; h = 0; d = 0;
        
        while (s > 60): m += 1; s -= 60
        while (m > 60): h += 1; m -= 60
        while (h > 24): d += 1; h -= 24

        # Build string
        d = f'{d} days  '  if d != 0 else ''
        h = f'{h} h  '     if h != 0 else ''
        m = f'{m} mins  '  if m != 0 else ''
        s = f'{s:.3f} sec' if s != 0 else ''

        # Print
        print(f"Time elapsed:\t{d}{h}{m}{s}\n");
    else:
        Die(f"Error: Timer '{timer}' wasn't set")
# Alias
Timestop = Stoptime

def Deletetime(timer=1):
    """Remove a timer without showing elapsed time"""
    global blang_timer

    if timer in blang_timer:
        del blang_timer[timer]
    # else:
    #     Die(f"Error: Timer '{timer}' wasn't set")

def Time(timer=1):
    """Start a new timer, or stop one and show elapsed time if it already exists"""
    global blang_timer

    # Create timer dictionary if it doesn't exist yet
    if "blang_timer" not in globals():
        blang_timer = {}

    if timer not in blang_timer:
        # Start timer
        blang_timer[timer] = time.perf_counter()
    else:
        # Stop timer and show elapsed time
        s = time.perf_counter() - blang_timer[timer]
        
        # Process
        m = 0; h = 0; d = 0;
        
        while (s > 60): m += 1; s -= 60
        while (m > 60): h += 1; m -= 60
        while (h > 24): d += 1; h -= 24

        # Build string
        d = f'{d} days  '  if d != 0 else ''
        h = f'{h} h  '     if h != 0 else ''
        m = f'{m} mins  '  if m != 0 else ''
        s = f'{s:.3f} sec' if s != 0 else ''

        # Print
        print(f"Time elapsed:\t{d}{h}{m}{s}\n");
        
        # Set to current time
        blang_timer[timer] = time.perf_counter()

def GetFasta(seqfile):
    seqs = {}
    with open(seqfile) as f:
        for title, seq in Fasta(f):
            # print(f">{title}\n{seq}")
            seqs[title] = seq
    # Return dictionary
    return seqs

# Characterise a distribution, rounding to 3 digits by default
# def Characterise(a, digits=3):
# Characterise a distribution, not rounding by default
def Characterise(a, digits=None):
    if len(a) > 0:
        if type(digits) == int:
            # With rounding
            print(f"NUM {Comma(len(a))}    SUM {round(sum(a), digits)}    MIN {round(min(a), digits)}    MED {round(median(a), digits)}    AVG {round(mean(a), digits)}    MAX {round(max(a), digits)}    STD {round(stdev(a), digits)}")
        else:
            # Without rounding
            print(f"NUM {Comma(len(a))}    SUM {sum(a)}    MIN {min(a)}    MED {median(a)}    AVG {mean(a)}    MAX {max(a)}    STD {stdev(a)}")
    else:
        print("NUM N/A    MIN N/A    MED N/A    AVG N/A    MAX N/A    STD N/A");
# Add alias
Summarise = Characterise

# # Pretty-print a pandas data frame
# def pp(df):
#     # print(df)
#     # using tabulate:
#     # print(tabulate(df, headers='keys', tablefmt='simple'))

# Count non-comment lines in a file and store in file.lines.txt
def Lines(f, ignore_comments=False):

    linesfile = f + ".lines.txt"

    if not Exists(f):
        Die(f"Error: File '{f}' doesn't exist")

    if not Exists(linesfile):
        # Check if compressed
        if rx(r'\.b?gz$', f):
            cat = 'zcat'
        else:
            cat = 'cat'

        # Note that this won't "re-count" if the ignore_comments setting is changed
        if ignore_comments == True:
            # Count lines and store in .lines.txt file (ignoring comment lines)
            # This actually has the problem that tq() usually progresses past it
            Run(f"Counting lines in file '{f}' (ignoring comment lines)", f"{cat} '{f}' | grep -v '^#' | wc -l > {linesfile}")
        else:
            # Count lines and store in .lines.txt file
            Run(f"Counting lines in file '{f}'", f"{cat} '{f}' | wc -l > {linesfile}")
    
    # Get line count from .lines.txt file and return it
    return Return(f"cat '{linesfile}'")

# Done!
def Done():
    print("\nDone!\n")
