#!/usr/bin/env python3
"""
alphauniprot.py:
Get UniProt sequences and annotation via UniProt's REST API and fill AlphaSync table 'alphauniprot'.
"""

# Initialize
import pandas as pd
# import numpy as np
# import networkx as nx
# import seaborn as sns
# import matplotlib.pyplot as mp
# import tarfile
# import gzip
# import io
import os
import sys
# from ftplib import FTP
# import re
import requests
from requests.adapters import HTTPAdapter, Retry
# import tempfile
import ijson
# from Bio import SeqIO
from blang_mysql import *
from blang import *

# Variables
alphauniprot = "alphauniprot"

# Start
Args(0, "[-reviewed] [-compara_species] [-all] [-debug] [-refresh]", " -reviewed: Get all reviewed proteins from UniProt (i.e. UniProtKB/Swiss-Prot), including isoforms, in JSON format, batch size 500 (default)\n -compara_species: Get all proteins for the 200 Ensembl Compara species from UniProt (i.e. UniProtKB/Swiss-Prot and UniProtKB/TrEMBL), including isoforms, in JSON format, batch size 500\n -all: Get all proteins from UniProt (i.e. UniProtKB/Swiss-Prot and UniProtKB/TrEMBL), including isoforms, in JSON format, batch size 500 (not really feasible)\n -refresh: Run even if the local UniProt release is still current\n -debug: Don't actually make any changes, just simulate", " -compara_species")

# Default to:
# - The 200 Ensembl Compara species,
# - plus the 48 "Model organisms" and "Global health proteomes" from the AlphaFold Protein Structure Database,
# - plus all reviewed proteins from any species (all of UniProtKB/Swiss-Prot),
# all including alternative isoforms.
proteins_to_get = "compara_species"

if Switch("compara_species"):
    proteins_to_get = "compara_species"
if Switch("reviewed"):
    proteins_to_get = "reviewed"
if Switch("all"):
    proteins_to_get = "all"

# Check if UniProt has been updated, and exit if not (unless -refresh is active)
check_uniprot_release()

# Clear table
if not Switch('debug'):
    Clear(alphauniprot)



# Define queries
# # url = 'https://rest.uniprot.org/uniprotkb/search?format=json'
# # url = 'https://rest.uniprot.org/uniprotkb/search?query=organism_id:9606&format=json&size=1000000'
# # url = 'https://rest.uniprot.org/uniprotkb/search?query=P04637&format=json&size=1000000'
# # url = 'https://rest.uniprot.org/uniprotkb/search?query=A0A0C5B5G6&format=json'
# # url = 'https://rest.uniprot.org/uniprotkb/search?query=Q7LBC6&format=json'
# # url = 'https://rest.uniprot.org/uniprotkb/search?query=(keyword:KW-1185)&format=tsv'
# # url = 'https://rest.uniprot.org/uniprotkb/search?query=P04637&format=json'
# # url = 'https://rest.uniprot.org/uniprotkb/search?query=P04637&includeIsoform=true&format=json'
# url = 'https://rest.uniprot.org/uniprotkb/search?query=reviewed:true%20and%20organism_id:9606&includeIsoform=true&format=json&size=500'
# Time(1)
# js = requests.get(url).text
# Time(1)
# data = json.loads(js)
# # df = pd.json_normalize(data['results'])
# d()

# Get all reviewed proteins from UniProt (i.e. UniProtKB/Swiss-Prot), including isoforms, in JSON format, batch size 500
if proteins_to_get == "reviewed":
    # url = 'https://rest.uniprot.org/uniprotkb/search?format=json&size=500&includeIsoform=true&query=%28%28reviewed%3Atrue%29%29'
    queries = ["(reviewed:true)"]

if proteins_to_get == "all":
    # Get all proteins from UniProt (i.e. UniProtKB/Swiss-Prot and UniProtKB/TrEMBL), including isoforms, in JSON format, batch size 500
    # Not really feasible since this is 250 million proteins
    # url = 'https://rest.uniprot.org/uniprotkb/search?format=json&size=500&includeIsoform=true&query=*'
    queries = ['*']

if proteins_to_get == "compara_species":
    
    # # Get all proteins for the 200 Ensembl Compara species from UniProt (i.e. UniProtKB/Swiss-Prot and UniProtKB/TrEMBL), including isoforms, in JSON format, batch size 500
    # # Hardcoded
    # url = 'https://rest.uniprot.org/uniprotkb/search?format=json&size=500&includeIsoform=true&query=%28taxonomy_id%3A6239%29+OR+%28taxonomy_id%3A7227%29+OR+%28taxonomy_id%3A7719%29+OR+%28taxonomy_id%3A7757%29+OR+%28taxonomy_id%3A7764%29+OR+%28taxonomy_id%3A7868%29+OR+%28taxonomy_id%3A7897%29+OR+%28taxonomy_id%3A7918%29+OR+%28taxonomy_id%3A7950%29+OR+%28taxonomy_id%3A7955%29+OR+%28taxonomy_id%3A7957%29+OR+%28taxonomy_id%3A7994%29+OR+%28taxonomy_id%3A7998%29+OR+%28taxonomy_id%3A8005%29+OR+%28taxonomy_id%3A8010%29+OR+%28taxonomy_id%3A8019%29+OR+%28taxonomy_id%3A8022%29+OR+%28taxonomy_id%3A8030%29+OR+%28taxonomy_id%3A8032%29+OR+%28taxonomy_id%3A8049%29+OR+%28taxonomy_id%3A8078%29+OR+%28taxonomy_id%3A8081%29+OR+%28taxonomy_id%3A8083%29+OR+%28taxonomy_id%3A8090%29+OR+%28taxonomy_id%3A8103%29+OR+%28taxonomy_id%3A8128%29+OR+%28taxonomy_id%3A8153%29+OR+%28taxonomy_id%3A8154%29+OR+%28taxonomy_id%3A8175%29+OR+%28taxonomy_id%3A8187%29+OR+%28taxonomy_id%3A8364%29+OR+%28taxonomy_id%3A8478%29+OR+%28taxonomy_id%3A8502%29+OR+%28taxonomy_id%3A8508%29+OR+%28taxonomy_id%3A8630%29+OR+%28taxonomy_id%3A8663%29+OR+%28taxonomy_id%3A8673%29+OR+%28taxonomy_id%3A8840%29+OR+%28taxonomy_id%3A9031%29+OR+%28taxonomy_id%3A9103%29+OR+%28taxonomy_id%3A9135%29+OR+%28taxonomy_id%3A9157%29+OR+%28taxonomy_id%3A9258%29+OR+%28taxonomy_id%3A9305%29+OR+%28taxonomy_id%3A9315%29+OR+%28taxonomy_id%3A9358%29+OR+%28taxonomy_id%3A9361%29+OR+%28taxonomy_id%3A9365%29+OR+%28taxonomy_id%3A9371%29+OR+%28taxonomy_id%3A9483%29+OR+%28taxonomy_id%3A9531%29+OR+%28taxonomy_id%3A9541%29+OR+%28taxonomy_id%3A9544%29+OR+%28taxonomy_id%3A9545%29+OR+%28taxonomy_id%3A9555%29+OR+%28taxonomy_id%3A9568%29+OR+%28taxonomy_id%3A9595%29+OR+%28taxonomy_id%3A9597%29+OR+%28taxonomy_id%3A9598%29+OR+%28taxonomy_id%3A9601%29+OR+%28taxonomy_id%3A9606%29+OR+%28taxonomy_id%3A9615%29+OR+%28taxonomy_id%3A9627%29+OR+%28taxonomy_id%3A9643%29+OR+%28taxonomy_id%3A9646%29+OR+%28taxonomy_id%3A9669%29+OR+%28taxonomy_id%3A9685%29+OR+%28taxonomy_id%3A9689%29+OR+%28taxonomy_id%3A9691%29+OR+%28taxonomy_id%3A9739%29+OR+%28taxonomy_id%3A9749%29+OR+%28taxonomy_id%3A9755%29+OR+%28taxonomy_id%3A9771%29+OR+%28taxonomy_id%3A9785%29+OR+%28taxonomy_id%3A9796%29+OR+%28taxonomy_id%3A9813%29+OR+%28taxonomy_id%3A9823%29+OR+%28taxonomy_id%3A9838%29+OR+%28taxonomy_id%3A9913%29+OR+%28taxonomy_id%3A9925%29+OR+%28taxonomy_id%3A9940%29+OR+%28taxonomy_id%3A9978%29+OR+%28taxonomy_id%3A9986%29+OR+%28taxonomy_id%3A9994%29+OR+%28taxonomy_id%3A9999%29+OR+%28taxonomy_id%3A10020%29+OR+%28taxonomy_id%3A10029%29+OR+%28taxonomy_id%3A10036%29+OR+%28taxonomy_id%3A10089%29+OR+%28taxonomy_id%3A10090%29+OR+%28taxonomy_id%3A10093%29+OR+%28taxonomy_id%3A10096%29+OR+%28taxonomy_id%3A10103%29+OR+%28taxonomy_id%3A10116%29+OR+%28taxonomy_id%3A10141%29+OR+%28taxonomy_id%3A10160%29+OR+%28taxonomy_id%3A10181%29+OR+%28taxonomy_id%3A13489%29+OR+%28taxonomy_id%3A13616%29+OR+%28taxonomy_id%3A13735%29+OR+%28taxonomy_id%3A27687%29+OR+%28taxonomy_id%3A28377%29+OR+%28taxonomy_id%3A28743%29+OR+%28taxonomy_id%3A29073%29+OR+%28taxonomy_id%3A29139%29+OR+%28taxonomy_id%3A30521%29+OR+%28taxonomy_id%3A30522%29+OR+%28taxonomy_id%3A30538%29+OR+%28taxonomy_id%3A30608%29+OR+%28taxonomy_id%3A30611%29+OR+%28taxonomy_id%3A30732%29+OR+%28taxonomy_id%3A31033%29+OR+%28taxonomy_id%3A32507%29+OR+%28taxonomy_id%3A34839%29+OR+%28taxonomy_id%3A35670%29+OR+%28taxonomy_id%3A37003%29+OR+%28taxonomy_id%3A37293%29+OR+%28taxonomy_id%3A37347%29+OR+%28taxonomy_id%3A38626%29+OR+%28taxonomy_id%3A39432%29+OR+%28taxonomy_id%3A40151%29+OR+%28taxonomy_id%3A41447%29+OR+%28taxonomy_id%3A42100%29+OR+%28taxonomy_id%3A42254%29+OR+%28taxonomy_id%3A42514%29+OR+%28taxonomy_id%3A43179%29+OR+%28taxonomy_id%3A43346%29+OR+%28taxonomy_id%3A48698%29+OR+%28taxonomy_id%3A48699%29+OR+%28taxonomy_id%3A48883%29+OR+%28taxonomy_id%3A51154%29+OR+%28taxonomy_id%3A51337%29+OR+%28taxonomy_id%3A51511%29+OR+%28taxonomy_id%3A52904%29+OR+%28taxonomy_id%3A55149%29+OR+%28taxonomy_id%3A56716%29+OR+%28taxonomy_id%3A56723%29+OR+%28taxonomy_id%3A59463%29+OR+%28taxonomy_id%3A59479%29+OR+%28taxonomy_id%3A59729%29+OR+%28taxonomy_id%3A59894%29+OR+%28taxonomy_id%3A60711%29+OR+%28taxonomy_id%3A61621%29+OR+%28taxonomy_id%3A61622%29+OR+%28taxonomy_id%3A61819%29+OR+%28taxonomy_id%3A61853%29+OR+%28taxonomy_id%3A62062%29+OR+%28taxonomy_id%3A64144%29+OR+%28taxonomy_id%3A64176%29+OR+%28taxonomy_id%3A68415%29+OR+%28taxonomy_id%3A69293%29+OR+%28taxonomy_id%3A72004%29+OR+%28taxonomy_id%3A74533%29+OR+%28taxonomy_id%3A74940%29+OR+%28taxonomy_id%3A75366%29+OR+%28taxonomy_id%3A79684%29+OR+%28taxonomy_id%3A80966%29+OR+%28taxonomy_id%3A80972%29+OR+%28taxonomy_id%3A83772%29+OR+%28taxonomy_id%3A84702%29+OR+%28taxonomy_id%3A93934%29+OR+%28taxonomy_id%3A96440%29+OR+%28taxonomy_id%3A99883%29+OR+%28taxonomy_id%3A105023%29+OR+%28taxonomy_id%3A106582%29+OR+%28taxonomy_id%3A106734%29+OR+%28taxonomy_id%3A109280%29+OR+%28taxonomy_id%3A113540%29+OR+%28taxonomy_id%3A123683%29+OR+%28taxonomy_id%3A132585%29+OR+%28taxonomy_id%3A132908%29+OR+%28taxonomy_id%3A144197%29+OR+%28taxonomy_id%3A158456%29+OR+%28taxonomy_id%3A161767%29+OR+%28taxonomy_id%3A183150%29+OR+%28taxonomy_id%3A205130%29+OR+%28taxonomy_id%3A215358%29+OR+%28taxonomy_id%3A223781%29+OR+%28taxonomy_id%3A230844%29+OR+%28taxonomy_id%3A244447%29+OR+%28taxonomy_id%3A283035%29+OR+%28taxonomy_id%3A286419%29+OR+%28taxonomy_id%3A299321%29+OR+%28taxonomy_id%3A303518%29+OR+%28taxonomy_id%3A379532%29+OR+%28taxonomy_id%3A441894%29+OR+%28taxonomy_id%3A445787%29+OR+%28taxonomy_id%3A452646%29+OR+%28taxonomy_id%3A559292%29+OR+%28taxonomy_id%3A586833%29+OR+%28taxonomy_id%3A630221%29+OR+%28taxonomy_id%3A1026970%29+OR+%28taxonomy_id%3A1328070%29+OR+%28taxonomy_id%3A1676925%29+OR+%28taxonomy_id%3A1825980%29+OR+%28taxonomy_id%3A1841481%29+OR+%28taxonomy_id%3A1868482%29+OR+%28taxonomy_id%3A2489341%29+OR+%28taxonomy_id%3A2587831%29+OR+%28taxonomy_id%3A2715852%29'
    # Build from query
    # taxa = {9606}
    # taxa = FetchSet(Query("SELECT tax FROM compara_species"))
    # taxa = ' OR '.join([f"(taxonomy_id:{tax})" for tax in taxa])
    # # urlencode
    # taxa = taxa.replace('(', '%28')
    # taxa = taxa.replace(')', '%29')
    # taxa = taxa.replace(':', '%3A')
    # taxa = taxa.replace(' ', '+')
    # # url = f'https://rest.uniprot.org/uniprotkb/stream?format=json&includeIsoform=true&query=%28{taxa}%29'
    # url = f'https://rest.uniprot.org/uniprotkb/search?format=json&size=500&includeIsoform=true&query=%28{taxa}%29'
    # taxa = FetchList(Query("SELECT tax FROM compara_species ORDER BY id"))
    # Exclude taxa that are already present in table 'alphauniprot' (in case this script is being re-run)
    taxa = FetchList(Query(f"SELECT tax FROM compara_species WHERE tax NOT IN (SELECT DISTINCT tax FROM {alphauniprot}) ORDER BY id"))
    
    # Also include the "Model organisms" and "Global health proteomes" from the AlphaFold Protein Structure Database, all of which start with UP...
    # taxa.append(FetchList(Query("SELECT DISTINCT tax FROM alphafrag WHERE source IN (SELECT DISTINCT source FROM alphafrag HAVING source LIKE 'UP%') AND tax NOT IN (SELECT DISTINCT tax FROM compara_species)")))
    # Exclude taxa that are already present in table 'alphauniprot' (and those that are in table 'compara_species')
    taxa.extend(FetchList(Query(f"SELECT DISTINCT tax FROM alphafrag WHERE tax NOT IN (SELECT DISTINCT tax FROM {alphauniprot}) AND source IN (SELECT DISTINCT source FROM alphafrag HAVING source LIKE 'UP%') AND tax NOT IN (SELECT DISTINCT tax FROM compara_species)")))
    
    # taxa = FetchList(Query("SELECT 3469 AS tax"))
    # taxa = FetchList(Query("SELECT 9358 AS tax"))
    # queries = [f"%28taxonomy_id%3A{tax}%29" for tax in taxa]  # Will also fetch all subspecies, e.g.: https://www.uniprot.org/uniprotkb?query=%28taxonomy_id%3A9606%29+NOT+%28organism_id%3A9606%29
    # queries = [f"%28organism_id%3A{tax}%29" for tax in taxa]  # Only those exact taxa
    queries = []
    # # First, get all reviewed proteins:
    # queries += ["%28reviewed%3Atrue%29"]
    # # Then, add unreviewed only, for the 200 compara_species taxa:
    # queries += [f"%28organism_id%3A{tax}%29+AND+%28reviewed%3Afalse%29" for tax in taxa]
    # First, get reviewed and unreviewed proteins for the 200 compara_species taxa:
    queries += [f"organism_id:{tax}" for tax in taxa]
    # Then, get reviewed proteins for all other (non-compara_species) taxa:
    # Add any taxa that are already fully in alphauniprot (in casethis script is being re-run)
    tmptaxa = taxa
    tmptaxa.extend(FetchList(Query(f"SELECT DISTINCT tax FROM {alphauniprot}")))
    # list(dict(...)): remove duplicates while maintaining order
    queries += ["(reviewed:true) AND NOT (" + " OR ".join([f"organism_id:{tax}" for tax in list(dict.fromkeys(tmptaxa))]) + ")"]
    
    taxa_set = set(taxa)

if Switch("debug"):
    # If debugging: only get p53 and its isoforms
    queries = ["accession:P04637"]



# Function to download data UniProt API to file, showing download progress
def download_uniprot_data(url, block_size=262144):  # 256 KB default
    # retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
    # session = requests.Session()
    # session.mount("https://", HTTPAdapter(max_retries=retries))
    # 
    # response = session.get(url, stream=True)
    # response.raise_for_status()
    # 
    # # # Display compression information
    # # content_encoding = response.headers.get('Content-Encoding', '')
    # # print(f" >> Content-Encoding: {content_encoding}")
    # # is_compressed = 'gzip' in content_encoding or 'deflate' in content_encoding
    # # print(f" >> Is compressed: {is_compressed}")
    # 
    # total_size = int(response.headers.get('content-length', 0))
    
    # # Create directory tmp_alphauniprot if it doesn't exist
    # os.makedirs(temp_dir, exist_ok=True)
    # Create a temporary file
    # with tempfile.NamedTemporaryFile(delete=False, mode='wb', dir=temp_dir) as temp_file:
    # with open(temp_file_name, 'wb') as temp_file:
    #     try:
    #         # with tqdm(total=total_size, unit='iB', unit_scale=True, desc="Downloading") as progress_bar:
    #         # with tqdm(total=total_size, unit='iB', unit_scale=True, desc=f"   >> Downloading to '{temp_file.name}'") as progress_bar:
    #         # with tqdm(total=total_size, file=sys.stdout, unit='B', unit_divisor=1024, unit_scale=True, desc=f"   >> Downloading", bar_format="{desc}: {n_fmt} >> {rate_fmt}{postfix}, {elapsed} elapsed") as progress_bar:
    #         # with tqdm(total=total_size, file=sys.stdout, unit='B', unit_divisor=1024, unit_scale=True, desc=f"   >> Downloading to '{temp_file.name}'", bar_format="{desc}: {n_fmt} >> {rate_fmt}{postfix}, {elapsed} elapsed") as progress_bar:
    #         with tqdm(total=total_size, file=sys.stdout, unit='B', unit_divisor=1024, unit_scale=True, desc=f"   >> Downloading", bar_format="{desc}: {n_fmt} >> {rate_fmt}{postfix}, {elapsed} elapsed") as progress_bar:
    #             for chunk in response.iter_content(chunk_size=block_size):
    #                 size = len(chunk)
    #                 temp_file.write(chunk)
    #                 progress_bar.update(size)
    #     except:
    #         # Clean up the temporary file if interrupted
    #         os.unlink(temp_file_name)
    #         # # Clean up the temporary directory if interrupted
    #         # os.rmdir(temp_dir)

    # Keep trying until successfully fully downloaded
    temp_file_name = None
    while temp_file_name is None:
        try:
            # retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
            retries = Retry(
                total=10,           # Total retries
                backoff_factor=0.5, # Delay between retries
                status_forcelist=[500, 502, 503, 504, 429, 408],  # Rate limit (429) and timeout (408)
                connect=5,          # Connection retries
                read=5,             # Read retries
                respect_retry_after_header=True  # Honor server's retry-after header
            )
            session = requests.Session()
            session.mount("https://", HTTPAdapter(max_retries=retries))
            
            response = session.get(url, stream=True)
            response.raise_for_status()

            total_size = int(response.headers.get('content-length', 0))

            # Create a temporary file
            temp_file_name = "tmp-alphauniprot.json"
            with open(temp_file_name, 'wb') as temp_file:
                with tqdm(total=total_size, file=sys.stdout, unit='B', unit_divisor=1024, unit_scale=True, desc=f"   >> Downloading", bar_format="{desc}: {n_fmt} >> {rate_fmt}{postfix}, {elapsed} elapsed") as progress_bar:
                    for chunk in response.iter_content(chunk_size=block_size):
                        size = len(chunk)
                        temp_file.write(chunk)
                        progress_bar.update(size)
            response.close()
            session.close()

        except:
            # Download failed
            print( "     >> Download failed, retrying...")
            response.close()
            session.close()
            # Clean up the temporary file if interrupted
            os.unlink(temp_file_name)
            temp_file_name = None

    return temp_file_name

def parse_uniprot_data(file_name):
    # Parse JSON from the downloaded file (streaming to avoid high memory usage)
    with open(file_name, 'rb') as json_file:
        try:
            # Use ijson.items to directly iterate over the 'results' array
            for item in ijson.items(json_file, 'results.item'):
            # for item in tqdm(ijson.items(json_file, 'results.item'), unit='', desc=f"   >> Parsing proteins"):
                yield item
        finally:    
            # Clean up the temporary file, including if interrupted
            os.unlink(file_name)
            # # Clean up the temporary directory, including if interrupted
            # os.rmdir(temp_dir)



# Query API

print(f"Getting all '{proteins_to_get}' proteins from UniProt API, including isoforms...")
# if proteins_to_get == "compara_species":
#     State(f"Note: Individual queries might appear to hang at first, especially for 'all reviewed' and 'human' (the first and second queries)")

# # Function to fetch from UniProt API, showing progress
# def fetch_uniprot_data(url):
#     retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
#     session = requests.Session()
#     session.mount("https://", HTTPAdapter(max_retries=retries))
#     
#     response = session.get(url, stream=True)
#     response.raise_for_status()
#     
#     total_size = int(response.headers.get('content-length', 0))
#     block_size = 1024  # 1 KB
#     
#     data = []
#     with tqdm(total=total_size, unit='iB', unit_scale=True, desc="Downloading") as progress_bar:
#         for chunk in response.iter_content(chunk_size=block_size):
#             size = len(chunk)
#             progress_bar.update(size)
#             data.append(chunk)
#     
#     json_data = json.loads(b''.join(data).decode('utf-8'))
#     return json_data


# Run individual queries (one per species)
for qi, query in enumerate(queries, 1):

#     # Standard UniProt API pagination from https://www.uniprot.org/help/api_queries (quite slow, though)
#     # url = f'https://rest.uniprot.org/uniprotkb/search?format=json&size=500&includeIsoform=true&query={query}'
#     print(f"Querying UniProt API ('{url}')...")
# 
#     re_next_link = re.compile(r'<(.+)>; rel="next"')
#     retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
#     session = requests.Session()
#     session.mount("https://", HTTPAdapter(max_retries=retries))
# 
#     def get_next_link(headers):
#         if "Link" in headers:
#             match = re_next_link.match(headers["Link"])
#             if match:
#                 return match.group(1)
# 
#     def get_batch(batch_url):
#         while batch_url:
#             response = session.get(batch_url)
#             response.raise_for_status()
#             total = response.headers["x-total-results"]
#             yield response, total
#             batch_url = get_next_link(response.headers)

    # First, get maximum number of items from API /search endpoint and response.headers["x-total-results"]
    url = f'https://rest.uniprot.org/uniprotkb/search?format=json&includeIsoform=true&size=500&query={query}'
    # retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
    retries = Retry(
        total=10,           # Total retries
        backoff_factor=0.5, # Delay between retries
        status_forcelist=[500, 502, 503, 504, 429, 408],  # Rate limit (429) and timeout (408)
        connect=5,          # Connection retries
        read=5,             # Read retries
        respect_retry_after_header=True  # Honor server's retry-after header
    )
    session = requests.Session()
    session.mount("https://", HTTPAdapter(max_retries=retries))
    response = session.get(url)
    response.raise_for_status()
    total = response.headers["x-total-results"]
    total = int(total)
    response.close()

    # Without pagination (much faster), and per species so the queries aren't too resource-intensive
    url = f'https://rest.uniprot.org/uniprotkb/stream?format=json&includeIsoform=true&query={query}'
    # print(f" >> Query {qi} / {len(queries)}: '{url}'")
    print(f" >> {qi} / {len(queries)} >> Query '{query}'")

    



    # retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
    # session = requests.Session()
    # session.mount("https://", HTTPAdapter(max_retries=retries))
    # response = session.get(url)
    # response.raise_for_status()
    # # total = response.headers["x-total-results"]
    # # total = int(total)
    # # data = requests.get(url).json()
    # batch = response
    # # print(f"Total proteins: {Comma(total)}")
    # data = fetch_uniprot_data(url)

    # Get UniProt data
    # i = 0
    # for batch, total in get_batch(url):
    # Get JSON from batch
    # data = batch.json()
    # total = int(total)
    # i += 1

    # Process each protein in batch
    # for e in tq(data['results'], initial=(i - 1) * 500, total=total, leave=False):
    # for e in tq(data['results'], initial=(i - 1) * 500, total=total):
    # for e in tq(data['results']):
    # for e in fetch_uniprot_data(url):
    temp_file = download_uniprot_data(url)
    # with tqdm(total=0, unit='', desc=f"   >> Parsing proteins") as progress_bar:
        # for qe, e in enumerate(parse_uniprot_data(temp_file), 1):
    # for e in tqdm(parse_uniprot_data(temp_file), total=0, unit='', desc=f"   >> Parsing proteins", smoothing=0, bar_format="{desc:<5.5}{percentage:3.0f}%[{bar:32}] {n:,} / {total:,} >> {rate_fmt}{postfix}, {elapsed} elapsed, {remaining} remaining"):
    # for e in tqdm(parse_uniprot_data(temp_file), total=total, unit='', desc=f"   >> Parsing proteins", smoothing=0, bar_format="{desc}: {n:,} >> {rate_fmt}{postfix}, {elapsed} elapsed"):
    print(f"   >> Parsing proteins:")
    # i = 0
    for e in tqd(parse_uniprot_data(temp_file), total=total):
        # i += 1
        if Switch("debug"):
            d()
        # data['results'][0]['primaryAccession']
        # data['results'][0]['sequence']['length']
        # data['results'][0]['sequence']['value']
        # data['results'][0]['organism']['taxonId']
        # data['results'][0]['organism']['commonName']
        # data['results'][0]['organism']['scientificName']
        # data['results'][0]['uniProtkbId']
        # data['results'][0]['genes'][0]['geneName']['value']
        # data['results'][0]['genes'][0]['synonyms'][0]['value']
        # data['results'][0]['proteinDescription']['recommendedName']['fullName']['value']
        # data['results'][0]['comments'][0]['texts'][0]['value']
        # data['results'][0]['comments'][0]['commentType']

        # Get protein annotation
        acc = e['primaryAccession']
        name = e['uniProtkbId']
        seq = e['sequence']['value']
        seqlen = e['sequence']['length']
        tax = e['organism']['taxonId']
        # "entryType":"UniProtKB reviewed (Swiss-Prot)"
        # "entryType":"UniProtKB unreviewed (TrEMBL)"
        reviewed = 0
        if 'entryType' in e and e['entryType'] == 'UniProtKB reviewed (Swiss-Prot)':
            reviewed = 1
        
        # Is this protein part of a UniProt Reference Proteome?
        # e.g. e['keywords'][1]['id'] == "KW-1185"
        # Note: only canonical entries (e.g. P04637, not P04637-2) have 'keywords', so only the canonical isoform will have refproteome=1.
        # This is actually incorrect, since UniProt lists e.g. A0A0G2KTI4-2 as being part of the zebrafish reference proteome: https://www.uniprot.org/uniprotkb?query=proteome%3AUP000000437+AND+accession%3AA0A0G2KTI4-2
        # Will need to retrieve reference proteome status via FTP.
        refproteome = 0
        if 'keywords' in e and any(k['id'] == "KW-1185" for k in e['keywords']):
            refproteome = 1

        # # Skip subspecies (e.g. taxon https://www.uniprot.org/taxonomy/69293, Gasterosteus aculeatus, will also return all of its subspecies such as https://www.uniprot.org/taxonomy/481459)
        # Actually keeping these now, since this is the official UniProt annotation for them. I won't actually be using these taxon IDs anywhere (definitely not for Compara).
        # if proteins_to_get == "compara_species":
        #     if tax not in taxa_set:
        #         Log("skipping undesired subspecies tax", tax)
        #         continue

        # Get full protein name
        fullname = ''
        # Recommended name, if available
        if 'recommendedName' in e['proteinDescription']:
            fullname = e['proteinDescription']['recommendedName']['fullName']['value']
        # Or: Try first submitted name from submissionNames ([0])
        elif 'submissionNames' in e['proteinDescription']:
            fullname = e['proteinDescription']['submissionNames'][0]['fullName']['value']

        # Get base canonical accession from acc if isoform (e.g. P04637-6 to P04637)
        m = rx(r'^(\w{6,10})(-\d+)?$', acc)
        canon = acc
        if m:
            if m[1]:
                canon = m[0]

        # Get UniProt species mnemonic from name (e.g. P53_HUMAN)
        m = rx(r'^\w+_(\w{3,5})$', name)
        if m:
            species = m[0]
        else:
            Die(f"Error: Could not extract species mnemonic from UniProt name '{name}'")

        # Get species names        
        species_common = ''
        if 'commonName' in e['organism']:
            species_common = e['organism']['commonName']
        else:
            Log("no common name for species (kept)", species)
        species_latin = ''
        if 'scientificName' in e['organism']:
            species_latin = e['organism']['scientificName']
        else:
            Log("no scientific name for species (kept)", species)

        # Unpack gene symbols, synonyms and function comments
        symbols = []
        synonyms = []
        comments = []
        if 'genes' in e and len(e['genes']) > 0:
            # for gene in e['genes']:
            #     if 'geneName' in gene:
            #         symbols.append(gene['geneName']['value'])
            #     else:
            #         # e.g. [{'orderedLocusNames': [{'value': 'At4g29180'}], 'orfNames': [{'value': 'F19B15.210'}]}]
            #         Log("no gene symbol for gene for acc (kept)", acc)
            #     if 'synonyms' in gene:
            #         for synonym in gene['synonyms']:
            #             synonyms.append(synonym['value'])
            symbols = [gene['geneName']['value'] for gene in e['genes'] if 'geneName' in gene]
            synonyms = [synonym['value'] for gene in e['genes'] if 'synonyms' in gene for synonym in gene['synonyms']]
        else:
            Log("no gene symbol for acc (kept)", acc)
        if 'comments' in e and len(e['comments']) > 0:
            comments = [text['value'] for comment in e['comments'] if 'comments' in e if comment['commentType'] == 'FUNCTION' for text in comment['texts']]
        else:
            Log("no function comment for acc (kept)", acc)

        # Check if any symbols already contain '|'
        for s in symbols:
            if '|' in s:
                Log("Warning: gene symbol contained '|' (kept)", s)
        for s in synonyms:
            if '|' in s:
                Log("Warning: gene synonym contained '|' (kept)", s)

        # Merge lists into strings (separated by |)
        symbols = '|'.join(symbols)
        synonyms = '|'.join(synonyms)

        # Merge function comments into a single string (separated by linebreak)
        # Make sure each comment ends in a period
        for c in range(len(comments)):
            if not comments[c].endswith('.'):
                comments[c] += '.'
        comments = '\n'.join(comments)

        # print(f" >> acc >> {acc}")
        # print(f" >> canon >> {canon}")
        # print(f" >> name >> {name}")
        # print(f" >> fullname >> {fullname}")
        # print(f" >> seqlen >> {seqlen}")
        # print(f" >> seq >> {seq}")
        # print(f" >> tax >> {tax}")
        # print(f" >> species >> {species}")
        # print(f" >> species_common >> {species_common}")
        # print(f" >> species_latin >> {species_latin}")
        # print(f" >> symbols >> {symbols}")
        # print(f" >> synonyms >> {synonyms}")
        # print(f" >> comments >> {comments}")

        # CREATE TABLE `alphauniprot` (
        #   `acc` char(13) NOT NULL,
        #   `canon` char(10) NOT NULL,
        #   `name` char(11) DEFAULT NULL,
        #   `fullname` varchar(250) DEFAULT NULL,
        #   `tax` mediumint DEFAULT NULL,
        #   `species` varchar(20) DEFAULT NULL,
        #   `species_common` varchar(50) DEFAULT NULL,
        #   `species_latin` varchar(50) DEFAULT NULL,
        #   `symbols` varchar(250) DEFAULT NULL,
        #   `synonyms` varchar(250) DEFAULT NULL,
        #   `func` varchar(15000) DEFAULT NULL,
        #   `seqlen` mediumint DEFAULT NULL,
        #   `seq` varchar(36000) DEFAULT NULL,
        #   PRIMARY KEY (`acc`),
        #   KEY `Canon` (`canon`),
        #   KEY `Name` (`name`),
        #   KEY `Species` (`species`),
        #   KEY `Seq` (`seq`(50)),
        #   KEY `Tax` (`tax`)
        # ) ENGINE=InnoDB DEFAULT CHARSET=latin1 COMMENT='UniProt annotation via API';

        # Insert into table
        q = f"INSERT INTO {alphauniprot} SET acc='{acc}', canon='{canon}', name='{name}', fullname='{Esc(fullname)}', tax='{tax}', species='{species}', species_common='{Esc(species_common)}', species_latin='{Esc(species_latin)}', reviewed='{reviewed}', refproteome='{refproteome}', symbols='{Esc(symbols)}', synonyms='{Esc(synonyms)}', func='{Esc(comments)}', seqlen='{seqlen}', seq='{seq}'"
        q = q.replace("=''", "=NULL")
        if not Switch('debug'):
            Query(q)
        else:
            State(q)

        Log("successfully inserted uniprot annotation for acc", acc)
        Log("successfully inserted uniprot annotation for canon", canon)
        Log("successfully inserted uniprot annotation for name", name)
        Log("successfully inserted uniprot annotation for tax", tax)
        Log("successfully inserted uniprot annotation for species", species)

        # # Stop after 100
        # if i >= 100:
        #     break


# print()
Show(lim=50, sort=True)

# Successfully finished: update locally recorded UniProt release version
# update_local_uniprot_release()
check_uniprot_release(update = 1)

if not Switch('debug'):
    Optimize(alphauniprot)

# SELECT COUNT(DISTINCT species), COUNT(DISTINCT tax), COUNT(DISTINCT species_common), COUNT(DISTINCT species_latin), COUNT(DISTINCT species, tax), COUNT(DISTINCT species_latin, tax), COUNT(DISTINCT species_latin, species, tax) FROM alphauniprot;

print("\nDone!")
