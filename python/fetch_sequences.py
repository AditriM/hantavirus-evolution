#!/usr/bin/env python3
"""
STEP 1: download the hantavirus sequences from NCBI.
 
  1. read our email, API key and search query from config.yaml (other people need their own email ad api key, rn this is using mines for local testing)
  2. skip everything if we already downloaded the data
  3. ask NCBI how many records match
  4. download them 200 at a time
  5. save it all to data/raw/hantavirus_raw.gb
 
Run with: python scripts/fetch_sequences.py
"""
 
import sys
import time
from pathlib import Path
 
import yaml
from Bio import Entrez
 
 
# ---- settings we might want to change ----
 
# True means download again even if we already have the file
FORCE_REDOWNLOAD = False
 
LIMIT = None
 
# How many to ask for at once. NCBI times out if we ask for too much.
BATCH_SIZE = 200

 
PROJECT_FOLDER = Path(__file__).resolve().parent.parent
CONFIG_FILE = PROJECT_FOLDER / "config.yaml"
OUTPUT_FILE = PROJECT_FOLDER / "data" / "raw" / "hantavirus_raw.gb"
 
 
# ---- step 1, read the config ----
 
# Th email and key live in a separate file so they don't reach GitHub
if not CONFIG_FILE.exists():
    print("Can't find config.yaml")
    print("Copy config.example.yaml to config.yaml and fill it in.")
    sys.exit()
 
config = yaml.safe_load(open(CONFIG_FILE))
 
if not config.get("email"):
    print("config.yaml needs an email. NCBI won't let us search without one.")
    sys.exit()
 
query = config.get("query")
if not query:
    print("config.yaml needs a query.")
    sys.exit()
 
# Biopython wants these set before we make any requests
Entrez.email = config["email"]
if config.get("api_key"):
    Entrez.api_key = config["api_key"]
else:
    print("No API key, so NCBI limits us to 3 requests a second.")
 
 
# ---- step 2, stop early if we already have the data ----
 
# Keeps re runs fast and means our demo survives bad wifi
if OUTPUT_FILE.exists() and not FORCE_REDOWNLOAD:
    size = OUTPUT_FILE.stat().st_size / 1000000
    print(f"We already have {OUTPUT_FILE.name} ({size:.1f} MB), skipping.")
    print("Set FORCE_REDOWNLOAD to True up top to get it again.")
    sys.exit()
 
 
# ---- step 3, ask NCBI what matches ----
 
print(f"Searching NCBI for: {query}")
 
# usehistory keeps the result list on NCBI's server and gives us a ticket, way easier than us passing 11000 ID numbers back and forth
search_handle = Entrez.esearch(db="nucleotide", term=query,
                               retmax=0, usehistory="y")
search_results = Entrez.read(search_handle)
search_handle.close()
 
total = int(search_results["Count"])
webenv = search_results["WebEnv"]
query_key = search_results["QueryKey"]
 
print(f"NCBI says there are {total} matching records.")
 
if total == 0:
    print("Got nothing back. Try the query on the NCBI website first.")
    sys.exit()
 
if LIMIT is not None and LIMIT < total:
    print(f"Only grabbing {LIMIT} because LIMIT is set.")
    total = LIMIT
 
 
# ---- step 4, download in batches ----
 
OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
 
# Writing to a temp file and renaming at the end, so a crashed download can't leave us with half a file that looks finished
temp_file = OUTPUT_FILE.with_suffix(".partial")
out = open(temp_file, "w")
 
for start in range(0, total, BATCH_SIZE):
    end = min(start + BATCH_SIZE, total)
    print(f"  downloading {start + 1} to {end} of {total}")
 
    # NCBI drops connections sometimes, so give each batch a few tries
    downloaded = False
    for attempt in [1, 2, 3]:
        try:
            fetch_handle = Entrez.efetch(db="nucleotide", rettype="gb",
                                         retmode="text", retstart=start,
                                         retmax=BATCH_SIZE, webenv=webenv,
                                         query_key=query_key)
            out.write(fetch_handle.read())
            fetch_handle.close()
            downloaded = True
            break
        except Exception as error:
            print(f"    try {attempt} failed: {error}")
            time.sleep(attempt * 5)
 
    if not downloaded:
        out.close()
        print("Couldn't get that batch after 3 tries, stopping here.")
        sys.exit()
 
    time.sleep(0.4)  # don't hammer their servers
 
out.close()
 
# Download finished, so rename it to the real filename
temp_file.rename(OUTPUT_FILE)
 
 
# ---- done ----
 
size = OUTPUT_FILE.stat().st_size / 1000000
print(f"\nDone. Saved {total} records to {OUTPUT_FILE} ({size:.1f} MB)")
print("Next step: python scripts/clean_sequences.py")
