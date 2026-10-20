#!/usr/bin/env python3
"""
STEP 2: clean up the GenBank records and sort them by segment.
 
  1. read every record out of data/raw/hantavirus_raw.gb
  2. throw out patents, duplicates and anything too short
  3. work out if each one is the S, M or L segment
  4. pull out the host animal and the collection year
  5. save one FASTA per segment plus a metadata table
  6. print a summary of what got dropped and why
 
Split by segment because the three segments swap between strains, so they
don't share one family tree.
 
Run with: python scripts/clean_sequences.py
"""
 
import csv
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
 
from Bio import SeqIO
 
 
# ---- settings ----
 
# Nothing under this is worth keeping at all
MIN_LENGTH = 400
 
# Proper minimums, used once we know the segment. A 400bp bit of the L
# gene is only about 6 percent of a real L segment, so lining it up
# against full genomes just makes a mess of gaps. Picked these off the
# length distribution in our own data.
MIN_LENGTH_BY_SEGMENT = {"S": 1000, "M": 2500, "L": 4000}
 
# Rough segment sizes, only used as a last ditch guess
SEGMENT_RANGES = {"S": (1200, 2400), "M": (3000, 4000), "L": (5500, 7000)}
 
# For the host check later. Need this many records before asking whether
# a virus is consistent about its host.
MIN_RECORDS_FOR_CHECK = 10
 
# Host genus under this share of a virus's records gets flagged
ODD_HOST_THRESHOLD = 0.10
 
 
# ---- where our files are ----
 
PROJECT_FOLDER = Path(__file__).resolve().parent.parent
RAW_FILE = PROJECT_FOLDER / "data" / "raw" / "hantavirus_raw.gb"
OUTPUT_FOLDER = PROJECT_FOLDER / "data" / "processed"
 
 
# ---- our host lookup table ----
 
# GenBank's host field is a blank box, so people type whatever. No way to
# automate this. The script prints what it couldn't match sorted by how
# often it shows up, so we read the top of that list and add the ones
# worth adding. One line here can fix a few hundred records at once.
HOST_TO_GENUS = {
    # Peromyscus, the Sin Nombre hosts in North America
    "peromyscus maniculatus": "Peromyscus",
    "peromyscus leucopus": "Peromyscus",
    "peromyscus boylii": "Peromyscus",
    "peromyscus levipes": "Peromyscus",
    "deer mouse": "Peromyscus",
    "white-footed mouse": "Peromyscus",
    # Apodemus, Hantaan and Dobrava over in Eurasia
    "apodemus agrarius": "Apodemus",
    "apodemus flavicollis": "Apodemus",
    "apodemus peninsulae": "Apodemus",
    "striped field mouse": "Apodemus",
    # Red backed voles, the Puumala hosts. People argue over whether this
    # genus is Myodes or Clethrionomys and NCBI can't decide either, it
    # files glareolus under one and rutilus under the other. Same animals,
    # so we dump them all in one bucket.
    "myodes glareolus": "Myodes",
    "clethrionomys glareolus": "Myodes",
    "myodes rufocanus": "Myodes",
    "clethrionomys rufocanus": "Myodes",
    "clethrionomys rutilus": "Myodes",
    "bank vole": "Myodes",
    # Microtus, Prospect Hill and Tula
    "microtus arvalis": "Microtus",
    "microtus pennsylvanicus": "Microtus",
    "microtus fortis": "Microtus",
    # Rattus, Seoul virus
    "rattus norvegicus": "Rattus",
    "rattus rattus": "Rattus",
    "norway rat": "Rattus",
    "brown rat": "Rattus",
    # Oligoryzomys, Andes virus in South America
    "oligoryzomys longicaudatus": "Oligoryzomys",
    "oligoryzomys flavescens": "Oligoryzomys",
    # Other rodents that turned up in the unmapped list
    "sigmodon hispidus": "Sigmodon",
    "eothenomys miletus": "Eothenomys",
    "eothenomys cachinus": "Eothenomys",
    "niviventer confucianus": "Niviventer",
    "niviventer niviventer": "Niviventer",
    "niviventer fulvescens": "Niviventer",
    "bandicota indica": "Bandicota",
    "mus musculus": "Mus",
    "mus booduga": "Mus",
    "akodon montensis": "Akodon",
    "arvicola amphibius": "Arvicola",
    "reithrodontomys megalotis": "Reithrodontomys",
    "reithrodontomys sumichrasti": "Reithrodontomys",
    "oryzomys couesi": "Oryzomys",
    "calomys laucha": "Calomys",
    "holochilus chacarius": "Holochilus",
    "hylomyscus simus": "Hylomyscus",
    "neodon fuscus": "Neodon",
    "leopoldamys edwardsi": "Leopoldamys",
    "zygodontomys brevicauda": "Zygodontomys",
    "lasiopodomys gregalis": "Lasiopodomys",
    "lagurus lagurus": "Lagurus",
    "tamias minimus": "Tamias",
    "stenocephalemys albipes": "Stenocephalemys",
    "eliurus majori": "Eliurus",
    "mesocricetus auratus": "Mesocricetus",
    "sundamys muelleri": "Sundamys",
    # Shrews, which aren't rodents. Part of why we kept the whole family.
    "sorex araneus": "Sorex",
    "sorex caecutiens": "Sorex",
    "sorex monticolus": "Sorex",
    "sorex vagrans": "Sorex",
    "sorex minutus": "Sorex",
    "crocidura lasiura": "Crocidura",
    "crocidura shantungensis": "Crocidura",
    "suncus murinus": "Suncus",
    "su. murinus": "Suncus",
    "anourosorex squamipes": "Anourosorex",
    "blarina brevicauda": "Blarina",
    "neomys fodiens": "Neomys",
    "chodsigoa hypsibia": "Chodsigoa",
    "chodsigoa smithii": "Chodsigoa",
    "myosorex zinki": "Myosorex",
    "myosorex geata": "Myosorex",
    "episoriculus caudatus": "Episoriculus",
    "episoriculus leucops": "Episoriculus",
    # Moles and a hedgehog
    "talpa europaea": "Talpa",
    "urotrichus talpoides": "Urotrichus",
    "scalopus aquaticus": "Scalopus",
    "neurotrichus gibbsii": "Neurotrichus",
    "erinaceus amurensis": "Erinaceus",
    # Bats
    "hipposideros pomona": "Hipposideros",
    "nyctalus noctula": "Nyctalus",
    "taphozous melanopogon": "Taphozous",
    "pipistrellus abramus": "Pipistrellus",
    "aselliscus stoliczkanus": "Aselliscus",
    "rhinolophus affinis": "Rhinolophus",
    "rousettus amplexicaudatus": "Rousettus",
    "hypsugo alaschanicus": "Hypsugo",
    "neoromicia tenuipinnis": "Neoromicia",
    "afronycteris nana": "Afronycteris",
    "pteropus alecto": "Pteropus",
    "mops condylurus": "Mops",
    # Fish, a gecko, a skink, two birds. Careful how we word this one:
    # these are in the FAMILY Hantaviridae but not the genus
    # Orthohantavirus. The fish and reptile ones sit in their own
    # subfamilies and everything that infects mammals sits in another.
    # So "hantaviruses infect fish" is wrong. The family just goes way
    # further than rodents.
    "perca fluviatilis": "Perca",
    "eviota zebrina": "Eviota",
    "eptatretus burgeri": "Eptatretus",
    "triacanthodes anomalus": "Triacanthodes",
    "lophius litulon": "Lophius",
    "halieutaea stellata": "Halieutaea",
    "retropinna retropinna": "Retropinna",
    "luposicya lupus": "Luposicya",
    "pogonophryne barsukovi": "Pogonophryne",
    "hemidactylus bowringii": "Hemidactylus",
    "scincella lateralis": "Scincella",
    "ardenna grisea": "Ardenna",
    "cyanocitta cristata": "Cyanocitta",
    # Patients. Real sequences but they came out of a person, so they say
    # nothing about what animal the virus actually lives in.
    "homo sapiens": "Homo",
    "human": "Homo",
}
 
# Not a real virus out of a real animal
EXCLUDE_WORDS = ["patent", "synthetic construct", "cloned vector",
                 "clone vector", "recombinant", "vaccine strain"]
 
 
# ---- small helper functions ----
 
def get_field(record, name):
    """Grab one field. All the sample info lives in the source feature."""
    for feature in record.features:
        if feature.type == "source":
            values = feature.qualifiers.get(name, [])
            if values:
                return values[0].strip()
    return ""
 
 
def should_exclude(record):
    """True if this looks like a patent or something made in a lab."""
    text = record.description.lower()
    text = text + " " + " ".join(record.annotations.get("keywords", [])).lower()
    for word in EXCLUDE_WORDS:
        if word in text:
            return True
    return False
 
 
def find_segment(record):
    """Figure out if this is S, M or L. None if we can't tell.
 
    Three tries: the segment label, then the protein name since each
    segment makes a different one, then the length as a straight guess.
    84 records out of 11660 failed all three, so it holds up fine.
    """
    # Easy case, somebody filled in the label
    for feature in record.features:
        if feature.type == "source":
            for value in feature.qualifiers.get("segment", []):
                label = value.strip().upper()
                if label in ["S", "M", "L"]:
                    return label
                if "SMALL" in label:
                    return "S"
                if "MEDIUM" in label:
                    return "M"
                if "LARGE" in label:
                    return "L"
 
    # Otherwise go by which protein it codes for
    text = record.description.lower()
    for feature in record.features:
        for key in ["gene", "product"]:
            text = text + " " + " ".join(feature.qualifiers.get(key, [])).lower()
 
    if "nucleocapsid" in text:
        return "S"
    if "glycoprotein" in text or "gpc" in text:
        return "M"
    if "polymerase" in text or "rdrp" in text:
        return "L"
 
    # Last resort, guess from how long it is
    length = len(record.seq)
    for segment in SEGMENT_RANGES:
        low, high = SEGMENT_RANGES[segment]
        if low <= length <= high:
            return segment
 
    return None
 
 
def find_genus(raw_host):
    """Messy host string in, genus name out. Blank if we can't match it."""
    if not raw_host:
        return ""
 
    text = raw_host.strip().lower()
    # Drop brackets, some say "Peromyscus maniculatus (deer mouse)"
    text = re.sub(r"\(.*?\)", "", text)
    # Cut anything after a semicolon, like "; specimen_voucher: TK150006"
    text = text.split(";")[0]
    text = re.sub(r"\s+", " ", text).strip()
 
    if text in HOST_TO_GENUS:
        return HOST_TO_GENUS[text]
 
    # Catches "Peromyscus sp." where they didn't know the species
    all_genera = list(HOST_TO_GENUS.values())
    first_word = text.split(" ")[0].capitalize()
    if first_word in all_genera:
        return first_word
 
    return ""
 
 
def find_year(raw_date):
    """Dates come in a million formats. We just want the year."""
    if not raw_date:
        return ""
    match = re.search(r"(19|20)\d\d", raw_date)
    if match:
        return match.group(0)
    return ""
 
 
# ---- step 1, check we have the raw data ----
 
if not RAW_FILE.exists():
    print(f"No raw data at {RAW_FILE}")
    print("Run this first: python scripts/fetch_sequences.py")
    sys.exit()
 
OUTPUT_FOLDER.mkdir(parents=True, exist_ok=True)
 
 
# ---- step 2, go through every record ----
 
counts = Counter()           # what got dropped where
sequences_by_segment = defaultdict(list)
already_seen = set()         # for catching duplicates
unmapped_hosts = Counter()   # host strings our table didn't know
table_rows = []
 
for record in SeqIO.parse(RAW_FILE, "genbank"):
    counts["total"] = counts["total"] + 1
 
    if should_exclude(record):
        counts["dropped_patent"] = counts["dropped_patent"] + 1
        continue
 
    sequence = str(record.seq).upper()
 
    if len(sequence) < MIN_LENGTH:
        counts["dropped_too_short"] = counts["dropped_too_short"] + 1
        continue
 
    # Those weird letters are IUPAC codes meaning "could be either base",
    # which is fine. Anything else means the record is broken.
    if set(sequence) - set("ACGTNRYKMSWBDHV-"):
        counts["dropped_bad_letters"] = counts["dropped_bad_letters"] + 1
        continue
 
    # Exact copies only. CD HIT handles the nearly identical outbreak
    # ones later.
    if sequence in already_seen:
        counts["dropped_duplicate"] = counts["dropped_duplicate"] + 1
        continue
    already_seen.add(sequence)
 
    segment = find_segment(record)
    if segment is None:
        counts["dropped_no_segment"] = counts["dropped_no_segment"] + 1
        continue
 
    # Now we know the segment we can use a real minimum length
    if len(sequence) < MIN_LENGTH_BY_SEGMENT[segment]:
        counts["dropped_short_for_segment"] += 1
        counts["short_" + segment] += 1
        continue
 
    raw_host = get_field(record, "host")
    genus = find_genus(raw_host)
    if raw_host and not genus:
        unmapped_hosts[raw_host] = unmapped_hosts[raw_host] + 1
 
    raw_date = get_field(record, "collection_date")
    year = find_year(raw_date)
 
    counts["kept"] = counts["kept"] + 1
    counts["kept_" + segment] += 1
    # Patients don't count as a useful host here
    if genus and genus != "Homo":
        counts["host_" + segment] += 1
    if year:
        counts["year_" + segment] += 1
 
    # Short FASTA headers keep MAFFT and IQ TREE happy later
    record.description = ""
    sequences_by_segment[segment].append(record)
 
    country = get_field(record, "country")
    if not country:
        # NCBI renamed this at some point so check both
        country = get_field(record, "geo_loc_name")
 
    table_rows.append({
        "accession": record.id,
        "segment": segment,
        "length": len(sequence),
        "organism": record.annotations.get("organism", ""),
        "strain": get_field(record, "strain"),
        "country": country,
        "host_raw": raw_host,
        "host_genus": genus,
        "collection_date_raw": raw_date,
        "year": year,
    })
 
 
# ---- step 3, save the files ----
 
for segment in sorted(sequences_by_segment):
    path = OUTPUT_FOLDER / (segment + ".fasta")
    SeqIO.write(sequences_by_segment[segment], path, "fasta")
    print(f"Wrote {len(sequences_by_segment[segment]):>5} records to {path}")
 
columns = ["accession", "segment", "length", "organism", "strain", "country",
           "host_raw", "host_genus", "collection_date_raw", "year"]
 
meta_path = OUTPUT_FOLDER / "metadata.csv"
meta_file = open(meta_path, "w", newline="")
writer = csv.DictWriter(meta_file, fieldnames=columns)
writer.writeheader()
writer.writerows(table_rows)
meta_file.close()
print(f"Wrote {len(table_rows):>5} rows to {meta_path}")
 
 
# ---- step 4, print the summary ----
 
print("\n" + "=" * 60)
print("FILTERING FUNNEL")
print("=" * 60)
print(f"  records retrieved          {counts['total']:>6}")
 
drop_labels = [("patent or lab made", "dropped_patent"),
               ("under 400bp", "dropped_too_short"),
               ("bad letters", "dropped_bad_letters"),
               ("exact duplicate", "dropped_duplicate"),
               ("segment unknown", "dropped_no_segment"),
               ("too short for segment", "dropped_short_for_segment")]
 
for label, key in drop_labels:
    print(f"    dropped, {label:<22} {counts[key]:>6}")
 
print(f"  records kept               {counts['kept']:>6}")
 
print("\n  the length filter, by segment:")
print(f"    {'segment':<10}{'minimum':>10}{'dropped':>10}")
for segment in ["S", "M", "L"]:
    print(f"    {segment:<10}{MIN_LENGTH_BY_SEGMENT[segment]:>10}"
          f"{counts['short_' + segment]:>10}")
 
print("\n  what survived:")
print(f"    {'segment':<10}{'kept':>8}{'w/ host':>10}{'w/ year':>10}")
for segment in ["S", "M", "L"]:
    print(f"    {segment:<10}{counts['kept_' + segment]:>8}"
          f"{counts['host_' + segment]:>10}{counts['year_' + segment]:>10}")
 
print("\n  (w/ host leaves out the patients since those say nothing about")
print("   the reservoir animal)")
 
 
# ---- step 5, sanity check the host mapping ----
 
# Each hantavirus usually sticks to one animal. So if a virus mostly
# shows up in one genus but a couple of records say otherwise, those are
# worth a look. Not saying they're wrong, just pointing at them.
 
hosts_per_virus = defaultdict(Counter)
 
for row in table_rows:
    virus = row["organism"]
    genus = row["host_genus"]
 
    if not virus or not genus or genus == "Homo":
        continue
 
    # Skip the catch all names. "Orthohantavirus sp." isn't one virus,
    # it's every unclassified one under the same label, so asking what
    # its usual host is means nothing.
    if "sp." in virus.lower() or "unclassified" in virus.lower():
        continue
 
    hosts_per_virus[virus][genus] += 1
 
flagged = []
for virus in hosts_per_virus:
    genus_counts = hosts_per_virus[virus]
    total_records = sum(genus_counts.values())
 
    if total_records < MIN_RECORDS_FOR_CHECK:
        continue
    if len(genus_counts) < 2:
        continue  # one host only, nothing to disagree with
 
    usual_host = genus_counts.most_common(1)[0][0]
 
    for genus in genus_counts:
        if genus == usual_host:
            continue
        share = genus_counts[genus] / total_records
        if share < ODD_HOST_THRESHOLD:
            flagged.append((genus_counts[genus], virus, usual_host, genus))
 
flagged.sort(reverse=True)
 
print("\n" + "=" * 60)
print("HOST MAPPING CHECK")
print("=" * 60)
print(f"Checked {len(hosts_per_virus)} virus names that had a mapped host.\n")
 
if not flagged:
    print("  Nothing looks off.")
else:
    print(f"  {len(flagged)} odd pairings to look at by hand:\n")
    print(f"    {'virus':<34}{'usual host':<16}{'odd host':<16}{'n':>4}")
    for count, virus, usual_host, odd_host in flagged[:20]:
        print(f"    {virus[:32]:<34}{usual_host:<16}{odd_host:<16}{count:>4}")
    if len(flagged) > 20:
        print(f"\n    ... and {len(flagged) - 20} more")
    print("\n  Could be real spillover, a mislabelled GenBank record, or a")
    print("  mistake in our lookup table.")
 
 
# ---- step 6, show what we couldn't map ----
 
if unmapped_hosts:
    print("\n" + "=" * 60)
    print(f"HOST STRINGS WE DIDN'T RECOGNISE ({len(unmapped_hosts)} distinct)")
    print("=" * 60)
    print("Add the common ones to HOST_TO_GENUS up top and run again.\n")
    for host, count in unmapped_hosts.most_common(25):
        print(f"  {count:>5}  {host}")