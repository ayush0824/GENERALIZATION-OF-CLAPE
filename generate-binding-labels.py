import requests
import re
import random
import time
from diskcache import Cache
from typing import List

session = requests.Session()
from diskcache import Cache
cache = Cache("./pdb_year_cache", size_limit=10**9)

def get_years_from_rcsb_batch(pdb_ids):
    url = "https://data.rcsb.org/graphql"
    ids_formatted = '"' + '","'.join(pdb_ids) + '"'
    query = f"""
    {{
      entries(entry_ids: [{ids_formatted}]) {{
        rcsb_id
        rcsb_accession_info {{
          initial_release_date
        }}
      }}
    }}
    """
    try:
        response = session.post(url, json={'query': query}, timeout=10)
        if response.status_code == 200:
            data = response.json()
            result_map = {}
            if "data" in data and "entries" in data["data"]:
                for entry in data["data"]["entries"]:
                    pdb = entry["rcsb_id"].lower()
                    date_str = entry["rcsb_accession_info"]["initial_release_date"]
                    if date_str:
                        result_map[pdb] = int(date_str.split("-")[0])
            return result_map
        else:
            print(f"Batch request failed: {response.status_code}")
            return {}
    except Exception as e:
        print(f"Batch request error: {e}")
        return {}

def get_pdb_release_year(pdb_id):
    pdb_id = pdb_id.lower()
    if pdb_id in cache:
        val = cache[pdb_id]
        if val != "NOT_FOUND":
            return val
        return None

    results = get_years_from_rcsb_batch([pdb_id])
    if pdb_id in results:
        cache[pdb_id] = results[pdb_id]
        return results[pdb_id]

    cache[pdb_id] = "NOT_FOUND"
    return None


def generate_binding_labels(sequence: str, binding_col9: str) -> str:
    """Generates 0/1 label string aligned with sequence"""
    labels = ["0"] * len(sequence)
    if not binding_col9.strip():
        return "".join(labels)

    tokens = binding_col9.split()
    pattern = re.compile(r"([A-Za-z])(\d+)")
    
    for token in tokens:
        match = pattern.match(token)
        if not match:
            continue
        expected_char = match.group(1)
        one_based_index = int(match.group(2))
        idx = one_based_index - 1
        if not (0 <= idx < len(sequence)):
            print(f"Index {idx} out of bounds. Token: {token}")
            continue
        if expected_char.upper() != 'X' and sequence[idx].upper() != expected_char.upper():
            raise ValueError(f"Mismatch: expected {expected_char} at {idx}, found {sequence[idx]}")
        labels[idx] = "1"
    return "".join(labels)


def parse_biolip_file(biolip_path, output_path, sample_size=15000):
    collected = []

    with open(biolip_path, "r") as f:
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) < 21:
                continue  # Incomplete line
            pdb_id = fields[0].strip()
            ligand = fields[4].strip()
            chain_id = fields[1].strip()
            binding_col9 = fields[8].strip()
            sequence = fields[20].strip()

            if ligand.upper() != "RNA":
                continue

            release_year = get_pdb_release_year(pdb_id)
            if release_year is None or release_year <= 2022:
                continue

            collected.append((pdb_id, chain_id, binding_col9, sequence))

    print(f"Collected {len(collected)} valid entries")
    sampled = random.sample(collected, min(sample_size, len(collected)))

    with open(output_path, "w") as fout:
        for pdb_id, chain_id, binding_col9, sequence in sampled:
            try:
                labels = generate_binding_labels(sequence, binding_col9)
            except ValueError as e:
                print(f"Skipping {pdb_id}{chain_id} due to label mismatch: {e}")
                continue
            # CLAPE-compatible FASTA-style block format
            fout.write(f">{pdb_id}_{chain_id}\n")
            fout.write(f"{sequence}\n")
            fout.write(f"{labels}\n")
            # If there was any previous tab-separated output, it is now removed/commented out

    print(f"Saved {len(sampled)} entries to {output_path}")


# Example usage
parse_biolip_file("BioLiP.txt", "RNA2022_1500_Test.txt", sample_size=1500)