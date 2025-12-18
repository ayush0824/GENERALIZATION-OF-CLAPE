import os
import subprocess

# ================= CONFIGURATION =================
TRAIN_FILE = "RNA-495_Train.txt"
TEST_FILE  = "T-2025.txt"
FINAL_OUTPUT = "T-2025_PDB30.txt"

# PDB-30 Standard Parameters
# BLAST returns both of these as percentages (0-100), not fractions.
IDENTITY_CUTOFF = 30.0  # 30%
COVERAGE_CUTOFF = 40.0  # 40%
# =================================================

def parse_original_to_fasta(input_path, output_fasta, file_type):
    """Converts original 4-line/3-line format to standard FASTA for BLAST."""
    print(f"--> Converting {input_path} to temporary FASTA...")
    step = 4 if file_type == 'train' else 3
    count = 0
    with open(input_path, 'r') as infile, open(output_fasta, 'w') as outfile:
        lines = [l.strip() for l in infile.readlines() if l.strip()]
        for i in range(0, len(lines), step):
            if i + 1 >= len(lines): break
            # Header cleanup: >101mA -> 101mA
            raw_header = lines[i].replace(">", "")
            seq = lines[i+1]
            outfile.write(f">{raw_header}\n{seq}\n")
            count += 1
    return count

def run_pipeline():
    # 1. Prepare FASTA files
    parse_original_to_fasta(TRAIN_FILE, "temp_train.fasta", 'train')
    parse_original_to_fasta(TEST_FILE, "temp_test.fasta", 'test')

    # 2. Build BLAST DB from Training Set
    print("--> Building BLAST database...")
    subprocess.run(
        "makeblastdb -in temp_train.fasta -dbtype prot -out train_db -parse_seqids",
        shell=True, check=True, stdout=subprocess.DEVNULL
    )

    # 3. Run BLASTP (Test vs Train)
    # -evalue 10: Ensures we see even weak matches (30% is weak)
    # -outfmt 6: Tabular output for easy parsing
    # qcovs: Returns INTEGER percentage (0-100)
    print(f"--> Running BLASTP to detect >{IDENTITY_CUTOFF}% homology...")
    cmd = (
        f"blastp -query temp_test.fasta -db train_db -out blast_results.tab "
        f"-evalue 10 -num_threads 4 -outfmt '6 qseqid sseqid pident qcovs'"
    )
    subprocess.run(cmd, shell=True, check=True)

    # 4. Filter (The Logic)
    print("--> Filtering results...")
    ids_to_remove = set()
    with open("blast_results.tab", 'r') as f:
        for line in f:
            cols = line.strip().split('\t')
            q_id = cols[0]
            pident = float(cols[2]) # e.g. 30.5
            qcovs = float(cols[3])  # e.g. 95.0 (Percentage, not fraction)
            
            # STRICT FILTER: Remove if identity > 30% AND it covers > 40% of the protein
            # Both metrics are now on the 0-100 scale.
            if pident >= IDENTITY_CUTOFF and qcovs >= COVERAGE_CUTOFF:
                ids_to_remove.add(q_id)

    print(f"    identified {len(ids_to_remove)} homologous sequences to remove.")

    # 5. Reconstruct Final Dataset
    kept = 0
    with open(TEST_FILE, 'r') as infile, open(FINAL_OUTPUT, 'w') as outfile:
        lines = [l.strip() for l in infile.readlines() if l.strip()]
        for i in range(0, len(lines), 3):
            if i + 2 >= len(lines): break
            header = lines[i]
            clean_id = header.replace(">", "")
            
            if clean_id not in ids_to_remove:
                # Write original 3 lines: Header, Seq, Label
                outfile.write(f"{lines[i]}\n{lines[i+1]}\n{lines[i+2]}\n")
                kept += 1

    print(f"--> SUCCESS. Final PDB-30 dataset: {FINAL_OUTPUT} ({kept} sequences).")
    
    # Cleanup
    subprocess.run("rm temp_* train_db* blast_results.tab", shell=True)

if __name__ == "__main__":
    run_pipeline()