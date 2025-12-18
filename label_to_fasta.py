import os

def biolip_to_fasta(input_path, output_path, dataset_type):
    """
    Converts BioLiP formats to standard FASTA for CD-HIT.
    
    Args:
        dataset_type (str): 'train' (4-line format) or 'test' (3-line format)
    """
    if not os.path.exists(input_path):
        print(f"Error: {input_path} not found.")
        return

    with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
        lines = [l.strip() for l in infile.readlines() if l.strip()]
        
        # Define stride based on file format
        # Train: >Header, Seq, Label1, Label2 (4 lines)
        # Test:  >Header, Seq, Label (3 lines)
        step = 4 if dataset_type == 'train' else 3
        
        count = 0
        for i in range(0, len(lines), step):
            if i + 1 >= len(lines): break
            
            # Line 0 is the Header (e.g. >101mA)
            header = lines[i]
            
            # Line 1 is the Sequence (e.g. MVLSE...)
            sequence = lines[i+1]
            
            # Write to standard FASTA
            outfile.write(f"{header}\n{sequence}\n")
            count += 1

    print(f"[{dataset_type.upper()}] Processed {count} sequences from {input_path} -> {output_path}")

# --- Execution ---
# 1. Process Training Set (4 lines)
biolip_to_fasta("RNA-495_Train.txt", "train.fasta", dataset_type='train')

# 2. Process Test Set (3 lines - derived from your snippet)
biolip_to_fasta("RNA_1500_Test.txt", "test.fasta", dataset_type='test')