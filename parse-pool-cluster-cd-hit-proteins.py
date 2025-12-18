import sys
import argparse

def parse_clusters(clstr_file, original_test_file, output_file):
    print(f"--> Parsing clusters from: {clstr_file}")
    print(f"--> Source Data: {original_test_file}")
    
    valid_test_ids = set()
    
    current_cluster = []
    has_train = False
    
    try:
        with open(clstr_file, 'r') as f:
            for line in f:
                if line.startswith(">Cluster"):
                    # Process the previous cluster
                    if current_cluster:
                        # If NO training sequence in cluster, these Test seqs are safe
                        if not has_train:
                            for seq_id in current_cluster:
                                valid_test_ids.add(seq_id)
                    
                    # Reset for new cluster
                    current_cluster = []
                    has_train = False
                
                elif line.strip(): # Sequence line
                    # Example line: 0	300aa, >TRAIN_1a2b... at 90.00%
                    parts = line.split('>')
                    if len(parts) > 1:
                        # Extract ID part (remove trailing '...' or statistics)
                        full_id = parts[1].split('...')[0].strip()
                        
                        if full_id.startswith("TRAIN_"):
                            has_train = True
                        elif full_id.startswith("TEST_"):
                            # Remove the prefix we added earlier to get the clean original ID
                            clean_id = full_id.replace("TEST_", "")
                            current_cluster.append(clean_id)

            # Process the very last cluster
            if current_cluster and not has_train:
                for seq_id in current_cluster:
                    valid_test_ids.add(seq_id)
                    
    except FileNotFoundError:
        print(f"Error: Cluster file '{clstr_file}' not found.")
        sys.exit(1)

    print(f"--> Found {len(valid_test_ids)} valid Test sequences (No homology to Train).")

    # --- Reconstruct the Final File ---
    print(f"--> Retrieving sequences and labels...")
    kept = 0
    try:
        with open(original_test_file, 'r') as infile, open(output_file, 'w') as outfile:
            lines = [l.strip() for l in infile.readlines() if l.strip()]
            
            # Iterate 3-line blocks
            for i in range(0, len(lines), 3):
                if i+2 >= len(lines): break
                
                raw_header = lines[i] # >1a2b
                # Clean header to match what we stored in valid_test_ids
                clean_id = raw_header.replace(">", "").strip()
                
                if clean_id in valid_test_ids:
                    outfile.write(f"{lines[i]}\n{lines[i+1]}\n{lines[i+2]}\n")
                    kept += 1
    except FileNotFoundError:
        print(f"Error: Original test file '{original_test_file}' not found.")
        sys.exit(1)
                
    print(f"--> Success! Final dataset saved to: {output_file}")
    print(f"    Total Sequences: {kept}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter Test set based on PSI-CD-HIT clusters.")
    
    # Define arguments
    parser.add_argument("--clstr", required=True, help="Input .clstr file from psi-cd-hit")
    parser.add_argument("--test", required=True, help="Original Test file (3-line format) containing labels")
    parser.add_argument("--out", required=True, help="Output filename for the clean dataset")
    
    args = parser.parse_args()
    
    # Run the function with arguments
    parse_clusters(args.clstr, args.test, args.out)