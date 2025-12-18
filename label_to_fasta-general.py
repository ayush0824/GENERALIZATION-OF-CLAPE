import argparse
import sys
import os

def convert_3line_to_fasta(input_path, output_path):
    """
    Converts a 3-line format file (Header, Sequence, Mask) to standard FASTA.
    
    Structure Expected:
    Line 1: >Header
    Line 2: Sequence
    Line 3: Mask (00101...) -> DISCARDED
    """
    if not os.path.exists(input_path):
        print(f"Error: Input file '{input_path}' not found.", file=sys.stderr)
        sys.exit(1)

    print(f"--> Converting {input_path} to FASTA...")
    
    count = 0
    with open(input_path, 'r') as f_in, open(output_path, 'w') as f_out:
        while True:
            # 1. Read Header
            header = f_in.readline().strip()
            if not header: break # End of file
            
            # Robustness: Skip empty lines or handle broken blocks
            if not header.startswith(">"):
                continue

            # 2. Read Sequence
            seq = f_in.readline().strip()
            if not seq:
                print(f"Warning: Unexpected end of file after header {header}", file=sys.stderr)
                break

            # 3. Read Mask (and discard it)
            mask = f_in.readline() # Read but don't strip/store to save time
            
            # 4. Write to Output
            f_out.write(f"{header}\n{seq}\n")
            count += 1
            
    print(f"--> Success! Wrote {count} sequences to {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="General Tool: Convert 3-line labeled files to FASTA")
    parser.add_argument("-i", "--input", required=True, help="Input text file (3-line format)")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA file")
    
    args = parser.parse_args()
    convert_3line_to_fasta(args.input, args.output)