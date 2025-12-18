import sys

def convert(input_txt, output_fasta):
    print(f"--> Converting {input_txt} to FASTA format...")
    count = 0
    with open(input_txt, 'r') as f_in, open(output_fasta, 'w') as f_out:
        lines = [l.strip() for l in f_in.readlines() if l.strip()]
        
        # Iterate 3-line blocks: Header, Seq, Mask
        for i in range(0, len(lines), 3):
            if i+1 >= len(lines): break
            
            header = lines[i]   # e.g. >7zod_A
            seq = lines[i+1]    # e.g. AVQ...
            # We ignore lines[i+2] (the mask) for the FASTA input
            
            f_out.write(f"{header}\n{seq}\n")
            count += 1
            
    print(f"--> Done. Created {output_fasta} with {count} sequences.")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python txt_to_fasta.py <input_3line.txt> <output.fasta>")
    else:
        convert(sys.argv[1], sys.argv[2])