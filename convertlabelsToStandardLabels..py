#!/usr/bin/env python3
import sys

if len(sys.argv) != 3:
    print("Usage: python3 convert_pm_to_binary.py <input_txt> <output_txt>")
    sys.exit(1)

infile = sys.argv[1]
outfile = sys.argv[2]

with open(infile, "r") as fin, open(outfile, "w") as fout:
    lines = [l.rstrip("\n") for l in fin]
    i = 0
    while i < len(lines):
        # Skip empty lines
        if not lines[i].strip():
            i += 1
            continue

        # Expect a header
        header = lines[i]
        seq = None
        mask = None

        # Some train files are 3‑line blocks: header / seq / mask
        if i + 2 < len(lines):
            seq = lines[i + 1]
            mask = lines[i + 2]
        else:
            print(f"Warning: incomplete record at line {i+1}")
            break

        # Convert + => 1, - => 0 on mask
        new_mask = mask.replace("+", "1").replace("-", "0")

        fout.write(header + "\n")
        fout.write(seq + "\n")
        fout.write(new_mask + "\n")

        i += 3

print(f"Converted binding masks and wrote to {outfile}")