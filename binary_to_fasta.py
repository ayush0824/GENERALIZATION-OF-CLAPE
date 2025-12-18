def txt_to_fasta(input_txt, output_fasta):
    with open(input_txt, 'r') as f, open(output_fasta, 'w') as out:
        lines = f.readlines()
        for i in range(0, len(lines), 4):
            pdb = lines[i].strip()
            chain = lines[i+1].strip()
            sequence = lines[i+3].strip()
            header = f">{pdb}_{chain}"
            out.write(f"{header}\n{sequence}\n")

# Convert both training and test sets
txt_to_fasta("RNA-495_Train.txt", "train.fasta")
txt_to_fasta("T-2025.txt", "test.fasta")

#this file is not being used in any pipeline, its de-functed