# predict_with_probs.py
import argparse
import os
import re
import warnings
import torch
from Bio import SeqIO
from transformers import BertModel, BertTokenizer
from clape.model import CNNOD  # Ensure this import works based on your folder structure

warnings.filterwarnings('ignore')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ligand', '-l', type=str, choices=['DNA', 'RNA', 'AB'], required=True)
    parser.add_argument('--input', '-i', required=True, help='Input FASTA')
    parser.add_argument('--output_probs', '-o', default='clape_probs.txt', help='Output file for probabilities')
    parser.add_argument('--model', '-p', required=True, help='Path to model folder')
    parser.add_argument('--cache', '-c', default='protbert', help='Cache dir')
    args = parser.parse_args()

    # 1. Load Data
    seq_ids, seqs = [], []
    for seq in SeqIO.parse(args.input, 'fasta'):
        seq_ids.append(seq.id)
        seqs.append(str(seq.seq))
    
    # 2. Load Models
    print("Loading ProtBert...")
    tokenizer = BertTokenizer.from_pretrained("Rostlab/prot_bert", do_lower_case=False, cache_dir=args.cache)
    pretrain_model = BertModel.from_pretrained("Rostlab/prot_bert", cache_dir=args.cache)
    
    print(f"Loading {args.ligand} Model...")
    predictor = CNNOD()
    model_path = os.path.join(args.model, args.ligand + ".pth")
    predictor.load_state_dict(torch.load(model_path))
    predictor.eval()

    # 3. Predict & Save Probabilities
    print(f"Predicting and saving probabilities to {args.output_probs}...")
    
    with open(args.output_probs, 'w') as f:
        for i, s in enumerate(seqs):
            # Generate Features
            sequence_Example = ' '.join(s)
            sequence_Example = re.sub(r"[UZOB]", "X", sequence_Example)
            encoded_input = tokenizer(sequence_Example, return_tensors='pt')
            features = pretrain_model(**encoded_input).last_hidden_state.squeeze(0)[1:-1, :]
            features = features.detach().unsqueeze(0)
            
            # Get Probabilities (The critical part)
            # The model outputs [Class 0 Prob, Class 1 Prob] -> We take Class 1
            probs = predictor(features).squeeze(0).detach().numpy()[:, 1]
            
            # Convert numpy floats to space-separated string
            prob_str = " ".join([f"{p:.4f}" for p in probs])
            
            # Write 3-line format: >ID, Sequence, Probabilities
            f.write(f">{seq_ids[i]}\n")
            f.write(f"{s}\n")
            f.write(f"{prob_str}\n")

    print("Done!")

if __name__ == "__main__":
    main()