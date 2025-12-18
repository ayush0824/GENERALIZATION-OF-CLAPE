#!/usr/bin/env python3
import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import (
    roc_curve, auc, precision_recall_curve, average_precision_score,
    roc_auc_score, f1_score
)

# --- 1. Data Loading Functions (Reused & Robust) ---
def load_true_labels(true_path):
    labels = {}
    try:
        with open(true_path) as f:
            while True:
                hdr = f.readline().strip()
                if not hdr: break
                if not hdr.startswith(">"): continue
                
                seq_id = hdr.lstrip(">").strip().upper().replace("_", "")
                f.readline() # Skip sequence
                mask = f.readline().strip()
                try:
                    labels[seq_id] = [int(c) for c in mask]
                except ValueError:
                    continue
    except FileNotFoundError:
        print(f"Error: File not found {true_path}", file=sys.stderr)
        sys.exit(1)
    return labels

def load_predictions(pred_path):
    preds = {}
    try:
        with open(pred_path) as f:
            while True:
                hdr = f.readline().strip()
                if not hdr: break
                if not hdr.startswith(">"): continue
                
                seq_id = hdr.lstrip(">").strip().upper().replace("_", "")
                f.readline() # Skip sequence
                prob_line = f.readline().strip()
                try:
                    preds[seq_id] = [float(x) for x in prob_line.split()]
                except ValueError:
                    continue
    except FileNotFoundError:
        print(f"Error: File not found {pred_path}", file=sys.stderr)
        sys.exit(1)
    return preds

def align_data(pred_file, true_file):
    """Aligns prediction and truth dictionaries into flat arrays."""
    true_dict = load_true_labels(true_file)
    pred_dict = load_predictions(pred_file)
    
    y_true_flat = []
    y_scores_flat = []
    
    # Only use intersecting IDs
    common_ids = set(true_dict.keys()).intersection(set(pred_dict.keys()))
    
    if not common_ids:
        print(f"Warning: No matching IDs found between {pred_file} and {true_file}!", file=sys.stderr)
        return np.array([]), np.array([])

    for seq_id in common_ids:
        t = true_dict[seq_id]
        p = pred_dict[seq_id]
        if len(t) == len(p):
            y_true_flat.extend(t)
            y_scores_flat.extend(p)
            
    return np.array(y_true_flat), np.array(y_scores_flat)

# --- 2. Main Logic ---
def main(datasets, output_plot):
    # Store results for plotting
    roc_data = [] # (fpr, tpr, auc, name)
    pr_data = []  # (recall, precision, auc, name)
    
    print(f"\n{'='*60}")
    print(f"{'DATASET':<20} | {'ROC-AUC':<10} | {'PR-AUC':<10} | {'Best F1':<10}")
    print(f"{'-'*60}")

    # Colors for plotting (auto-cycle)
    colors = plt.cm.tab10.colors 

    for idx, dataset_str in enumerate(datasets):
        # Parse "Name:PredFile:TruthFile"
        try:
            name, pred_file, true_file = dataset_str.split(":")
        except ValueError:
            print(f"Error: Dataset arg '{dataset_str}' format invalid. Use Name:Pred:True")
            continue

        # Get Data
        y_true, y_scores = align_data(pred_file, true_file)
        
        if len(y_true) == 0:
            continue

        # Calculate Metrics
        roc_auc = roc_auc_score(y_true, y_scores)
        pr_auc = average_precision_score(y_true, y_scores)
        
        # Calculate Curves
        fpr, tpr, _ = roc_curve(y_true, y_scores)
        precision, recall, thresholds = precision_recall_curve(y_true, y_scores)
        
        # Calculate Max F1
        p_trim, r_trim = precision[:-1], recall[:-1]
        f1_scores = 2 * (p_trim * r_trim) / (p_trim + r_trim + 1e-10)
        best_f1 = np.max(f1_scores) if len(f1_scores) > 0 else 0.0

        # Print Row
        print(f"{name:<20} | {roc_auc:.4f}     | {pr_auc:.4f}     | {best_f1:.4f}")

        # Store for Plotting
        color = colors[idx % len(colors)]
        roc_data.append((fpr, tpr, roc_auc, name, color))
        pr_data.append((recall, precision, pr_auc, name, color))

    print(f"{'='*60}\n")

    # --- 3. Plotting ---
    plt.figure(figsize=(14, 6))

    # Subplot 1: ROC Curve
    plt.subplot(1, 2, 1)
    for fpr, tpr, auc_score, name, color in roc_data:
        plt.plot(fpr, tpr, color=color, lw=2, label=f'{name} (AUC={auc_score:.3f})')
    plt.plot([0, 1], [0, 1], 'k--', lw=1, alpha=0.5)
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC Curve Comparison')
    plt.legend(loc="lower right")
    plt.grid(alpha=0.3)

    # Subplot 2: Precision-Recall Curve
    plt.subplot(1, 2, 2)
    for recall, precision, auc_score, name, color in pr_data:
        plt.plot(recall, precision, color=color, lw=2, label=f'{name} (AUC={auc_score:.3f})')
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('Precision-Recall Comparison')
    plt.legend(loc="upper right") # PR legend usually better top-right
    plt.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_plot, dpi=300)
    print(f"--> Comparison plot saved to: {output_plot}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compare ROC/PR curves for multiple datasets.",
        usage="%(prog)s --inputs Name1:Pred1:True1 Name2:Pred2:True2 ..."
    )
    parser.add_argument(
        "--inputs", 
        nargs="+", 
        required=True, 
        help="List of datasets in format: 'Name:PredFile:TrueFile' (space separated)"
    )
    parser.add_argument("--out", default="Comparison_ROC.png", help="Output filename")
    
    args = parser.parse_args()
    main(args.inputs, args.out)