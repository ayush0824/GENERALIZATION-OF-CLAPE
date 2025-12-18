#!/usr/bin/env python3
import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import (
    precision_score, recall_score, f1_score,
    matthews_corrcoef, roc_auc_score, roc_curve, auc,
    average_precision_score, precision_recall_curve, confusion_matrix
)

def load_true_labels(true_path):
    labels = {}
    print(f"--> Loading Ground Truth from: {true_path}")
    with open(true_path) as f:
        while True:
            hdr = f.readline().strip()
            if not hdr: break
            if not hdr.startswith(">"): continue
            
            # Normalize ID: >7zod_A -> 7ZODA
            seq_id = hdr.lstrip(">").strip().upper().replace("_", "")
            f.readline() # Skip sequence
            mask = f.readline().strip()
            try:
                labels[seq_id] = [int(c) for c in mask]
            except ValueError:
                print(f"Warning: Non-binary mask for {seq_id}", file=sys.stderr)
    return labels

def load_predictions(pred_path):
    preds = {}
    print(f"--> Loading Predictions from: {pred_path}")
    with open(pred_path) as f:
        while True:
            hdr = f.readline().strip()
            if not hdr: break
            if not hdr.startswith(">"): continue
            
            # Normalize ID exactly same way as above
            seq_id = hdr.lstrip(">").strip().upper().replace("_", "")
            f.readline() # Skip sequence
            prob_line = f.readline().strip()
            try:
                # Handle space-separated floats
                preds[seq_id] = [float(x) for x in prob_line.split()]
            except ValueError:
                print(f"Warning: Non-float score for {seq_id}", file=sys.stderr)
    return preds

def plot_curves(y_true, y_scores, output_file):
    fpr, tpr, _ = roc_curve(y_true, y_scores)
    roc_auc = auc(fpr, tpr)
    precision, recall, _ = precision_recall_curve(y_true, y_scores)
    pr_auc = average_precision_score(y_true, y_scores)

    plt.figure(figsize=(12, 5))
    
    # ROC Plot
    plt.subplot(1, 2, 1)
    plt.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC (AUC={roc_auc:.3f})')
    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC Curve')
    plt.legend(loc="lower right")
    plt.grid(alpha=0.3)

    # Precision-Recall Plot
    plt.subplot(1, 2, 2)
    plt.plot(recall, precision, color='blue', lw=2, label=f'PR (AUC={pr_auc:.3f})')
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('Precision-Recall')
    plt.legend(loc="lower left")
    plt.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    print(f"--> Plot saved to {output_file}")

def main(pred_file, true_file, output_plot):
    true_dict = load_true_labels(true_file)
    pred_dict = load_predictions(pred_file)

    y_true_all = []
    y_scores_all = []
    
    # COUNTERS
    total_preds = len(pred_dict)
    matched = 0
    missing_truth = 0
    mismatched_len = 0
    
    # Align Data
    for seq_id, scores in pred_dict.items():
        if seq_id not in true_dict:
            missing_truth += 1
            if missing_truth <= 3: # Debug first few missing IDs
                print(f"Debug: Pred ID '{seq_id}' not found in Truth file.", file=sys.stderr)
            continue
            
        true_labels = true_dict[seq_id]
        
        if len(scores) != len(true_labels):
            mismatched_len += 1
            if mismatched_len <= 3:
                print(f"Debug: Length mismatch {seq_id} (Pred {len(scores)} vs True {len(true_labels)})", file=sys.stderr)
            continue
            
        matched += 1
        y_true_all.extend(true_labels)
        y_scores_all.extend(scores)

    # --- DATA HEALTH REPORT ---
    print("\n" + "="*40)
    print("      DATA HEALTH SUMMARY")
    print("="*40)
    print(f"Total Predictions:      {total_preds}")
    print(f"Matched Successfully:   {matched}")
    print(f"Missing Ground Truth:   {missing_truth}")
    print(f"Length Mismatches:      {mismatched_len}")
    print("="*40)

    if matched == 0:
        print("\nCRITICAL ERROR: No sequences matched. Check ID formats.", file=sys.stderr)
        return

    # Convert to numpy
    y_true = np.array(y_true_all)
    y_scores = np.array(y_scores_all)

    # 1. Threshold-Independent Metrics
    roc_auc = roc_auc_score(y_true, y_scores)
    pr_auc = average_precision_score(y_true, y_scores)

    # 2. Find Optimal Threshold (Maximize F1)
    precisions, recalls, thresholds = precision_recall_curve(y_true, y_scores)
    
    # [FIX] Slice to exclude the last point (Recall=0) to match 'thresholds' length
    p_trim = precisions[:-1]
    r_trim = recalls[:-1]
    
    # Calculate F1 for every threshold
    f1_scores = 2 * (p_trim * r_trim) / (p_trim + r_trim + 1e-10)
    
    best_idx = np.argmax(f1_scores)
    best_thresh = thresholds[best_idx]
    best_f1 = f1_scores[best_idx]

    # 3. Binary Metrics (Default Threshold 0.5)
    y_pred_bin = (y_scores >= 0.5).astype(int)
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred_bin).ravel()

    # REPORT
    print("\n" + "="*40)
    print("      MODEL PERFORMANCE")
    print("="*40)
    print(f"Total Residues:    {len(y_true)}")
    print(f"Positives:         {sum(y_true)} ({sum(y_true)/len(y_true):.2%})")
    print("-" * 40)
    print(f"ROC-AUC:           {roc_auc:.4f}")
    print(f"PR-AUC:            {pr_auc:.4f}")
    print("-" * 40)
    print(f"Optimal F1:        {best_f1:.4f} (at Threshold {best_thresh:.4f})")
    print("-" * 40)
    print("Default Threshold (0.5):")
    print(f"  F1 Score:        {f1_score(y_true, y_pred_bin):.4f}")
    print(f"  MCC:             {matthews_corrcoef(y_true, y_pred_bin):.4f}")
    print(f"  Precision:       {precision_score(y_true, y_pred_bin):.4f}")
    print(f"  Recall:          {recall_score(y_true, y_pred_bin):.4f}")
    print(f"  Confusion:       TP={tp} FP={fp} TN={tn} FN={fn}")
    print("="*40)

    # Generate Plot
    plot_curves(y_true, y_scores, output_plot)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--pred", required=True, help="Prediction Probabilities File")
    parser.add_argument("--true", required=True, help="Ground Truth Binary File")
    parser.add_argument("--plot", default="ROC_Curve.png", help="Output Image")
    args = parser.parse_args()
    
    main(args.pred, args.true, args.plot)