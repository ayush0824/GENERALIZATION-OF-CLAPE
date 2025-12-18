#!/usr/bin/env python3
import argparse
import sys
from sklearn.metrics import (
    precision_score, recall_score, f1_score,
    matthews_corrcoef, roc_auc_score,
    average_precision_score, confusion_matrix
)

def load_true_labels(true_path):
    """
    Load true labels from test file with 3‑line blocks:
    >ID
    sequence
    0/1 mask
    Returns dict ID -> [0,1,...]
    """
    labels = {}
    with open(true_path) as f:
        while True:
            hdr = f.readline()
            if not hdr:
                break
            hdr = hdr.strip()
            if not hdr:
                continue
            if not hdr.startswith(">"):
                print("Warning: expected header line starting with '>' but got:", hdr, file=sys.stderr)
                continue
            raw = hdr.lstrip(">").strip()
            seq_id = raw.upper().replace("_", "")
            seq = f.readline()
            if not seq:
                print("Warning: missing sequence line after header", hdr, file=sys.stderr)
                break
            seq = seq.strip()
            mask = f.readline()
            if not mask:
                print("Warning: missing mask line after sequence for", hdr, file=sys.stderr)
                break
            mask = mask.strip()
            # optional: check mask length == len(seq)
            if len(mask) != len(seq):
                print(f"Warning: mask length {len(mask)} != sequence length {len(seq)} for {seq_id}", file=sys.stderr)
            try:
                labels[seq_id] = [int(c) for c in mask]
            except ValueError:
                print(f"Warning: non-binary character in mask for {seq_id}: {mask}", file=sys.stderr)
    return labels

def load_predictions(pred_path):
    """
    Load prediction file (3 lines per record: >ID, sequence, then 0/1 mask).
    Returns dict ID -> [0,1,...]
    """
    preds = {}
    with open(pred_path) as f:
        while True:
            hdr = f.readline()
            if not hdr:
                break
            hdr = hdr.strip()
            if not hdr:
                continue
            if not hdr.startswith(">"):
                print("Warning: expected header line starting with '>' but got:", hdr, file=sys.stderr)
                continue
            raw = hdr.lstrip(">").strip()
            seq_id = raw.upper().replace("_", "")
            seq_line = f.readline()   # ignore the sequence itself — but read to advance
            if not seq_line:
                print("Warning: missing sequence line after header", hdr, file=sys.stderr)
                break
            mask = f.readline()
            if not mask:
                print("Warning: missing mask line after sequence for", hdr, file=sys.stderr)
                break
            mask = mask.strip()
            try:
                preds[seq_id] = [int(c) for c in mask]
            except ValueError:
                print(f"Warning: non-binary character in prediction mask for {seq_id}: {mask}", file=sys.stderr)
    return preds

def main(pred_file, true_file):
    true_dict = load_true_labels(true_file)
    pred_dict = load_predictions(pred_file)

    y_true, y_pred = [], []
    missing = 0
    mismatched = 0
    matched = 0

    for seq_id, pred_mask in pred_dict.items():
        if seq_id not in true_dict:
            missing += 1
            print(f"Missing ground truth for prediction {seq_id}", file=sys.stderr)
            continue
        true_mask = true_dict[seq_id]
        if len(pred_mask) != len(true_mask):
            mismatched += 1
            print(f"Length mismatch for {seq_id}: pred {len(pred_mask)} vs true {len(true_mask)}", file=sys.stderr)
            continue
        matched += 1
        y_true.extend(true_mask)
        y_pred.extend(pred_mask)

    if matched == 0:
        print("No matching records to evaluate.", file=sys.stderr)
        return

    # Metrics
    precision = precision_score(y_true, y_pred, zero_division=0)
    recall    = recall_score(y_true, y_pred, zero_division=0)
    f1        = f1_score(y_true, y_pred, zero_division=0)
    mcc       = matthews_corrcoef(y_true, y_pred)
    try:
        roc_auc = roc_auc_score(y_true, y_pred)
        pr_auc  = average_precision_score(y_true, y_pred)
    except ValueError:
        roc_auc = pr_auc = float('nan')

    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
    total = len(y_true)
    positives = sum(y_true)
    pos_pct = positives / total * 100 if total else 0

    print(f"\nMatched sequences: {matched}")
    print(f"Missing sequences: {missing}")
    print(f"Mismatched lengths: {mismatched}")
    print(f"Total residues: {total}")
    print(f"Positive residues: {positives} ({pos_pct:.2f}%)\n")

    print(f"Precision:    {precision:.4f}")
    print(f"Recall:       {recall:.4f}")
    print(f"F1 Score:     {f1:.4f}")
    print(f"MCC:          {mcc:.4f}")
    print(f"ROC-AUC:      {roc_auc:.4f}")
    print(f"PR-AUC:       {pr_auc:.4f}\n")
    print("Confusion matrix:")
    print(f" TP: {tp}, FP: {fp}")
    print(f" TN: {tn}, FN: {fn}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--pred", required=True, help="Prediction file (3-line blocks: header, sequence, mask)")
    parser.add_argument("--true", required=True, help="True labels file (3-line blocks: header, sequence, mask)")
    args = parser.parse_args()
    main(args.pred, args.true)