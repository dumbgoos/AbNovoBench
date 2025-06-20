import re
import argparse
import logging
import sys
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd
from spectrum_utils.utils import mass_diff
from sklearn.metrics import auc

# --- Mass constants ---
mass_H = 1.0078
mass_H2O = 18.0106
mass_NH3 = 17.0265
mass_N_terminus = 1.0078
mass_C_terminus = 17.0027
mass_CO = 27.9949
mass_Phosphorylation = 79.96633

# Monoisotopic masses for amino acids and modifications
mass_AA = {
    "_PAD": 0.0,
    "_GO": mass_N_terminus - mass_H,
    "_EOS": mass_C_terminus + mass_H,
    "A": 71.03711,
    "R": 156.10111,
    "N": 114.04293,
    "n": 115.02695,  # modified N
    "D": 115.02694,
    "C": 160.03065,
    "E": 129.04259,
    "Q": 128.05858,
    "q": 129.0426,   # modified Q
    "G": 57.02146,
    "H": 137.05891,
    "I": 113.08406,
    "L": 113.08406,
    "K": 128.09496,
    "M": 131.04049,
    "m": 147.0354,  # oxidized M
    "F": 147.06841,
    "P": 97.05276,
    "S": 87.03203,
    "T": 101.04768,
    "W": 186.07931,
    "Y": 163.06333,
    "V": 99.06841,
    "d": 25.980265,  # deamidated D
    "e": -17.026549, # loss of NH3 from E
    "f": 43.005814,
    "g": 42.010565,
    "p": 111.032,
}

# List of PTM tokens to track
ptm_list = ['n', 'q', 'C', 'm', 'd', 'e', 'f', 'g', 'p']

# Tools to evaluate
tools_list = [
    "PepNet", "pNovo3", "InstaNovo", "CasanovoV1", "ContraNovo",
    "PGPointNovo", "AdaNovo", "SMSNet", "DeepNovo", "PointNovo",
    "π-PrimeNovo", "CasanovoV2", "π-HelixNovo"
]


def split_peptide(peptide: str, aa_dict: Dict[str, float]) -> List[str]:
    """
    Split a peptide string into amino acid tokens based on aa_dict keys.
    Raises ValueError if unmatched characters are found.
    """
    tokens = sorted(aa_dict.keys(), key=len, reverse=True)
    pattern = '|'.join(map(re.escape, tokens))
    parts = re.findall(pattern, peptide)
    if ''.join(parts) != peptide:
        raise ValueError(f"Sequence contains invalid tokens: '{peptide}'")
    return parts


def aa_match_prefix(
    p1: List[str], p2: List[str], aa_dict: Dict[str,float], ptm_list: List[str],
    cum_thresh: float = 0.5, ind_thresh: float = 0.1
) -> Tuple[np.ndarray, bool, np.ndarray, np.ndarray]:
    """
    Match prefix of two peptides by cumulative mass. Return arrays of matches and full-match flag.
    """
    max_len = max(len(p1), len(p2))
    matches = np.zeros(max_len, bool)
    ptm1 = np.zeros(max_len, bool)
    ptm2 = np.zeros(max_len, bool)
    i1 = i2 = 0
    cum1 = cum2 = 0.0
    while i1 < len(p1) and i2 < len(p2):
        m1 = aa_dict.get(p1[i1], 0)
        m2 = aa_dict.get(p2[i2], 0)
        if abs(mass_diff(cum1 + m1, cum2 + m2, True)) < cum_thresh:
            idx = max(i1, i2)
            if abs(mass_diff(m1, m2, True)) < ind_thresh:
                matches[idx] = True
                ptm1[idx] = p1[i1] in ptm_list
                ptm2[idx] = p2[i2] in ptm_list
            cum1 += m1; cum2 += m2; i1 += 1; i2 += 1
        elif cum2 + m2 > cum1 + m1:
            cum1 += m1; i1 += 1
        else:
            cum2 += m2; i2 += 1
    full = matches.all()
    return matches, full, ptm1, ptm2


def aa_match(
    p1: List[str], p2: List[str], aa_dict: Dict[str,float], ptm_list: List[str],
    cum_thresh: float = 0.5, ind_thresh: float = 0.1, mode: str = "best"
) -> Tuple[np.ndarray, bool, np.ndarray, np.ndarray]:
    """
    Compute AA-level matches for two peptides. Uses prefix and then suffix matching if needed.
    """
    matches, full, ptm1, ptm2 = aa_match_prefix(p1, p2, aa_dict, ptm_list, cum_thresh, ind_thresh)
    if full:
        return matches, True, ptm1, ptm2
    i1, i2 = len(p1)-1, len(p2)-1
    stop = np.argwhere(~matches)[0][0]
    cum1 = cum2 = 0.0
    while i1 >= stop and i2 >= stop:
        m1 = aa_dict.get(p1[i1], 0); m2 = aa_dict.get(p2[i2], 0)
        if abs(mass_diff(cum1 + m1, cum2 + m2, True)) < cum_thresh:
            idx = max(i1, i2)
            if abs(mass_diff(m1, m2, True)) < ind_thresh:
                matches[idx] = True
                ptm1[idx] = p1[i1] in ptm_list
                ptm2[idx] = p2[i2] in ptm_list
            cum1 += m1; cum2 += m2; i1 -= 1; i2 -= 1
        elif cum2 + m2 > cum1 + m1:
            cum1 += m1; i1 -= 1
        else:
            cum2 += m2; i2 -= 1
    return matches, matches.all(), ptm1, ptm2


def aa_match_batch(
    gts: Iterable[str], preds: Iterable[str], aa_dict: Dict[str,float],
    ptm_list: List[str], cum_thresh: float, ind_thresh: float, mode: str
) -> Tuple[List[Tuple[np.ndarray,bool,np.ndarray,np.ndarray]], int, int, int, int, int]:
    """
    Batch comparison of ground-truth and predicted peptide lists.

    Returns:
    - batch of (matches, full_match, ptm1, ptm2)
    - number of peptides
    - total true AA count
    - total predicted AA count
    - total true PTM count
    - total predicted PTM count
    """
    batch = []
    n_pep = len(list(gts))
    n_aa_true = sum(len(split_peptide(x, aa_dict)) for x in gts)
    n_aa_pred = 0
    n_ptm_true = n_ptm_pred = 0
    for gt, pr in zip(gts, preds):
        seq_gt = split_peptide(gt, aa_dict)
        seq_pr = split_peptide(pr, aa_dict) if isinstance(pr, str) else []
        n_aa_pred += len(seq_pr)
        n_ptm_true += sum(aa in ptm_list for aa in seq_gt)
        n_ptm_pred += sum(aa in ptm_list for aa in seq_pr)
        if not seq_pr:
            batch.append((np.zeros(len(seq_gt), bool), False,
                          np.zeros(len(seq_gt), bool), np.zeros(len(seq_gt), bool)))
        else:
            batch.append(aa_match(seq_gt, seq_pr, aa_dict, ptm_list, cum_thresh, ind_thresh, mode))
    return batch, n_pep, n_aa_true, n_aa_pred, n_ptm_true, n_ptm_pred


def aa_match_metrics(
    batch: List[Tuple[np.ndarray,bool,np.ndarray,np.ndarray]],
    n_pep_true: int, n_aa_true: int, n_aa_pred: int,
    n_ptm_true: int, n_ptm_pred: int, scores: List[float]
) -> Dict[str, float]:
    """
    Compute evaluation metrics:
    - aa_precision, aa_recall
    - pep_precision, pep_recall
    - ptm_precision, ptm_recall
    - peptide-level PR AUC
    """
    # AA
    n_aa_corr = sum(m[0].sum() for m in batch)
    aa_precision = n_aa_corr / (n_aa_pred + 1e-8)
    aa_recall = n_aa_corr / (n_aa_true + 1e-8)
    # Peptide
    n_pep_corr = sum(m[1] for m in batch)
    pep_precision = n_pep_corr / (len(batch) + 1e-8)
    pep_recall = n_pep_corr / (n_pep_true + 1e-8)
    # PTM
    ptm_recall = sum(m[2].sum() for m in batch) / (n_ptm_true + 1e-8)
    ptm_precision = sum(m[3].sum() for m in batch) / (n_ptm_pred + 1e-8)
    # PR AUC
    bools = [m[1] for m in batch]
    combined = sorted(zip(scores, bools), key=lambda x: x[0], reverse=True)
    sorted_bools = [b for _, b in combined]
    prec_curve = np.cumsum(sorted_bools) / np.arange(1, len(sorted_bools) + 1)
    rec_curve = np.cumsum(sorted_bools) / n_pep_true
    curve_auc = auc(rec_curve, prec_curve)
    return {
        'aa_precision': aa_precision,
        'aa_recall': aa_recall,
        'pep_precision': pep_precision,
        'pep_recall': pep_recall,
        'ptm_precision': ptm_precision,
        'ptm_recall': ptm_recall,
        'curve_auc': curve_auc
    }


def setup_logger():
    """Initialize logging format with timestamps and levels."""
    fmt = "%Y-%m-%d %H:%M:%S" + " [%(levelname)s] %(message)s"
    logging.basicConfig(level=logging.INFO, format="%(asctime)s" + " [%(levelname)s] %(message)s")


if __name__ == '__main__':
    setup_logger()
    parser = argparse.ArgumentParser(description='Evaluate peptide predictions from summary CSV.')
    parser.add_argument('-i', '--input', required=True, help='Path to summary CSV file')
    parser.add_argument('-o', '--output', required=True, help='Path to output file (with .csv or .json extension)')
    parser.add_argument('--prefix-mass-threshold', type=float, default=0.5, help='Cumulative mass match threshold')
    parser.add_argument('--ind-mass-threshold', type=float, default=0.1, help='Individual mass match threshold')
    args = parser.parse_args()

    logging.info('Reading summary file: %s', args.input)
    df = pd.read_csv(args.input)
    gts = df['Modified Sequence'].tolist()
    n_pep_true = len(gts)

    results = []
    for tool in tools_list:
        pep_col = f"{tool} Peptide"
        score_col = f"{tool} Score"
        if pep_col not in df or score_col not in df:
            logging.warning('Skipping %s: missing columns', tool)
            continue
        logging.info('Processing tool: %s', tool)
        preds = df[pep_col].tolist()
        scores = df[score_col].tolist()
        batch, _, n_aa_true, n_aa_pred, n_ptm_true, n_ptm_pred = aa_match_batch(
            gts, preds, mass_AA, ptm_list,
            args.prefix_mass_threshold, args.ind_mass_threshold, 'best'
        )
        metrics = aa_match_metrics(
            batch, n_pep_true, n_aa_true, n_aa_pred, n_ptm_true, n_ptm_pred, scores
        )
        metrics['tool'] = tool
        results.append(metrics)

    out_df = pd.DataFrame(results)
    ext = args.output.rsplit('.',1)[-1].lower()
    if ext == 'csv':
        out_df.to_csv(args.output, index=False)
    elif ext == 'json':
        out_df.to_json(args.output, orient='records', lines=True)
    else:
        logging.error('Unsupported output extension: %s', ext)
        sys.exit(1)

    logging.info('Results saved to %s', args.output)
