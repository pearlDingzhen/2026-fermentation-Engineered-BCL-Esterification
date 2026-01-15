"""
Step 4: Aggregate Scores

This script combines correlation data with LigandMPNN scores to generate
a final mutation ranking.

Usage:
    python 04_aggregate_scores.py --correlation_dir <dir> --processed_dir <dir> --output_file <file>
    
Example:
    python 04_aggregate_scores.py \
        --correlation_dir data/outputs/correlation \
        --processed_dir data/outputs/processed \
        --output_file data/outputs/correlation/final_ranking.out
"""

import argparse
import pandas as pd
import numpy as np
import os
from pathlib import Path
from typing import Dict, Tuple


# Amino acid index for mutation lookup
AMINO_ACID_INDEX = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 
                    'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
AA_TO_INDEX = {aa: idx for idx, aa in enumerate(AMINO_ACID_INDEX)}


def cal_value_score(value: float, max_value: float = 5.0) -> float:
    """
    Calculate value score as percentage of maximum.
    
    Formula: Score = (value / max_value) * 100
    """
    return (value / max_value) * 100


def cal_correlation_score(corr: float, max_corr: float = 0.6, 
                          min_corr: float = -0.7) -> float:
    """
    Calculate correlation score.
    
    Lower correlation is better (indicates mutation sensitivity).
    
    Formula: Score = (max - corr) / (max - min)^0.7 * 100
    """
    return (max_corr - corr) / ((max_corr - min_corr) ** 0.7) * 100


def process_correlation_file(file_path: str) -> Tuple[str, str, float, float, str]:
    """
    Parse a correlation output file.
    
    Returns:
        Tuple of (wt_aa, position, correlation, mpnn_score, mutant_aa)
    """
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    # Find data lines (after header)
    for i, line in enumerate(lines):
        if line.startswith('Wild AA'):
            data_start = i + 1
            break
    else:
        return None
    
    results = []
    for line in lines[data_start:]:
        if not line.strip():
            continue
        parts = line.strip().split(',')
        if len(parts) >= 5:
            wt_aa = parts[0]
            position = int(parts[1])
            correlation = float(parts[2])
            mpnn_score = float(parts[3])
            mutant_aa = parts[4]
            results.append((wt_aa, position, correlation, mpnn_score, mutant_aa))
    
    return results


def load_score_matrices(processed_dir: str) -> Dict[str, pd.DataFrame]:
    """
    Load all processed score matrices.
    
    Returns:
        Dictionary mapping filename to DataFrame
    """
    matrices = {}
    for root, dirs, files in os.walk(processed_dir):
        for f in files:
            if f.endswith('_default.csv'):
                path = os.path.join(root, f)
                df = pd.read_csv(path, header=None, skiprows=1, nrows=20,
                                usecols=range(1, 285))
                matrices[f.replace('_default.csv', '')] = df
    return matrices


def aggregate_scores(correlation_dir: str, score_matrices: Dict[str, pd.DataFrame]) -> Dict:
    """
    Aggregate scores across all comparison pairs.
    
    Returns:
        Dictionary mapping mutation key to (total_score, details_list)
    """
    final_score_list = {}
    
    # Process each correlation file
    for f in os.listdir(correlation_dir):
        if not f.startswith('output_') or not f.endswith('.csv'):
            continue
        
        file_path = os.path.join(correlation_dir, f)
        results = process_correlation_file(file_path)
        
        if not results:
            continue
        
        # Extract system names from filename
        name_parts = f.replace('output_', '').replace('.csv', '')
        
        for wt_aa, position, correlation, mpnn_score, mutant_aa in results:
            # Calculate individual scores
            corr_score = cal_correlation_score(correlation)
            value_score = cal_value_score(mpnn_score)
            
            # Combined score (weighted)
            individual_score = (corr_score * 5 + value_score) / 6
            
            # Skip if no mutation
            if wt_aa == mutant_aa:
                continue
            
            # Create unique mutation key
            key = f"{wt_aa}{position}{mutant_aa}"
            
            # Format details
            details = f"{name_parts}: {wt_aa} {position} {correlation:.3f} {mpnn_score:.2f} {mutant_aa}"
            
            if key not in final_score_list:
                final_score_list[key] = [individual_score, [details]]
            else:
                final_score_list[key][0] += individual_score
                final_score_list[key][1].append(details)
    
    return final_score_list


def save_final_ranking(final_scores: Dict, output_file: str) -> None:
    """
    Save final mutation ranking to file.
    """
    # Sort by score descending
    sorted_items = sorted(final_scores.items(), 
                         key=lambda x: x[1][0], reverse=True)
    
    with open(output_file, 'w') as f:
        f.write("=" * 60 + "\n")
        f.write("FINAL MUTATION RANKING\n")
        f.write("=" * 60 + "\n\n")
        f.write("Scoring Formula:\n")
        f.write("  Final Score = (CorrScore × 5 + ValueScore) / 6\n")
        f.write("  CorrScore = (0.6 - corr) / (0.6 - (-0.7))^0.7 × 100\n")
        f.write("  ValueScore = (MPNN_score / 5.0) × 100\n\n")
        f.write("-" * 60 + "\n\n")
        
        for rank, (key, (total_score, details_list)) in enumerate(sorted_items[:50], 1):
            f.write(f"Rank {rank}: {key}, Total Score={total_score:.2f}\n")
            f.write("  Supporting evidence:\n")
            for detail in details_list:
                f.write(f"    - {detail}\n")
            f.write("\n")
        
        f.write("-" * 60 + "\n")
        f.write(f"Total mutations analyzed: {len(sorted_items)}\n")
        f.write("=" * 60 + "\n")
    
    print(f"Final ranking saved to: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Aggregate mutation scores and generate final ranking"
    )
    parser.add_argument("--correlation_dir", type=str, required=True,
                        help="Directory containing correlation analysis results")
    parser.add_argument("--processed_dir", type=str, required=True,
                        help="Directory containing processed score matrices")
    parser.add_argument("--output_file", type=str, required=True,
                        help="Output file for final ranking")
    
    args = parser.parse_args()
    
    print("Step 4: Aggregating mutation scores")
    print("-" * 50)
    
    # Load score matrices
    print("Loading score matrices...")
    score_matrices = load_score_matrices(args.processed_dir)
    print(f"  Loaded {len(score_matrices)} score matrices")
    
    # Aggregate scores
    print("\nAggregating scores across comparisons...")
    final_scores = aggregate_scores(args.correlation_dir, score_matrices)
    print(f"  Total unique mutations: {len(final_scores)}")
    
    # Save final ranking
    print("\nSaving final ranking...")
    save_final_ranking(final_scores, args.output_file)
    
    print("\nDone!")


if __name__ == "__main__":
    main()
