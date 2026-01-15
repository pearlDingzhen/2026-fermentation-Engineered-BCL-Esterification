"""
Step 2: Process LigandMPNN Scores

This script converts raw LigandMPNN output (.pt) to CSV format with
calculated mutation scores.

Usage:
    python 02_process_scores.py --input <input.pt> --output <output.csv>
    
Example:
    python 02_process_scores.py \
        --input data/outputs/raw/ethanol/ethanol_cluster0/ethanol_cluster0.pt \
        --output data/outputs/processed/ethanol/ethanol_cluster0_default.csv
"""

import argparse
import torch
import re
import pandas as pd
import numpy as np
from pathlib import Path


# Standard amino acid index (20 amino acids)
AMINO_ACID_INDEX = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 
                    'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']


def load_pt_file(file_path: str) -> dict:
    """
    Load LigandMPNN output file in .pt format.
    
    Args:
        file_path: Path to the .pt file
        
    Returns:
        Dictionary with 'sequence' and 'probs' keys
    """
    data = torch.load(file_path, weights_only=False)
    return {
        "sequence": data["sequence"],
        "probs": data["mean_of_probs"]
    }


def calculate_mutation_scores(probs_df: pd.DataFrame, sequence: str) -> pd.DataFrame:
    """
    Calculate mutation scores using logarithmic transformation.
    
    Formula: Score = log(P_mutant) - log(P_wild_type)
    
    Args:
        probs_df: DataFrame of amino acid probabilities
        sequence: Wild-type amino acid sequence
        
    Returns:
        DataFrame with mutation scores
    """
    score_df = probs_df.copy()
    
    # Calculate wild-type reference values
    wt_values = []
    for i, pos in enumerate(probs_df.columns):
        wt_aa = sequence[i]
        wt_val = probs_df.loc[wt_aa, pos] if wt_aa in probs_df.index else 0.001
        wt_values.append(wt_val)
    
    # Calculate mutation scores
    for col_idx, col in enumerate(score_df.columns):
        score_df[col] = np.log(score_df[col].astype(float)) - np.log(wt_values[col_idx])
        score_df[col] = round(score_df[col], 2)
    
    return score_df


def identify_best_mutations(score_df: pd.DataFrame, sequence: str, 
                           threshold: float = 0.4) -> tuple:
    """
    Identify the best mutation for each residue position.
    
    Args:
        score_df: DataFrame with mutation scores
        sequence: Wild-type amino acid sequence
        threshold: Minimum score for reporting (default: 0.4)
        
    Returns:
        Tuple of (best_mutations, best_aa_list)
    """
    best_mutations = []
    best_aa_list = []
    
    for i, col in enumerate(score_df.columns):
        max_idx = np.argmax(score_df[col].values)
        max_aa = AMINO_ACID_INDEX[max_idx]
        max_value = score_df[col].max()
        
        # Extract residue number
        match = re.search(r'\d+', col)
        res_num = match.group() if match else str(i + 1)
        
        wt_aa = sequence[i]
        mutation_label = f"{wt_aa}{res_num}{max_aa}:{max_value}"
        best_mutations.append(mutation_label)
        best_aa_list.append(max_aa)
        
        if float(max_value) >= threshold:
            print(f"  Significant: {mutation_label}")
    
    return best_mutations, best_aa_list


def main():
    parser = argparse.ArgumentParser(
        description="Convert LigandMPNN .pt output to CSV with mutation scores"
    )
    parser.add_argument("--input", type=str, required=True,
                        help="Input .pt file path")
    parser.add_argument("--output", type=str, required=True,
                        help="Output CSV file path")
    parser.add_argument("--threshold", type=float, default=0.4,
                        help="Mutation significance threshold (default: 0.4)")
    
    args = parser.parse_args()
    
    # Ensure output directory exists
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    
    print(f"Processing: {args.input}")
    
    # Step 1: Load data
    data = load_pt_file(args.input)
    sequence = data["sequence"]
    probs = data["probs"]
    print(f"  Sequence length: {len(sequence)}")
    
    # Step 2: Create probability DataFrame
    probs_df = pd.DataFrame(probs, index=AMINO_ACID_INDEX)
    
    # Step 3: Calculate mutation scores
    score_df = calculate_mutation_scores(probs_df, sequence)
    
    # Step 4: Identify best mutations
    best_mutations, best_aa_list = identify_best_mutations(score_df, sequence, args.threshold)
    
    # Step 5: Save results with metadata
    result_df = score_df.copy()
    result_df.loc['Sequence'] = list(sequence)
    result_df.loc['Max_Values'] = best_mutations
    result_df.loc['Mutations'] = best_aa_list
    result_df.to_csv(args.output)
    
    print(f"  Results saved to: {args.output}")


if __name__ == "__main__":
    main()
