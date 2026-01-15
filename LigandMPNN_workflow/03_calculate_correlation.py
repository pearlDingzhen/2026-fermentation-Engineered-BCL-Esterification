"""
Step 3: Calculate Inter-Conformational Correlations

This script calculates Pearson correlation coefficients between mutation
scores from different conformational states.

Usage:
    python 03_calculate_correlation.py --input_dir <dir> --output_dir <dir>
    
Example:
    python 03_calculate_correlation.py \
        --input_dir data/outputs/processed \
        --output_dir data/outputs/correlation
"""

import argparse
import pandas as pd
import numpy as np
import os
from pathlib import Path
from typing import List, Dict


# Protein sequence (284 residues)
SEQUENCE = list("HPVFVLVHGAWHGAWCYAHVAAALAERGYLSIARDLPAHGINARFPASYLERPLDKDAFGAEPSPVANTTLDDYATQVMEAVDDAYALGHGKVVLVGHSMGGLAITAAAERAPEKIAKIVYLAAFMPASGVPGLDYVRAPENKGEMLAPLMLASPRVAGALRIDPRSGDAAYRALAKRALYDDAAQADFEAMANLMTCDVPAAPFATAIPTTAARWGAIDRHYIKCLADRVILPALQQRFIDEADAFVPGNPTHVHQLDSSHSPFVSQPGVLAGVLVDIAKSIA")


def get_csv_files(input_dir: str) -> Dict[str, List[str]]:
    """
    Organize CSV files by system.
    
    Returns:
        Dictionary with structure:
        {
            'ethanol': ['ethanol_cluster0_default.csv', ...],
            'isoamyl_alcohol': [...],
            'isoamyl_acetate': [...]
        }
    """
    systems = ['ethanol', 'isoamyl_alcohol', 'isoamyl_acetate']
    csv_files = {system: [] for system in systems}
    
    for system in systems:
        system_dir = Path(input_dir) / system
        if system_dir.exists():
            for f in system_dir.glob("*_default.csv"):
                csv_files[system].append(f.name)
    
    return csv_files


def calculate_correlation(df1: pd.DataFrame, df2: pd.DataFrame) -> List[float]:
    """
    Calculate Pearson correlation for each residue position.
    
    Args:
        df1: First score matrix
        df2: Second score matrix
        
    Returns:
        List of correlation coefficients for each position
    """
    correlations = []
    for i in range(len(SEQUENCE)):
        corr = df1.iloc[:, i].corr(df2.iloc[:, i])
        correlations.append(corr)
    return correlations


def calculate_and_save_correlations(csv_files: Dict[str, List[str]], 
                                    input_dir: str, 
                                    output_dir: str) -> None:
    """
    Calculate correlations between all system pairs and save results.
    """
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Define comparison pairs (ethanol vs isoamyl_alcohol)
    comparisons = []
    for ethanol_cluster in range(3):
        for ia_cluster in range(3):
            comparisons.append((
                'ethanol', f'ethanol_cluster{ethanol_cluster}_default.csv',
                'isoamyl_alcohol', f'isoamyl_alcohol_cluster{ia_cluster}_default.csv',
                f'ethanol_vs_isoamyl_alcohol_{ethanol_cluster}_{ia_cluster}'
            ))
    
    # Calculate and save each comparison
    for sys1, file1, sys2, file2, name in comparisons:
        path1 = Path(input_dir) / sys1 / file1
        path2 = Path(input_dir) / sys2 / file2
        
        if not path1.exists() or not path2.exists():
            print(f"  Warning: Skipping {name} (missing files)")
            continue
        
        print(f"  Calculating: {name}")
        
        # Load data (skip metadata rows)
        df1 = pd.read_csv(path1, header=None, skiprows=1, nrows=20, 
                         usecols=range(1, 285))
        df2 = pd.read_csv(path2, header=None, skiprows=1, nrows=20,
                         usecols=range(1, 285))
        
        # Calculate correlations
        corr_list = calculate_correlation(df1, df2)
        mean_corr = np.mean(corr_list)
        print(f"    Mean correlation: {mean_corr:.3f}")
        
        # Save correlation values
        txt_file = Path(output_dir) / f"corr_{name}.txt"
        with open(txt_file, 'w') as f:
            f.write(f"# {name}\n")
            f.write(f"# Mean correlation: {mean_corr:.4f}\n")
            for i, corr in enumerate(corr_list):
                f.write(f"{i+1} {corr:.4f}\n")
        
        # Save CSV with metadata
        csv_file = Path(output_dir) / f"output_{name}.csv"
        with open(csv_file, 'w') as f:
            f.write(f"Name,{name}\n")
            f.write(f"Average Correlation,{mean_corr:.2f}\n")
            f.write("Wild AA,Residue Number,Correlation,MPNN Score,MPNN AA\n")
            for i, corr in enumerate(corr_list):
                wt_aa = SEQUENCE[i]
                df_subset = pd.read_csv(path1, header=None, skiprows=1, nrows=20,
                                       usecols=range(1, 285))
                mpnn_score = df_subset.iloc[:, i].max()
                max_aa_idx = df_subset.iloc[:, i].idxmax()
                f.write(f"{wt_aa},{i+1},{corr:.4f},{mpnn_score:.2f},{max_aa_idx}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Calculate inter-conformational correlations"
    )
    parser.add_argument("--input_dir", type=str, required=True,
                        help="Directory containing processed CSV files")
    parser.add_argument("--output_dir", type=str, required=True,
                        help="Directory for correlation output files")
    
    args = parser.parse_args()
    
    print("Step 3: Calculating inter-conformational correlations")
    print("-" * 50)
    
    # Get organized CSV files
    csv_files = get_csv_files(args.input_dir)
    
    # Calculate and save correlations
    calculate_and_save_correlations(csv_files, args.input_dir, args.output_dir)
    
    print(f"\nCorrelation results saved to: {args.output_dir}")


if __name__ == "__main__":
    main()
