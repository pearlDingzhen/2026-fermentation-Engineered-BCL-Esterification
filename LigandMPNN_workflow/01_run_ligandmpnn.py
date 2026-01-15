"""
Step 1: Run LigandMPNN Scoring

This script executes LigandMPNN scoring on PDB structures.

Usage:
    python 01_run_ligandmpnn.py --pdb_path <input.pdb> --output_folder <output_dir>
    
Example:
    python 01_run_ligandmpnn.py \
        --pdb_path data/inputs/ethanol/ethanol_cluster0.pdb \
        --output_folder data/outputs/raw/ethanol/ethanol_cluster0
"""

import argparse
import subprocess
import os
from pathlib import Path


# LigandMPNN configuration - UPDATE THIS PATH to your LigandMPNN installation
LIGANDMPNN_PATH = "/mnt/hdd2/data/enzyme/ligandMPNN/LigandMPNN/score.py"
MODEL_PATH = "/mnt/hdd2/data/enzyme/ligandMPNN/LigandMPNN/model_params/ligandmpnn_v_32_010_25.pt"


def run_ligandmpnn(pdb_path: str, output_folder: str, seed: int = 111) -> None:
    """
    Execute LigandMPNN scoring on a PDB file.
    
    Args:
        pdb_path: Path to input PDB file
        output_folder: Directory for output files
        seed: Random seed for reproducibility (default: 111)
    """
    # Create output directory
    Path(output_folder).mkdir(parents=True, exist_ok=True)
    
    # Build command
    cmd = [
        "python", LIGANDMPNN_PATH,
        "--checkpoint_ligand_mpnn", MODEL_PATH,
        "--model_type", "ligand_mpnn",
        "--seed", str(seed),
        "--single_aa_score", "1",
        "--pdb_path", pdb_path,
        "--out_folder", output_folder,
        "--use_sequence", "1",
        "--batch_size", "1",
        "--number_of_batches", "10",
        "--homo_oligomer", "1"
    ]
    
    # Execute
    print(f"  Running LigandMPNN on: {os.path.basename(pdb_path)}")
    subprocess.run(cmd, check=True)
    print(f"  Output saved to: {output_folder}")


def main():
    parser = argparse.ArgumentParser(
        description="Run LigandMPNN scoring on PDB structures"
    )
    parser.add_argument("--pdb_path", type=str, required=True,
                        help="Path to input PDB file")
    parser.add_argument("--output_folder", type=str, required=True,
                        help="Directory for output files")
    parser.add_argument("--seed", type=int, default=111,
                        help="Random seed (default: 111)")
    
    args = parser.parse_args()
    
    run_ligandmpnn(args.pdb_path, args.output_folder, args.seed)


if __name__ == "__main__":
    main()
