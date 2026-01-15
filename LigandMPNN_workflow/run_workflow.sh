#!/bin/bash
# =============================================================================
# Lipase Mutational Sensitivity Analysis - Complete Workflow
# =============================================================================
# This script executes the complete pipeline for mutation analysis
# 
# Steps:
#   1. Run LigandMPNN scoring on all PDB structures
#   2. Process .pt files into CSV format
#   3. Calculate inter-conformational correlations
#   4. Generate final mutation rankings
# =============================================================================

set -e  # Exit on error

echo "=============================================="
echo "Lipase Mutational Sensitivity Analysis"
echo "=============================================="
echo ""

# Define paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${SCRIPT_DIR}/data"
INPUT_DIR="${DATA_DIR}/inputs"
OUTPUT_DIR="${DATA_DIR}/outputs"

# Define substrate systems
SYSTEMS=("ethanol" "isoamyl_alcohol" "isoamyl_acetate")
NUM_CLUSTERS=3

# =============================================================================
# Step 1: Run LigandMPNN Scoring
# =============================================================================
echo "[Step 1/4] Running LigandMPNN scoring..."
echo "----------------------------------------------"

# Create output directories
mkdir -p "${OUTPUT_DIR}/raw"
mkdir -p "${OUTPUT_DIR}/processed"
mkdir -p "${OUTPUT_DIR}/correlation"

# Run LigandMPNN for each system and cluster
for system in "${SYSTEMS[@]}"; do
    for cluster in $(seq 0 $((NUM_CLUSTERS - 1))); do
        pdb_file="${INPUT_DIR}/${system}/${system}_cluster${cluster}.pdb"
        out_folder="${OUTPUT_DIR}/raw/${system}/${system}_cluster${cluster}"
        
        # Skip if input file doesn't exist
        if [[ ! -f "$pdb_file" ]]; then
            echo "  Warning: ${pdb_file} not found, skipping..."
            continue
        fi
        
        echo "  Processing: ${system}_cluster${cluster}"
        
        # Run LigandMPNN (adjust path to score.py as needed)
        python "${SCRIPT_DIR}/01_run_ligandmpnn.py" \
            --pdb_path "$pdb_file" \
            --output_folder "$out_folder" \
            --seed 111
    done
done

echo ""
echo "[Step 1/4] LigandMPNN scoring completed."
echo ""

# =============================================================================
# Step 2: Process Score Files
# =============================================================================
echo "[Step 2/4] Processing score files (.pt → .csv)..."
echo "----------------------------------------------"

for system in "${SYSTEMS[@]}"; do
    for cluster in $(seq 0 $((NUM_CLUSTERS - 1))); do
        pt_file="${OUTPUT_DIR}/raw/${system}/${system}_cluster${cluster}/${system}_cluster${cluster}.pt"
        csv_file="${OUTPUT_DIR}/processed/${system}/${system}_cluster${cluster}_default.csv"
        
        # Skip if PT file doesn't exist
        if [[ ! -f "$pt_file" ]]; then
            echo "  Warning: ${pt_file} not found, skipping..."
            continue
        fi
        
        echo "  Processing: ${system}_cluster${cluster}"
        
        # Process PT to CSV
        python "${SCRIPT_DIR}/02_process_scores.py" \
            --input "$pt_file" \
            --output "$csv_file"
    done
done

echo ""
echo "[Step 2/4] Score processing completed."
echo ""

# =============================================================================
# Step 3: Calculate Correlations
# =============================================================================
echo "[Step 3/4] Calculating inter-conformational correlations..."
echo "----------------------------------------------"

python "${SCRIPT_DIR}/03_calculate_correlation.py" \
    --input_dir "${OUTPUT_DIR}/processed" \
    --output_dir "${OUTPUT_DIR}/correlation"

echo ""
echo "[Step 3/4] Correlation analysis completed."
echo ""

# =============================================================================
# Step 4: Aggregate Scores
# =============================================================================
echo "[Step 4/4] Aggregating final mutation rankings..."
echo "----------------------------------------------"

python "${SCRIPT_DIR}/04_aggregate_scores.py" \
    --correlation_dir "${OUTPUT_DIR}/correlation" \
    --processed_dir "${OUTPUT_DIR}/processed" \
    --output_file "${OUTPUT_DIR}/correlation/final_ranking.out"

echo ""
echo "[Step 4/4] Score aggregation completed."
echo ""

# =============================================================================
# Summary
# =============================================================================
echo "=============================================="
echo "Pipeline completed successfully!"
echo "=============================================="
echo ""
echo "Output files:"
echo "  - Raw scores: ${OUTPUT_DIR}/raw/"
echo "  - CSV matrices: ${OUTPUT_DIR}/processed/"
echo "  - Correlation results: ${OUTPUT_DIR}/correlation/"
echo ""
echo "Final ranking: ${OUTPUT_DIR}/correlation/final_ranking.out"
echo ""
