# Lipase Mutational Sensitivity Analysis Pipeline

## Overview

This pipeline uses **LigandMPNN** deep learning model to analyze mutational sensitivity of lipase enzymes under different substrate conditions, identifying key mutation sites to enhance enzymatic activity.

## Substrate Systems

| Abbreviation | Full Name | Reaction Type |
|--------------|-----------|---------------|
| **ethanol** | Ethanol (CH₃CH₂OH) | Esterification |
| **isoamyl_alcohol** | Isoamyl Alcohol ((CH₃)₂CHCH₂CH₂OH) | Esterification |
| **isoamyl_acetate** | Isoamyl Acetate (CH₃COOCH₂CH₂CH(CH₃)₂) | Hydrolysis |

## Pipeline Workflow

```
1. 01_run_ligandmpnn.py      # Run LigandMPNN on PDB structures
   ↓
2. 02_process_scores.py      # Convert .pt to .csv with mutation scores
   ↓
3. 03_calculate_correlation.py  # Calculate inter-conformational correlations
   ↓
4. 04_aggregate_scores.py    # Generate final mutation ranking
```

## Quick Start

```bash
# Execute complete pipeline
bash run_workflow.sh

# Or step by step:
python 01_run_ligandmpnn.py
python 02_process_scores.py
python 03_calculate_correlation.py
python 04_aggregate_scores.py
```

## Input/Output Structure

### Input Files
```
data/inputs/
├── ethanol/
│   ├── ethanol_cluster0.pdb
│   ├── ethanol_cluster1.pdb
│   └── ethanol_cluster2.pdb
├── isoamyl_alcohol/
│   ├── isoamyl_alcohol_cluster0.pdb
│   ├── isoamyl_alcohol_cluster1.pdb
│   └── isoamyl_alcohol_cluster2.pdb
└── isoamyl_acetate/
    ├── isoamyl_acetate_cluster0.pdb
    ├── isoamyl_acetate_cluster1.pdb
    └── isoamyl_acetate_cluster2.pdb
```

### Output Files
```
data/outputs/
├── raw/                  # LigandMPNN .pt output files
├── processed/            # CSV score matrices
└── correlation/          # Correlation analysis results
    ├── correlation_matrix.csv
    └── final_ranking.out
```

## Key Methods

### Mutation Score Calculation

For each residue position \( i \), mutation to amino acid \( aa \):

\[
\text{Score}_{i \to aa} = \log(P_{aa}) - \log(P_{wt})
\]

### Aggregate Score Formula

\[
\text{Final Score} = (\text{CorrScore} \times 5 + \text{ValueScore}) \times W
\]

Where:
- **CorrScore**: Correlation score based on cross-conformational consistency
- **ValueScore**: Normalized LigandMPNN mutation score
- **W**: Weight based on cluster population

## Output Interpretation

The final ranking file `final_ranking.out` lists mutations in descending order:

```
Rank 1: A23F, Score=45.32
  - ethanol_cluster0 vs isoamyl_alcohol_cluster0: A 23 0.65 0.8 F
  - ethanol_cluster0 vs isoamyl_alcohol_cluster1: A 23 0.72 0.85 F
```

**Mutation format:** `{WildType}{Position}{Mutant}` (e.g., A23F)

## Citation

If using this pipeline, please cite:
- LigandMPNN: Protein sequence design with differentiable models
