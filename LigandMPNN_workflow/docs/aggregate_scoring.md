# Aggregate Scoring Methodology

## Overview

This document describes the aggregate scoring system that combines LigandMPNN mutation scores with cross-conformational correlation analysis to rank protein mutations.

## Theoretical Framework

### Motivation

Individual mutation scores from different conformational states may vary significantly. The aggregate scoring system aims to:

1. **Identify consistent mutations** that perform well across multiple conformational states
2. **Account for conformational diversity** in the protein's functional cycle
3. **Provide a unified ranking** for mutation prioritization

### Mathematical Formulation

The aggregate score \( S \) for a mutation \( m \) is defined as:

\[
S(m) = \sum_{i,j} w_{ij} \cdot s_{ij}(m)
\]

Where:
- \( s_{ij}(m) \) is the individual score for mutation \( m \) comparing conformations \( i \) and \( j \)
- \( w_{ij} \) is the weight based on cluster population

#### Individual Score Components

For each conformation pair comparison, the individual score is:

\[
s_{ij}(m) = (\text{CorrScore}_{ij} \times 5 + \text{ValueScore}_{ij}) \times W_{ij}
\]

##### Correlation Score

\[
\text{CorrScore}_{ij} = \frac{\max_{corr} - \rho_{ij}}{(\max_{corr} - \min_{corr})^{0.7}} \times 100
\]

Where:
- \( \rho_{ij} \) = Pearson correlation coefficient between conformation \( i \) and \( j \)
- \( \max_{corr} = 0.6 \) (upper bound for correlation)
- \( \min_{corr} = -0.7 \) (lower bound for correlation)

**Interpretation:** Lower correlation (more negative) yields higher correlation score, indicating mutations that behave differently across conformations.

##### Value Score

\[
\text{ValueScore}_{ij} = \frac{V_{ij}}{\max_{value}} \times 100
\]

Where:
- \( V_{ij} \) = LigandMPNN mutation score for the specific amino acid change
- \( \max_{value} = 5.0 \) (normalization constant)

##### Weight Factor

\[
W_{ij} = \frac{n_i \times n_j}{n_{total}} / 6
\]

Where:
- \( n_i, n_j \) = Cluster populations for conformations \( i \) and \( j \)
- \( n_{total} \) = Total combined population of all compared conformations

## Algorithm Implementation

### Step 1: Data Collection

```python
# Load all mutation score files
score_files = [
    'output_ywc0_yc0.csv',
    'output_ywc0_yc1.csv',
    # ... all 18 comparison pairs
]

for file in score_files:
    df = pd.read_csv(file)
    # Extract mutation scores and correlations
```

### Step 2: Score Aggregation

```python
# Dictionary to accumulate scores for each unique mutation
final_score_list = {}

for comparison in comparisons:
    # Process each comparison file
    results = process_csv(comparison_file, system1, system2)
    
    for result in results:
        # Unique mutation key: wild_aa + position + mutant_aa
        key = result[0] + str(result[1]) + result[4]
        
        if key not in final_score_list:
            final_score_list[key] = [result[7], [information]]
        else:
            final_score_list[key][0] += result[7]
            final_score_list[key][1].append(information)
```

### Step 3: Ranking

```python
# Sort by total score in descending order
sorted_items = sorted(final_score_list.items(), 
                      key=lambda x: x[1][0], 
                      reverse=True)

# Output ranked mutations
for rank, (mutation, (total_score, details)) in enumerate(sorted_items, 1):
    print(f"Rank {rank}: {mutation}, Score={total_score:.2f}")
```

## Output Interpretation

### Score Components

| Component | Range | Interpretation |
|-----------|-------|----------------|
| Correlation Score | 0-100 | Higher = more variable across conformations |
| Value Score | 0-100 | Higher = more beneficial mutation |
| Weight | Variable | Reflects conformational sampling quality |
| Final Score | Variable | Weighted combination of above |

### Ranking Criteria

Mutations are ranked by the aggregate score, which considers:

1. **Magnitude** - Higher raw mutation scores are preferred
2. **Consistency** - Mutations that score well across multiple comparisons are reinforced
3. **Conformational Coverage** - Mutations supported by well-sampled conformations are weighted more heavily

## Limitations and Considerations

### 1. Correlation Threshold Selection

The correlation score formula uses fixed thresholds (0.6 and -0.7). These values were empirically determined and may need adjustment for different protein systems.

### 2. Weight Normalization

The division by 6 in the weight formula is a normalization constant derived from the 18 comparison pairs in this study. Other systems may require different normalization.

### 3. Log Transformation

Mutation scores use logarithmic transformation to ensure:
- Additivity of probabilities
- Symmetry between stabilizing/destabilizing mutations
- Numerical stability across probability ranges

## Validation

The aggregate scoring system was validated by:
1. Comparing with known functional mutations
2. Analyzing conservation patterns at high-scoring positions
3. Cross-validation with separate MD trajectory sets

## References

1. Pearson, K. (1909). Mathematical contributions to the theory of evolution. - XV. Regression, heredity and panmixia. *Philosophical Transactions of the Royal Society A*.
2. LigandMPNN: Learning protein sequence design with differentiable models. *GitHub Repository*.
