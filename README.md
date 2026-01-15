# Engineered BCL for Esterification Reactions

This repository contains computational and experimental data for the engineering of BCL (butyrolactone biosynthesis-related enzyme) for improved esterification activity, supporting the 2026 fermentation research project.

## Project Overview

This dataset accompanies research on enzyme engineering of Bacillus subtilis lipase A (BSLA) variants with enhanced esterification capabilities for biotechnological applications. The study combines protein structure prediction, molecular dynamics simulations, and LigandMPNN workflow for enzyme optimization.

## Data Description

The repository includes the following directories and their contents:

### 1. LigandMPNN_workflow
Contains inputs and outputs from the LigandMPNN protein design pipeline used for enzyme optimization. This includes:
- Ligand structure files used for binding site analysis
- Protein sequence designs and mutations
- Scoring and ranking results

### 2. Lipase_Esterification_Ethanol
Data for esterification experiments with ethanol as the alcohol substrate:
- Experimental design and parameters
- GC-MS analysis results
- Conversion yield data
- Time-course measurements

### 3. Lipase_Esterification_IsoamylAlcohol
Data for esterification experiments with isoamyl alcohol:
- Experimental protocols
- Product analysis data
- Activity assays
- Statistical replicates

### 4. Lipase_Hydrolysis_IsoamylAcetate
Control hydrolysis experiments using isoamyl acetate substrate:
- Hydrolysis rate measurements
- Substrate specificity data
- Kinetic parameters

### 5. Molecule_Dynamic_Parameters
Molecular dynamics simulation parameters and results:
- Force field parameters
- Simulation input files
- Trajectory analysis results
- Binding free energy calculations

### 6. Structure_prediction
Protein structure prediction outputs:
- AlphaFold/ColabFold prediction results
- Structural models and confidence scores
- Active site analysis
- Mutant structure predictions

## Usage Instructions

### Reproducibility

1. Clone this repository:
```bash
git clone https://github.com/[username]/[repository-name].git
cd [repository-name]
```

2. Required software and versions:
- Python 3.8+
- GROMACS 2023+ (for molecular dynamics)
- AlphaFold/ColabFold
- LigandMPNN
- GC-MS analysis software

3. Directory structure preservation:
   - Do not modify original data filenames
   - Use relative paths for reproducibility

### Data Analysis

Each subdirectory contains a `README.md` with specific analysis protocols used for that dataset.

## Citation

If you use this data in your research, please cite:

**APA Format:**
[Author Name(s)]. (2026). Engineered BCL for Esterification Reactions [Dataset]. GitHub. https://github.com/[username]/[repository-name]

**BibTeX:**
```bibtex
@misc{author2026engineered,
  author = {[Author Names]},
  title = {Engineered BCL for Esterification Reactions},
  year = {2026},
  publisher = {GitHub},
  url = {https://github.com/[username]/[repository-name]}
}
```

## License

This work is licensed under a [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/).

You are free to:
- Share: Copy and redistribute the material in any medium or format
- Adapt: Remix, transform, and build upon the material for any purpose, even commercially

Under the following terms:
- Attribution: You must give appropriate credit to the original authors

## Contact Information

**Corresponding Author:**
- Name: [Your Name]
- Email: [your.email@institution.edu]
- Institution: [Your Institution]

**For Questions:**
- Open an issue on GitHub for technical questions
- Contact the corresponding author for research inquiries

## Acknowledgments

This research was supported by:
- [Funding agency or grant number]
- [Institution name]
- [Any other acknowledgments]

## Related Publications

- [Paper Title 1] - [Journal, Year]
- [Paper Title 2] - [Journal, Year]

## Last Updated

January 2026

## Version

v1.0.0 - Initial release
