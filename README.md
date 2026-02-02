# Yeast Growth Prediction

Predicting and analysing yeast growth phenotypes from genotype data across multiple Bloom lab crosses. This repository brings together three analysis pipelines that cover machine learning, metabolic modelling, and statistical/biological interpretation.

## Repository Structure

| Directory | Description |
|---|---|
| [`Boosting/`](Boosting/) | ML model training and evaluation (LightGBM, SVM, RF, Elastic Net, etc.) with Optuna hyperparameter tuning and SHAP analysis |
| [`MATLAB/`](MATLAB/) | Genome-scale metabolic modelling — strain-specific GIMME models, pFBA, and flux sampling |
| [`R_analysis/`](R_analysis/) | Statistical analysis, visualisation, and biological interpretation of results from the other pipelines |

Each directory is self-contained with its own README, dependencies, and configuration. All scripts within a directory assume that directory is the working directory.

## Getting Started

Clone the repo and navigate into the pipeline you want to run:

```bash
git clone https://github.com/Rajeeva-IITM/yeast_growth_analysis.git
cd yeast_growth_analysis

# Example: run the ML pipeline
cd Boosting
# follow Boosting/README.md

# Example: run metabolic modelling
cd MATLAB
# follow MATLAB/README.md

# Example: run R analysis
cd R_analysis
# follow R_analysis/README.md
```

Refer to the README inside each subdirectory for environment setup, dependencies, and usage instructions. Each pipeline must be run from its own subdirectory (`Boosting/`, `MATLAB/`, `R_analysis/`).

> **Note:** Some scripts contain hardcoded values (e.g., number of parallel cores, LP solver choice, file indices) that were configured for the original development environment. Please review the configuration notes in each subdirectory's README before running.

## Data

<!-- TODO: Add Zenodo DOI link once published -->
All data required to reproduce the analyses will be available on [Zenodo](https://zenodo.org/). Download and extract the data archive into the repository root.

## Data Flow

```
Boosting (ML models, SHAP values)
    |
    v
R_analysis (statistical analysis, visualisation)
    ^
    |
MATLAB (metabolic models, flux sampling)
```

The R analysis pipeline consumes outputs from both the Boosting and MATLAB pipelines. Run those first if you need to regenerate results from scratch.
