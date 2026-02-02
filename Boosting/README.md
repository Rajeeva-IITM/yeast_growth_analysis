# Boosting - Yeast Growth Prediction

Hyperparameter tuning, training, and evaluation of machine learning models for predicting yeast growth phenotypes from genotype and chemical environment features. Uses data from Bloom et al. yeast crosses (Bloom2013, Bloom2015, Bloom2019 BYxRM).

Models include LightGBM (gradient boosting), Logistic Regression, SVM/SVR, Random Forest, Elastic Net, LDA/QDA, KNN, and Naive Bayes. Hyperparameter optimization is done via Optuna with TPE sampling and cross-validation.

## Setup

### Environment

Create the conda environment from the provided specification:

```bash
conda env create -f environment.yaml
conda activate boost_env
```

### Configuration

All scripts use [Hydra](https://hydra.cc/) for configuration and read paths from a `.env` file in the project root. Edit `.env` to match your machine:

```
ROOT_DIR="/path/to/your/project/root/"
PROJECT_DIR="/path/to/Boosting-yeast_growth_pred"
RUN_DIR="/path/to/runs"
DATA_DIR="/path/to/data"
NUM_THREADS=12
```

Hydra config files are in `configs/`. The main config is `configs/conf.yaml`, which references model parameters (`model_params.yaml`), metrics (`metrics/clf.yaml` or `metrics/reg.yaml`), and utility settings (`utilities.yaml`).

## Data Format

Input data should be in Feather (`.feather`), Parquet (`.parquet`), or CSV (`.csv`) format with the following columns:

- `Strain` — strain identifier
- `Condition` — experimental condition
- `Phenotype` — target variable (binary for classification, continuous for regression)
- Columns starting with `Y` — genotype features (yeast systematic gene names)
- Columns containing `latent` — chemical environment features

## Running

All scripts are run from the project root directory.

### LightGBM Tuning and Training

```bash
python src/tune_model.py --multirun 'seed=1,2,3,4,5' \
  'data.path=${oc.env:DATA_DIR}/your_data.feather' \
  data.savename=Your_Run_Name
```

See `scripts/run_script.sh` for regression examples and `scripts/regression.sh` for additional configurations.

### Simple Models (LogReg, SVM, RF, Elastic Net, SVR)

```bash
python src/tune_and_train_simple_models.py --multirun 'seed=1,2,3,4,5' \
  'data.path=${oc.env:DATA_DIR}/your_data.feather' \
  data.savename=Your_Run_Name
```

### Additional Classifiers (LDA, QDA, KNN, Naive Bayes)

```bash
python src/tune_and_train_simple_models_2.py --multirun 'seed=1,2,3,4,5' \
  'data.path=${oc.env:DATA_DIR}/your_data.feather' \
  data.savename=Your_Run_Name
```

See `scripts/run_simple.sh` for batch examples.

### Model Evaluation

```bash
python src/compare_performance.py \
  'run_path=${oc.env:RUN_DIR}/your_runs/' \
  'out_path=${oc.env:ROOT_DIR}/performance/output' \
  'data_path=${oc.env:DATA_DIR}/your_data_dir'
```

See `scripts/evaluate.sh` for a full example.

### Data Sampling Analysis

```bash
python src/data_sampling_analysis.py \
  data.train_data='${oc.env:DATA_DIR}/your_train.feather' \
  data.test_data='${oc.env:DATA_DIR}/your_test.feather' \
  data.savename=Your_Run_Name
```

See `scripts/run_tuning.sh` for batch examples across datasets.

### SHAP Analysis

```bash
python src/shap_analysis.py
```

Configured via `configs/interpret.yaml`.

## Project Structure

```
src/
  tune_model.py                    — LightGBM hyperparameter tuning + training
  tune_and_train_simple_models.py  — LogReg, SVM, RF, Elastic Net, SVR tuning
  tune_and_train_simple_models_2.py — LDA, QDA, KNN, NB tuning
  train.py                         — Standalone LightGBM training from saved params
  compare_performance.py           — Model evaluation and comparison
  data_sampling_analysis.py        — Data sampling robustness analysis
  shap_analysis.py                 — SHAP feature importance analysis
  utils.py                         — Shared utilities (data loading, model loading, metrics)
  utils/
    data.py                        — Data statistics utility

configs/                           — Hydra YAML configuration files
scripts/                           — Shell scripts for batch runs
notebooks/                         — Jupyter notebooks for visualization and evaluation
```
