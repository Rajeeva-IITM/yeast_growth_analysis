# SHAP Analysis and Latent Interpretation

Processing of ML feature importance (SHAP values) and interpretation of latent variable representations.

## Files

### `result_shap.R`
Processes SHAP values from multiple ML models (Boosting, LogReg, RF, SVM) across datasets. Generates top-feature barplots and SHAP distribution histograms. Identifies consensus important features (genes appearing in the top 1% across models and datasets) and saves common/unique feature sets as RDS.

### `result_shap_enrichment.R`
Functional enrichment of SHAP-identified genes using GO (clusterProfiler) and pathway databases (KEGG/TF via gProfiler2). Runs enrichment at multiple levels: genes common across all models, full per-model gene sets, and model-unique genes. Outputs tree plots and heatmaps.

### `result_latent_interpretability.R`
Interprets latent representations from autoencoder models via gradient-based attribution methods (Integrated Gradients, Noise Tunnel, Saliency, GradientSHAP, InputxGradient). For each of 256 latent dimensions, identifies top contributing features and correlates them with chemical functional groups. Generates per-latent variable attribution plots and feature-count summaries.

## Dependencies

- `src/utils/theme_set.R`, `src/utils/enrichment_analysis.R`
- Upstream: SHAP CSVs and latent attribution data from the Python ML pipeline
