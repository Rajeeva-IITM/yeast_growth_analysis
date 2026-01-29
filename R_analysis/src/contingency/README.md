# Contingency Analysis

Statistical association between genotype mutations and growth phenotypes using Fisher's exact test, with integration of SHAP values and QTL data.

## Files

### `contingency_table_analysis.R`
Performs Fisher's exact test on 2x2 contingency tables (mutated/unmutated vs high/low growth) for each gene across conditions. Outputs odds ratios with FDR-corrected p-values as CSVs and condition-wise volcano plots.

### `prediction_contingency.R`
Applies contingency analysis to ML model predictions rather than observed phenotypes. Uses parallelization across cross-validation folds and combines p-values via Fisher's method. Outputs parquet files for cross-dataset prediction contingency.

### `shap_contingency_analysis.R`
Directly compares SHAP feature importance against Fisher's exact test results. Produces scatter plots of SHAP vs odds ratio, colored by significance and annotated with QTL regions.

### `shap_qtl_analysis.R`
Examines overlap between SHAP-identified genes, known QTLs, and contingency-significant genes. Generates Venn diagrams across traits and GO enrichment (tree plots) for contingency genes grouped by stress type (carbon, oxidative, genotoxic).

### `qtl_shap_ranking_analysis.R`
Ranks genes per condition by combining contingency odds ratios and SHAP values. Overlays SHAP magnitude on genomic position plots with QTL annotations. Includes locus-specific visualizations for known QTLs (4NQO, sorbitol, congo red, neomycin).

### `bloom2019_shap_qtl_contingency.R`
Bloom2019-specific analysis integrating SHAP values, previously detected causal genes, and contingency results. Outputs Venn diagrams comparing causal vs contingency-identified vs SHAP-identified gene sets.

## Dependencies

- `src/utils/theme_set.R`, `src/utils/enrichment_analysis.R`
- Upstream: SHAP values (from Python pipeline), QTL data, contingency tables
