# Data Exploration

Initial data quality assessment, phenotype distribution visualization, and exploratory GO enrichment by mutation class.

## Files

### `result_fig_gen.R`
Generates dataset overview figures: growth rate histograms per condition, correlation heatmaps (using ComplexHeatmap), and a UMAP embedding of phenotype space. Covers all five Bloom datasets.

### `preliminary_go_analysis.R`
GO enrichment analysis of genes grouped by mutation class (unmutated, missense, nonsense, frameshift). Tests whether different mutation types affect distinct biological processes. Also generates tile plots comparing continuous vs binarized growth phenotypes.

## Dependencies

- `src/utils/theme_set.R`, `src/utils/enrichment_analysis.R`
