# Gene Regulatory Network Analysis

Construction and analysis of gene regulatory networks (GRN), integrating SHAP feature importance, metabolic data, and differential expression.

## Files

### `grn_building.R`
Constructs the GRN using the BioNERO package. Filters expression data (non-expressed genes, top 90% variance), applies PC correction, and infers regulator-target edges using a curated transcription factor list. Must be run before `result_grn.R`.

### `result_grn.R`
Comprehensive downstream analysis of the GRN:
- GO enrichment of regulator and target gene sets (tree/heatmap plots)
- Overlay of SHAP values, contingency results, and pFBA-optimal genes onto the network
- Network visualizations with node sizes scaled by SHAP importance
- Per-regulator subgraph analysis with GO annotations
- Focused analysis of novel regulators (MDS3, PDR8) including differential expression enrichment

## Run Order

`grn_building.R` must run first to produce the GRN object (`data/grn/GRN_results_gimme_0.5_sigma_reimand.sav`), which `result_grn.R` then loads.

## Dependencies

- `src/utils/theme_set.R`, `src/utils/enrichment_analysis.R`
- Upstream: SHAP values, contingency tables, pFBA results, DGE results
