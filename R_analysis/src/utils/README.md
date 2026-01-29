# Utilities

Shared functions used across all analysis modules.

## Files

### `theme_set.R`
Sets a consistent ggplot2 theme for all figures (12pt font, centered titles, facet styling). Defines `save_plot()` (a `ggsave` wrapper with 300 DPI and transparent background) and `cbPalette` (colorblind-friendly palette).

### `enrichment_analysis.R`
Wrapper functions for gene set enrichment:
- `enrichGO_yeast()` -- standardized GO enrichment using clusterProfiler and org.Sc.sgd.db
- `parse_enrich_results()` -- computes fold change for enrichment results
- `custom_treeplot()`, `go_heatplot()` -- publication-ready GO visualizations
- `yeast_gost()`, `plot_gost_result()` -- gProfiler2 wrappers for KEGG and TF enrichment
- `get_gene_annotation()` -- extracts gene lists per GO term

Loads `all_genes` (gene universe) and `scGO` (semantic similarity data) from `data/miscellaneous/`.

### `flux_enrichment_analysis.R`
Defines `FEA()`, a hypergeometric enrichment function for metabolic reaction subsystems. Tests whether sampled reactions are over-represented in specific pathway groups. Returns fold-change and FDR-adjusted p-values.
