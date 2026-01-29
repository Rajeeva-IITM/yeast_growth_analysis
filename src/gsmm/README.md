# Genome-Scale Metabolic Model Analysis

Flux enrichment analysis (FEA) on metabolic reaction data from constraint-based modeling.

## Files

### `result_gsmm.R`
Performs hypergeometric enrichment to identify metabolic subsystems (pathway groups) over-represented among reactions that differ between high and low growth strains. Analyzes unique reactions, fold-change-sorted reactions, and active reactions. Outputs enrichment barplots with FDR-adjusted p-values per subsystem.

## Dependencies

- `src/utils/flux_enrichment_analysis.R`, `src/utils/theme_set.R`
- Upstream: flux sampling results from the MATLAB pipeline (`data/gsmm/`)
