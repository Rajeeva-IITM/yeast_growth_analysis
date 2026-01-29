# Expression Analysis

Differential gene expression (DGE) analysis using DESeq2 on RNA-seq data from extreme-growth strains.

## Files

### `dge_analysis.R`
Two-pronged DGE analysis:
1. **Growth phenotype**: High vs Low growth strains (DESeq2 with ashr LFC shrinkage)
2. **Per-gene mutation status**: Mutated vs Unmutated for each gene individually

Inputs RNA-seq counts (Albert2018) and strain classifications. Outputs DESeq2 results (RDS), normalized count plots for significant genes, and per-gene fold-change vs significance plots.

## Dependencies

- `src/utils/theme_set.R`
- `data/expression/counts_Albert2018.RData`
