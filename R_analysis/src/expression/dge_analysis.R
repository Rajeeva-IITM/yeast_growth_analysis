library(DESeq2)
library(data.table)
library(dplyr)
library(purrr)
library(tibble)
library(ggplot2)
library(ggpubr)
library(readr)

set.seed(42)

source("./src/utils/theme_set.R")

# Get extreme strains
# concordant_strains <- fread("./data/phenotype/Experimentally_extreme_strains.csv")
# concordant_strains <-
#   concordant_strains[, Growth := factor(Growth, levels = c("Low", "High"), ordered = FALSE)]
# strains <- concordant_strains$Strain


concordant_strains <- read_csv(
  "./data/training_data/strains_0.5_ynb.csv",
  show_col_types = F,
  num_threads = 4
)

concordant_strains <- concordant_strains %>% 
  mutate(Growth = ifelse(Phenotype==0, 'Low','High')) %>% 
  mutate(Growth = factor(Growth, levels = c("Low", "High"), ordered = FALSE))

strains <- concordant_strains$Strain

# Preparing Data ----
load("./data/expression/counts_Albert2018.RData")

## Clean up strain names ----
strain_names <- colnames(counts$pheno) %>%
  map_chr( ~ (strsplit(.x, split = "-", fixed = T) %>% pluck(1) %>% pluck(1)))
count_data <- counts$pheno

colnames(count_data) <- strain_names
storage.mode(count_data) <- "integer"

common_strains <- intersect(strain_names, strains)


concordant_strains <- concordant_strains %>% filter(Strain %in% common_strains)
count_data <-
  count_data[, common_strains]  # Select only extreme strains

# Differential expression w.r.t Growth Phenotype ----

dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = concordant_strains,
                              design = ~ Growth)
# dds$Growth <- relevel(dds$Growth, ref = 'Low')

## Prefiltering ----

smallestGroupSize <- round(ncol(dds) * 0.1) # At least 10% of samples
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep, ]

## Differential Expression Analysis ----

dds <- DESeq(dds)
res <-
  results(dds,
          contrast = c("Growth", "High", "Low"),
          alpha = 0.05)
res <- res[order(res$pvalue), ]

resLFC <- lfcShrink(dds, "Growth_High_vs_Low", type = "ashr")

resLFC <- resLFC[order(resLFC$pvalue), ]


## Table of Results ----

plot_gene_counts <- function(Gene, dds, intgroup="Growth") {
  data = plotCounts(
    dds,
    gene = Gene,
    normalized = TRUE,
    intgroup = intgroup,
    returnData = T
  )
  columns <- colnames(data) # count Factor
  # print(columns)
  # ggplot(data, mapping = aes(.data[[columns[2]]], .data[[columns[1]]])) +
  #   geom_boxplot(fill = 'cornflowerblue') +
  #   labs(title = Gene, x = columns[2], y = 'Normalized Counts')
  
  ggboxplot(
    data,
    columns[2], columns[1],
    title = Gene, xlab =  columns[2], ylab = 'Normalized Count'
  )
}

result_df <- resLFC %>%
  as_tibble(rownames='Gene')

normalized_counts <- counts(dds, normalized = T) %>%
  as_tibble(rownames = 'Gene')


## Saving results ----

DESeq_results <- list(
  result_df = result_df,
  normalized_counts = normalized_counts,
  plot_gene_counts = plot_gene_counts
)

saveRDS(DESeq_results, file = "./results/result_dge/0.5/DESeq_resutls_0.5_sigma.RDS") # Save variable

## Plotting for all genes ----

# Get only significantly differing genes
genes <- result_df %>% 
  filter(padj < 0.05) %>% 
  pluck('Gene')

for(gene in genes){
  
  savename <- paste0("results/result_dge/0.5_redo/gene_count_vs_growth/", gene, ".svg")
  
  if(file.exists(savename)) { # If already present, skip
    next()
  }
  
  p <- plot_gene_counts(gene, dds)
  
  save_plot(
    savename,
    height = 5,
    width = 5,
    units = 'in'
  )
}


# Differential expression w.r.t mutational status ----

geno_df <- read_csv(
  "./data/training_data/bloom2013_clf_sigma0.5_ynb.csv",
  show_col_types = F,
  num_threads = 4
) %>% 
  filter(Strain %in% common_strains) %>%  # only concordant strains
  filter(Condition == 'ynb_glucose')

geno_df <- geno_df %>%
  select(-Strain) %>% 
  mutate(across(starts_with('Y'), ~if_else(.x==0, 0, 1)))

genes <- colnames(geno_df)

res_list <- list()

for(gene in genes){
  
  # Filter out non-varying genes
  if(dim(table(geno_df[[gene]])) != 2){
    next
  }
  # Skip if already done
  savename <- paste0("./results/result_dge/0.5_redo/genewise_anlaysis/", gene, ".svg")
  if(file.exists(savename)){
    next 
  }
  
  col_data <- tibble(
    Strain = common_strains,
  )
  col_data[[gene]] <- ifelse(geno_df[[gene]]==0, "Unmutated", "Mutated") %>% 
    as.factor()
  
  formula_gene <- paste("~", gene) %>%
    as.formula()
  
  dds_gene <- DESeqDataSetFromMatrix(
    countData = count_data,
    colData = col_data,
    design = formula_gene
  )
  
  ## Pre-filtering 
  
  smallestGroupSize <- round(ncol(dds_gene) * 0.1) # At least 10% of samples
  keep <- rowSums(counts(dds_gene) >= 10) >= smallestGroupSize
  dds_gene <- dds_gene[keep, ]
  
  dds_gene <- DESeq(dds_gene)
  coef_name <- paste(gene, 'Unmutated', 'vs', 'Mutated', sep = '_')
  res_gene <- lfcShrink(dds_gene, coef = coef_name, type = "ashr")
  res_gene <- res_gene[order(res_gene$pvalue), ]
  
  # dds_list[[gene]] <- dds_gene
  res_list[[gene]] <- res_gene
  
  
  
  p <- ggmaplot(
    res_gene,
    main = paste0(gene, " : ","Mutated vs Unmutated"),
    genenames = rownames(res_gene),
    size = 0.4,
    top = 10,
    fc = 2,
    font.label = c("bold", 11),
    label.rectangle = TRUE
  )
  
  save_plot(
    savename,
    plot = p,
    height = 5,
    width = 9,
    units = "in"
  )
}

saveRDS(res_list, file = "./results/result_dge/0.5_redo/DESeq_genenwise_results.RDS")
