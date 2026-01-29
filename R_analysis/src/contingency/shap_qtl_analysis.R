# SHAP and previously identified QTLs

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(purrr)
# library(VennDiagram)
library(ggvenn)


source("./src/utils/theme_set.R")
source("./src/utils/enrichment_analysis.R")

shap_df <- read_csv(
  "./data/shap/sigmas/shap_classification_0.5/shap_Bloom2013_Boosting.csv"
) %>%
  filter(Value > 1e-5)

shap_genes <- keep(shap_df$Feature, ~grepl("^Y", .x))

qtl_df <- read_csv(
  "./data/qtl/detected_qtl_bloom2013.csv"
)

all_qtl_df <- qtl_df %>% 
  select(Trait, Pheno_fraction_explained, Peak_pos, Chromosome, Genes) %>% 
  mutate(Genes = strsplit(as.character(Genes), split = "|", fixed=T)) %>% 
  unnest_longer(Genes) 

all_qtl_gene_df <- all_qtl_df %>% 
  group_by(Trait) %>% 
  reframe(all_genes = list(Genes))

for(i in 1:dim(all_qtl_gene_df)[1]){
  
  trait <- all_qtl_gene_df[i,1]$Trait
  x <- list()
  x[[trait]] <- all_qtl_gene_df[i,2]$all_genes %>% unlist
  x[['SHAP']] <- shap_genes
  
  p <- ggvenn(x)
  save_plot(paste0("./results/result_contingency/venn/", trait, ".png"), 
            plot=p)
  # venn.diagram(
  #   x =x,
  #   filename = paste0("./results/result_contingency/venn/", trait, ".png"),
  #   main = trait,
  #   imagetype = 'png',
  #   disable.logging = T,
  #   
  #   # Circle customization
  #   lwd = 2,
  #   lty = 'blank',
  #   fill = myCol,
  #   
  #   # Numbers
  #   cex = .6,
  #   fontface = "bold",
  #   fontfamily = "sans",
  #   
  #   # # Set names
  #   # cat.cex = 0.6,
  #   # cat.fontface = "bold",
  #   # cat.default.pos = "outer",
  #   # cat.fontfamily = "sans",
  #   # rotation = 1
  # )
}

gene_counts <- all_qtl_df %>% group_by(Trait) %>% summarise(Genes = length(Genes))


shap_qtl_df <- all_qtl_df %>% 
  filter(Genes %in% shap_genes)

shap_gene_counts <- shap_qtl_df %>% group_by(Trait) %>% summarise(Genes = length(Genes))

# Filtering ----


qtl_gene_counts_df <- qtl_df %>% 
   select(Trait, Pheno_fraction_explained, Peak_pos, Chromosome,
          Genes, `1.5 LOD Confidence Interval (Left) (bp)`, `1.5 LOD Confidence Interval (Right) (bp)`) %>% 
     mutate(
       Genes = strsplit(as.character(Genes), split = "|", fixed=T),
       interval_length = `1.5 LOD Confidence Interval (Right) (bp)` - `1.5 LOD Confidence Interval (Left) (bp)`
     ) %>%
  mutate(num_genes = map_int(Genes, length)) %>%
  arrange(-num_genes) %>%
  select(Trait, Pheno_fraction_explained, Peak_pos, interval_length, num_genes, Genes)

filtered_qtl_gene_counts_df <- qtl_gene_counts_df %>%
  filter(Pheno_fraction_explained >0.05) %>% 
  group_by(Trait) %>% 
  reframe(all_genes = list(Genes), Total_variance = sum(Pheno_fraction_explained), num_genes=sum(num_genes))

for(i in 1:dim(filtered_qtl_gene_counts_df)[1]){
  
  trait <- filtered_qtl_gene_counts_df[i,1]$Trait
  x <- list()
  x[[trait]] <- filtered_qtl_gene_counts_df[i,2]$all_genes %>% unlist
  x[['SHAP']] <- shap_genes
  
  p <- ggvenn(x)
  save_plot(paste0("./results/result_contingency/venn/filtered_0.05/", trait, ".png"), 
            plot=p)
}

filtered_qtl_df <- filtered_qtl_gene_counts_df %>% 
  unnest_longer(all_genes) 

# QTL and contingency analysis ----

files <- list.files("./results/result_contingency/0.5_condition-wise/Bloom2013/", 
                    pattern = "csv$", full.names = TRUE)
filenames   <- list.files("./results/result_contingency/0.5_condition-wise/Bloom2013/", 
                          pattern = "csv$") %>% map_chr(
                            ~strsplit(.x, '.', fixed = T) %>% unlist %>% pluck(1)
                          )
dfs <- list()
i=1
for(file in files){
  df <- read_csv(file, )
  df <- df %>% mutate(Condition=filenames[i])
  dfs[[i]] <- df
  i = i+1
}
final_qtl_df <- bind_rows(dfs)

all_contingency_df <- final_qtl_df %>% 
  filter(adj.pval < 0.05) %>% 
  group_by(Condition) %>% 
  reframe(Genes = list(Gene))

for(i in 1:dim(all_contingency_df)[1]){
  
  trait <- all_contingency_df[i,1]$Condition
  x <- list()
  x[[paste0(trait,"_QTL")]] <- filter(all_qtl_gene_df, Trait==trait) %>%
    pluck('all_genes') %>% unlist
  x[[paste0(trait, '_odds')]] <- filter(all_contingency_df, Condition==trait) %>%
    pluck('Genes') %>%  unlist
  
  p <- ggvenn(x)
  save_plot(paste0("./results/result_contingency/0.5_condition-wise/Bloom2013/venn/",
                   trait, "-small.svg"), 
            plot=p)
}
## Filtered gene counts ----
for(i in 1:dim(all_contingency_df)[1]){
  
  trait <- all_contingency_df[i,1]$Condition
  x <- list()
  x[[paste0(trait,"_QTL")]] <- filter(filtered_qtl_df, Trait==trait) %>%
    pluck('all_genes') %>% unlist
  x[[paste0(trait, '_odds')]] <- filter(all_contingency_df, Condition==trait) %>%
    pluck('Genes') %>%  unlist
  
  p <- ggvenn(x, auto_scale = T, show_percentage = F)
  save_plot(paste0(
    "./results/result_contingency/0.5_condition-wise/Bloom2013/venn/filtered_0.05/",
                   trait, "_small.svg"), width=7, height=5, units='cm',
            plot=p, )
}

# Carbon source commonality ----

c_sources <- c('galactose', 'ethanol', 'lactose','maltose','mannose','raffinose',
               'sorbitol', 'trehalose', 'xylose')

top_genes_carbons <- final_qtl_df %>% 
  mutate(logodds = -log10(odds_ratio)) %>% 
  filter(adj.pval < 0.05) %>% 
  filter(Condition %in% c_sources) %>% 
  slice_max(n=50, logodds, by = Condition) %>% 
  pluck('Gene') %>% 
  unique()

enrichGO_yeast(top_genes_carbons) 

bottom_genes_carbons <- final_qtl_df %>% 
  mutate(logodds = -log10(odds_ratio)) %>% 
  filter(adj.pval < 0.05) %>% 
  filter(Condition %in% c_sources) %>% 
  slice_min(n=50, logodds, by = Condition) %>% 
  pluck('Gene') %>% 
  unique()

enrichGO_yeast(bottom_genes_carbons) 


top_genes_carbons <- final_qtl_df %>% 
  mutate(logodds = -log10(odds_ratio)) %>% 
  filter(adj.pval < 0.05) %>% 
  filter(Condition %in% c_sources) %>% 
  filter(logodds > 0) %>% 
  pluck('Gene') %>% 
  unique()

enrichGO_yeast(top_genes_carbons) 

oxidative_stresses <- c('formamide', 'paraquat', 'hydroquinone', 'menadione',
                        'CdCl2', 'ethanol', '4-hydroxybenzaldehyde', 'H2O2',
                        'MnSO4')

top_genes_oxidative <- final_qtl_df %>% 
  mutate(logodds = -log10(odds_ratio)) %>% 
  filter(adj.pval < 0.05) %>% 
  filter(Condition %in% oxidative_stresses) %>% 
  slice_max(n=50, logodds, by = Condition) %>% 
  pluck('Gene') %>% 
  unique()

enrichGO_yeast(top_genes_oxidative) %>% custom_treeplot()
save_plot(
  "./results/result_contingency/0.5_condition-wise/Bloom2013/oxidative_stress_tree_top.svg"
)

bottom_genes_oxidative <- final_qtl_df %>% 
  mutate(logodds = -log10(odds_ratio)) %>% 
  filter(adj.pval < 0.05) %>% 
  filter(Condition %in% oxidative_stresses) %>% 
  slice_min(n=50, logodds, by = Condition) %>% 
  pluck('Gene') %>% 
  unique()

enrichGO_yeast(bottom_genes_oxidative) 

genotoxic_stresses <- c('CoCl2', 'hydroxyurea','6-azauracil','fluorocytosine',
                        '4NQO', 'caffeine', 'fluorouracil', 'zeocin', 'cisplatin',
                        'methotrexate')
top_genes_genotoxic <- final_qtl_df %>% 
  mutate(logodds = -log10(odds_ratio)) %>% 
  filter(adj.pval < 0.05) %>% 
  filter(Condition %in% genotoxic_stresses) %>% 
  slice_max(n=50, logodds, by = Condition) %>% 
  pluck('Gene') %>% 
  unique()

enrichGO_yeast(top_genes_genotoxic) %>% custom_treeplot()

bottom_genes_genotoxic <- final_qtl_df %>% 
  mutate(logodds = -log10(odds_ratio)) %>% 
  filter(adj.pval < 0.05) %>% 
  filter(Condition %in% genotoxic_stresses) %>% 
  slice_min(n=50, logodds, by = Condition) %>% 
  pluck('Gene') %>% 
  unique()

enrichGO_yeast(bottom_genes_oxidative) 


# Pleiotropic genes

bloom2013_qtl_full = read_csv("./results/result_contingency/contingency_bloom2013_0.5sigma_full.csv") %>% 
  mutate(adj_pval = p.adjust(pval, 'fdr')) 
pleiotropic_qtl_df <- all_qtl_df %>%
  filter(grepl("^Y", Genes), Pheno_fraction_explained>0.05) %>%
  group_by(Genes) %>%
  reframe(traits=unique(list(Trait)), Count=n()) %>% 
  filter(Count>1) %>% 
  left_join((rename(shap_df, Genes=Feature)), 'Genes') %>% 
  left_join(rename(bloom2013_qtl_full, Genes=Gene), 'Genes') %>% 
  mutate(Significance = ifelse(adj_pval > 0.05, 'Not Significant', 'Significant'))

pleiotropic_qtl_df %>%
  ggplot(aes(x=log(Value), y=-log10(odds_ratio))) +
  geom_point(aes(color=Significance, size=Count)) +
  ggrepel::geom_text_repel(mapping = aes(label=Genes)) +
  labs(y="-log2(Odds Ratio)", x='log(SHAP value)')

save_plot('results/result_contingency/bloom2013_pleiotropic_odds_vs_shap.svg',
          width=10, height=6)
