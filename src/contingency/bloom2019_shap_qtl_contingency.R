library(dplyr)
library(data.table)
library(readr)
library(purrr)
library(ggplot2)
library(ggpubr)

source("./src/utils/theme_set.R")
source("./src/utils/enrichment_analysis.R")

shap_df_2019 <- read_csv("data/shap/sigmas/shap_classification_0.5/shap_Bloom2019_BYxRM_Boosting.csv") %>%
  filter(Value > 1e-5) %>% rename(Gene=Feature)

shap_genes <- keep(shap_df_2019$Gene, ~grepl("^Y", .x)) 

causal_genes_2019 <- read_csv("data/papers/causal_genes_bloom2019.csv")
causal_genes_2019 <- causal_genes_2019 %>% filter(!grepl('^Y', trait)) %>% rename(Gene=ORF)

contingency_2019 <- read_csv("results/result_contingency/contingency_bloom2019_0.5sigma_full.csv") %>% 
  mutate(adj_pval = p.adjust(pval, 'fdr'))

final_df <- causal_genes_2019 %>% select(trait, Gene, jointmaxPPC)  %>% 
  left_join(shap_df_2019) %>% left_join(contingency_2019) %>% 
  mutate(signifance = ifelse(adj_pval < 0.05, 'Significant', 'Not-significant'))

plot_df <- final_df %>% filter(!is.na(Value)&!is.na(odds_ratio)) %>% group_by(Gene) %>% 
  summarise(jointmaxPPC = mean(jointmaxPPC),
            Value = mean(Value, na.rm=T),
            odds_ratio = mean(odds_ratio, na.rm=T),
            adj_pval = mean(adj_pval, na.rm=T),
            signifance = unique(signifance))



list(Causal = causal_genes_2019$Gene,
     Bloom2019_contingency =(filter(contingency_2019, adj_pval<0.05) %>% pluck('Gene'))) %>%
  ggvenn::ggvenn(auto_scale = TRUE)

save_plot("./results/result_contingency/venn/Bloom2019/causal_vs_contingency.svg",
          width=5, height = 4, units='in')

list(Causal = causal_genes_2019$Gene,
     Bloom2019_SHAP =shap_genes) %>%
  ggvenn::ggvenn(auto_scale = TRUE)

save_plot("./results/result_contingency/venn/Bloom2019/causal_vs_SHAP.svg",
          width=5, height = 4, units='in')

pleiotropic_genes <- causal_genes_2019 %>%
  group_by(Gene) %>%
  summarise(count=n()) %>%
  filter(count>1) %>% 
  pluck('Gene')

list(Causal_pleiotropic = pleiotropic_genes,
     Bloom2019_SHAP =shap_genes) %>%
  ggvenn::ggvenn(auto_scale = TRUE)

save_plot("./results/result_contingency/venn/Bloom2019/causal_pleiotropic_vs_SHAP.svg",
          width=5, height = 4, units='in')

list(Causal_pleiotropic = pleiotropic_genes,
     Bloom2019_contingency =(filter(contingency_2019, adj_pval<0.05) %>% pluck('Gene'))) %>%
  ggvenn::ggvenn(auto_scale = TRUE)

save_plot("./results/result_contingency/venn/Bloom2019/causal_pleiotropic_vs_contingency.svg",
          width=5, height = 4, units='in')


# Identify genes from SHAP not identified by contingency analysis

common_shap_pleiotropic <- intersect(pleiotropic_genes, shap_genes)
common_contingency_pleiotropic <- intersect(
  filter(contingency_2019, adj_pval<0.05) %>% pluck('Gene'),
  pleiotropic_genes
)
unique_shap_genes <- setdiff(common_shap_pleiotropic, common_contingency_pleiotropic)
