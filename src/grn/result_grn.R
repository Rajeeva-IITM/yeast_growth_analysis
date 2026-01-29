# Doing GRN analysis

library(gprofiler2)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(readr)
library(BioNERO)
library(doParallel)
library(plotly)
library(purrr)
library(tidygraph)
library(ggraph)

source("./src/utils/enrichment_analysis.R") # Load enrichment functions
source("./src/utils/theme_set.R") # Load theme

load("./data/grn/GRN_results_gimme_0.5_sigma_reimand.sav")  # Load grn

# Changing filenames appopriately
query_high <- readLines(
  "./data/gsmm/Intersection-gimme_pfba_0.5_sigma/intersection_high_strains-pFBAoptima.txt"
)
query_low <- readLines(
  "./data/gsmm/Intersection-gimme_pfba_0.5_sigma/intersection_low_strains-pFBAoptima.txt"
)
high_genes <- setdiff(query_high, query_low)

# geno_df <- read_csv("./data/genotype/Bloom2013_ynb_geno.csv", show_col_types = F, num_threads = 4)
geno_df <- read_csv(
  "./data/training_data/bloom2013_clf_sigma0.5_ynb.csv",
  show_col_types = F,
  num_threads = 4
)

# Binarizing mutations to mutated and unmutated
geno_df <- geno_df %>%
  mutate(across(starts_with('Y'), ~ if_else(.x == 0, 0, 1)))

# contingency_table <- read_csv("./results/contingency_bloom2013_ynb.csv", 
#show_col_types = F) %>%
#   mutate(adj_pval = p.adjust(pval,method='BH')) %>%
#   filter(adj_pval < 0.05)

contingency_table <- read_csv(
  "./results/result_contingency/contingency_bloom2013_ynb_0.5sigma.csv",
  show_col_types = F
) %>%
  mutate(adj_pval = p.adjust(pval, method = 'BH')) %>%
  filter(adj_pval < 0.05)

## Plotting the regulators and targets ----

# grn$Regulator %>%
#   unique() %>%
#   enrichGO_yeast() %>%
#   custom_treeplot(cluster.params.n = 4) %>%
#   save_plot("results/result_grn/Regulators_all_heatplot.svg",
#plot=., width=10, height=6, units='in')

grn$Target %>%
  unique() %>%
  enrichGO_yeast() %>%
  custom_treeplot(cluster.params.n = 4) %>%
  save_plot(
    "results/result_grn/0.5_reimand/Targets_all_heatplot_0.5_sigma.svg",
    plot = .,
    width = 12,
    height = 6,
    units = 'in'
  )

# high_grn <- filter(grn, Regulator %in% tfs_high) # Filtering those that are high genes
# high_grn2 <- filter(grn, Target %in% tfs_high)
high_grn3 <- filter(grn, Target %in% high_genes)
# plot_high_grn <- plot_grn(high_grn)
# plot_high_grn2 <- plot_grn(high_grn2)
plot_high_grn3 <- plot_grn(high_grn3, show_labels = "all")


# Check all regulators  ----
regulators_grn <- unique(grn$Regulator)
grn_go <- enrichGO_yeast(regulators_grn)

p1 <- custom_treeplot(grn_go, cluster.params.n = 4)
save_plot(
  "./results/result_grn/0.5_reimand/GRN_regulators_GO_full_treeplot.svg",
  height = 6,
  width = 10.5,
  units = 'in'
)

p2 <- go_heatplot(grn_go)
save_plot(
  "./results/result_grn/0.5_reimand/GRN_regulators_GO_full_heatplot.svg",
  height = 6,
  width = 18,
  units = 'in'
)

# Check for filtered for pFBA optimal high unique genes ----

regulators_grn_high <- unique(high_grn3$Regulator)
grn_go <- enrichGO_yeast(regulators_grn_high)

p1 <- custom_treeplot(grn_go, cluster.params.n = 4)
save_plot(
  "./results/result_grn/0.5_reimand/GRN_regulators_GO_pfba-high_treeplot.svg",
  height = 6,
  width = 10.5,
  units = 'in'
)

p2 <- go_heatplot(grn_go)
save_plot(
  "./results/result_grn/0.5_reimand/GRN_regulators_GO_pfba-high_heatplot.svg",
  height = 6,
  width = 13,
  units = 'in'
)

# SHAP commonality ----

shap_genes <- read_csv(
  "./data/shap/sigmas/shap_classification_0.5/shap_Bloom2013_Boosting.csv") %>%
  filter(Value > 1e-5)

## SHAP in regulator and target ----
grn_shap_regulators <- filter(grn, Regulator %in% shap_genes$Feature)
grn_shap_targets <- filter(grn, Target %in% shap_genes$Feature)

# Count the number of edges
grn_shap_regulator_genes <-  grn_shap_regulators$Regulator %>%
  table() %>%
  sort(decreasing = T)

shap_genes_regulators <- filter(shap_genes, Feature %in% names(grn_shap_regulator_genes)) %>%
  arrange(-Value) %>%
  mutate(Degree = as.integer(grn_shap_regulator_genes[Feature]))

# Similarly for Targets
grn_shap_target_genes <- grn_shap_targets$Target %>% table() %>% sort(decreasing = T)
shap_genes_targets <- filter(shap_genes, Feature %in% names(grn_shap_target_genes)) %>%
  arrange(-Value) %>%
  mutate(Degree = as.integer(grn_shap_target_genes[Feature]))

regulators_go <- grn_shap_regulator_genes %>% names %>% enrichGO_yeast()
targets_go <- grn_shap_target_genes %>% names %>% enrichGO_yeast()
regulator_shap_targets_go <- grn_shap_regulators$Target %>% enrichGO_yeast()

regulators_go %>%
  custom_treeplot(cluster.params.n = 3) %>%
  save_plot(
    "./results/result_grn/0.5_reimand/SHAP_regulators.svg",
    plot = .,
    height = 6,
    width = 10,
    units = 'in'
  )

targets_go %>%
  custom_treeplot(cluster.params.n = 4) %>%
  save_plot(
    "./results/result_grn/0.5_reimand/SHAP_targets.svg",
    plot = .,
    height = 6,
    width = 10,
    units = 'in'
  )

regulator_shap_targets_go %>%
  custom_treeplot(cluster.params.n = 3) %>%
  save_plot(
    "./results/result_grn/0.5_reimand/SHAP_regulators_targets.svg",
    plot = .,
    height = 6,
    width = 10,
    units = 'in'
  )

# SHAP and pFBA on same side ----

c_grn <- high_grn3 %>% filter(Target %in% shap_genes$Feature)

# SHAP and pFBA ----

connection_grn <- high_grn3 %>% filter(Regulator %in% shap_genes$Feature)
connection_go_target <- connection_grn$Target %>% enrichGO_yeast()
connection_go_regulator <- connection_grn$Regulator %>% enrichGO_yeast()

connection_go_target %>% custom_treeplot(20, cluster.params.n = 4) %>%
  save_plot(
    "./results/result_grn/0.5_reimand/SHAP_connection_targets.svg",
    plot = .,
    height = 6,
    width = 10,
    units = 'in'
  )
connection_go_regulator %>% custom_treeplot(20, cluster.params.n = 4) %>%
  save_plot(
    "./results/result_grn/0.5_reimand/SHAP_connection_regulators.svg",
    plot = .,
    height = 6,
    width = 10,
    units = 'in'
  )
p <- connection_go_regulator %>%
  custom_emapplot(20, 3, 'fr') + scale_fill_viridis_c(direction = -1, option =
                                                        'rocket')

save_plot(
  "./results/result_grn/0.5_reimand/SHAP_connection_regulators_map.svg",
  plot = p,
  height = 7,
  width = 9,
  units = 'in',
  scale = 1.25
)

# Plotting the graph ----

graph <- connection_grn %>%
  # mutate(`SHAP value` = connection_regulator_shap) %>%
  as_tbl_graph(directed = T) %>%
  mutate(Centrality = centrality_degree()) %>%
  mutate(Class = factor(ifelse(
    name %in% grn$Regulator, 'Regulator', 'Target'
  ))) %>%
  mutate(component = group_components())

get_shap <- function(gene) {
  value <- shap_genes %>%
    filter(Feature == gene) %>%
    pluck('Value')
  if (length(value) == 0) {
    return(0)
  }
  return(value)
}

get_odds_ratio <- function(gene) {
  value <- contingency_table %>%
    filter(Gene == gene) %>%
    pluck('odds_ratio')
  if (length(value) == 0) {
    return(-1)
  }
  return(value)
}

graph <- graph %>%
  activate(nodes) %>%
  mutate(SHAP = map_dbl(name, get_shap)) %>%
  mutate(xpos = ifelse(Class == 'Regulator', 0, 1)) %>%
  group_by(Class) %>%
  mutate(ypos = seq(0, 1, length.out = n()))

ggraph(graph, layout = 'kk') +
  geom_edge_fan(
    arrow = arrow(length = unit(1, 'mm')),
    start_cap = circle(1, 'mm'),
    end_cap = circle(1, 'mm')
  ) +
  geom_node_point(aes(size = SHAP, color = Class), fill = 'white') +
  geom_node_label(aes(label = name, color = Class), repel = T)

save_plot(
  "results/result_grn/0.5_reimand/connection_grn.svg",
  height = 14,
  width = 20,
  units = 'in'
)

# Stronger Connection GRN ----

str_cnxn_grn <- grn %>%
  filter(Regulator %in% contingency_table$Gene) %>%
  filter(Target %in% high_genes)

str_cnxn_grn %>% write_csv("./results/result_grn/0.5_reimand/str_connection_grn_df_edgelist.csv")

str_cnxn_grn_reg_go <- str_cnxn_grn$Regulator %>% enrichGO_yeast()
str_cnxn_grn_target_go <- str_cnxn_grn$Target %>% enrichGO_yeast()

reg_annot <- str_cnxn_grn_reg_go %>%
  simplify(select_fun = max) %>%
  get_gene_annotation() %>%
  mutate(Class = 'Regulator')
target_annot <- str_cnxn_grn_target_go %>%
  simplify(select_fun = max) %>%
  get_gene_annotation() %>%
  mutate(Class = 'Target')

bind_rows(reg_annot, target_annot) %>%
  write_csv("./results/result_grn/0.5_reimand/Go_annotations_str_grn.csv")

str_cnxn_graph <- str_cnxn_grn %>%
  # mutate(`SHAP value` = connection_regulator_shap) %>%
  as_tbl_graph(directed = T) %>%
  mutate(Centrality = centrality_degree()) %>%
  mutate(Class = factor(ifelse(
    name %in% grn$Regulator, 'Regulator', 'Target'
  ))) %>%
  mutate(component = group_components()) %>%
  activate(nodes) %>%
  mutate(SHAP = map_dbl(name, get_shap)) %>%
  mutate(odds_ratio = map_dbl(name, get_odds_ratio)) %>%
  mutate(xpos = ifelse(Class == 'Regulator', 0, 1)) %>%
  group_by(Class) %>%
  mutate(ypos = seq(0, 1, length.out = n()))

as_data_frame(str_cnxn_graph) %>%
  write_csv("./results/result_grn/0.5_reimand/str_connection_grn_df.csv")


## HIGHLY CUSTOMIZED FIGURE WARNING ----- Cannot be used generally

ggraph(str_cnxn_graph,
       layout = 'manual',
       x = xpos,
       y = ypos) +
  geom_edge_fan(
    arrow = arrow(length = unit(1, 'mm')),
    start_cap = circle(1, 'mm'),
    end_cap = circle(1, 'mm')
  ) +
  geom_node_point(aes(
    size = SHAP,
    color = -log10(odds_ratio),
    shape = Class
  ), fill = 'white') +
  geom_node_label(aes(label = name), repel = T, ) +
  scale_color_continuous(type = "viridis") +
  scale_size(range = c(3, 8))

save_plot(
  "results/result_grn/0.5_reimand/connection_grn_strong.svg",
  height = 20,
  width = 15,
  units = 'in'
)


# MIG specific networks ----

mig_genes <- c("YGL209W", "YGL035C", "YER028C")

mig_graph <- grn %>% filter(Regulator %in% mig_genes) %>% as_tbl_graph(directed = T) %>%
  mutate(Class = factor(ifelse(
    name %in% grn$Regulator, 'Regulator', 'Target'
  ))) %>%
  activate(nodes) %>%
  mutate(SHAP = map_dbl(name, get_shap))

ggraph(mig_graph, layout = 'kk') +
  geom_edge_fan(
    arrow = arrow(length = unit(1, 'mm')),
    start_cap = circle(1, 'mm'),
    end_cap = circle(1, 'mm')
  ) +
  geom_node_point(aes(size = SHAP, color = Class), fill = 'white') +
  geom_node_label(aes(label = name, ), repel = T)  +
  scale_color_manual(values = cbPalette)

# All regulator specific network in strong_connection_grn ----

str_regulators <- str_cnxn_grn$Regulator %>% unique()

for (regulator in str_regulators) {
  save_loc <- "./results/result_grn/0.5_reimand/subgraphs/"
  graph_filename <- paste0(save_loc, regulator, ".svg")
  enrich_filename <- paste0(save_loc, regulator, '_treeplot', ".svg")
  barplot_filename <- paste0(save_loc, regulator, '_barplot', ".svg")
  
  subgrn <- str_cnxn_grn %>%  filter(Regulator == {{regulator}})
  subgraph <- subgrn  %>%
    as_tbl_graph(directed = T) %>%
    mutate(Class = factor(ifelse(
      name %in% grn$Regulator, 'Regulator', 'Target'
    ))) %>%
    mutate(SHAP = map_dbl(name, get_shap)) %>%
    mutate(odds_ratio = map_dbl(name, get_odds_ratio))
  
  p1 <- ggraph(subgraph, layout = "kk") +
    geom_edge_fan(
      arrow = arrow(length = unit(1, 'mm')),
      start_cap = circle(1, 'mm'),
      end_cap = circle(1, 'mm')
    ) +
    geom_node_point(aes(size = SHAP, color = Class), fill = 'white') +
    geom_node_label(aes(label = name, ), repel = T)  +
    scale_color_manual(values = cbPalette)
  
  save_plot(
    graph_filename,
    plot = p1,
    height = 7,
    width = 7,
    units = 'in'
  )
  
  enrichment <- subgrn$Target %>% enrichGO_yeast()
  
  if (dim(enrichment)[1] == 0) {
    print(paste0("No enrichment for ", regulator))
    next()
  }
  
  p2 <- enrichment %>% parse_enrich_results() %>%
    pluck('result') %>%
    ggplot(aes(
      x = Fold.Change,
      y = reorder(Description, Fold.Change),
      fill = -log10(p.adjust)
    )) +
    geom_bar(stat = 'identity') +
    labs(x = 'Fold Change', y = 'GO Process', fill = '-log10(Adj. p-value)')
  
  save_plot(
    barplot_filename,
    plot = p2,
    height = 10,
    width = 15,
    units = 'in'
  )
  
  if (dim(enrichment)[1] == 1) {
    print(paste0("only 1 enrichment for ", regulator))
    next()
  }
  
  enrichment %>% custom_treeplot() %>%
    save_plot(
      enrich_filename,
      plot = .,
      height = 8,
      width = 14,
      units = 'in'
    )
  
}

# Genes mentioned in result ----

mentioned_genes <- list(
  YER045C = 'ACA1',
  YJL110C = 'GZF3',
  YGL254W = 'FZF1',
  YNL167C = 'SKO1',
  YNL068C = 'FKH1',
  YMR019W = 'STB4',
  YGL209W = 'MIG2'
)

custom_grn <- connection_grn %>%
  filter(Regulator %in% names(mentioned_genes))

# Two genes MSD3 and PDR8 ----

mds3 <- readRDS('results/result_dge/0.5_redo/MDS3.rds')
pdr8 <- readRDS('results/result_dge/0.5_redo/PDR8.rds')
result_path <- 'results/result_grn/0.5_reimand/novel-genes/'

## MDS3 first ----

regulator <- 'YGL197W'
coef_name <- paste(regulator, 'Unmutated', 'vs', 'Mutated', sep = '_')

res_mds3 <- DESeq2::lfcShrink(mds3, coef = coef_name, type = "ashr")

targets <- filter(str_cnxn_grn, Regulator == regulator) %>% pluck('Target')
for (target in targets) {
  count_data <- plotCounts(mds3, target, regulator, returnData = TRUE)
  
  pval_df <- tibble(
    group1 = 'Mutated',
    group2 = 'Unmutated',
    padj = formatC(res_mds3[target, 'padj'], format = 'e', digits = 2),
    y.position = max(count_data$count) + 50
  )
  
  count_plot <- count_data %>%
    ggboxplot(regulator,
              'count',
              fill = regulator,
              ylab = paste(target, 'expression')) +
    stat_pvalue_manual(pval_df, label = 'Adjusted p-val = {padj}') + 
    theme(legend.position = 'none')
  
  save_plot(
    file.path(result_path, paste0('mds3_',target,'.svg')),
    height=95, width=95, units='mm'
  )
}

## PDR8 first ----

regulator <- 'YLR266C'
coef_name <- paste(regulator, 'Unmutated', 'vs', 'Mutated', sep = '_')

res_pdr8 <- DESeq2::lfcShrink(pdr8, coef = coef_name, type = "ashr")

targets <- filter(str_cnxn_grn, Regulator == regulator) %>% pluck('Target')
for (target in targets) {
  count_data <- plotCounts(pdr8, target, regulator, returnData = TRUE)
  
  pval_df <- tibble(
    group1 = 'Mutated',
    group2 = 'Unmutated',
    padj = formatC(res_pdr8[target, 'padj'], format = 'e', digits = 2),
    y.position = max(count_data$count) + 50
  )
  
  count_plot <- count_data %>%
    ggboxplot(regulator,
              'count',
              fill = regulator,
              ylab = paste(target, 'expression')) +
    stat_pvalue_manual(pval_df, label = 'Adjusted p-val = {padj}')+ 
    theme(legend.position = 'none')
  
  save_plot(
    file.path(result_path, paste0('pdr8_',target,'.svg')),
    height=95, width=95, units='mm'
  )
}


# Differentially expressed PDR8 gene

 count_plot <- count_data %>%
    ggboxplot(regulator,
              'count',
              fill = regulator,
              ylab = paste(target, 'expression')) +
    stat_pvalue_manual(pval_df, label = 'Adjusted p-val = {padj}')
 
 save_plot(
    file.path(result_path, paste0('pdr8_','dummy','.svg')),
    height=95, width=95, units='mm'
  )

res <- c('YDL095W','YOR321W','YLR028C','YGR010W','YDL093W') %>%
  enrichGO_yeast() %>%
  parse_enrich_results() %>% 
  pluck('result')

p <- res %>%
  filter(p.adjust < 0.05) %>% 
  slice_min(p.adjust, n=10) %>%
  ggplot(
    .,
    aes(
      x=FoldEnrichment,
      y=reorder(Description, FoldEnrichment),
      fill=-log10(p.adjust)
    )
  ) + 
  geom_bar(stat='identity',) +
  scale_fill_continuous(guide=guide_colorbar(), type='viridis') +
  labs(x = 'Fold Change', y= 'GO Biological Process', fill='-log10(Adj. p-value)')


save_plot(
  "./results/result_grn/0.5_reimand/novel-genes/pdr8_diff_expression_enrichment.svg",
  width=147.374, height=105.267, units='mm', plot=p, scale=1.5)
