library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(purrr)
# library(VennDiagram)
library(ggvenn)
library(gt)
library(dbscan)

source("./src/utils/theme_set.R")
source("./src/utils/enrichment_analysis.R")

# We've generated SHAP values from global model for each condition by 
# changing the dataset used for prediction

# Reading Datasets -----

bloom2013_shap_conditionwise <- arrow::read_parquet(
  "data/shap/sigmas/shap_classification_0.5/Bloom2013_conditionwise_shap.parquet"
) %>% 
  filter(Value > 0, grepl('^Y', Feature)) %>% 
  rename(Gene = Feature, SHAP=Value)

# bloom2015_shap_conditionwise <- arrow::read_parquet(
#   "data/shap/sigmas/shap_classification_0.5/Bloom2015_conditionwise_shap.parquet"
# ) %>% 
#   filter(Value > 0, grepl('^Y', Feature)) %>% 
#   rename(Gene = Feature, SHAP=Value)

bloom2019_shap_conditionwise <- arrow::read_parquet(
  "data/shap/sigmas/shap_classification_0.5/Bloom2019_conditionwise_shap.parquet"
) %>% 
  filter(Value > 0, grepl('^Y', Feature)) %>% 
  rename(Gene = Feature, SHAP=Value)

## QTL data ----

qtl_df <- read_csv(
  "./data/qtl/detected_qtl_bloom2013.csv"
)

all_qtl_df <- qtl_df %>% 
  select(Trait, Pheno_fraction_explained, Peak_pos, Chromosome, Genes) %>% 
  mutate(Genes = strsplit(as.character(Genes), split = "|", fixed=T)) %>% 
  unnest_longer(Genes) 

all_qtl_gene_df <- all_qtl_df %>% 
  group_by(Trait, Peak_pos) %>% 
  reframe(all_genes = list(Genes))

## Contingency data ---- 
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

final_contingency_df <- bind_rows(dfs)

all_contingency_df <- final_contingency_df %>% 
  filter(adj.pval < 0.05) %>% 
  group_by(Condition) %>% 
  reframe(Genes = list(Gene))

## Helping functions ----

get_genomic_position <- function(genename) {
  genename_vec <- strsplit(genename, "", T)[[1]] # get every character
  chrom_num <- match(genename_vec[2], LETTERS) # get chromosome number
  pos_on_chrom_arm <- purrr::reduce(genename_vec[4:6], paste0) %>% as.integer
  arm_factor <- ifelse(genename_vec[3]=='R', 0.3, -0.3)
  
  # Constructing final genomic position
  final_loc <- (chrom_num + pos_on_chrom_arm * 1e-4 + arm_factor) 
  return(final_loc)
}

get_final_df <- function(qtl_df, shap_df, contingency_df, condition) {
  qtl <- qtl_df %>% filter(Trait == condition) %>% rename(Gene=Genes)
  shap <- shap_df %>% filter(condition == {{condition}}) 
  contingency <- contingency_df %>% filter(Condition == condition, adj.pval < 0.05) 
  all <- contingency %>% 
    left_join(shap, 'Gene', suffix = c('_full', '_conditionwise')) %>%
    left_join(qtl, 'Gene') %>%
    replace_na(list(Pheno_fraction_explained=0)) %>% 
    select(-Trait, -count, -condition, -pval) %>%
    mutate(abs_log_odds = abs(log2(odds_ratio)), genomic_pos= map_dbl(Gene, get_genomic_position))
  
  # dbscan_results <- all %>% select(odds_ratio) %>% dbscan(0.1) %>% pluck('cluster') 
  # all$dbscan_cluster <- dbscan_results
  return(all)
}

get_final_df_bloom2013 <- partial(
  get_final_df, 
  qtl_df=all_qtl_df,
  contingency_df=final_contingency_df,
  shap_df = bloom2013_shap_conditionwise
)

get_final_df_bloom2019 <- partial(
  get_final_df, 
  qtl_df=all_qtl_df,
  contingency_df=final_contingency_df,
  shap_df = bloom2019_shap_conditionwise
)

plot_position_odds <- function(df, title, shap_size=T) {
  p <- df %>% arrange(Gene) %>% 
    mutate(abs_log_odds = abs(log2(odds_ratio)), index=1:n()) %>% 
    filter(adj.pval < 0.05) %>% 
    ggplot(aes(x=index, y=abs_log_odds)) +
    ggrepel::geom_text_repel(aes(label=Gene))
  
  if(shap_size) {
     p <- p + geom_point(mapping=aes(size=SHAP_conditionwise), alpha=0.7,
                    color='cornflowerblue') +
      labs(x='Genomic Position (relative)',
           y='Absolute Log-Odds',
           title=title,
           size='SHAP')
  } else {
    p <- p + geom_point( alpha=0.7, color='cornflowerblue') +
      labs(x='Genomic Position (relative)',
           y='Absolute Log-Odds',
           title=title)
  }
  return(p)
}

### 4NQO ----

nqo <- all_qtl_df %>% filter(Trait == '4NQO') %>% rename(Gene=Genes)
nqo_shap <- bloom2013_shap_conditionwise %>% filter(condition == '4NQO') %>% 
  mutate(rank_shap_conditionwise = rank(SHAP), shap_scaled = scale(SHAP))
nqo_contingency <- final_contingency_df %>%
  filter(adj.pval < 0.05, Condition == '4NQO')
nqo_all <- nqo_contingency %>%
  left_join(nqo_shap, 'Gene', suffix = c('_full', '_conditionwise')) %>%
  left_join(nqo, 'Gene') %>%
  replace_na(list(Pheno_fraction_explained=0)) %>% 
  select(-Trait, -count, -condition, -pval)

nqo_all <- get_final_df_bloom2013('4NQO')  

### maltose ----

mal <- all_qtl_df %>% filter(Trait == 'maltose') %>% rename(Gene = Genes)
mal_shap <- bloom2013_shap_conditionwise %>% filter(condition == 'maltose')
mal_contingency <- final_contingency_df %>% 
  filter(Condition == 'maltose') 
mal_all <- mal_contingency %>%
  left_join(mal_shap, 'Gene', suffix = c('_full', '_conditionwise')) %>%
  left_join(mal, 'Gene') %>%
  replace_na(list(Pheno_fraction_explained=0)) %>% 
  select(-Trait, -count, -condition, -pval)

### sorbitol ----

sorbitol_all <- get_final_df_bloom2013('sorbitol')

### cdcl2 
cdcl2 <- get_final_df_bloom2013('CdCl2')

### mannose

mannose_all <- get_final_df_bloom2013('mannose')

### neomycin

neomycin_all <- get_final_df_bloom2013('neomycin')

### congo_red

congo_all <- get_final_df_bloom2013('congo_red')

# Plotting everything ----

save_loc <- "results/result_contingency/ranking/Bloom2013/"

for(condition in unique(bloom2013_shap_conditionwise$condition)){
  print(condition)
  savename_shap <- paste0(save_loc, condition, '_shap.svg')
  savename_noshap <- paste0(save_loc, condition, '_noshap.svg')
  condition_df <- get_final_df_bloom2013(condition)
  
  p_shap <- plot_position_odds(condition_df, condition)
  p_noshap <- plot_position_odds(condition_df, condition, F)
  
  write_csv(condition_df, paste0(save_loc, condition, '.csv'))
  save_plot(filename = savename_noshap, plot=p_noshap)
  save_plot(filename = savename_shap, plot=p_shap)
}

save_loc <- "results/result_contingency/ranking/Bloom2019/"

for(condition in unique(bloom2019_shap_conditionwise$condition)){
  print(condition)
  savename_shap <- paste0(save_loc, condition, '_shap.svg')
  savename_noshap <- paste0(save_loc, condition, '_noshap.svg')
  condition_df <- get_final_df_bloom2019(condition)
  
  p_shap <- plot_position_odds(condition_df, condition)
  p_noshap <- plot_position_odds(condition_df, condition, F)
  
  write_csv(condition_df, paste0(save_loc, condition, '.csv'))
  save_plot(filename = savename_noshap, plot=p_noshap)
  save_plot(filename = savename_shap, plot=p_shap)
}


# Specific examples ----

plot_shap_odds <- function(df, title_val){
  df %>% 
    ggplot(., aes(x=abs_log_odds, y=SHAP_conditionwise)) + 
    geom_point(aes(color=targets)) +
    ggrepel::geom_text_repel(aes(label=Gene), max.overlaps = 7) +
    labs(title=title_val,x="Abs. Log Odds", y='SHAP Value', color='Targets with Evidence') %>% 
    return
}
save_loc <- "results/result_contingency/0.5_condition-wise/specific_qtls-colored/"

## 4NQO  - two QTLs 
nqo_all %>%
  filter(grepl('^YLR0[2|3|4|5]', Gene)) %>%
  mutate(targets = ifelse(Gene %in% c('YLR034C','YLR035C'), 'True', 'False')) %>% 
  plot_shap_odds(title_val = "4NQO - YLR023-YLR053 locus") %>% 
save_plot(
  paste0(save_loc, '4nqo-qtl1.svg'),
  plot=.,
  height=7, width=15, units='cm'
)

nqo_all %>%
  filter(grepl('^YNL0[7-9]', Gene)) %>%
  mutate(targets = ifelse(Gene %in% c('YNL085W','YLR035C'), 'True', 'False')) %>% 
  plot_shap_odds(title_val = "4NQO - YNL070-YNL090 locus") %>% 
save_plot(
  paste0(save_loc, '4nqo-qtl2.svg'),
  plot=.,
  height=7, width=15, units='cm'
)

## Sorbitol - 3 QTLs

sorbitol_all %>% 
  filter(grepl('^YNL1[0-4]', Gene)) %>% 
  mutate(targets = ifelse(Gene %in% c('YNL127W'), 'True', 'False')) %>% 
  plot_shap_odds("Sorbitol - YNL102 - YNL142 locus") %>% 
  save_plot(
    paste0(save_loc, 'sorbitol-qtl1.svg'),
    plot=.,
    height=7, width=15, units='cm'
  )
sorbitol_all %>% 
  filter(grepl('^YOL[0|1][6-9]', Gene)) %>% 
  mutate(targets = ifelse(Gene %in% c('YOL081W'), 'True', 'False')) %>%
  plot_shap_odds("Sorbitol - YOL060 - YOL106 locus") %>% 
  save_plot(
    paste0(save_loc, 'sorbitol-qtl2.svg'),
    plot=.,
    height=7, width=15, units='cm'
  )
sorbitol_all %>% 
  filter(grepl('^YOR0[3-5]', Gene)) %>% 
  mutate(targets = ifelse(Gene %in% c('YOR054C','YOR043W'), 'True', 'False')) %>%
  plot_shap_odds("Sorbitol - YOR035-YOR054 locus") %>% 
  save_plot(
    paste0(save_loc, 'sorbitol-qtl3.svg'),
    plot=.,
    height=7, width=15, units='cm'
  )

## Congo-red - 2 QTLs

congo_all %>%
  filter(grepl('^YB[L|R]0[1|0]', Gene)) %>% 
  mutate(targets = ifelse(Gene %in% c('YOR054C','YOR043W'), 'True', 'False')) %>%
  plot_shap_odds("CongoRed - YBL010-YBR010 locus") %>% 
  save_plot(
    paste0(save_loc, 'congo_red-qtl1.svg'),
    plot=.,
    height=7, width=15, units='cm'
  )
congo_all %>%
  filter(grepl('^YML[0|1][0|9|1|2]', Gene)) %>% 
  mutate(targets = ifelse(Gene %in% c('YML111W','YOR043W'), 'True', 'False')) %>%
  plot_shap_odds("CongoRed - YML070-YML100 locus") %>% 
  save_plot(
    paste0(save_loc, 'congo_red-qtl2.svg'),
    plot=.,
    height=7, width=15, units='cm'
  )

## Neomycin - one qtl

neomycin_all %>%
  filter(grepl('^YM[L|R]0[0-1]', Gene)) %>% 
  mutate(targets = ifelse(Gene %in% c('YOR054C','YOR043W'), 'True', 'False')) %>%
  plot_shap_odds("Neomycin - YML070-YML100 locus") %>% 
  save_plot(
    paste0(save_loc, 'neomycin-qtl1.svg'),
    plot=.,
    height=8, width=15, units='cm'
  )
