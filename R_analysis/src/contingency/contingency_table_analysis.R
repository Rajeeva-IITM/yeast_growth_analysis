# Contingency table analysis - Only for YNB Media

library(dplyr)
library(data.table)
library(readr)
library(purrr)
library(ggplot2)
library(ggpubr)

source("./src/utils/theme_set.R")
source("./src/utils/enrichment_analysis.R")


# concordant_strains <- fread("./data/phenotype/Experimentally_extreme_strains.csv")
# concordant_strains <-
#   concordant_strains[, Growth := factor(Growth, levels = c("Low", "High"), ordered = FALSE)]

# Other sigma
# concordant_strains <- read_csv(
#   "./data/training_data/bloom2013_clf_sigma0.25_ynb.csv",
#   col_select =  c("Strain", "Phenotype")
# )

# strains <- concordant_strains$Strain
# 
# geno_df <- read_csv("./data/genotype/Bloom2013_ynb_geno.csv", show_col_types = F, num_threads = 4) 
# 
# geno_df <- geno_df %>% 
#   mutate(across(starts_with('Y'), ~if_else(.x==0, 0, 1))) %>% 
#   filter(Strain %in% strains)

geno_df <- read_csv(
  "./data/training_data/bloom2013_clf_sigma0.5_ynb.csv",
  num_threads = 4,
  show_col_types = F
) # Varying sigma

geno_df <- geno_df %>% 
  select(!starts_with("l")) %>% 
  filter(Phenotype != -1, Condition=='ynb_glucose') %>% 
  mutate(Growth = ifelse(Phenotype==0, "Low", "High")) %>% 
  mutate(Growth = factor(Growth, levels = c("Low", "High"), ordered = FALSE))

genenames <- colnames(geno_df)[grepl("^Y", colnames((geno_df)))]

result <- list(
  odds_ratio = rep(-1, 6014),
  pval = rep(-1, 6014)
)

for(i in 1:6014){
  
  mut_table <- table( # Contingency table
    geno_df$Growth, geno_df[[genenames[i]]]
  ) 
  
  if(sum(dim(mut_table))!=4){
    next()
  }
  
  fisher_result <- fisher.test(mut_table)
  
  result$odds_ratio[i] <- fisher_result %>% pluck("estimate") %>%  pluck("odds ratio")
  result$pval[i] <- fisher_result %>% pluck("p.value")
}

contingency_df <- result %>% 
  as_tibble() %>% 
  mutate(Gene = genenames)


contingency_df %>%
  filter(pval != -1, odds_ratio!=-1) %>% 
  write_csv(., "./results/result_contingency/contingency_bloom2013_ynb_0.5sigma.csv")


# Across the entirety of the dataset ----

geno_df <- arrow::read_feather("./data/training_data/0.5_sigma/bloom2015_clf.feather") %>% 
  mutate(Growth = ifelse(Phenotype==0, "Low", "High")) %>% 
  mutate(Growth = factor(Growth, levels = c("Low", "High"), ordered = FALSE))

genenames <- colnames(geno_df)[grepl("^Y", colnames((geno_df)))]

result_global <- list(
  odds_ratio = rep(-1, length(genenames)),
  pval = rep(-1,  length(genenames))
)


pb1 <- progress::progress_bar$new(format = "[:bar],[Estimated time left - :eta]",
                                  total = length(genenames))
for(i in 1:length(genenames)){
  
  mut_table <- table( # Contingency table
    geno_df$Phenotype, geno_df[[genenames[i]]]
  ) 
  
  if(sum(dim(mut_table))!=4){
    next()
    pb1$tick()
  }
  
  fisher_result <- fisher.test(mut_table)
  
  result_global$odds_ratio[i] <- fisher_result %>% pluck("estimate") %>%  pluck("odds ratio")
  result_global$pval[i] <- fisher_result %>% pluck("p.value")
  pb1$tick()
}


contingency_df_global <- result_global %>% 
  as_tibble() %>% 
  mutate(Gene = genenames)


contingency_df_global %>%
  filter(pval != -1, odds_ratio!=-1) %>% 
  write_csv(., "./results/result_contingency/contingency_bloom2015_0.5sigma_full.csv")


# Condition-wise -----

plot_contingency <- function(contingency_df, bounds=c(0.005,0.995)){
  contingency_df <- contingency_df %>% 
    mutate(Significance = ifelse(adj.pval < 0.05, "Significant", "Not Significant"))
  # print(contingency_df)
  q <- quantile(-log10(contingency_df[['odds_ratio']]), bounds)
  # print(q)
  annotation_df <- contingency_df %>% 
    filter(!log10(odds_ratio) %between% q)
  # print(annotation_df)
  
  contingency_df %>% 
    ggplot(aes(y = -log10(adj.pval), x = -log10(odds_ratio) )) +
    geom_point(mapping=aes(color=Significance,size=SHAP)) +
  labs(x = "-log10(Odds value)", y = "-log10(Adj. P-value)",
       title = "Odds Ratio Volcano")  +
    ggrepel::geom_text_repel(
      mapping = aes(label = Gene),
      data = annotation_df,
      color = "black",
      nudge_x = -0.0,
      nudge_y = -0.05,
      size = 4
    ) %>% 
    return
    
}


conditions <- geno_df$Condition %>% unique()
pb2 <- progress::progress_bar$new(total = length(conditions))
shap_df <- read_csv("data/shap/sigmas/shap_classification_0.5/shap_Bloom2015_Boosting.csv") %>% 
  rename(Gene = Feature, SHAP=Value)
for(condition in conditions){
  condition_df <- geno_df %>% filter(Condition == condition)
  result_condition <- list(
    odds_ratio = rep(-1, length(genenames)),
    pval = rep(1,  length(genenames))
  )
  for(i in 1:length(genenames)){
    
    mut_table <- table( # Contingency table
      condition_df$Phenotype, condition_df[[genenames[i]]]
    ) 
    
    if(sum(dim(mut_table))!=4){
      next()
    }
    
    fisher_result <- fisher.test(mut_table)
    
    result_condition$odds_ratio[i] <- fisher_result %>% pluck("estimate") %>%  pluck("odds ratio")
    result_condition$pval[i] <- fisher_result %>% pluck("p.value")
    
    result_condition_df <- result_condition %>% 
      as_tibble() %>% 
      mutate(Gene=genenames, adj.pval = p.adjust(pval, 'fdr')) %>% 
      filter(pval != 1, odds_ratio!=-1) %>% 
      left_join(shap_df, 'Gene')
      
  }
  result_condition_df %>% write_csv(
    paste0("./results/result_contingency/0.5_condition-wise/Bloom2015/", condition, '.csv')
  )
  p <- plot_contingency(result_condition_df)
  save_plot(paste0("./results/result_contingency/0.5_condition-wise/Bloom2015/", condition, '.svg'),
            # width=5, height = 4, units='in'
            )
  print(condition)
  # pb2$tick()
}

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

final_df_max <- bind_rows(dfs) %>% 
  filter(adj.pval<0.05) %>%
  mutate(logodds = -log10(odds_ratio)) %>%
  slice_max(n=10, logodds, by=Condition)

final_df_min <- bind_rows(dfs) %>% 
  filter(adj.pval<0.05) %>%
  mutate(logodds = -log10(odds_ratio)) %>%
  slice_min(n=10, logodds, by=Condition)


# Carbon source commonality ----

c_sources <- c('galactose', 'ethanol', 'lactose','maltose','mannose','raffinose',
               'sorbitol', 'trehalose', 'xylose')

oxidative_stresses <- c('formamide', 'paraquat', 'hydroquinone', 'menadione',
                        'CdCl2', 'ethanol', '4-hydroxybenzaldehyde', 'H2O2',
                        'MnSO4')
genotoxic_stresses <- c('CoCl2', 'hydroxyurea','6-azauracil','fluorocytosine',
                        '4NQO', 'caffeine', 'fluorouracil', 'zeocin', 'cisplatin',
                        'methotrexate')

