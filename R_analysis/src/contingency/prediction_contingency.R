library(dplyr)
library(data.table)
library(readr)
library(purrr)
library(ggplot2)
library(ggpubr)
library(arrow)
library(progress)
library(foreach)
library(doParallel)
library(metap)

cl <- makeCluster(4)
registerDoParallel(cl)

source("./src/utils/theme_set.R")
source("./src/utils/enrichment_analysis.R")

prediction_dir <- list.files("data/predictions/", full.names = T)
geno_files <- list.files("data/training_data/0.5_sigma/", pattern = "b.{8}_clf.feather", full.names = T)

df <- read_parquet(prediction_dir[2])
geno_df <- read_feather(geno_files[2])
genenames  <- read_csv(
  "./data/training_data/bloom2013_clf_sigma0.5_ynb.csv",
  num_threads = 4,
  show_col_types = F
) %>% 
  colnames %>% 
  keep(~ grepl('^Y',.x))

single_contingency <- function(pred_df_fold, geno_df_condition, gene) {
  if(length(pred_df_fold[['Preds']]) != length( geno_df_condition[[gene]])) {
    return(NULL)
  }
  result <- list()
  mut_table <- table( # Contingency table
    pred_df_fold[['Preds']], geno_df_condition[[gene]] )
  if(sum(dim(mut_table))!=4){ return(NULL) }
  fisher_result <- fisher.test(mut_table)
  result$gene <- gene
  result$odds_ratio <- fisher_result %>% pluck("estimate") %>% pluck("odds ratio")
  result$pval <- fisher_result %>% pluck("p.value")
  
  return(result)
}

perform_contingency_analysis <- function(geno_df, pred_df, genes=genenames) {
  # Ensure geno_df matches with tested dataset
  fold <- unique(pred_df[['Fold']]) # Get each fold 
  conditions <- intersect(pred_df[['Condition']], geno_df[['Condition']]) %>% unique()
  genes = genes
  
  cl <- makeCluster(12)
  registerDoParallel(cl)
  
  results_fold <- list()
  
  
  for(i in fold){
    pb1 <- progress_bar$new(total = length(conditions),
                            format = "Fold :fold Condition (:condition) - [:bar], Elapsed - :percent  ,ETA - :eta")
    result_condition <- list()
    for(condition in conditions){
      # pb2 <- progress_bar$new(total = 6014, format = "Genes - [:bar], Elapsed - :elapsed  ,ETA - :eta")
      geno_df_condition <- filter(geno_df, Condition==condition)
      pred_df_fold <- filter(pred_df, Fold == i, Condition==condition)
      
      result <- foreach(
        gene = genes,
        # .final = function(x) setNames(x, setNames(x, genes)),
        .packages = c('dplyr', 'purrr'),
        # .verbose = T,
        .export = c('single_contingency')
      ) %dopar% {
        single_contingency(pred_df_fold, geno_df_condition, gene)
      }

      pb1$tick(tokens = list(condition = condition, fold=i))
      
      result_df <- map_dfr(
        result,
        ~as_tibble(.x)
      )
      
      # return(result_df)
      
      result_condition[[condition]] <- result_df %>% 
        mutate(Condition = condition, Fold=i)   
      # print(result_condition)
    }
    results_fold[[i+1]] <- bind_rows(result_condition)
    
    }
  contingency_df <- bind_rows(results_fold)
  stopCluster(cl)
  return(contingency_df)
}

combine_odds <- function(odds){
  map_dbl(odds, log) %>% 
    mean %>% 
    exp %>% 
    return
}

combine_pvals <- function(pvals){
  fisher(pvals) %>% 
    pluck('p') %>% 
    return
}


bloom2013_2015 <- perform_contingency_analysis(
  read_feather(geno_files[2]),
  read_parquet(prediction_dir[2])
)

write_parquet(bloom2013_2015, "results/result_contingency/more/bloom2013_bloom2015.parquet")


bloom2013_2019 <- perform_contingency_analysis(
  read_feather(geno_files[3]),
  read_parquet(prediction_dir[3])
)

bloom2015_2013 <- perform_contingency_analysis(
  read_feather(geno_files[1]),
  read_parquet(prediction_dir[4])
)

write_parquet(bloom2015_2013, "results/result_contingency/more/bloom2015_bloom2013.parquet")

bloom2015_2019 <- perform_contingency_analysis(
  read_feather(geno_files[3]),
  read_parquet(prediction_dir[6])
)

bloom2019_2013 <- perform_contingency_analysis(
  read_feather(geno_files[1]),
  read_parquet(prediction_dir[7])
)

bloom2015_2019 <- perform_contingency_analysis(
  read_feather(geno_files[2]),
  read_parquet(prediction_dir[8])
)



# Condition-wise ----

shap_df <- read_parquet("data/shap/sigmas/shap_classification_0.5/shap_bloom2015_folds.parquet") %>%
  rename(Gene=Feature) %>%
  group_by(Gene) %>%
  reframe(SHAP=mean(Value)) 

process_df <- function(df, shap_df) {
  df %>%
    group_by(Condition, Fold) %>%
    reframe(Gene =gene, odds_ratio, adj_pval = p.adjust(pval, 'fdr')) %>%
    # filter(adj_pval < 0.05) %>%
    group_by(Condition, Gene) %>%
    reframe(
            odds_ratio = combine_odds(odds_ratio),
            adj.pval = combine_pvals(adj_pval)) %>% 
    left_join(shap_df, 'Gene') %>% 
    return
}

plot_contingency <- function(contingency_df, extreme=c(0.9995,0.0005)){
  contingency_df <- contingency_df %>% 
    mutate(Significance = ifelse(adj.pval < 0.05, "Significant", "Not Significant"))
  # print(contingency_df)
  q <- quantile(-log10(contingency_df[['odds_ratio']]), extreme)
  # print(q)
  annotation_df <- contingency_df %>% 
    # slice_max(SHAP, n= sum(extreme))
    filter(!(-log10(odds_ratio) %between% q))
  # print(annotation_df)
  
  contingency_df %>% 
    ggplot(aes(y = -log10(adj.pval), x = -log10(odds_ratio) )) +
    geom_point(mapping=aes(size=SHAP, color=Significance)) +
    labs(x = "-log10(Odds value)", y = "-log10(Adj. P-value)",
         title = "Odds Ratio Volcano")  +
    ggrepel::geom_text_repel(
      mapping = aes(label = Gene),
      data = annotation_df,
      color = "black",
      nudge_x = -0.0,
      nudge_y = -0.05,
      size = 4,
      max.overlaps = 20
    ) %>% 
    return
  
}

processed_2015_2013 <- process_df(bloom2015_2013, shap_df = shap_df)
