library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggcorrplot)
library(purrr)

source("./src/utils/theme_set.R")

# Read genotype dataset
geno_df <- arrow::read_feather("./data/training_data/0.5_sigma/bloom2013_clf.feather")

geno_df <- geno_df %>% select(Strain, starts_with('Y')) %>% distinct()
varying_cols <- geno_df %>% # Select only varying columns
  select(starts_with('Y')) %>% 
  summarise(across(everything(), sd)) %>%
  as.list %>% 
  (\(x) x!= 0)
varying_cols <- colnames(select(geno_df, starts_with('Y')))[varying_cols]

chromosomes <- map(LETTERS[1:16], ~ paste0('Y', .x))
dfs <- map(chromosomes, ~select(geno_df, starts_with(.x)))


dfs <- map(dfs, ~select(.x, any_of(varying_cols)))
cor_chrom <- map(dfs, cor)
plots <- map(cor_chrom, ~ggcorrplot(.x, tl.cex = 5))

walk(1:16, ~save_plot(paste0("results/result_genotype/bloom2013/", .x, '.svg'),
                      width=10, height=10, units='in', plot=plots[[.x]]))
