library(readr)
library(arrow)

source("./src/utils/enrichment_analysis.R")
source("./src/utils/theme_set.R")

get_unmutated_genes <- function(df) {
  genes <- df %>% 
    select(starts_with('Y')) %>% 
    select(\(x) sum(x)==0) %>% 
    colnames 
  print(paste("Number of unmutated genes: ", length(genes)))
  print(paste("Number of mutated genes:", 6014 - length(genes)))
  return(genes)
}

get_missense_genes <- function(df) {
  genes <- df %>% 
    select(starts_with('Y')) %>% 
    select(\(x) 1 %in% x) %>% 
    colnames
  print(paste0("Number of missense genes: ", length(genes)))
  return(genes)
}

get_nonsense_genes <- function(df) {
  genes <- df %>% 
    select(starts_with('Y')) %>% 
    select(\(x) 2 %in% x) %>% 
    colnames
  print(paste("Number of nonsense genes: ", length(genes)))
  return(genes)
}

get_frameshift_genes <- function(df) {
  genes <- df %>% 
    select(starts_with('Y')) %>% 
    select(\(x) 3 %in% x) %>% 
    colnames 
    print(paste("Number of frameshift genes: ", length(genes)))
  return(genes)
}

files <- list(
  Bloom2013 = "./data/training_data/0.5_sigma/bloom2013_clf.feather",
  Bloom2015 = "./data/training_data/0.5_sigma/bloom2015_clf.feather",
  Bloom2019 = "./data/training_data/0.5_sigma/bloom2019_clf.feather",
  Bloom2019_BYxM22 = "./data/training_data/0.5_sigma/bloom2019_BYxM22_clf.feather",
  Bloom2019_RMxYPS163 = "./data/training_data/0.5_sigma/bloom2019_RMxYPS163_clf.feather"
)

sigma <- 0.5

# read_functions <- list(
#   Bloom2013 = read_feather,
#   Bloom2015 = read_feather,
#   Bloom2019 = read_feather,
#   Bloom2019_BYxM22 = read_feather,
#   Bloom2019_RMxYPS163 = read_feather
# )

save_loc <- "./results/result_1/enrichment/"

for(file in names(files)) {
  
  print(file)
  
  df <- read_feather(files[[file]])
  
  unmutated_go <- get_unmutated_genes(df) %>% enrichGO_yeast()
  nonsense_go <- get_nonsense_genes(df) %>% enrichGO_yeast()
  frameshift_go <- get_frameshift_genes(df) %>% enrichGO_yeast()
  
  if(nrow(unmutated_go)!=0) {
    unmutated_go %>% custom_treeplot(20) %>%
      save_plot(paste0(save_loc, file, "_", sigma, "_", "unmutated_tree.svg"), plot=., 
                 height = 8,
                width = 12,
                units = 'in',)
    unmutated_go %>% go_heatplot() %>%
      save_plot(paste0(save_loc, file, "_", sigma, "_", "unmutated_heatplot.svg"), plot=.,
                 height = 8,
                width = 12,
                units = 'in',)
  }
  
  if(nrow(nonsense_go)!=0) {
    nonsense_go %>% custom_treeplot(20) %>%
      save_plot(paste0(save_loc, file,"_", sigma, "_", "unmutated_tree.svg"), plot=.,
                 height = 8,
                width = 12,
                units = 'in',)
    nonsense_go %>% go_heatplot() %>%
      save_plot(paste0(save_loc, file,"_", sigma, "_", "unmutated_heatplot.svg"), plot=.,
                 height = 8,
                width = 12,
                units = 'in',)
  }
  
  if(nrow(frameshift_go)!=0) {
    frameshift_go %>% custom_treeplot(20) %>% 
      save_plot(paste0(save_loc, file,"_", sigma, "_", "unmutated_tree.svg"), plot=.,
                 height = 8,
                width = 12,
                units = 'in',)
    frameshift_go %>% go_heatplot() %>%
      save_plot(paste0(save_loc, file,"_", sigma, "_", "unmutated_heatplot.svg"), plot=.,
                 height = 8,
                width = 12,
                units = 'in',)
  }
 
  
}

# Extra analysis ----

## Common missense

common_missense <- list(
  bloom2013 = read_feather(files$Bloom2013),
  bloom2015 = read_feather(files$Bloom2015),
  bloom2019 = read_feather(files$Bloom2019)
) %>% 
  map(get_missense_genes) %>%
  purrr::reduce(intersect)

p <- common_missense %>% parse_enrich_results(calc_pairwise = T) %>% slice_ego(n=1:10) %>% emapplot(
  showCategory = 10,
  force = 2,
  nodel_label = 'category',
  layout.params = list(layout = 'kk'),
  edge.params = list(show=T),
  cluster.params = list(
    cluster=T,  # Cluster the points
    legend=T,   # Put it in the legend
    n=4 # Number of clusters
  ),
  cex.params=list(category_node=0.7, category_label=0.7),
  color = 'Fold.Change'
) + 
  scale_fill_viridis_c(direction = -1, option = 'rocket') +
  scale_color_discrete()
  labs(fill = 'Fold Change', size='Number of Genes')
  
  
  
# Regression data ----
  
df1 <- read_feather(
  "./data/training_data/regression_data/bloom2013_regression_std.feather",
  c('Strain','Condition','Phenotype')
)

sampled_df1 <- df1 %>%
  group_by(Condition) %>%
  reframe(Growth = ifelse(Phenotype < mean(Phenotype), 'Low', 'High'), Strain, Phenotype) %>% 
  group_by(Condition) %>% 
  slice_head(n=10,)

p1 <- ggplot(sampled_df1, aes(x=Strain, y=Condition, fill=Phenotype)) +
  geom_tile() +
  scale_fill_continuous(type='viridis') + theme_void()

save_plot(
  "./results/Composite figures/Components/Continuous_growth.svg",
  width= 5,
  height = 7,
  units ='in',
  plot=p1
)

p2 <- ggplot(sampled_df1, aes(x=Strain, y=Condition, fill=Growth)) +
  geom_tile()  + 
  scale_fill_viridis_d() + theme_void()

save_plot(
  "./results/Composite figures/Components/Binarized_growth.svg",
  width= 5,
  height = 7,
  units ='in',
  plot=p2
)
