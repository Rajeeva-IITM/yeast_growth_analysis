# Code for plotting shap analysis data. Data is in the form of a csv file with feature name,
# shap value 

library(dplyr)
library(ggplot2)
library(readr)
library(purrr)

source('./src/utils/theme_set.R')

# Modify this section as per the data ----
# Data folder
data_path <-  "./data/shap/sigmas/shap_classification_0.5/" # Need to change for classification
data_files <-  list.files(path=data_path, pattern = ".*\\.csv",)
data_files <- keep(data_files, ~!grepl('RandomForest', .x))  # RandomForest the erroneous models

# Output dir
out_dir <-  "./results/result_4_shap/sigmas/0.5/"

# Don't modify this section per the data ----
get_plotting_shap_data <- function (shap_df,  top=10){
  # read the dataset
  # Take the top 9 and sum the SHAP values of the rest of the features
  # Final dataset should have 9 observations
  
  # shap_df <- read_csv(shap_df, col_types = 'cd')
  
  if(is.null(top)){
    return(shap_df)
  }
  
  top_10_df <- shap_df %>% arrange(-Value) %>% slice_head(n=top) # get first top
  
  # rest_SHAP <- shap_df %>% # calculate rest of the SHAP values
  #   arrange(-Value) %>%
  #   slice((top+1):n()) %>% 
  #   pluck('Value') %>% 
  #   sum()
  # 
  # top_10_df <- top_10_df %>% 
  #   add_row(Feature = 'Sum of the other 6261 features', Value = rest_SHAP)
  
  return(top_10_df)
}

plot_shap_data <- function(shap_df) {
  p1 <- ggplot(shap_df, aes(x=Value, y=reorder(Feature, desc(Value)))) +
    geom_bar(stat = "identity") +
    xlab("mean(|SHAP Value|)") +
    ylab("Features") + 
    theme(axis.title.y = element_text(angle = 90))
  
  return(p1)
}

generate_shap_plots <- function(data_files, latent_only=F) {  #Run this after all above
  for(file in data_files) {
    
    shap_df <- paste0(data_path, file) %>% read_csv(., col_types = "cd")
    
    savename <- strsplit(file, split='.', fixed=TRUE) %>%  # get file root
      unlist() %>%                             # flatten the result
      pluck(1)                              # The first element
    
    if(latent_only){
      shap_df <- filter(shap_df, grepl('latent', Feature))
      savename <- paste0(savename, "_latent")
    }
    
    p <- shap_df %>% 
      get_plotting_shap_data() %>% 
      plot_shap_data() 
    
    savename %>% 
      paste0(out_dir, ., '.svg') %>% 
      save_plot(, height=5, width=7, units="in", scale=1.2)
    
    p_hist <- shap_df %>% 
      filter(Value>1e-5) %>%
      ggplot(., aes(x=Value)) +
      geom_histogram(bins=50, fill='cornflowerblue')
    
    savename %>% 
      paste0(out_dir, 'shap_histograms/', ., '.svg') %>%  # Save in a subfolder
      save_plot(, height=5, width=7, units='in', scale=1.2)
  }
}

get_unique_features <- function(model_shap) {
  final_result <- list()
  for(model in names(model_shap)){ # 4 Models
    relevant_shap <- model_shap[[model]]
    irrlevant_shap <- discard_at(model_shap, model)
    final_result[[model]] <- map(irrlevant_shap, ~setdiff(relevant_shap, .x)) %>% 
      reduce(., intersect)
  }
  return(final_result)
}

all_datasets <- c('Bloom2013', 'Bloom2015', 'Bloom2019_BYxRM', 'Bloom2019_BYxM22', 'Bloom2019_RMxYPS163')
all_model_names <- c('Boosting', 'LogReg', 'RF', 'SVM')


# Comparing across different types of models 
generate_consensus_unique_shap <- function(savedir, datasets=all_datasets, model_names = all_model_names) {
  
  model_shap <- list()
  
  for(dataset in datasets){
    small_data_files <- keep(data_files, ~grepl(dataset, .x)) # Get the relevant csv files

    for(i in 1:length(small_data_files)){
      print(small_data_files[i])
      shap_df <- paste0(data_path, small_data_files[i]) %>% read_csv(., col_types = "cd")
      shap_features <- shap_df %>% filter(Value>1e-5) %>%  pluck('Feature')
      model_shap[[model_names[i]]] <- shap_features # First add all features and then reduce
    }
    
    if(length(small_data_files)>1){
    common_shap <- purrr::reduce(model_shap, intersect)  # Get the common SHAP features
    model_unique_shap <- get_unique_features(model_shap)  # Getting the unique genes for each model
    } else {
      return(0)
    }
    # Write the results:
    paste0(savedir, dataset, '_common.RDS') %>% saveRDS(common_shap, file=.)
    paste0(savedir, dataset, '_full.RDS') %>% saveRDS(model_shap, file=.)
    paste0(savedir, dataset, '_model_unique.RDS') %>% saveRDS(model_unique_shap, file=.)
  }
}

# generate_consensus_unique_shap("./data/shap_features/")

# Recurring features across Bloom2013, 2015 and 2019 that are in top 1 percentile

relevant_datasets <- c('Bloom2013', 'Bloom2015', 'Bloom2019_BYxRM')
relevant_model_names <- c('Boosting')
relevant_files <- keep(data_files, ~grepl(relevant_model_names, .x)) %>% 
  keep(., ~grepl('Bloom2013|Bloom2015|Bloom2019_BYxRM', .x))
data_shap <- list()
for(i in 1:length(relevant_files)){
  shap_df <- paste0(data_path, relevant_files[i]) %>% read_csv(., col_types = "cd")
  top_features <- shap_df %>% filter(Value > quantile(Value, 0.99)) %>% pluck('Feature')
  data_shap[[relevant_datasets[i]]] <- top_features
}

common_shap <- reduce(data_shap, intersect)
