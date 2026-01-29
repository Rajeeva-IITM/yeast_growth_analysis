# latent interpretability results

# To do:
# 1. Identify chemical groups for each latent variable
# 2. Correlate with important features from GenoChem models

library(dplyr)
library(arrow)
library(ggplot2)
library(readr)
library(purrr)

source('./src/utils/theme_set.R')



data_path <- "./data/latent_interpretability/"
savedirs <- list.files("./results/result_latent/", full.names = T)
ig_files <- list.files(data_path, pattern= "__IG", full.names = T) # Integrated gradients
nt_files <- list.files(data_path, pattern="__NT", full.names = T)  # Noise Tunnel
sl_files <- list.files(data_path, pattern='__Saliency', full.names = T) # Saliency
gs_files <- list.files(data_path, pattern='__GradientShap', full.names = T) # Gradient SHAP
ixg_files <- list.files(data_path, pattern='__Inpu', full.names = T) # InputXGradient

# Raw feature counts
feature_count_files <- list.files("./data/latent_interpretability/feature_count_chemicals/", full.names = T)

# Defining a function to filter extreme values and plot the attribution plots
plot_all_latent <- function(attr_df, savedir, savename, title, top_val=10) {
  
  for(i in 0:255) {
    
    feature_name <- paste0('latent_', i)
    final_save_loc <- paste0(savedir, '/', savename, '_', feature_name, '.svg') 
    
    # Postive Features
    top <- attr_df %>%
      select('feature_names', 'feature_types', feature_name) %>%
      slice_max( .data[[feature_name]], n=top_val)
    
    # Negative Features
    bottom <- attr_df %>% 
      select('feature_names', 'feature_types', feature_name) %>%
      slice_min( .data[[feature_name]], n=top_val)
    
    plot_df <- bind_rows(top, bottom)
    
    
    p <- plot_df %>% 
      mutate(
        feature_types = as.factor(
          case_match(
            feature_types,
            'fg'~'Functional Groups',
            'mfg'~'Mined Functional Groups'
          )
        ),
        feature_names = make.unique(feature_names, '_')
      ) %>% 
      ggplot(aes(y=reorder(feature_names, .data[[feature_name]]), x = .data[[feature_name]], fill=feature_types)) +
      geom_bar(stat = 'identity') + 
      labs(y='Features', x='Attributions', title = paste(feature_name, title, sep = " - ")) +
      guides(fill=guide_legend(title = 'Feature Type'))
    
    save_plot(final_save_loc, p, height=6, width=10, units='in')
  }
  
}

for(i in 1:6){
  
  savedir <- savedirs[i]
  print(savedir)
  
  read_parquet(ig_files[i]) %>% 
    plot_all_latent(., savedir, 'IG', 'Integrated Gradients')
  print('Integrated Gradients Done')
  
  read_parquet(nt_files[i]) %>% 
    plot_all_latent(., savedir, 'NT', 'Noise Tunnel - Integrated Gradients')
  print('Noise Tunnel Done')
  
  read_parquet(sl_files[i]) %>% 
    plot_all_latent(., savedir, 'SL', 'Saliency')
  print('Saliency Done')
  
  read_parquet(gs_files[i]) %>% 
    plot_all_latent(., savedir, 'GS', "Gradient Shap")
  print('GradientShap Done')
  
  read_parquet(ixg_files[i]) %>% 
    plot_all_latent(., savedir, 'IxG', 'InputxGradients')
  print('InputxGradients Done')
}

# Get feature counts at the dataset level ----

fc_files <- list.files("./data/latent_interpretability/feature_count_chemicals/", )
for(i in 1:6){
  
  savename <- strsplit(fc_files[i], split = '\\.') %>% pluck(1) %>% pluck(1)
  p <- read_parquet(
    paste0("./data/latent_interpretability/feature_count_chemicals/", fc_files[i]),
    col_select = c('feature_names', 'feature_types', 'feature_count')
  ) %>% 
    mutate(
      feature_types = as.factor(
        case_match(
          feature_types,
          'fg'~'Functional Groups',
          'mfg'~'Mined Functional Groups'
        )
      ),
      feature_names = make.unique(feature_names, '_')
    ) %>% 
    slice_max(feature_count, n=20) %>% 
    ggplot(., aes(x = feature_count, y = reorder(feature_names, feature_count), fill=feature_types)) +
    geom_bar(stat = 'identity') +
    labs(y='Features', x='Counts', title = savename) +
    guides(fill=guide_legend(title = 'Feature Type'))
  
  paste0("./results/result_latent/", savename, '.svg') %>% 
  save_plot(., p, height=6, width=10, units='in')
}

# Function to plot attributions across different methods for a given latent_variable

plot_latent <- function(
    latent_name,                 # Name of the variable considered
    path_to_data = data_path,       # Path to the data
    data_name = "bloom2013",     # Name of the dataset 
    plot_top=5,                  # The top groups to be plotted
    percentile_considered = 0.9  # Top features to be considered in all
                                 # datasets to identify the `plot_top` features
  )  {
  
  top_dataframes <- list()
  
  data_files <- list.files(path_to_data, pattern = data_name, full.names = T)
  
  i = 0
  for(data_file in data_files){
    
    if(grepl("GradientShap", data_file)) {
      model_name = "GradientShap"
    } else if(grepl("NT", data_file)) {
      next()
    } else if(grepl("IG", data_file)) {
      model_name = "Integrated Gradients"
    } else if(grepl("InputxGradient", data_file)) {
      model_name = "InputxGradient"
    } else {
      model_name = "Saliency"
      next()
    }
    
    i = i+1
    
    latent <- read_parquet(data_file) %>%  # Read file
      select(feature_names, feature_types, all_of(latent_name)) %>%  # Choose only the necessary columns
      # Create two columns with one indicating method used and the attribution
      mutate(Attribution_method = model_name, attribution = .data[[latent_name]]) %>% 
      # Modify attribution such that total attribution is 1
      # mutate(
      #   attribution = attribution/sum(attribution)
      # ) %>%
      # Filter the top `percentile_considered` attribution
      filter(attribution > quantile(attribution, percentile_considered)) 
      
    
    top_dataframes[[i]] <- latent
    
    
  }
  # print(top_dataframes)
  
  final_df <- bind_rows(top_dataframes) %>% # Concatenate all the dataframes
    filter(n()>2, .by = feature_names) %>%  # Filter out features that agreed upon by all models
    mutate(avg_attribution = mean(attribution), .by = feature_names) %>%  
    # arrange(-avg_attribution) %>% 
    slice_max(avg_attribution, n=5, by=Attribution_method) 
  
  plot <- ggplot(
    final_df,
    aes(x=feature_names, y=attribution, fill=Attribution_method)
  ) +
    geom_bar(stat='identity', position='dodge', color='black') +
    scale_fill_manual(values = cbPalette) + 
    theme(
      axis.text.x = element_text(angle = 15, hjust = 1)
    )
  
  return(plot)
  
}
plot_latent('latent_206') %>% save_plot("./results/result_latent/Bloom2013/Combined results/Latent_206.svg", plot=., width=8, height=6, units='in')
  
  

