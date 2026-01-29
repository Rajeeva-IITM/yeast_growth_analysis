library(dplyr)
library(ggplot2)
library(purrr)

source('./src/utils/theme_set.R') # Plotting consistency
source('./src/utils/enrichment_analysis.R')  # Enrichment codes

data_path <- "./data/shap_features/classification/0.5_sigma/"
data_files <- list.files(data_path)

out_path_plot <- "./results/result_4_shap/sigmas/0.5/enrichment/"
out_path_enrichment <- './data/enrichment/classification/0.5/'

get_suffix <- function(filename) {
  if (grepl('common', filename)) {
    return('common')
  } else if (grepl('full', filename)) {
    return('full')
  } else {
    return('model_unique')
  }
}

check_empty <- function(enrichResult) {
  if (is.null(enrichResult)) {
    return(TRUE)
  }
  enrichResult %>%
    pluck('result') %>%
    filter(p.adjust < 0.05) %>%
    isEmpty() %>%
    return
}

for (file in data_files) {
  print(file)
  
  suffix <- get_suffix(file)
  
  # If the suffix is common, we have a list of genes that we can directly use
  if (suffix == 'common') {
    genes <- paste0(data_path, file) %>% readRDS()
    
    enrichment_go <- genes %>% enrichGO_yeast()
    
    if (check_empty(enrichment_go)) {
      print(paste0(file, " : Pretty empty enrichment"))
      next
    }
    
    paste0(out_path_enrichment, file) %>% saveRDS(enrichment_go, file =
                                                    .) # Save the file
    
    p <- enrichment_go %>% parse_enrich_results(calc_pairwise = T) %>% custom_treeplot(top_category = 20)
    
    gsub('.RDS', '.svg', file) %>%
      paste0(out_path_plot, .) %>%
      save_plot(
        .,
        height = 5,
        width = 8,
        units = 'in',
        scale = 1.5
      )
    
    # Save heatplot of the data
    p_heatplot <- enrichment_go %>% go_heatplot(10)
    
    gsub('.RDS', '_heatplot.svg', file) %>%
      paste0(out_path_plot, .) %>%
      save_plot(
        .,
        height = 5,
        width = 8,
        units = 'in',
        scale = 1.5
      )
    
  } else {
    gene_list <-  paste0(data_path, file) %>% readRDS()
    
    for (model in names(gene_list)) {
      #Need to do for each model
      
      genes <- gene_list[[model]]
      
      if (is_empty(genes)) {
        print(paste0(file, " : ", model, " : Pretty empty gene set"))
        next
      }
      
      enrichment_go <- genes %>% enrichGO_yeast()
      gost_df <- get_gost_df(genes)
      
      if (!is.null(gost_df)) {
        # If enrichment present
        
        p <- plot_gost_result(gost_df, 'KEGG') # Get plot for kegg
        
        if (!is.null(p)) {
          # If enrichment for Kegg present
          filename <- gsub('.RDS', '.csv', file) %>%
            paste0(out_path_enrichment, 'KEGG_', .)
          
          readr::write_csv(gost_df, filename)
          
          #saving figure
          filename %>%
            gsub(".csv", ".svg", .) %>%
            save_plot(
              .,
              height = 5,
              width = 7,
              units = 'in',
              scale = 1.5
            )
        }
        
        p <- plot_gost_result(gost_df, 'TF') # Get plot for TF
        
        if (!is.null(p)) {
          # If enrichment for TF present
          filename <- gsub('.RDS', '.csv', file) %>%
            paste0(out_path_enrichment, 'TF_', .)
          
          readr::write_csv(gost_df, filename)
          
          filename %>%
            gsub(".csv", ".svg", .) %>%
            save_plot(
              .,
              height = 5,
              width = 7,
              units = 'in',
              scale = 1.5
            )
        }
      }
      
      if (check_empty(enrichment_go)) {
        print(paste0(file, " : ", model, " : Pretty empty enrichment"))
        next
      }
      
      paste0(out_path_enrichment, file) %>% saveRDS(enrichment_go, file =
                                                      .)
      p <- enrichment_go %>%
        parse_enrich_results(calc_pairwise = T) %>%
        (\(x) tryCatch(
          custom_treeplot(x, top_category = 20),
          error = function(e) {
            paste0("ERROR :", conditionMessage(e), "\n")
          }
        ))
      
      paste0("_", model, ".svg") %>%
        gsub('.RDS', ., file) %>%
        paste0(out_path_plot, .) %>%
        save_plot(
          .,
          height = 5,
          width = 8,
          units = 'in',
          scale = 1.5
        )
      
      p_heatplot <- enrichment_go %>% go_heatplot(10)
      
      paste0("_", model, "_heatplot.svg") %>%
        gsub('.RDS', ., file) %>%
        paste0(out_path_plot, .) %>%
        save_plot(
          .,
          height = 5,
          width = 8,
          units = 'in',
          scale = 1.5
        )
    }
    
  }
}