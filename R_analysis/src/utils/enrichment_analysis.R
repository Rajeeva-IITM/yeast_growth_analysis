# Functions for enrichment analysis  

library(purrr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(gprofiler2)
library(GOSemSim)
library(DESeq2)
library(clusterProfiler)
library(data.table)
library(enrichplot)

# load("./data/miscellaneous/plotting_data_functions.Rdata")
load("./data/miscellaneous/scGO.Rdata") # Database file

load("./data/miscellaneous/all_genes.Rdata")

# GO enrichment function customized for yeast genes
enrichGO_yeast <- partial(
  enrichGO,
  universe = all_genes,
  keyType = 'ENSEMBL',
  OrgDb = "org.Sc.sgd.db",
  ont = 'BP',
  pAdjustMethod = 'BH',
  pvalueCutoff = 0.05,
  minGSSize = 10,
  maxGSSize = Inf
)


#' Parse enrichment Results
#'
#' This function changes the GeneRatio, BgRatio to a numeric values and
#' calculates the Fold Change of enrichment
#'
#' @param eresult  - `enrichResult` instance
#' @param calc_pairwise - parameter for calculating pairwise distances between
#'   GO terms
#'
#' @return
#'
#' @examples
parse_enrich_results <-  function(eresult, calc_pairwise=FALSE) {
  eresult <- mutate(  # Step 1: Converting strings to values and calculating fold change
    eresult,
    GeneRatio = map_dbl(GeneRatio, ~eval(parse(text=.x))),
    BgRatio = map_dbl(BgRatio, ~eval(parse(text=.x))),
    Fold.Change = GeneRatio/BgRatio
  ) %>% 
    filter(., p.adjust<0.05) %>%  # Step 2: filter out insignificant results
    arrange(desc(Fold.Change))
  if (calc_pairwise) {
    eresult <- pairwise_termsim(eresult)
  }
  return(eresult)
}

slice_ego <- clusterProfiler::slice  # To avoid confusion

list2df <- function(inputList) {
  # ldf <- lapply(1:length(inputList), function(i) {
  ldf <- lapply(seq_len(length(inputList)), function(i) {
    data.frame(categoryID=rep(names(inputList[i]),
                              length(inputList[[i]])),
               Gene=inputList[[i]])
  })
  
  do.call('rbind', ldf)
}



#' Obtain annotations for  genes
#'
#' @param eresult - `enrichResult` instance
#' @param showCategory - Number of top categories/GO terms to consider. Can be
#'   number or `n()`
#'
#' @return `tibble` with Genes and their annotations
#'
#' @examples
get_gene_annotation <- function(eresult, showCategory=20) {
  
  go_result <- eresult %>%
    parse_enrich_results() %>%
    slice_ego( n=c(1:showCategory)) 
  
  result <- as.data.frame(go_result) %>%
    select(ID, Description, geneID, Fold.Change) %>% 
    mutate(Genes = map(geneID, ~as_vector(strsplit(.x, split="/", fixed = T)))) %>% 
    select(-geneID) %>% 
    tidyr::unnest_longer(Genes)  
  
  return(result)
}



#' Plot a heatmap of Gene Annotation
#' 
#' This code plots a table of genes and their respective GO annotations
#'
#' @param eresult - `enrichResult` instance
#' @param showCategory - Number of top categories/GO terms to consider. Can be
#'   number or `n()` 
#'
#' @return `ggplot2` object
#'
#' @examples
go_heatplot <- function(eresult, showCategory=20) {
  
  go_result <- get_gene_annotation(eresult, showCategory = showCategory)
  
  plot <- go_result %>% 
    ggplot(aes(Genes, reorder(Description, Fold.Change))) + 
    geom_raster() + 
    theme(
      panel.grid.major = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) + 
    labs(y='GO Term')
  
  return(plot)
}

custom_treeplot <- function(
    go_result,
    top_category = 10,
    cluster.params.n = 5,
    cluster.params.label_words_n = 4,
    cluster.params.label_format = 2
    # = list(n=5, label_words_n=4, label_format=2)
) {
  go_result %>% 
    parse_enrich_results(calc_pairwise = T) %>%
    slice_ego(., n=c(1:top_category)) %>% 
    treeplot(
      cex_category=0.7,
      cluster.params = list(
        n = cluster.params.n, 
        label_words_n = cluster.params.label_words_n,
        label_format = cluster.params.label_format
      ),
      color='Fold.Change',
    ) +
    scale_color_continuous(guide=guide_colorbar(), type='viridis') +
    labs(color='Fold Change') %>% 
    return()
}

custom_emapplot <- function(
  go_result,
  top_category = 10,
  cluster.params.n = 4,
  layout = 'kk'
  # ...
){
  go_result %>% 
    parse_enrich_results(calc_pairwise = T) %>% 
    slice_ego(n=c(1:top_category)) %>% 
    emapplot(
      showCategory = top_category,
      force = 2,
      nodel_label = 'category',
      layout.params = list(layout = layout),
      edge.params = list(show=T),
      cluster.params = list(
        cluster=T,  # Cluster the points
        legend=T,   # Put it in the legend
        n=cluster.params.n # Number of clusters
      ),
      cex.params=list(category_node=0.7, category_label=0.7),
      color = 'Fold.Change'
    ) + 
    scale_fill_continuous(guide=guide_colorbar(), type='viridis') + 
    labs(fill = 'Fold Change', size='Number of Genes')
}


# gProfiler codes ----

sources <- c('KEGG', 'TF')

#' Perform GSEA for yeast for KEGG and Transcription factors
#'
#' @param x - vector of gene names 
#' @param ...  - extra params to pass to `gprofiler2::gost()`
#'
#' @return Enrichment result
#'
#' @examples
yeast_gost <- partial(
  gost,
  organism='scerevisiae',
  correction_method='fdr',
  sources = sources
)

#' Get dataframe from GOST result
#'
#' @param x - vector of gene names
#'
#' @return `tibble`
#'
#' @examples
get_gost_df <- function(x, ...){
  
  query_size <- length(x)
  universe_size <- length(all_genes)
  gost_result <- yeast_gost(x, ...)
  
  if(is.null(gost_result)) {
    print('No Significant Enrichment')
    return(NULL)
  }
  
  gost_df <- tibble(
    source = gost_result$result$source,
    term_id = gost_result$result$term_id,
    term_name = gost_result$result$term_name,
    term_size = gost_result$result$term_size,
    intersection_size = gost_result$result$intersection_size,
    adj_pval = gost_result$result$p_value
  ) %>% 
    mutate(
      `Fold Change` = (intersection_size/query_size)/(term_size/universe_size)
    )
  return(gost_df)
}

#' Function to clean up the transcription factor result from GOST
#'
#' @param gost_df_tf - GOST Dataframe containing only transcription factors, 
#' ideally after a filter operation
#'
#' @return
#' @export
#'
#' @examples
clean_tf <- function(gost_df_tf) {
  gost_df_tf %>% 
    filter(!grepl("_1", term_id)) %>% 
    # mutate(term_name = gsub("; .*", "", term_name)) %>%  # Removing the motif
    mutate(term_name = gsub("Factor: ", "", term_name)) %>% 
    return
}

plot_gost_result <- function(gost_df, data_source, top_n=10) {
  
  if(is.null(gost_df)) {
    return(NULL)
  }
  
  if(data_source=='KEGG') {
    gost_df <- filter(gost_df, `source`==data_source) 
    
    if(isEmpty(gost_df)) {
      print('No KEGG Enrichment')
      return(NULL)
    }
    
    gost_df %>% 
      arrange(desc(adj_pval)) %>% 
      slice_head(n=top_n) %>% 
      ggplot(aes(x=`Fold Change`, y =reorder(term_name, `Fold Change`), fill = -log10(adj_pval))) +
      geom_bar(stat="identity")  +
      labs(
        x = "Fold Change",
        y = "KEGG Terms",
      ) + 
      scale_fill_continuous(guide=guide_colorbar(title="-log10(Adj. p value)"),
                            type = "viridis") %>% 
      return
  } else if(data_source=='TF') {
   
   gost_df <-  filter(gost_df, `source`==data_source) 
   
   if(isEmpty(gost_df)) {
     print('No Transfac Enrichment')
     return(NULL)
   }
   
   gost_df %>% 
      clean_tf %>% 
      arrange(desc(adj_pval)) %>% 
      slice_head(n=top_n) %>% 
      ggplot(aes(x=`Fold Change`, y =reorder(term_name, `Fold Change`), fill = -log10(adj_pval))) +
      geom_bar(stat="identity")  +
      labs(
        x = "Fold Change",
        y = "Transcription Factors",
      ) + 
      scale_fill_continuous(guide=guide_colorbar(title="-log10(Adj. p value)"),
                            type = "viridis") %>% 
      return
  }
}
  