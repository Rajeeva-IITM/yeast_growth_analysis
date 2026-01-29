# GRN preprocessing and construction

library(readr)
library(BioNERO)
library(dplyr)
set.seed(42)
# Get extreme strains
concordant_strains <- read_csv(
  "./data/training_data/strains_0.5_ynb.csv",
  show_col_types = F,
  num_threads = 4
)

ge_matrix <- read_csv(
  "./data/expression/expressionValues_imputed_median_antilog.csv",
  show_col_types = F,
  num_threads = 4
)

ge_df <- as.data.frame(ge_matrix)
row.names(ge_df) <- ge_df$Genes
ge_df <- select(ge_df, -Genes)
ge_df <- select(ge_df, any_of(concordant_strains$Strain)) # Take only those strains common between these two

colData <- concordant_strains %>%
  mutate(Phenotype = ifelse(Phenotype==1, "High", "Low")) %>% 
  filter(Strain %in% colnames(ge_df)) %>%
  select(Phenotype) 
  
high_strains <- concordant_strains %>%
  filter(Strain %in% colnames(ge_matrix), Phenotype==1) %>% 
  select(Strain)
low_strains <- concordant_strains %>% 
  filter(Strain %in% colnames(ge_matrix), Phenotype==0) %>% 
  select(Strain)

# Preprocessing ----

# Filtering the data
ge_filtered <- remove_nonexp(
  ge_df, 
  method = 'median',
  min_exp = 10
)
ge_filtered <- filter_by_variance(ge_filtered, percentile = 90)
ge_filtered_final <- PC_correction(ge_filtered)

final_data <- SummarizedExperiment::SummarizedExperiment(list(expression = ge_filtered_final), 
                                                         colData = colData)

# Exploratory analysis ----

p1 <- plot_heatmap(final_data, type = "samplecor", show_rownames=FALSE)
p2 <- plot_heatmap(final_data, type = "expr", show_colnames=FALSE, 
                   show_rownames=FALSE)
plot_PCA(final_data)


# GRN building ----

## Loading data ----

# transcription_factors <- read_csv("./data/grn/tfs_metadata.csv", show_col_types = F)
transcription_factors <- read_csv("./data/grn/transcription_factors_Reimand_systemic.csv",show_col_types = F)
regulators <- transcription_factors$Gene.secondaryIdentifier  # Get systemic names


grn <- exp2grn(
  exp = final_data,
  regulators = regulators,
  nCores=6 # For GENIE
)

# Saving 
save(grn, file = "./data/grn/GRN_results_gimme_0.5_sigma_reimand.sav")
