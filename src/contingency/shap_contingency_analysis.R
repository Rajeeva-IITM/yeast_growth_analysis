# Code to compare SHAP and contingency table obtained for YNB media

library(dplyr)
library(data.table)
library(readr)
library(purrr)
library(ggplot2)
library(ggpubr)
library(ggforce)

source("./src/utils/theme_set.R")

qtl_info <- read_csv("./data/qtl/YNB_qtls.csv") %>%
  mutate(gene_list = strsplit(`Genes Under Confidence Interval (Standard Name)`, "|", fixed = TRUE))

qtl_annotate <- function(gene) {
  presence <- map(qtl_info$gene_list, ~ (gene %in% .x)) %>%
    unlist() %>%
    qtl_info$qtl_id[.] %>%
    ifelse(length(.) == 0, "Non QTL", .)
  return(presence)
}


contingency_data <- read_csv(
  "./results/result_contingency/contingency_bloom2013_ynb_0.5sigma.csv" # Variable
) %>%
  mutate(adj_pval = p.adjust(pval, "fdr")) %>%
  mutate(Significance = ifelse(adj_pval < 0.05, "Significant", "Not Significant"))

shap_df <- read_csv("./data/shap/sigmas/shap_classification_0.5/shap_Bloom2013_Boosting.csv") %>%
  filter(grepl("^Y", Feature)) %>%
  rename(Gene = Feature)

# Relation between shap value and contingency_data

contingency_data <- left_join(contingency_data, shap_df) %>%
  mutate(qtl = map_chr(Gene, qtl_annotate)) %>%
  mutate(qtl = as.factor(gsub('_', " ", qtl, fixed = TRUE)))

annotation_df <- contingency_data %>% filter(!log10(odds_ratio) %between% c(-0.4, 0.4))

fig <- contingency_data %>%
  ggplot(aes(x = log10(Value), y = -log10(odds_ratio), )) +
  geom_point(data = annotation_df,
             inherit.aes = T,
             size = 2) +
  geom_point(mapping = aes(color = Significance)) +
  geom_abline(slope = 0,
              intercept = 0,
              color = "blue") +
  labs(x = "log10(SHAP value)", y = "-log10(Odds Ratio)", title = "Comparison of SHAP analysis and Fisher Test")  +
  ggrepel::geom_text_repel(
    mapping = aes(label = Gene),
    data = annotation_df,
    color = "black",
    nudge_x = -0.0,
    nudge_y = -0.05,
    size = 4
  )

save_plot(
  "./results/result_contingency/bloom2013_ynb_0.5sigma_shapcomparison.svg",
  height = 8,
  width=12,
  units='in'
) # Variable

# QTLs

fig2 <- ggplot(contingency_data, aes(log10(Value), -log10(odds_ratio), )) +
  ggforce::geom_mark_ellipse(aes(fill = qtl, filter = qtl != "Non QTL")) +
  geom_point(data = annotation_df,
             inherit.aes = T,
             size = 2) +
  geom_point(mapping = aes(color = Significance)) +
  geom_abline(slope = 0,
              intercept = 0,
              color = "blue") +
  labs(x = "log10(mean(|SHAP value|)",
       y = "-log10(Odds Ratio)",
       title = "Comparison of SHAP analysis and Fisher Test") +
  geom_text(
    mapping = aes(label = Gene),
    data = annotation_df,
    color = "black",
    nudge_x = -0.0,
    nudge_y = -0.1,
    size = 4
  ) +
  facet_wrap( ~ qtl)

save_plot(
  "./results/result_contingency/bloom2013_ynb_0.5sigma_shapcomparison_qtl.svg",
  height=12,
  width=16,
  units='in'
) # Variable

# One QTL and One Non-QTL

subset_contingency <- contingency_data %>%
  filter(qtl %in% c("qtl chr 14-3", "Non QTL")) # Variable

annotation_df <- filter(subset_contingency, !(log10(odds_ratio) %between% c(-0.4, 0.4)))

fig3 <- ggplot(subset_contingency, aes(x = log10(Value), y = -log10(odds_ratio))) +
  geom_point(
    data = annotation_df,
    color = "black",
    inherit.aes = TRUE,
    size = 2
  ) +
  ggrepel::geom_text_repel(
    mapping = aes(label = Gene),
    data = annotation_df,
    color = "black",
    nudge_x = -0.0,
    nudge_y = -0.05,
    size = 4
  ) +
  geom_point(aes(color = Significance)) + geom_hline(yintercept = 0) +
  labs(# title = "Comparing SHAP identified genes and their relationship to growth",
    x = "log10(mean(|SHAP value|))", y = "-log10(Odds Ratio)") +
  guides(color = guide_legend(title = "Fisher Test\nSignificance")) +
  facet_wrap( ~ qtl, labeller = as_labeller(c(
    `qtl chr 14-3` = "QTL - Chromosome 14", # `qtl_chr_15_2` = "QTL - Chromosome 15", # Variable
    `Non QTL` = "Non QTL"
  )))

save_plot(
  "./results/result_contingency/bloom2013_ynb_0.5sigma_shapcomparison_qtl_2.svg", # Variable
  height = 8,
  width = 12,
  units='in'
)

# Global Bloom2013 ----

contingency_data_full <- read_csv(
  "./results/result_contingency/contingency_bloom2013_0.5sigma_full.csv" # Variable
) %>%
  mutate(adj_pval = p.adjust(pval, "fdr")) %>%
  mutate(Significance = ifelse(adj_pval < 0.05, "Significant", "Not Significant"))

contingency_data_full <- left_join(contingency_data_full, shap_df) %>%
  mutate(qtl = map_chr(Gene, qtl_annotate)) %>%
  mutate(qtl = as.factor(gsub('_', " ", qtl, fixed = TRUE)))

annotation_df_full <- contingency_data_full %>% filter(!log10(odds_ratio) %between% c(-0.2, 0.2))

# fig 1
contingency_data_full %>%
  ggplot(aes(x = log10(Value), y = -log10(odds_ratio), )) +
  geom_point(data = annotation_df_full,
             inherit.aes = T,
             size = 2) +
  geom_point(mapping = aes(color = Significance)) +
  geom_abline(slope = 0,
              intercept = 0,
              color = "blue") +
  labs(x = "log10(SHAP value)", y = "-log10(Odds Ratio)", title = "Comparison of SHAP analysis and Fisher Test")  +
  ggrepel::geom_text_repel(
    mapping = aes(label = Gene),
    data = annotation_df_full,
    color = "black",
    nudge_x = -0.0,
    nudge_y = -0.05,
    size = 4
  )

# fig 2
ggplot(contingency_data_full, aes(log10(Value), -log10(odds_ratio), )) +
  ggforce::geom_mark_ellipse(aes(fill = qtl, filter = qtl != "Non QTL")) +
  geom_point(data = annotation_df_full,
             inherit.aes = T,
             size = 2) +
  geom_point(mapping = aes(color = Significance)) +
  geom_abline(slope = 0,
              intercept = 0,
              color = "blue") +
  labs(x = "log10(mean(|SHAP value|)",
       y = "-log10(Odds Ratio)",
       title = "Comparison of SHAP analysis and Fisher Test") +
  geom_text(
    mapping = aes(label = Gene),
    data = annotation_df_full,
    color = "black",
    nudge_x = -0.0,
    nudge_y = -0.1,
    size = 4
  ) +
  facet_wrap( ~ qtl)

# fig 3
subset_contingency <- contingency_data_full %>%
  filter(qtl %in% c("qtl chr 14-3", "Non QTL")) # Variable

annotation_df_full <- filter(subset_contingency, !(log10(odds_ratio) %between% c(-0.4, 0.4)))

ggplot(subset_contingency, aes(x = log10(Value), y = -log10(odds_ratio))) +
  geom_point(
    data = annotation_df_full,
    color = "black",
    inherit.aes = TRUE,
    size = 2
  ) +
  ggrepel::geom_text_repel(
    mapping = aes(label = Gene),
    data = annotation_df_full,
    color = "black",
    nudge_x = -0.0,
    nudge_y = -0.05,
    size = 4
  ) +
  geom_point(aes(color = Significance)) + geom_hline(yintercept = 0) +
  labs(# title = "Comparing SHAP identified genes and their relationship to growth",
    x = "log10(mean(|SHAP value|))", y = "-log10(Odds Ratio)") +
  guides(color = guide_legend(title = "Fisher Test\nSignificance")) +
  facet_wrap( ~ qtl, labeller = as_labeller(c(
    `qtl chr 14-3` = "QTL - Chromosome 14", # `qtl_chr_15_2` = "QTL - Chromosome 15", # Variable
    `Non QTL` = "Non QTL"
  )))
