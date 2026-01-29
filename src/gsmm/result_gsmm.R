
source("./src/utils/flux_enrichment_analysis.R")
source("./src/utils/theme_set.R")


fea_barplot <- function(fea_df) {
  p <- fea_df %>% 
    filter(Adj_pval <0.05) %>% 
  ggplot(
    .,
    aes(x=Fold_change, y=reorder(Group, Fold_change), fill=-log10(Adj_pval))
  ) + 
    geom_bar(stat='identity') +
    labs(x="Fold Change", y="Reaction Subsystem", fill='-log10(Adj. p-value)') +
    scale_fill_continuous(type = "viridis")
    return(p)
}

all_reactions <- read_csv("./data/gsmm/reaction_data.csv", show_col_types = F)

fea <- partial(
  FEA, total_reactions=all_reactions$reaction, groups=all_reactions$subsystem
)

# High unique reactions

s1 <- read_csv("./data/gsmm/Sampling-0.5_sigma/Unique_High_Reactions.csv", show_col_types = F)
res1 <- fea(s1$Reaction)
res1 %>%
  fea_barplot() %>%
  save_plot(
    "./results/result_5_gsmm/sigma_0.5/FEA_Unique_High.svg",
    plot=., width=8, height=6, units='in'
  )

# Low unique reactions 

res2 <- read_csv(
  "./data/gsmm/Sampling-0.5_sigma/Unique_Low_Reactions.csv", show_col_types = F
) %>% 
  pluck('Reaction') %>% 
  fea()

res2 %>% 
  fea_barplot() %>%
  save_plot(
    "./results/result_5_gsmm/sigma_0.5/FEA_Unique_Low.svg",
    plot=., width=8, height=6, units='in'
  )

# High FC reaction

s3 <- read_csv(
  "./data/gsmm/Sampling-0.5_sigma/high_FC_reactions.csv", show_col_types = F
) 

res3 <- s3 %>% 
  pluck('reaction') %>% 
  fea()

res3 %>% 
  fea_barplot() %>%
  save_plot(
    "./results/result_5_gsmm/sigma_0.5/FEA_High_FC.svg",
    plot=., width=8, height=6, units='in'
  )

# Low FC reaction

s4 <- read_csv(
  "./data/gsmm/Sampling-0.5_sigma/low_FC_reactions.csv", show_col_types = F
)
res4 <- fea(s4$reaction)

res4 %>% 
  fea_barplot() %>%
  save_plot(
    "./results/result_5_gsmm/sigma_0.5/FEA_Low_FC.svg",
    plot=., width=8, height=6, units='in'
  )

# Active low unique reactions

s5 <- read_csv(
  "./data/gsmm/Sampling-0.5_sigma/active_low_unique_reactions.csv", show_col_types = F
)
res5 <- fea(s5$Reaction)

res5 %>% 
  fea_barplot() %>%
  save_plot(
    "./results/result_5_gsmm/sigma_0.5/FEA_active_low_unique.svg",
    plot=., width=8, height=6, units='in'
  )

# Active High unique reactions

s6 <- read_csv(
  "./data/gsmm/Sampling-0.5_sigma/active_high_unique_reactions.csv", show_col_types = F
)
res6 <- fea(s6$Reaction)

res6 %>%
  fea_barplot() %>%
  save_plot(
    "./results/result_5_gsmm/sigma_0.5/FEA_active_Unique_High.svg",
    plot=., width=8, height=6, units='in'
  )

# All active reactions high

res7 <- union(s6$Reaction, s3$reaction) %>% fea()

res7 %>%
  fea_barplot() %>%
  save_plot(
    "./results/result_5_gsmm/sigma_0.5/FEA_active_all_high.svg",
    plot=., width=8, height=6, units='in'
  )

# All active reactions low

res8 <- union(s5$Reaction, s4$reaction) %>% fea()

res8 %>%
  fea_barplot() %>%
  save_plot(
    "./results/result_5_gsmm/sigma_0.5/FEA_active_all_low.svg",
    plot=., width=8, height=6, units='in'
  )

# pFBA high genes
