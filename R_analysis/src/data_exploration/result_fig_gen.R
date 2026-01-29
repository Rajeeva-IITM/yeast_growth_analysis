# Generate initial data figures

library(dplyr)
library(tidyr)
library(data.table)
library(arrow)
library(ComplexHeatmap)
library(viridisLite)

source("./src/utils/theme_set.R")

create_heatmap <- function(df, ...) {
  df %>% 
    cor() %>% 
    Heatmap(
      name = "Correlation",
      col = viridis(256),
      rect_gp = gpar(col="black", lwd=0.5),
      column_names_gp = gpar(fontsize=4),
      column_names_rot = 45,
      row_names_gp = gpar(fontsize=4),
      heatmap_legend_param = list(
        fontsize=6,
        title_gp = gpar(fontsize=6),
        labels_gp = gpar(fontsize=5)
      ),
      ...
    ) %>% 
    return()
}

# Peter Data ----
df <- fread("../../Original Paper/pheno_matrix.txt", sep = "\t")
df %>% 
  mutate(across(!V1, scale)) %>% 
  pivot_longer(
    cols = starts_with("Y"),
    names_to = "Condition",
    values_to = "Growth Rate"
  ) %>% 
  mutate(Condition = gsub("(YPD|YP)", "", Condition)) %>% 
  filter(!(Condition %in% c("14", "40", "42"))) %>% 
  ggplot(., aes(x=`Growth Rate`)) +
  geom_histogram(bins = 50, fill="cornflowerblue", color="black") +
  facet_wrap(~Condition, ncol=4)

save_plot("./results/result_1/Peter_hist.png", width=5, height=6, units="in", scale=2)

colnames(df) <- gsub("(YPD|YP)", "", colnames(df))

png(
  "./results/result_1/Peter_correlation.png",
  width = 7, height = 6, units = "in",
  res=300, bg="transparent"
)
df %>% 
  select(-c(V1, `14`, `40`, `42`)) %>% 
  create_heatmap(
    row_split=5, column_split=5,
    heatmap_height=unit(5, "in"), heatmap_width=unit(5.5, "in"
                                                                  )) %>% 
  draw()
dev.off()

# Bloom2013 ----

df_2013 <- read_feather(
  "../../Original Paper/Encoding Studies 2/regression_data/bloom2013_regression.feather",
  col_select = c("Strain", "Phenotype", "Condition")
) 

df_2013 %>% 
  rename(`Colony size` = Phenotype) %>% 
  mutate(Condition = stringr::str_to_upper(Condition)) %>% 
  ggplot(., aes(x=`Colony size`)) +
  geom_histogram(bins=50, fill='cornflowerblue', color="black") + 
  facet_wrap(~Condition, ncol=5)
save_plot("results/result_1/Bloom2013_hist.png", width=5, height=6, units="in", scale =2)

png(
  "./results/result_1/Bloom2013_correlation.png",
  width = 6, height = 5, units = "in",
  res=300, bg="transparent"
)
df_2013 %>% 
  mutate(Condition = stringr::str_to_upper(Condition)) %>% 
  pivot_wider(names_from = Condition, values_from = Phenotype) %>% 
  select(!Strain) %>% 
  create_heatmap(
    row_split=4, column_split=4,
    heatmap_height=unit(4, "in"), heatmap_width=unit(4.5, "in")
  )
dev.off()

# Bloom2015 ----

df_2015 <- read_feather(
  "../../Original Paper/Encoding Studies 2/regression_data/bloom2015_regression.feather",
  col_select = c("Strain", "Phenotype", "Condition")
) 

df_2015 %>% 
  rename(`Colony size` = Phenotype) %>% 
  mutate(Condition = stringr::str_to_upper(Condition)) %>% 
  ggplot(., aes(x=`Colony size`)) +
  geom_histogram(bins=50, fill='cornflowerblue', color="black") + 
  facet_wrap(~Condition, ncol=5)
save_plot("results/result_1/Bloom2015_hist.png", width=5, height=6, units="in", scale =2)

png(
  "./results/result_1/Bloom2015_correlation.png",
  width = 6, height = 5, units = "in",
  res=300, bg="transparent"
)
df_2015 %>% 
  mutate(Condition = stringr::str_to_upper(Condition)) %>% 
  group_by(Strain, Condition) %>% 
  summarise(Phenotype = mean(Phenotype)) %>% 
  pivot_wider(names_from = Condition, values_from = Phenotype) %>% 
  ungroup(Strain) %>% 
  select(!Strain) %>% 
  create_heatmap(
    row_split=3, column_split=3,
    heatmap_height=unit(4, "in"), heatmap_width=unit(4.5, "in")
  )
dev.off()

# Bloom2019 ----

df_2019 <- read_parquet(
  "./data/phenotype/Bloom2019_phenotype.parquet"
)

colnames(df_2019) <- stringr::str_to_upper(colnames(df_2019))

df_2019 %>% 
  mutate(across(!ID, scale)) %>%
  pivot_longer(
    cols = !starts_with("i"),
    names_to = "Condition",
    values_to = "Colony Size"
  ) %>% 
  # mutate(Condition = stringr::str_to_upper(Condition)) %>% 
  ggplot(., aes(x=`Colony Size`)) +
  geom_histogram(bins = 50, fill="cornflowerblue", color="black") +
  facet_wrap(~Condition, ncol=4)
save_plot("results/result_1/Bloom2019_hist.png", width=5, height=6, units="in", scale =2)


png(
  "./results/result_1/Bloom2019_correlation.png",
  width = 6, height = 5, units = "in",
  res=300, bg="transparent"
)
df_2019 %>% 
  # mutate(Condition = stringr::str_to_upper(Condition)) %>% 
  # pivot_wider(names_from = Condition, values_from = Phenotype) %>% 
  select(!ID) %>% 
  create_heatmap(
    row_split=4, column_split=4,
    heatmap_height=unit(4, "in"), heatmap_width=unit(4.5, "in")
  )
dev.off()

# Bloom2019_BYxM22 ----

df_2019_1 <- read_parquet(
  "./data/phenotype/Bloom2019_BYxM22_phenotype.parquet"
)
colnames(df_2019_1) <- stringr::str_to_upper(colnames(df_2019_1))

df_2019_1 %>% 
  mutate(across(!ID, scale)) %>%
  pivot_longer(
    cols = !starts_with("i"),
    names_to = "Condition",
    values_to = "Colony Size"
  ) %>% 
  # mutate(Condition = stringr::str_to_upper(Condition)) %>% 
  ggplot(., aes(x=`Colony Size`)) +
  geom_histogram(bins = 50, fill="cornflowerblue", color="black") +
  facet_wrap(~Condition, ncol=4)
save_plot("results/result_1/Bloom2019_BYxM22_hist.png", width=5, height=6, units="in", scale =2)

png(
  "./results/result_1/Bloom2019_BYxM22_correlation.png",
  width = 6, height = 5, units = "in",
  res=300, bg="transparent"
)
df_2019_1 %>% 
  # mutate(Condition = stringr::str_to_upper(Condition)) %>% 
  # pivot_wider(names_from = Condition, values_from = Phenotype) %>% 
  select(!ID) %>% 
  create_heatmap(
    row_split=4, column_split=4,
    heatmap_height=unit(4, "in"), heatmap_width=unit(4.5, "in")
  )
dev.off()

# Bloom2019_RMxYPS163 ----

df_2019_2 <- read_parquet(
  "./data/phenotype/Bloom2019_RMxYPS163_phenotype.parquet"
)

colnames(df_2019_2) <- stringr::str_to_upper(colnames(df_2019_2))

df_2019_2 %>% 
  mutate(across(!ID, scale)) %>%
  pivot_longer(
    cols = !starts_with("i"),
    names_to = "Condition",
    values_to = "Colony Size"
  ) %>% 
  # mutate(Condition = stringr::str_to_upper(Condition)) %>% 
  ggplot(., aes(x=`Colony Size`)) +
  geom_histogram(bins = 50, fill="cornflowerblue", color="black") +
  facet_wrap(~Condition, ncol=4)
save_plot("results/result_1/Bloom2019_RMxYPS163_hist.png", width=5, height=6, units="in", scale =2)

png(
  "./results/result_1/Bloom2019_RMxYPS163_correlation.png",
  width = 6, height = 5, units = "in",
  res=300, bg="transparent"
)
df_2019_2 %>% 
  # mutate(Condition = stringr::str_to_upper(Condition)) %>% 
  # pivot_wider(names_from = Condition, values_from = Phenotype) %>% 
  select(!ID) %>% 
  create_heatmap(
    row_split=4, column_split=4,
    heatmap_height=unit(4, "in"), heatmap_width=unit(4.5, "in")
  )
dev.off()

# UMAP results ----

umap_df <- read_csv("data/phenotype/umap_results.csv", col_types = "ddf")
umap_df %>% 
  ggplot(., aes(`UMAP 1`, `UMAP 2`, color=Dataset)) +
  geom_point() 

save_plot("results/result_1/UMAP.png")
