library(dplyr)
library(purrr)
library(progress)
library(glmnet)
library(randomForest)

geno_df <- arrow::read_feather(
  "./data/training_data/0.5_sigma/bloom2013_clf.feather"
) %>%
  mutate(Growth = ifelse(Phenotype == 0, "Low", "High")) %>%
  mutate(Growth = factor(Growth, levels = c("Low", "High"), ordered = FALSE))

# Create a subset
small <- geno_df %>%
  filter(Condition == '4NQO') %>%
  mutate(across(starts_with('Y'), ~ map_int(.x, ~ ifelse(.x >= 1, 1, 0))))

geno_df <- geno_df %>%
  select(Strain, starts_with('Y')) %>%
  distinct() %>%
  mutate(across(starts_with('Y'), ~ map_int(.x, ~ ifelse(.x >= 1, 1, 0))))

varying_cols <- geno_df %>% # Select only varying columns
  select(starts_with('Y')) %>%
  summarise(across(everything(), sd)) %>%
  as.list %>%
  (\(x) x != 0)
varying_cols <- colnames(select(geno_df, starts_with('Y')))[varying_cols]

small <- small %>% select(all_of(varying_cols), Growth)

calculate_uniqueness <- function(geno, window_size = 6) {
  max_length <- length(geno)
  half_w <- window_size / 2
  uniqueness <- numeric(max_length)

  for (i in 1:max_length) {
    # 1. Middle cases: Full window available
    if ((i > half_w) & (i <= (max_length - half_w))) {
      neighbor_indices <- c((i - half_w):(i - 1), (i + 1):(i + half_w))

      # 2. Start cases: Not enough neighbors to the left
    } else if (i <= half_w) {
      # Use as many as possible from the right to keep window size consistent
      # We take the first 'window_size' elements and exclude 'i'
      neighbor_indices <- setdiff(1:window_size, i)

      # 3. End cases: Not enough neighbors to the right
    } else {
      # Take the last 'window_size' elements and exclude 'i'
      neighbor_indices <- setdiff((max_length - window_size + 1):max_length, i)
    }

    # Calculate average using the fixed denominator (5)
    average_genotype <- sum(geno[neighbor_indices]) / (window_size - 1)
    uniqueness[i] <- geno[i] - average_genotype
  }

  return(uniqueness)
}

# Work on small

uniqueness_matrix <- small %>%
  select(all_of(varying_cols)) %>%
  as.matrix() %>%
  array_branch(1) %>%
  map(calculate_uniqueness, .progress = T) %>%
  do.call(rbind, .)
colnames(uniqueness_matrix) <- varying_cols
modified_geno_df <- uniqueness_matrix %>% as_tibble()
modified_geno_df['Growth'] <- small$Growth

res1 <- glmnet(
  as.matrix(modified_geno_df[, varying_cols]),
  modified_geno_df$Growth,
  family = 'binomial',
  alpha = 1,
  intercept = F
)
res2 <- glmnet(
  as.matrix(small[, varying_cols]),
  small$Growth,
  family = 'binomial',
  alpha = 1,
  intercept = F
)

cv_res1 <- cv.glmnet(
  as.matrix(modified_geno_df[, varying_cols]),
  modified_geno_df$Growth,
  family = 'binomial',
  alpha = 1,
  intercept = F,
  type.measure = 'class'
)
cv_res2 <- cv.glmnet(
  as.matrix(small[, varying_cols]),
  small$Growth,
  family = 'binomial',
  alpha = 1,
  intercept = F,
  type.measure = 'class'
)


rf_mod <- randomForest(
  as.factor(modified_geno_df$Growth) ~ .,
  data = modified_geno_df[, varying_cols],
  ntree = 500
)
print(rf_mod)

rf_mod2 <- randomForest(
  as.factor(modified_geno_df$Growth) ~ .,
  data = small[, varying_cols],
  ntree = 500
)

rf_mod_full <- randomForest(
  as.factor(modified_geno_df$Growth) ~ .,
  data = modified_geno_df[, varying_cols],
  mtry = 1000,
  ntree = 500
)
print(rf_mod_full)
rf_mod_full2 <- randomForest(
  as.factor(modified_geno_df$Growth) ~ .,
  data = small[, varying_cols],
  mtry = 1000,
  ntree = 500
)
print(rf_mod_full2)

rf_mod_all <- randomForest(
  as.factor(modified_geno_df$Growth) ~ .,
  data = modified_geno_df[, varying_cols],
  mtry = 2700,
  ntree = 500
)
print(rf_mod_full)
rf_mod_all2 <- randomForest(
  as.factor(modified_geno_df$Growth) ~ .,
  data = small[, varying_cols],
  mtry = 2700,
  ntree = 500
)
print(rf_mod_full2)
