# Does chemical proximity imply phenotypic similarity?
#
# The method assumes that conditions close in the learned chemical embedding
# elicit similar growth responses across strains. This script tests that
# assumption directly, as requested by review: for every pair of conditions
# within a dataset it computes the distance between their embedding vectors and
# the correlation of their strain-level growth phenotypes, and reports the pairs
# that break the assumption - chemically near, phenotypically unrelated.
#



library(arrow)
library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(purrr)
library(ggplot2)
library(gt)


source("./src/utils/theme_set.R")

set.seed(20260910) # random number 

N_PERM <- 9999 # number of permutations
MIN_OVERLAP <- 30
OUT_DIR <- "results/result_2/embedding_cka"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Embedding


emb_tbl <- read_parquet("data/latent_conditions.parquet", as_data_frame = FALSE)
chr_col <- function(tb, nm) as.vector(tb[[nm]]$cast(arrow::string()))

emb_meta <- tibble(
  Key = chr_col(emb_tbl, "Key"),
  Name = chr_col(emb_tbl, "Name"),
  Class = chr_col(emb_tbl, "Class"),
  RawLabels = chr_col(emb_tbl, "RawLabels"),
  SourceLabel = chr_col(emb_tbl, "SourceLabel"),
  VectorID = as.integer(chr_col(emb_tbl, "VectorID")),
  Bloom2013 = as.integer(chr_col(emb_tbl, "Bloom2013")),
  Bloom2015 = as.integer(chr_col(emb_tbl, "Bloom2015")),
  Bloom2019_BYxRM = as.integer(chr_col(emb_tbl, "Bloom2019_BYxRM"))
)

latent_cols <- grep("^latent_", names(emb_tbl), value = TRUE)
latent_cols <- latent_cols[order(as.integer(sub("^latent_", "", latent_cols)))]
LATENT <- do.call(cbind, lapply(latent_cols, function(cn) as.vector(emb_tbl[[cn]])))
dimnames(LATENT) <- list(emb_meta$Key, latent_cols)

# "etoh" is the one condition label absent from RawLabels; "sds"/"SDS" and
# "copper"/"CuSO4" resolve through RawLabels plus tolower().
EXTRA_ALIASES <- c(etoh = "ethanol")
alias_lookup <- emb_meta %>%
  transmute(Key, label = strsplit(paste(Key, SourceLabel, RawLabels, sep = ","), ",")) %>%
  unnest(label) %>%
  mutate(label = tolower(trimws(label))) %>%
  distinct(label, .keep_all = TRUE) %>%
  { setNames(.$Key, .$label) }
alias_lookup[names(EXTRA_ALIASES)] <- unname(EXTRA_ALIASES)


# Phenotype


DATASETS <- tribble(
  ~dataset,          ~file,           ~flag,
  "Bloom2013",       "bloom2013_clf", "Bloom2013",
  "Bloom2015",       "bloom2015_clf", "Bloom2015",
  "Bloom2019_BYxRM", "bloom2019_clf", "Bloom2019_BYxRM"
)

load_matrix <- function(file) {
  raw <- read_feather(
    sprintf("data/training_data/0.5_sigma/%s.feather", file),
    col_select = c("Strain", "Condition", "Phenotype")
  ) %>%
    mutate(
      Strain = as.character(Strain),
      Condition = as.character(Condition),
      Phenotype = as.numeric(Phenotype)
    )

  # Drop duplicates
  kept <- raw %>%
    group_by(Strain, Condition) %>%
    summarise(Phenotype = mean(Phenotype), .groups = "drop") %>%
    filter(Phenotype == 0 | Phenotype == 1) %>%
    mutate(Key = unname(alias_lookup[tolower(Condition)]))
  stopifnot(!any(is.na(kept$Key)))

  wide <- kept %>%
    select(Key, Strain, Phenotype) %>%
    pivot_wider(names_from = Strain, values_from = Phenotype)
  M <- as.matrix(wide[, -1, drop = FALSE])
  rownames(M) <- wide$Key
  M
}


# Pair table: chemical distance against phenotype correlation


pair_table <- function(dataset, flag, P) {
  keys <- rownames(P)
  stopifnot(setequal(keys, emb_meta$Key[emb_meta[[flag]] == 1]))

  X <- LATENT[keys, , drop = FALSE]
  chem_d <- as.matrix(dist(X))
  # Pearson correlation between condition vectors. The embedding is column-
  # centered over the conditions of this dataset first: the raw vectors are
  # near-collinear so correlating them uncentered returns 0.988-1.000 for every
  # pair and carries no usable contrast. Centering removes that common component
  # and leaves what actually distinguishes the conditions
  Xc <- sweep(X, 2, colMeans(X), "-")
  emb_r <- stats::cor(t(Xc))

  phen_r <- stats::cor(t(P), use = "pairwise.complete.obs")
  observed <- !is.na(P)
  overlap <- tcrossprod(observed * 1)

  cls <- setNames(emb_meta$Class[match(keys, emb_meta$Key)], keys)
  nm <- setNames(emb_meta$Name[match(keys, emb_meta$Key)], keys)

  idx <- which(upper.tri(diag(length(keys))), arr.ind = TRUE)
  tibble(
    dataset = dataset,
    condition_a = keys[idx[, 1]],
    condition_b = keys[idx[, 2]],
    name_a = nm[keys[idx[, 1]]],
    name_b = nm[keys[idx[, 2]]],
    class_a = cls[keys[idx[, 1]]],
    class_b = cls[keys[idx[, 2]]],
    chem_distance = chem_d[idx],
    embedding_correlation = emb_r[idx],
    phenotype_r = phen_r[idx],
    shared_strains = overlap[idx]
  ) %>%
    filter(shared_strains >= MIN_OVERLAP, is.finite(phenotype_r)) %>%
    mutate(
      phenotype_r2 = phenotype_r^2,
      same_class = ifelse(class_a == class_b, "Same chemical class", "Different class")
    )
}

message("Building condition pairs")
pairs <- pmap(
  list(DATASETS$dataset, DATASETS$flag, map(DATASETS$file, load_matrix)),
  pair_table
) %>% list_rbind()

for (ds in unique(pairs$dataset)) {
  d <- pairs %>% filter(dataset == ds)
  message(sprintf("  %s: %d pairs, shared strains %d-%d, embedding correlation %.3f-%.3f",
                  ds, nrow(d), min(d$shared_strains), max(d$shared_strains),
                  min(d$embedding_correlation), max(d$embedding_correlation)))
}

write_csv(pairs, file.path(OUT_DIR, "chemical_distance_vs_phenotype-pairs.csv"))


# Is the assumption supported? Mantel-style permutation test


# Pairs are not independent (each condition appears in many), so significance
# comes from permuting condition labels rather than from a naive correlation
# test on the pair list.
mantel_test <- function(d, distance_col) {
  keys <- union(d$condition_a, d$condition_b)
  n <- length(keys)
  D <- matrix(NA_real_, n, n, dimnames = list(keys, keys))
  R <- D
  ia <- match(d$condition_a, keys)
  ib <- match(d$condition_b, keys)
  D[cbind(ia, ib)] <- d[[distance_col]]
  D[cbind(ib, ia)] <- d[[distance_col]]
  R[cbind(ia, ib)] <- d$phenotype_r2
  R[cbind(ib, ia)] <- d$phenotype_r2

  stat <- function(perm) {
    Dp <- D[perm, perm]
    ok <- upper.tri(D) & is.finite(Dp) & is.finite(R)
    suppressWarnings(stats::cor(Dp[ok], R[ok], method = "spearman"))
  }
  obs <- stat(seq_len(n))
  null <- vapply(seq_len(N_PERM), function(i) stat(sample.int(n)), numeric(1))
  null <- null[is.finite(null)]
  tibble(
    n_conditions = n, n_pairs = nrow(d), distance = distance_col,
    rho = obs,
    perm_p = (1 + sum(abs(null) >= abs(obs))) / (1 + length(null))
  )
}

message("Mantel permutation tests")
mantel <- map(unique(pairs$dataset), function(ds) {
  d <- pairs %>% filter(dataset == ds)
  bind_rows(
    mantel_test(d, "chem_distance"),
    mantel_test(d, "embedding_correlation")
  ) %>% mutate(dataset = ds, .before = 1)
}) %>% list_rbind()

print(as.data.frame(mantel), digits = 3)
write_csv(mantel, file.path(OUT_DIR, "chemical_distance_vs_phenotype-mantel.csv"))


# Counter-examples: chemically near, phenotypically unrelated


# The embedding assigns identical vectors to two isomer pairs, so those are the
# strongest possible test of the assumption: zero chemical distance by
# construction.
identical_pairs <- pairs %>%
  filter(chem_distance < 1e-9) %>%
  arrange(phenotype_r2)

counter <- pairs %>%
  group_by(dataset) %>%
  filter(embedding_correlation >= quantile(embedding_correlation, 0.90)) %>%
  slice_min(phenotype_r2, n = 6) %>%
  ungroup() %>%
  arrange(dataset, phenotype_r2)

message("\nCondition pairs with identical embedding vectors:")
if (nrow(identical_pairs)) {
  print(as.data.frame(identical_pairs %>%
    select(dataset, condition_a, condition_b, embedding_correlation, phenotype_r, phenotype_r2,
           shared_strains)), digits = 3)
}

write_csv(counter, file.path(OUT_DIR, "chemical_distance_vs_phenotype-counterexamples.csv"))

# Figure


plot_df <- pairs %>%
  mutate(dataset = sub("Bloom2019_BYxRM", "Bloom2019", dataset))

label_df <- counter %>%
  mutate(
    dataset = sub("Bloom2019_BYxRM", "Bloom2019", dataset),
    pair_label = paste(condition_a, condition_b, sep = " / ")
  ) %>%
  group_by(dataset) %>%
  slice_min(phenotype_r2, n = 3) %>%
  ungroup()

scatter <- ggplot(plot_df, aes(x = embedding_correlation, y = phenotype_r)) +
  geom_point(colour = cbPalette[6], alpha = 0.6, size = 1.4) +
  geom_smooth(method = "lm", formula = y ~ x, colour = "black", linewidth = 0.6, se = TRUE) +
  ggrepel::geom_text_repel(
    data = label_df, aes(label = pair_label),
    size = 2.6, min.segment.length = 0, max.overlaps = 20
  ) +
  facet_wrap(~dataset) +
  labs(
    x = "Correlation between condition embeddings, r",
    y = expression(paste("Phenotype correlation across strains, ", r)),
    title = "Chemical similarity against phenotypic similarity, all condition pairs"
  )

save_plot(
  file.path(OUT_DIR, "chemical_distance_vs_phenotype.svg"),
  scatter, width = 11, height = 4.4
)

# Counter-example table

counter_gt <- counter %>%
  mutate(
    dataset = sub("Bloom2019_BYxRM", "Bloom2019", dataset),
    pair = paste(name_a, "/", name_b),
    classes = ifelse(class_a == class_b, class_a, paste(class_a, "/", class_b))
  ) %>%
  select(dataset, pair, classes, embedding_correlation, phenotype_r, phenotype_r2, shared_strains) %>%
  gt(rowname_col = "pair", groupname_col = "dataset") %>%
  fmt_number(columns = c(embedding_correlation, phenotype_r, phenotype_r2), decimals = 3) %>%
  fmt_integer(columns = shared_strains) %>%
  cols_label(
    classes = "Chemical class", embedding_correlation = "Embedding correlation",
    phenotype_r = "r", phenotype_r2 = "RSQUARED", shared_strains = "Strains"
  ) %>%
  tab_stubhead("Condition pair") %>%
  cols_align(align = "center", columns = -pair) %>%
  tab_footnote(
    "Pearson correlation of binary growth phenotypes across the strains measured in both conditions.",
    locations = cells_column_labels(columns = phenotype_r)
  ) %>%
  opt_footnote_marks(marks = "letters") %>%
  tab_options(
    table.border.bottom.color = "black",
    table.border.bottom.width = px(1.5),
    column_labels.font.weight = "bold",
    column_labels.border.top.width = px(2),
    column_labels.border.bottom.width = px(2),
    column_labels.border.top.color = "black",
    column_labels.border.bottom.color = "black",
    row_group.font.weight = "bold",
    stub.font.weight = "bold",
    data_row.padding = px(4),
    table.font.size = px(12)
  ) %>%
  opt_table_font(font = "CMU Serif")

as_tex <- function(gt_tbl, caption) {
  lines <- strsplit(as.character(as_latex(gt_tbl)), "\n", fixed = TRUE)[[1]]
  at <- grep("\\begin{table}", lines, fixed = TRUE)[1]
  stopifnot(!is.na(at))
  lines <- append(lines, paste0("\\caption{", caption, "}"), after = at)
  x <- paste(lines, collapse = "\n")
  x <- gsub("±", "$\\pm$", x, fixed = TRUE)
  x <- gsub("Δ", "$\\Delta$", x, fixed = TRUE)
  gsub("RSQUARED", "$r^{2}$", x, fixed = TRUE)
}

writeLines(
  as_tex(counter_gt, paste(
    "Condition pairs that are close in chemical embedding space but elicit",
    "uncorrelated growth responses. For each dataset the six pairs with the",
    "lowest phenotype correlation among the decile of pairs with the highest",
    "embedding correlation",
    "are shown. Phenotype correlation is computed across the strains measured in",
    "both conditions at 0.5 sigma."
  )),
  file.path(OUT_DIR, "chemical_distance_vs_phenotype-counterexamples.tex")
)

message("\nDone. Outputs in ", OUT_DIR)
