# Centered Kernel Alignment (CKA) between the chemical-condition embedding and
# the strain-level growth-phenotype geometry.
#
#
# Why CKA and not CCA: n is 13-39 conditions against 256 embedding dimensions
# and ~1,000-4,300 strains. Canonical correlation returns exactly 1.000 in that
# regime and means nothing. CKA compares the two Gram matrices directly, needs
# no dimension reduction, is invariant to orthogonal transforms and isotropic
# scaling, and has a clean permutation null.
#

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(purrr)
  library(ggplot2)
  library(gt)
})

source("./src/utils/theme_set.R")

set.seed(20260910)

N_PERM <- 9999
N_BOOT <- 2000
N_RANDOM <- 200
MIN_OVERLAP <- 30 # minimum strains shared by a pair of conditions
OUT_DIR <- "results/result_2/embedding_cka"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# HSIC / CKA

# Unbiased HSIC (Song et al. 2012). 
hsic_unbiased <- function(K, L) {
  n <- nrow(K)
  stopifnot(nrow(L) == n, ncol(K) == n, ncol(L) == n, n > 3)
  diag(K) <- 0
  diag(L) <- 0
  term1 <- sum(K * L)
  term2 <- sum(K) * sum(L) / ((n - 1) * (n - 2))
  term3 <- 2 * sum(colSums(K) * colSums(L)) / (n - 2)
  (term1 + term2 - term3) / (n * (n - 3))
}

hsic_biased <- function(K, L) {
  n <- nrow(K)
  H <- diag(n) - 1 / n
  sum((H %*% K %*% H) * L) / (n - 1)^2
}

cka <- function(K, L, unbiased = TRUE) {
  f <- if (unbiased) hsic_unbiased else hsic_biased
  denom <- f(K, K) * f(L, L)
  if (!is.finite(denom) || denom <= 0) {
    return(NA_real_)
  }
  f(K, L) / sqrt(denom)
}

# --- kernels -----------------------------------------------------------------

kernel_linear <- function(X) {
  Xc <- sweep(X, 2, colMeans(X), "-")
  tcrossprod(Xc)
}

kernel_cor <- function(X) stats::cor(t(X))

kernel_rbf_from_dist <- function(D) {
  bw <- stats::median(D[upper.tri(D)])
  if (!is.finite(bw) || bw <= 0) bw <- 1
  exp(-(D^2) / (2 * bw^2))
}

cor_to_dist <- function(C) sqrt(pmax(2 * (1 - C), 0))

# --- diagnostics -------------------------------------------------------------

# A pairwise-complete similarity estimate is not guaranteed positive
# semi-definite. 
psd_report <- function(L) {
  ev <- eigen((L + t(L)) / 2, symmetric = TRUE, only.values = TRUE)$values
  list(neg_mass = sum(pmax(-ev, 0)) / sum(abs(ev)), min_eig = min(ev))
}

# --- inference ---------------------------------------------------------------

perm_test_cka <- function(K, L, B = N_PERM, unbiased = TRUE) {
  obs <- cka(K, L, unbiased)
  n <- nrow(K)
  null <- vapply(seq_len(B), function(i) {
    p <- sample.int(n)
    cka(K, L[p, p, drop = FALSE], unbiased)
  }, numeric(1))
  null <- null[is.finite(null)]
  list(
    cka = obs,
    p = (1 + sum(null >= obs)) / (1 + length(null)),
    null_mean = mean(null),
    null_q95 = unname(stats::quantile(null, 0.95)),
    null = null
  )
}

# Resampling conditions with replacement duplicates rows of both kernels, which
# the unbiased estimator was not designed for. Reported as indicative only.
boot_ci_cka <- function(K, L, B = N_BOOT, unbiased = TRUE) {
  n <- nrow(K)
  vals <- vapply(seq_len(B), function(i) {
    idx <- sample.int(n, replace = TRUE)
    if (length(unique(idx)) < 6) {
      return(NA_real_)
    }
    cka(K[idx, idx, drop = FALSE], L[idx, idx, drop = FALSE], unbiased)
  }, numeric(1))
  vals <- vals[is.finite(vals)]
  unname(stats::quantile(vals, c(0.025, 0.975)))
}

# Condition embedding

emb_tbl <- read_parquet("data/latent_conditions.parquet", as_data_frame = FALSE)

# Class, and several other columns, are dictionary-encoded; arrow refuses to
# convert those to R without an explicit cast.
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

# "etoh" (Bloom2019) is the one condition label absent from RawLabels. "sds" vs
# "SDS" and "copper" vs "CuSO4" both resolve through RawLabels + tolower().
EXTRA_ALIASES <- c(etoh = "ethanol")
stopifnot(all(EXTRA_ALIASES %in% emb_meta$Key))

alias_lookup <- emb_meta %>%
  transmute(Key, label = strsplit(paste(Key, SourceLabel, RawLabels, sep = ","), ",")) %>%
  unnest(label) %>%
  mutate(label = tolower(trimws(label))) %>%
  distinct(label, .keep_all = TRUE) %>%
  { setNames(.$Key, .$label) }
alias_lookup[names(EXTRA_ALIASES)] <- unname(EXTRA_ALIASES)

# Conditions whose embedding vector is not unique. lactose/maltose and
# galactose/mannose are isomer pairs (C12H22O11, C6H12O6) that the embedding
# maps to identical vectors, so no model can tell them apart from chemistry.
dup_vec <- emb_meta$VectorID[duplicated(emb_meta$VectorID)]
ISOMER_DROP <- emb_meta %>%
  filter(VectorID %in% dup_vec) %>%
  group_by(VectorID) %>%
  slice(2) %>%
  pull(Key)

# Phenotype matrices

DATASETS <- tribble(
  ~dataset,          ~file,           ~flag,
  "Bloom2013",       "bloom2013_clf", "Bloom2013",
  "Bloom2015",       "bloom2015_clf", "Bloom2015",
  "Bloom2019_BYxRM", "bloom2019_clf", "Bloom2019_BYxRM"
)

# Drop duplicates
load_phenotype <- function(file) {
  raw <- read_feather(
    sprintf("data/training_data/0.5_sigma/%s.feather", file),
    col_select = c("Strain", "Condition", "Phenotype")
  ) %>%
    mutate(
      Strain = as.character(Strain),
      Condition = as.character(Condition),
      Phenotype = as.numeric(Phenotype)
    )

  collapsed <- raw %>%
    group_by(Strain, Condition) %>%
    summarise(Phenotype = mean(Phenotype), reps = n(), .groups = "drop")

  conflicting <- collapsed %>% filter(Phenotype != 0, Phenotype != 1)
  kept <- collapsed %>% filter(Phenotype == 0 | Phenotype == 1)

  message(sprintf(
    "  %s: %d rows -> %d (Strain, Condition) pairs; %d replicated, %d conflicting dropped, %d kept",
    file, nrow(raw), nrow(collapsed), sum(collapsed$reps > 1), nrow(conflicting), nrow(kept)
  ))
  stopifnot(nrow(conflicting) + nrow(kept) == nrow(collapsed))
  kept
}

phenotype_matrix <- function(d) {
  d <- d %>% mutate(Key = unname(alias_lookup[tolower(Condition)]))
  unmatched <- unique(d$Condition[is.na(d$Key)])
  if (length(unmatched)) {
    stop("Conditions with no embedding row: ", paste(unmatched, collapse = ", "))
  }
  wide <- d %>%
    select(Key, Strain, Phenotype) %>%
    pivot_wider(names_from = Strain, values_from = Phenotype)
  M <- as.matrix(wide[, -1, drop = FALSE])
  rownames(M) <- wide$Key
  M
}

# Condition x condition kernel from a conditions x strains 0/1 matrix with
# missing cells. 
phenotype_kernel <- function(P, min_overlap = MIN_OVERLAP) {
  Pc <- sweep(P, 2, colMeans(P, na.rm = TRUE), "-")
  observed <- !is.na(Pc)
  Pc[!observed] <- 0
  overlap <- tcrossprod(observed * 1)
  L <- tcrossprod(Pc) / pmax(overlap, 1)
  L[overlap < min_overlap] <- NA_real_
  dimnames(L) <- list(rownames(P), rownames(P))
  list(L = L, overlap = overlap)
}

# Drop the condition responsible for the most unusable pairs until the kernel is
# complete. 
drop_sparse_conditions <- function(L, overlap) {
  dropped <- character(0)
  while (anyNA(L)) {
    worst <- which.max(rowSums(is.na(L)))
    dropped <- c(dropped, rownames(L)[worst])
    L <- L[-worst, -worst, drop = FALSE]
    overlap <- overlap[-worst, -worst, drop = FALSE]
  }
  list(L = L, overlap = overlap, dropped = dropped)
}

message("Loading phenotype data")
pheno <- DATASETS %>%
  mutate(
    data = map(file, load_phenotype),
    P = map(data, phenotype_matrix)
  )

# Assemble per-dataset kernels

build_dataset <- function(dataset, flag, P) {
  keys_flagged <- emb_meta$Key[emb_meta[[flag]] == 1]
  keys_pheno <- rownames(P)
  stopifnot(setequal(keys_flagged, keys_pheno))

  pk <- phenotype_kernel(P)
  cleaned <- drop_sparse_conditions(pk$L, pk$overlap)
  keys <- rownames(cleaned$L)
  if (length(cleaned$dropped)) {
    message(sprintf(
      "  %s: dropped %s (fewer than %d shared strains with some other condition)",
      dataset, paste(cleaned$dropped, collapse = ", "), MIN_OVERLAP
    ))
  }

  P <- P[keys, , drop = FALSE]
  X <- LATENT[keys, , drop = FALSE]
  L_cor <- stats::cor(t(P), use = "pairwise.complete.obs")[keys, keys, drop = FALSE]
  ov <- cleaned$overlap
  cls <- emb_meta$Class[match(keys, emb_meta$Key)]

  list(
    dataset = dataset,
    keys = keys,
    n = length(keys),
    Class = cls,
    X = X,
    P = P,
    L_lin = cleaned$L,
    L_cor = L_cor,
    K_lin = kernel_linear(X),
    K_cor = kernel_cor(X),
    K_rbf = kernel_rbf_from_dist(as.matrix(stats::dist(X))),
    L_rbf = kernel_rbf_from_dist(cor_to_dist(L_cor)),
    min_overlap = min(ov[upper.tri(ov)]),
    psd = psd_report(cleaned$L)
  )
}

message("Building kernels")
DS <- pmap(list(pheno$dataset, pheno$flag, pheno$P), build_dataset)
names(DS) <- pheno$dataset

for (d in DS) {
  message(sprintf(
    "  %s: n = %d conditions, min shared strains = %d, negative eigenvalue mass = %.4f",
    d$dataset, d$n, d$min_overlap, d$psd$neg_mass
  ))
}

# Estimator sanity checks - run before any result is produced

message("Estimator sanity checks")
local({
  X <- DS[[1]]$X
  K <- kernel_linear(X)
  stopifnot(abs(cka(K, K) - 1) < 1e-10)
  stopifnot(abs(cka(K, K, unbiased = FALSE) - 1) < 1e-10)

  # CKA must be blind to an orthogonal rotation of the feature space.
  set.seed(1)
  Q <- qr.Q(qr(matrix(rnorm(ncol(X)^2), ncol(X))))
  stopifnot(abs(cka(K, kernel_linear(X %*% Q)) - 1) < 1e-8)
  # ... and to isotropic rescaling.
  stopifnot(abs(cka(K, kernel_linear(X * 7.3)) - 1) < 1e-8)
  message("  CKA(K,K) = 1, rotation- and scale-invariant: OK")
})

# Null calibration: if the alignment is destroyed before the test begins, the
# permutation p must be roughly uniform. This is the check that decides whether
# any p-value below is worth reporting.
message("Null calibration (condition labels shuffled before testing)")
local({
  d <- DS[[1]]
  ps <- vapply(1:40, function(s) {
    set.seed(1000 + s)
    p <- sample.int(d$n)
    perm_test_cka(d$K_lin, d$L_lin[p, p, drop = FALSE], B = 499)$p
  }, numeric(1))
  message(sprintf(
    "  %d shuffled runs: mean p = %.3f (want ~0.5), share below 0.05 = %.3f (want ~0.05)",
    length(ps), mean(ps), mean(ps < 0.05)
  ))
})

# Analyses

run_alignment <- function(analysis, comparison, K, L, n, kernel, extra = list()) {
  pt <- perm_test_cka(K, L)
  ci <- boot_ci_cka(K, L)
  tibble(
    analysis = analysis,
    comparison = comparison,
    kernel = kernel,
    n_conditions = n,
    cka = pt$cka,
    perm_p = pt$p,
    null_mean = pt$null_mean,
    null_q95 = pt$null_q95,
    boot_lo = ci[1],
    boot_hi = ci[2],
    !!!extra
  )
}

# --- A. Embedding vs phenotype ----------------------------------------------

message("A. Embedding vs phenotype alignment")
res_A <- map(DS, function(d) {
  bind_rows(
    run_alignment("A. Embedding vs phenotype", d$dataset, d$K_lin, d$L_lin, d$n, "linear"),
    run_alignment("A. Embedding vs phenotype", d$dataset, d$K_cor, d$L_cor, d$n, "correlation"),
    run_alignment("A. Embedding vs phenotype", d$dataset, d$K_rbf, d$L_rbf, d$n, "RBF")
  )
}) %>% list_rbind()

null_curves <- map(DS, function(d) {
  pt <- perm_test_cka(d$K_lin, d$L_lin)
  tibble(dataset = d$dataset, observed = pt$cka, null = pt$null)
}) %>% list_rbind()

# Robustness: the isomer pairs share an embedding vector, so drop one of each
# and confirm the alignment is not an artefact of duplicated rows.
message("A'. Isomer-collapsed robustness")
res_A_iso <- map(DS, function(d) {
  keep <- setdiff(d$keys, ISOMER_DROP)
  if (length(keep) == d$n) {
    return(NULL)
  }
  i <- match(keep, d$keys)
  run_alignment(
    "A'. Embedding vs phenotype, isomers collapsed", d$dataset,
    d$K_lin[i, i, drop = FALSE], d$L_lin[i, i, drop = FALSE], length(keep), "linear"
  )
}) %>% compact() %>% list_rbind()

# --- B. Phenotype geometry across panels ------------------------------------

# Does the condition-level phenotype structure reproduce across independently
# generated segregant panels? 
message("B. Cross-panel phenotype reproducibility")
pairs_B <- combn(names(DS), 2, simplify = FALSE)
res_B <- map(pairs_B, function(pr) {
  a <- DS[[pr[1]]]
  b <- DS[[pr[2]]]
  shared <- intersect(a$keys, b$keys)
  if (length(shared) < 8) {
    return(NULL)
  }
  ia <- match(shared, a$keys)
  ib <- match(shared, b$keys)
  bind_rows(
    run_alignment(
      "B. Phenotype vs phenotype", paste(pr[1], "vs", pr[2]),
      a$L_lin[ia, ia, drop = FALSE], b$L_lin[ib, ib, drop = FALSE],
      length(shared), "linear"
    ),
    run_alignment(
      "B. Phenotype vs phenotype", paste(pr[1], "vs", pr[2]),
      a$L_cor[ia, ia, drop = FALSE], b$L_cor[ib, ib, drop = FALSE],
      length(shared), "correlation"
    )
  )
}) %>% compact() %>% list_rbind()

# --- C. Baselines ------------------------------------------------------------

# Two reference points, so the headline number has a scale: a random embedding
# of the same shape (chance), and a three-level chemical class indicator (does
# the learned 256-dim space beat a trivial categorisation?).
message("C. Baselines")
res_C <- map(DS, function(d) {
  rand <- vapply(seq_len(N_RANDOM), function(i) {
    Xr <- matrix(rnorm(d$n * ncol(d$X)), nrow = d$n)
    cka(kernel_linear(Xr), d$L_lin)
  }, numeric(1))

  cls <- model.matrix(~ 0 + factor(d$Class))
  class_row <- run_alignment(
    "C. Baseline", paste(d$dataset, "- chemical class indicator"),
    kernel_linear(cls), d$L_lin, d$n, "linear"
  )

  bind_rows(
    class_row,
    tibble(
      analysis = "C. Baseline",
      comparison = paste(d$dataset, "- random 256-dim embedding"),
      kernel = "linear",
      n_conditions = d$n,
      cka = mean(rand),
      perm_p = NA_real_,
      null_mean = mean(rand),
      null_q95 = unname(stats::quantile(rand, 0.95)),
      boot_lo = unname(stats::quantile(rand, 0.025)),
      boot_hi = unname(stats::quantile(rand, 0.975))
    )
  )
}) %>% list_rbind()

# --- D. Does the embedding add anything beyond chemical class? --------------

# The class indicator in C aligns with phenotype geometry better than the
# learned embedding does, which raises the obvious question: is the embedding's
# alignment just the coarse class signal in disguise? Partial correlation
# between the off-diagonal entries of the centred kernels answers it directly.
upper_vec <- function(M) {
  n <- nrow(M)
  H <- diag(n) - 1 / n
  Mc <- H %*% M %*% H
  Mc[upper.tri(Mc)]
}

partial_r <- function(k, l, c) {
  r_kl <- stats::cor(k, l)
  r_kc <- stats::cor(k, c)
  r_lc <- stats::cor(l, c)
  (r_kl - r_kc * r_lc) / sqrt((1 - r_kc^2) * (1 - r_lc^2))
}

message("D. Embedding vs phenotype, controlling for chemical class")
res_D <- map(DS, function(d) {
  C <- kernel_linear(model.matrix(~ 0 + factor(d$Class)))
  kv <- upper_vec(d$K_lin)
  lv <- upper_vec(d$L_lin)
  cv <- upper_vec(C)

  stat <- function(K) partial_r(upper_vec(K), lv, cv)
  obs <- stat(d$K_lin)
  null <- vapply(seq_len(N_PERM), function(i) {
    p <- sample.int(d$n)
    stat(d$K_lin[p, p, drop = FALSE])
  }, numeric(1))

  bind_rows(
    tibble(
      analysis = "D. Embedding vs phenotype, class controlled",
      comparison = d$dataset, kernel = "linear", n_conditions = d$n,
      statistic = "Pearson r of kernel entries",
      cka = stats::cor(kv, lv), perm_p = NA_real_,
      null_mean = NA_real_, null_q95 = NA_real_,
      boot_lo = NA_real_, boot_hi = NA_real_
    ),
    tibble(
      analysis = "D. Embedding vs phenotype, class controlled",
      comparison = d$dataset, kernel = "linear", n_conditions = d$n,
      statistic = "Partial r, chemical class removed",
      cka = obs,
      perm_p = (1 + sum(null >= obs)) / (1 + length(null)),
      null_mean = mean(null),
      null_q95 = unname(stats::quantile(null, 0.95)),
      boot_lo = NA_real_, boot_hi = NA_real_
    )
  )
}) %>% list_rbind()

results <- bind_rows(res_A, res_A_iso, res_B, res_C, res_D) %>%
  mutate(statistic = coalesce(statistic, "CKA")) %>%
  relocate(statistic, .after = kernel) %>%
  left_join(
    map(DS, function(d) {
      tibble(comparison = d$dataset, min_shared_strains = d$min_overlap,
             neg_eigen_mass = d$psd$neg_mass)
    }) %>% list_rbind(),
    by = "comparison"
  )

write_csv(results, file.path(OUT_DIR, "condition_embedding_cka-0.5_sigma.csv"))
message("Wrote ", file.path(OUT_DIR, "condition_embedding_cka-0.5_sigma.csv"))

print(as.data.frame(results %>% select(analysis, comparison, kernel, statistic, n_conditions, cka, perm_p, null_q95)), digits = 3)

# Table

cka_gt <- results %>%
  mutate(comparison = sub("^Bloom2019_BYxRM", "Bloom2019", comparison),
         comparison = gsub("Bloom2019_BYxRM", "Bloom2019", comparison)) %>%
  select(analysis, comparison, kernel, statistic, n_conditions, cka, boot_lo, boot_hi, null_q95, perm_p) %>%
  gt(rowname_col = "comparison", groupname_col = "analysis") %>%
  fmt_number(columns = c(cka, boot_lo, boot_hi, null_q95), decimals = 3) %>%
  fmt_number(columns = perm_p, decimals = 4) %>%
  sub_missing(missing_text = "--") %>%
  cols_merge(columns = c(boot_lo, boot_hi), pattern = "[{1}, {2}]") %>%
  cols_label(
    kernel = "Kernel", statistic = "Statistic", n_conditions = "n", cka = "Value",
    boot_lo = "95% CI", null_q95 = "Null 95th", perm_p = "Perm. p"
  ) %>%
  tab_stubhead("Comparison") %>%
  cols_align(align = "center", columns = -comparison) %>%
  tab_footnote(
    "Unbiased HSIC estimator (Song et al. 2012); a null alignment sits at approximately 0, not at a positive offset.",
    locations = cells_column_labels(columns = cka)
  ) %>%
  tab_footnote(
    "Percentile bootstrap over conditions, 2,000 resamples. Indicative only - resampling duplicates kernel rows.",
    locations = cells_column_labels(columns = boot_lo)
  ) %>%
  tab_footnote(
    "95th percentile of the permutation null, shown so the observed value can be read against its own scale.",
    locations = cells_column_labels(columns = null_q95)
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
    table.font.size = px(13)
  ) %>%
  opt_table_font(font = "CMU Serif")

# gt emits no caption of its own; insert a \caption directly after the opening
# environment, and swap the unicode gt writes for math macros so the table also
# compiles under pdfLaTeX. Same helper as the significance tables.
as_tex <- function(gt_tbl, caption) {
  lines <- strsplit(as.character(as_latex(gt_tbl)), "\n", fixed = TRUE)[[1]]
  at <- grep("\\begin{table}", lines, fixed = TRUE)[1]
  stopifnot(!is.na(at))
  lines <- append(lines, paste0("\\caption{", caption, "}"), after = at)
  x <- paste(lines, collapse = "\n")
  x <- gsub("±", "$\\pm$", x, fixed = TRUE)
  gsub("Δ", "$\\Delta$", x, fixed = TRUE)
}

writeLines(
  as_tex(cka_gt, paste(
    "Centered kernel alignment between the chemical-condition embedding and the strain-level",
    "growth-phenotype geometry, at 0.5 sigma. Permutation p from 9,999 permutations of the",
    "condition labels."
  )),
  file.path(OUT_DIR, "condition_embedding_cka-0.5_sigma.tex")
)

# Figures

null_plot <- null_curves %>%
  mutate(dataset = sub("Bloom2019_BYxRM", "Bloom2019", dataset)) %>%
  ggplot(aes(x = null)) +
  geom_histogram(bins = 60, fill = cbPalette[3], colour = NA) +
  geom_vline(aes(xintercept = observed), colour = cbPalette[7], linewidth = 0.9) +
  facet_wrap(~dataset, scales = "free_y") +
  labs(
    x = "CKA under permuted condition labels",
    y = "Permutations",
    title = "Observed embedding-phenotype alignment against its permutation null"
  )

save_plot(
  file.path(OUT_DIR, "condition_embedding_cka-permutation_null.svg"),
  null_plot, width = 10, height = 3.6
)

mat_long <- function(M, keys, side, dataset) {
  n <- length(keys)
  tibble(
    row = factor(rep(keys, times = n), levels = keys),
    col = factor(rep(keys, each = n), levels = keys),
    value = as.vector(M),
    side = side,
    dataset = dataset
  )
}

kernel_long <- map(DS, function(d) {
  ord <- order(d$Class, d$keys)
  keys <- d$keys[ord]
  label <- sub("Bloom2019_BYxRM", "Bloom2019", d$dataset)
  scale01 <- function(M) {
    M <- M[ord, ord, drop = FALSE]
    (M - min(M)) / (max(M) - min(M))
  }
  bind_rows(
    mat_long(scale01(d$K_lin), keys, "Chemical embedding", label),
    mat_long(scale01(d$L_lin), keys, "Growth phenotype", label)
  )
}) %>% list_rbind()

heatmap_plot <- ggplot(kernel_long, aes(row, col, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Similarity\n(scaled)") +
  facet_grid(side ~ dataset, scales = "free", space = "free") +
  labs(x = NULL, y = NULL, title = "Condition-level similarity, chemistry versus phenotype") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
    axis.text.y = element_text(size = 6),
    panel.grid = element_blank()
  )

save_plot(
  file.path(OUT_DIR, "condition_similarity_heatmaps.svg"),
  heatmap_plot, width = 13, height = 9
)

message("Done. Outputs in ", OUT_DIR)
