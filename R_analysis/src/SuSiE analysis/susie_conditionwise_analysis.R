# Condition-wise SuSiE fine-mapping and comparison with conditionwise SHAP
# ---------------------------------------------------------------------------
# The training feather is a long/stacked table (one row per Strain x Condition),
# so each strain's genotype is repeated across every condition. Pooling all rows
# pseudo-replicates the genotypes and inflates n / the
# z-scores. Here we instead fine-map ONE condition at a time: subsetting to a
# condition leaves exactly one row per strain, giving a clean, chemical-specific
# fine-mapping whose credible sets line up 1:1 with the conditionwise SHAP file.
#
# Method (Zou et al. 2022, SuSiE-RSS): per condition we compute univariate
# logistic-regression z-scores (sanctioned for binary traits), an in-sample LD
# matrix, and pass them to susie_rss(). Design matrix is genotype-only (Y*),
# matching the '^Y' filter used on the SHAP output.

library(susieR)
library(fastglm)
library(Rfast)
library(arrow)
library(purrr)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ggvenn)
library(gt)

source("./src/utils/theme_set.R")
source("./src/utils/genomic_position.R")

# Fine-mapping helpers

#' Marginal (single-marker) logistic-regression z-scores.
#'
#' @param X numeric matrix (samples x markers), genotype-only, all varying.
#' @param pheno binary 0/1 response.
#' @return numeric vector of z-scores (beta / se), one per column of X.
uni_logreg_z <- function(X, pheno) {
  int <- rep(1, nrow(X))
  vapply(
    seq_len(ncol(X)),
    function(j) {
      fit <- fastglm(cbind(int, X[, j]), pheno, family = binomial(), method = 2)
      fit$coefficients[2] / fit$se[2]
    },
    numeric(1)
  )
}

#' Turn a fitted susie object into tidy per-marker and per-credible-set tables.
extract_susie <- function(res, feature_names, condition) {
  pip <- res$pip
  cs_list <- res$sets$cs                       # named list L1, L2, ... or NULL

  cs_id <- rep(NA_character_, length(pip))
  if (!is.null(cs_list) && length(cs_list) > 0) {
    for (nm in names(cs_list)) cs_id[cs_list[[nm]]] <- nm
  }

  pip_df <- tibble(
    Condition = condition,
    Gene      = feature_names,
    PIP       = pip,
    cs_id     = cs_id
  )

  if (!is.null(cs_list) && length(cs_list) > 0) {
    cs_names <- names(cs_list)
    cs_df <- map_dfr(seq_along(cs_list), function(k) { # map_dfr combines outputs
      idx <- cs_list[[k]]
      tibble(
        Condition    = condition,
        cs_id        = cs_names[k],
        Gene         = feature_names[idx],
        PIP          = pip[idx],
        cs_size      = length(idx),
        min_abs_corr = res$sets$purity[k, "min.abs.corr"],
        coverage     = res$sets$coverage[k]
      )
    }) %>% arrange(cs_id, desc(PIP))
  } else {
    cs_df <- tibble(
      Condition = character(), cs_id = character(), Gene = character(),
      PIP = numeric(), cs_size = integer(),
      min_abs_corr = numeric(), coverage = numeric()
    )
  }

  list(pip_df = pip_df, cs_df = cs_df)
}

#' Fine-map a single condition with SuSiE-RSS.
#'
#' @param df_cond rows of the training frame for ONE condition.
#' @param condition condition name (for labelling).
#' @param pheno_col name of the binary phenotype column.
#' @param L max number of single effects (fewer signals per single condition).
#' @param coverage target credible-set coverage.
fit_susie_condition <- function(df_cond, condition,
                                pheno_col = "Phenotype",
                                L = 10, coverage = 0.95) {
  ycols <- grep("^Y", colnames(df_cond), value = TRUE)
  Xall <- as.matrix(df_cond[, ycols])
  storage.mode(Xall) <- "double"

  # Keep only markers that actually vary among THIS condition's strains
  keep <- which(Rfast::colVars(Xall) > 0)
  X <- Xall[, keep, drop = FALSE]
  feats <- ycols[keep]
  rm(Xall)

  pheno <- df_cond[[pheno_col]]
  n <- nrow(X)

  z <- uni_logreg_z(X, pheno)

  # Guard against fit failures / extreme separation producing non-finite z.
  ok <- is.finite(z)
  if (any(!ok)) {
    X <- X[, ok, drop = FALSE]; feats <- feats[ok]; z <- z[ok]
  }

  R <- Rfast::cora(X, large = TRUE)
  res <- susie_rss(z = z, R = R, n = n, L = L,
                   var_y = var(pheno), coverage = coverage)

  out <- extract_susie(res, feats, condition)
  out$meta <- tibble(
    Condition = condition, n = n,
    n_markers = length(feats), n_cs = length(res$sets$cs)
  )
  rm(X, R, res); gc(verbose = FALSE)
  out
}

# fine-map every condition for a dataset

#' @param dataname feather stem, e.g. "bloom2013".
#' @param shap_name capitalised dataset name used in results/SHAP paths.
#' @param conditions optional subset of conditions (NULL = all).
run_conditionwise_susie <- function(dataname, shap_name,
                                     conditions = NULL, L = 10) {
  df <- read_feather(sprintf(
    "data/training_data/0.5_sigma/%s_clf.feather", dataname
  )) %>%
    select(Strain, Condition, Phenotype, starts_with("Y"))   # drop latents

  if (is.null(conditions)) conditions <- sort(unique(df$Condition))

  out_dir <- sprintf("results/result_susie/0.5_condition-wise/%s", shap_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  pip_all <- vector("list", length(conditions))
  cs_all  <- vector("list", length(conditions))
  meta_all <- vector("list", length(conditions))

  for (i in seq_along(conditions)) {
    cond <- conditions[i]
    message(sprintf("[%s] %2d/%d  %s", dataname, i, length(conditions), cond))
    fit <- fit_susie_condition(filter(df, Condition == cond), cond, L = L)

    write_csv(fit$cs_df, file.path(out_dir, paste0(cond, "_credible_sets.csv")))
    pip_all[[i]] <- fit$pip_df
    cs_all[[i]]  <- fit$cs_df
    meta_all[[i]] <- fit$meta
  }
  rm(df); gc(verbose = FALSE)

  combined_pip <- bind_rows(pip_all)
  write_csv(combined_pip,
            sprintf("results/result_susie/%s_conditionwise_pip.csv", shap_name))
  write_csv(bind_rows(meta_all),
            sprintf("results/result_susie/%s_conditionwise_meta.csv", shap_name))

  list(pip = combined_pip, cs = bind_rows(cs_all), meta = bind_rows(meta_all))
}

# Comparison with conditionwise SHAP

#' Load conditionwise SHAP, keeping only positive-SHAP gene (Y*) features and
#' harmonising the lowercase `condition` column to `Condition`.
load_conditionwise_shap <- function(shap_name) {
  read_parquet(sprintf(
    "data/shap/sigmas/shap_classification_0.5/%s_conditionwise_shap.parquet",
    shap_name
  )) %>%
    filter(Value > 0, grepl("^Y", Feature)) %>%
    rename(Gene = Feature, SHAP = Value, Condition = condition)
}

#' Compare SuSiE credible sets / PIPs against conditionwise SHAP.
#'
#' Produces, per dataset: (i) a Spearman(PIP, SHAP) correlation table + gt/plot,
#' (ii) per-condition PIP-vs-SHAP scatter coloured by credible-set membership,
#' (iii) per-condition Venn of credible-set genes vs top-SHAP genes, and
#' (iv) an overlap summary (where do credible-set genes rank in the SHAP list).
#'
#' @param susie_res result list from run_conditionwise_susie().
#' @param shap_name capitalised dataset name (e.g. "Bloom2013").
#' @param top_shap_n number of top-SHAP genes per condition for the Venn.
compare_susie_shap <- function(susie_res, shap_name, top_shap_n = 50) {
  shap_file <- sprintf(
    "data/shap/sigmas/shap_classification_0.5/%s_conditionwise_shap.parquet",
    shap_name
  )
  if (!file.exists(shap_file)) {
    message(sprintf("No conditionwise SHAP for %s - skipping comparison.",
                    shap_name))
    return(invisible(NULL))
  }

  shap <- load_conditionwise_shap(shap_name)
  cmp_dir <- sprintf(
    "results/result_susie/0.5_condition-wise/%s/comparison", shap_name
  )
  dir.create(file.path(cmp_dir, "scatter"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(cmp_dir, "venn"),    recursive = TRUE, showWarnings = FALSE)

  # Join PIP and SHAP on Gene + Condition (markers scored by both methods).
  joined <- susie_res$pip %>%
    inner_join(shap, by = c("Gene", "Condition")) %>%
    mutate(in_cs = !is.na(cs_id))

  conditions <- sort(unique(joined$Condition))
  cor_rows <- list(); overlap_rows <- list()

  for (cond in conditions) {
    jc <- filter(joined, Condition == cond)
    cs_genes <- susie_res$cs %>% filter(Condition == cond) %>% pull(Gene) %>% unique()

    # (i) Spearman correlation between PIP and SHAP. This genuinely needs paired
    # (PIP, SHAP) per gene, so the inner join `jc` is the right universe here:
    # imputing SHAP = 0 for markers the model ignored would inject a large block
    # of concordant (PIP~0, SHAP=0) points and inflate the rank correlation.
    rho <- if (nrow(jc) >= 3) {
      suppressWarnings(cor(jc$PIP, jc$SHAP, method = "spearman"))
    } else NA_real_
    cor_rows[[cond]] <- tibble(
      Condition = cond, n_common = nrow(jc), spearman = rho
    )

    # (ii) PIP-vs-SHAP scatter, credible-set members highlighted.
    p_sc <- jc %>%
      ggplot(aes(x = PIP, y = SHAP, color = in_cs)) +
      geom_point(alpha = 0.6) +
      ggrepel::geom_text_repel(
        data = filter(jc, in_cs), aes(label = Gene),
        size = 2.5, max.overlaps = 12, show.legend = FALSE
      ) +
      scale_color_manual(values = c(`FALSE` = "grey70", `TRUE` = "#D55E00"),
                         labels = c("Not in CS", "Credible set")) +
      labs(title = paste(shap_name, "-", cond),
           x = "SuSiE PIP", y = "SHAP value", color = NULL)
    save_plot(file.path(cmp_dir, "scatter", paste0(cond, "_pip_vs_shap.svg")),
              plot = p_sc, width = 15, height = 10, units = "cm")

    # (iii) Venn: credible-set genes vs top-SHAP genes.
    top_shap <- shap %>% filter(Condition == cond) %>%
      slice_max(SHAP, n = top_shap_n) %>% pull(Gene)
    if (length(cs_genes) > 0) {
      venn_sets <- list(`Credible sets` = cs_genes,
                        `Top SHAP` = top_shap)
      p_v <- ggvenn(venn_sets, auto_scale = TRUE, show_percentage = FALSE)
      save_plot(file.path(cmp_dir, "venn", paste0(cond, "_venn.svg")),
                plot = p_v, width = 12, height = 8, units = "cm")
    }

    # (iv) Overlap metrics: where do credible-set genes sit in the SHAP ranking?
    # Rank each CS gene's SHAP against the FULL per-condition SHAP distribution
    # (every positive-SHAP gene the model scored this condition)
    shap_c <- shap %>% filter(Condition == cond)
    cs_scored <- intersect(cs_genes, shap_c$Gene)
    shap_pct <- if (nrow(shap_c) >= 3 && length(cs_scored) > 0) {
      shap_c <- shap_c %>% mutate(shap_pctile = rank(SHAP) / n())
      median(shap_c$shap_pctile[shap_c$Gene %in% cs_genes])
    } else NA_real_
    overlap_rows[[cond]] <- tibble(
      Condition            = cond,
      n_cs_genes           = length(cs_genes),
      n_cs_genes_scored    = length(cs_scored),
      n_top_shap_hit       = length(intersect(cs_genes, top_shap)),
      median_shap_pctile   = shap_pct
    )
  }

  cor_df <- bind_rows(cor_rows)
  overlap_df <- bind_rows(overlap_rows)
  write_csv(cor_df, file.path(cmp_dir, "susie_shap_spearman.csv"))
  write_csv(overlap_df, file.path(cmp_dir, "susie_shap_overlap_summary.csv"))

  # Spearman summary: gt table  + bar plot.
  tryCatch({
    gt_tbl <- cor_df %>% arrange(desc(spearman)) %>%
      gt() %>% fmt_number(columns = spearman, decimals = 3) %>%
      tab_header(title = paste(shap_name, "- SuSiE PIP vs SHAP (Spearman)"))
    gtsave(gt_tbl, file.path(cmp_dir, "susie_shap_spearman.html"))
  }, error = function(e) message("gt table skipped: ", conditionMessage(e)))

  p_bar <- cor_df %>% filter(!is.na(spearman)) %>%
    ggplot(aes(x = reorder(Condition, spearman), y = spearman)) +
    geom_col(fill = "#0072B2") + coord_flip() +
    labs(x = NULL, y = "Spearman(PIP, SHAP)",
         title = paste(shap_name, "- per-condition PIP/SHAP rank correlation"))
  save_plot(file.path(cmp_dir, "susie_shap_spearman_bar.svg"),
            plot = p_bar, width = 16, height = 20, units = "cm")

  list(spearman = cor_df, overlap = overlap_df, joined = joined)
}

#' Cross-reference credible-set genes with previously detected QTL genes.
#' Only available for Bloom2013 (data/qtl/detected_qtl_bloom2013.csv).
compare_susie_qtl_bloom2013 <- function(susie_res, shap_name = "Bloom2013") {
  qtl_file <- "data/qtl/detected_qtl_bloom2013.csv"
  if (!file.exists(qtl_file)) {
    message("No QTL file - skipping QTL cross-reference."); return(invisible(NULL))
  }
  qtl_genes <- read_csv(qtl_file, show_col_types = FALSE) %>%
    select(Trait, Genes) %>%
    mutate(Genes = strsplit(as.character(Genes), split = "|", fixed = TRUE)) %>%
    unnest_longer(Genes) %>%
    filter(grepl("^Y", Genes)) %>%
    group_by(Trait) %>% summarise(qtl_genes = list(unique(Genes)), .groups = "drop")

  cs_by_cond <- susie_res$cs %>%
    group_by(Condition) %>% summarise(cs_genes = list(unique(Gene)), .groups = "drop")

  summ <- cs_by_cond %>%
    left_join(qtl_genes, by = c("Condition" = "Trait")) %>%
    mutate(
      n_cs_genes  = map_int(cs_genes, length),
      n_qtl_genes = map_int(qtl_genes, ~ if(is.null(.x)) 0L else length(.x)),
      n_overlap   = map2_int(cs_genes, qtl_genes,
                             ~ if (is.null(.y)) 0L else length(intersect(.x, .y)))
    ) %>%
    select(Condition, n_cs_genes, n_qtl_genes, n_overlap)

  out_dir <- sprintf("results/result_susie/0.5_condition-wise/%s/comparison", shap_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write_csv(summ, file.path(out_dir, "susie_qtl_overlap_summary.csv"))
  summ
}

#  (credible-set members across evidence lines)

marginal_stats_for_genes <- function(df_cond, genes, pheno_col = "Phenotype") {
  genes <- intersect(genes, colnames(df_cond))
  if (length(genes) == 0) {
    return(tibble(Gene = character(), beta = numeric(), se = numeric(),
                  z = numeric(), converged = logical(), maf = numeric()))
  }
  pheno <- df_cond[[pheno_col]]
  int <- rep(1, nrow(df_cond))
  map_dfr(genes, function(g) {
    x <- as.double(df_cond[[g]] >= 1)              # collapse 0/1/2/3 -> 0/1
    fit <- fastglm(cbind(int, x), pheno, family = binomial(), method = 2)
    p <- mean(x)
    tibble(Gene = g, beta = fit$coefficients[2], se = fit$se[2],
           z = fit$coefficients[2] / fit$se[2],
           converged = fit$converged, maf = min(p, 1 - p))
  })
}

#' Build the triangulation table for one condition's credible-set members.
#'
#' @param condition condition name.
#' @param cs_df credible-set table (Condition, cs_id, Gene, PIP, ...).
#' @param shap conditionwise SHAP (Gene, SHAP, Condition); may be empty.
#' @param qtl_all tibble(Trait, Gene, Pheno_fraction_explained) or NULL (non-Bloom2013).
#' @param df_cond genotype rows for THIS condition (one row per strain).
triangulate_condition <- function(condition, cs_df, shap, qtl_all, df_cond) {
  cs <- cs_df %>% filter(Condition == condition)
  if (nrow(cs) == 0) return(NULL)

  marg <- marginal_stats_for_genes(df_cond, unique(cs$Gene))
  shap_ranked <- shap %>% filter(Condition == condition) %>%
    arrange(desc(SHAP)) %>% mutate(shap_rank = row_number())

  tri <- cs %>%
    left_join(marg, by = "Gene") %>%
    left_join(select(shap_ranked, Gene, SHAP, shap_rank), by = "Gene")

  if (!is.null(qtl_all)) {
    qc <- qtl_all %>% filter(Trait == condition) %>%
      group_by(Gene) %>%
      summarise(Pheno_fraction_explained = max(Pheno_fraction_explained),
                .groups = "drop") %>%
      mutate(in_QTL = TRUE)
    tri <- tri %>% left_join(qc, by = "Gene") %>%
      mutate(in_QTL = coalesce(in_QTL, FALSE))
  } else {
    tri <- tri %>% mutate(in_QTL = NA, Pheno_fraction_explained = NA_real_)
  }

  tri %>%
    mutate(genomic_pos = map_dbl(Gene, get_genomic_position)) %>%
    group_by(cs_id) %>% mutate(lead = PIP == max(PIP)) %>% ungroup() %>%
    arrange(cs_id, desc(PIP))
}

#' Genomic-position vs PIP plot (analog of plot_position_odds): point size = SHAP,
#' colour = credible set, top-PIP gene per CS labelled.
plot_position_pip <- function(tri_df, title) {
  tri_df <- tri_df %>% mutate(SHAP_size = coalesce(SHAP, 0))
  ggplot(tri_df, aes(x = genomic_pos, y = PIP)) +
    geom_point(aes(size = SHAP_size, color = cs_id), alpha = 0.75) +
    ggrepel::geom_text_repel(data = filter(tri_df, lead),
                             aes(label = Gene), size = 2.8, max.overlaps = 20) +
    scale_size(range = c(1, 6)) +
    labs(x = "Genomic position (relative)", y = "SuSiE PIP",
         color = "Credible set", size = "SHAP", title = title)
}

#' Triangulate every condition for a dataset. Reads the per-condition
#' `*_credible_sets.csv` (written by run_conditionwise_susie) + the feather
#' (once) + SHAP + QTL, and writes `<condition>_triangulation.csv` +
#' `<condition>_pip_position.svg` under .../<shap_name>/triangulation/.
run_triangulation <- function(dataname, shap_name) {
  out_dir <- sprintf("results/result_susie/0.5_condition-wise/%s", shap_name)
  cs_files <- list.files(out_dir, pattern = "_credible_sets\\.csv$",
                         full.names = TRUE)
  if (length(cs_files) == 0) {
    message("No credible-set CSVs for ", shap_name,
            " - run run_conditionwise_susie() first.")
    return(invisible(NULL))
  }
  cs_df <- map_dfr(cs_files, read_csv, show_col_types = FALSE)

  shap_file <- sprintf(
    "data/shap/sigmas/shap_classification_0.5/%s_conditionwise_shap.parquet",
    shap_name
  )
  shap <- if (file.exists(shap_file)) load_conditionwise_shap(shap_name)
          else tibble(Gene = character(), SHAP = numeric(), Condition = character())

  qtl_file <- "data/qtl/detected_qtl_bloom2013.csv"
  qtl_all <- if (shap_name == "Bloom2013" && file.exists(qtl_file)) {
    read_csv(qtl_file, show_col_types = FALSE) %>%
      select(Trait, Genes, Pheno_fraction_explained) %>%
      mutate(Genes = strsplit(as.character(Genes), split = "|", fixed = TRUE)) %>%
      unnest_longer(Genes) %>% filter(grepl("^Y", Genes)) %>%
      rename(Gene = Genes)
  } else NULL

  df <- read_feather(sprintf(
    "data/training_data/0.5_sigma/%s_clf.feather", dataname
  )) %>% select(Condition, Phenotype, starts_with("Y"))

  tri_dir <- file.path(out_dir, "triangulation")
  dir.create(tri_dir, recursive = TRUE, showWarnings = FALSE)

  conditions <- sort(unique(cs_df$Condition))
  tri_all <- vector("list", length(conditions))
  for (i in seq_along(conditions)) {
    cond <- conditions[i]
    tri <- triangulate_condition(cond, cs_df, shap, qtl_all,
                                 filter(df, Condition == cond))
    if (is.null(tri) || nrow(tri) == 0) next
    write_csv(tri, file.path(tri_dir, paste0(cond, "_triangulation.csv")))
    save_plot(file.path(tri_dir, paste0(cond, "_pip_position.svg")),
              plot = plot_position_pip(tri, paste(shap_name, "-", cond)),
              width = 18, height = 10, units = "cm")
    tri_all[[i]] <- tri
  }
  rm(df); gc(verbose = FALSE)

  combined <- bind_rows(tri_all)
  write_csv(combined,
            sprintf("results/result_susie/%s_triangulation.csv", shap_name))
  invisible(combined)
}

# ---------------------------------------------------------------------------
# Full pipeline runner
# ---------------------------------------------------------------------------

# Bloom2019 = the base BY x RM cross (bloom2019_clf.feather), NOT the _BYxM22 /
# _RMxYPS163 crosses. Bloom2015 has no conditionwise SHAP file, so it is
# fine-mapped but its SHAP comparison is skipped.
DATASETS <- tibble::tribble(
  ~dataname,     ~shap_name,
  "bloom2013",   "Bloom2013",
  "bloom2015",   "Bloom2015",
  "bloom2019",   "Bloom2019"
)

#' Run the full condition-wise SuSiE + SHAP comparison over all datasets.
#' @param conditions optional condition subset (for quick tests); NULL = all.
main_conditionwise <- function(datasets = DATASETS, conditions = NULL) {
  results <- list()
  for (k in seq_len(nrow(datasets))) {
    dn <- datasets$dataname[k]; sn <- datasets$shap_name[k]
    susie_res <- run_conditionwise_susie(dn, sn, conditions = conditions)
    shap_cmp  <- compare_susie_shap(susie_res, sn)
    qtl_cmp   <- if (sn == "Bloom2013")
      compare_susie_qtl_bloom2013(susie_res, sn) else NULL
    tri_cmp   <- run_triangulation(dn, sn)
    results[[sn]] <- list(susie = susie_res, shap = shap_cmp,
                          qtl = qtl_cmp, triangulation = tri_cmp)
  }
  invisible(results)
}

# To run the full pipeline (all datasets, all conditions; ~30-40 min):
#   results <- main_conditionwise()
#
# To run a quick subset on one dataset:
  # res <- run_conditionwise_susie("bloom2013", "Bloom2013",
  #                                conditions = c("4NQO", "maltose"))
  # cmp <- compare_susie_shap(res, "Bloom2013")
  # qtl <- compare_susie_qtl_bloom2013(res)
