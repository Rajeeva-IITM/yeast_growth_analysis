# Condition-wise FINEMAP / CAVIAR fine-mapping and comparison with SHAP
# ---------------------------------------------------------------------------
# Companion to `src/SuSiE analysis/susie_conditionwise_analysis.R`. Same data,
# same phenotype coding, same conditionwise SHAP / QTL comparisons - but the
# fine-mapping is done by FINEMAP and CAVIAR (driven through the `finemapr`
# package) instead of SuSiE, so the three methods can be contrasted.
#
# Prerequisite: run `src/finemapping/install_finemap_caviar.sh` inside WSL once.
# It fetches FINEMAP v1.1 + builds CAVIAR into `tools/` (gitignored) and writes
# the Windows .cmd shims that let native Windows R call the Linux binaries.
# See `src/finemapping/SETUP_LOG.md` for the full setup record.
#
# ---------------------------------------------------------------------------
# Why this script is region-based while the SuSiE one is not
# ---------------------------------------------------------------------------
# susie_conditionwise_analysis.R fine-maps all ~2.7k-3k varying markers of a
# condition in a single susie_rss() call. That is fine for SuSiE, but neither
# FINEMAP nor CAVIAR can be used that way here:
#
#   * n < m. Each condition has ~600 strains but ~2,700 markers, so the
#     in-sample LD matrix has rank <= 600 and is badly singular. CAVIAR
#     factorises the LD matrix directly and degrades sharply when it is
#     rank-deficient; FINEMAP is more tolerant but still assumes a
#     well-behaved correlation matrix.
#   * CAVIAR enumerates causal configurations (O(m^c) for `-c` causals). With
#     m = 2,700 and c = 2 that is ~3.6M configurations per condition, each
#     involving linear algebra on the full matrix - not tractable. CAVIAR is
#     designed for locus-scale input (its own manual uses ~100 SNPs).
#
# Both tools are locus fine-mappers, so markers are split into regions and each
# region is fine-mapped independently, then results are stitched back into one
# genome-wide table per condition.
#
# Region definition: fine-map ASSOCIATION PEAKS, not the whole genome
# ---------------------------------------------------------------------------
# FINEMAP and CAVIAR are locus fine-mappers: they assume the region handed to
# them already contains a signal (a genome-wide-significant peak) and answer
# "which marker in this locus is causal?". Tiling the whole genome into fixed
# marker-count windows (the old region_by = "window") violates that assumption:
# most windows in a condition carry no signal, yet both tools' priors condition
# on >= 1 causal variant, so a null window is forced to spread ~1 unit of
# posterior mass diffusely and the 95%-coverage credible set swallows that
# whole diffuse tail. That, not any localization failure, is what produced the
# genome-sized (~1,300-1,700 gene) credible sets: leads were already sharp
# (maxPIP ~ 0.99), the sets were just padded by every signal-free window.
#
# We therefore define loci the way a GWAS -> fine-map pipeline does. Per
# condition we threshold the marginal association scan at genome-wide
# significance (region_by = "signal": two-sided Bonferroni |z| >=
# qnorm(1 - alpha/(2m)), alpha = 0.05, ~ LOD 4 for m ~ 2.7k markers), and
# fine-map a window of +/- `flank` position-ordered markers around each
# surviving peak (overlapping windows merged, each locus confined to its lead
# marker's chromosome). Signal-free regions are never fine-mapped, so the
# diffuse-tail bloat cannot arise - no post-hoc purity filter is needed. A
# condition with no genome-wide-significant peak simply produces no locus (see
# the *_no_signal_conditions.csv record). n_causal = 1 keeps sum(PIP) ~ 1 per
# locus. The old genome-tiling behaviour (and the "why it bloated" narrative)
# stays reproducible via region_by = "window" / "chromosome".
#
# Two further practical steps, both defaults, both recorded in the meta table:
#   * markers in near-perfect LD (|r| >= `prune_r2`) are collapsed to a single
#     tag marker (highest |z| wins); the collapsed partners are carried along
#     so credible sets can be read back in terms of all implicated genes.
#     Without this, perfectly collinear markers make CAVIAR singular.
#   * a region is capped at `max_markers`; if it is larger, only the top-|z|
#     markers are kept. Anything dropped is counted in the meta table rather
#     than silently discarded.
#
# Method: as in the SuSiE script, per condition we compute univariate logistic
# regression z-scores (sanctioned for binary traits) and an in-sample LD matrix
# from the genotype-only (Y*) design, matching the '^Y' filter used on SHAP.

library(finemapr)
library(fastglm)
library(Rfast)
library(arrow)
library(purrr)
library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(ggplot2)
library(ggvenn)
library(gt)

source("./src/utils/theme_set.R")
source("./src/utils/genomic_position.R")
source("./src/finemapping/finemapr_setup.R")

# Output roots (results/ is gitignored, as is tools/).
FM_RESULT_ROOT <- "results/result_finemapping/0.5_condition-wise"
FM_WORK_ROOT   <- "tools/work"        # scratch dirs for the tools' own I/O

# z-scores

#' Marginal (single-marker) logistic-regression z-scores.
#'
#' Intentionally identical to `uni_logreg_z()` in
#' `src/SuSiE analysis/susie_conditionwise_analysis.R`: the two pipelines must
#' be fed exactly the same summary statistics for their results to be
#' comparable. Duplicated rather than shared so either script can be run
#' standalone; keep the two in sync if you change one.
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

# Region construction

#' Chromosome number from a systematic ORF name (letter 2 -> 1..16).
#'
#' Same encoding `get_genomic_position()` uses for its chromosome component.
gene_chromosome <- function(gene) {
  match(substr(gene, 2, 2), LETTERS)
}

#' Split markers into fine-mapping regions.
#'
#' @param feats character vector of marker (gene) names.
#' @param by "chromosome" (default) or "window".
#' @param window_size markers per region when `by = "window"`; markers are
#'   ordered by genomic position first so windows are contiguous.
#' @return integer/character vector of region labels, one per marker.
make_regions <- function(feats, by = c("chromosome", "window"),
                         window_size = 150) {
  by <- match.arg(by)
  if (by == "chromosome") {
    chr <- gene_chromosome(feats)
    # Markers whose name does not decode (non-standard ORF) go to their own
    # bucket rather than being dropped.
    return(ifelse(is.na(chr), "unplaced", paste0("chr", chr)))
  }

  pos <- map_dbl(feats, get_genomic_position)
  ord <- order(pos)
  lab <- rep(NA_character_, length(feats))
  lab[ord] <- paste0("win", ((seq_along(ord) - 1) %/% window_size) + 1)
  lab
}

#' Call genome-wide-significant association peaks and turn each into a locus.
#'
#' The signal-anchored region rule (region_by = "signal"). FINEMAP/CAVIAR are
#' locus fine-mappers, so rather than tiling the genome we hand them only the
#' regions that carry a signal: threshold the marginal scan at genome-wide
#' significance and build a window around each surviving peak. Needs only the
#' z-scores already computed by the caller - no genome-wide LD matrix.
#'
#' @param feats character vector of marker (ORF) names, aligned to `z`.
#' @param z marginal logistic-regression z-scores (from uni_logreg_z()).
#' @param alpha genome-wide significance level; two-sided Bonferroni threshold
#'   |z| >= qnorm(1 - alpha/(2m)), m = number of markers.
#' @param flank markers to keep on each side of a significant marker (in
#'   genome-wide position-rank space) when forming its locus window.
#' @return character vector aligned to `feats`: "locus1".."locusk" for markers
#'   inside a called locus, NA otherwise (NA markers are skipped by the region
#'   loop). All-NA when no marker reaches significance.
call_signal_loci <- function(feats, z, alpha = 0.05, flank = 50) {
  m <- length(z)
  labels <- rep(NA_character_, m)
  z_thr <- stats::qnorm(1 - alpha / (2 * m))
  sig <- which(abs(z) >= z_thr)
  if (length(sig) == 0) return(labels)

  pos  <- map_dbl(feats, get_genomic_position)
  rnk  <- order(order(pos))            # genome-wide position rank, 1..m
  chr  <- gene_chromosome(feats)

  # Candidate window (in position-rank space) around each significant marker,
  # then merge overlapping windows into loci (interval union over sorted ranks).
  sig_rank <- sort(rnk[sig])
  starts <- sig_rank - flank
  ends   <- sig_rank + flank
  cur_s <- starts[1]; cur_e <- ends[1]
  merged <- list()
  for (i in seq_along(sig_rank)[-1]) {
    if (starts[i] <= cur_e) {
      cur_e <- max(cur_e, ends[i])
    } else {
      merged[[length(merged) + 1]] <- c(cur_s, cur_e)
      cur_s <- starts[i]; cur_e <- ends[i]
    }
  }
  merged[[length(merged) + 1]] <- c(cur_s, cur_e)

  # Assign markers to loci; confine each locus to its lead marker's chromosome
  # so a +/- flank window near a chromosome end cannot spill across the boundary.
  k <- 0L
  for (iv in merged) {
    in_win <- which(rnk >= iv[1] & rnk <= iv[2])
    sig_in <- in_win[abs(z[in_win]) >= z_thr]
    if (length(sig_in) == 0) next
    lead_chr <- chr[sig_in[which.max(abs(z[sig_in]))]]
    keep <- in_win[which(chr[in_win] == lead_chr)]   # which() drops NA chr
    if (length(keep) < 2) next
    k <- k + 1L
    labels[keep] <- paste0("locus", k)
  }
  labels
}

#' Collapse markers in near-perfect LD to a single tag marker.
#'
#' Greedy: markers are visited in decreasing |z|; the first unassigned marker
#' becomes a tag and absorbs every still-unassigned marker with |r| >= `thr`.
#' This is what keeps CAVIAR's LD matrix non-singular in a cross, where
#' neighbouring markers are frequently perfectly correlated.
#'
#' @param R correlation matrix for the region (markers x markers, named).
#' @param z z-scores aligned to `colnames(R)`.
#' @param thr absolute-correlation threshold.
#' @return tibble(tag, member) mapping each tag to the markers it represents
#'   (a tag is always a member of itself).
ld_tag_groups <- function(R, z, thr = 0.99) {
  feats <- colnames(R)
  ord <- order(abs(z), decreasing = TRUE)
  assigned <- rep(NA_character_, length(feats))
  names(assigned) <- feats

  for (i in ord) {
    if (!is.na(assigned[i])) next
    tag <- feats[i]
    partners <- which(is.na(assigned) & abs(R[i, ]) >= thr)
    assigned[partners] <- tag
    assigned[i] <- tag
  }
  tibble(tag = unname(assigned), member = feats)
}

# Tool runners (one region at a time)

#' Run FINEMAP on one region.
#'
#' @param tab tibble(snp, zscore) for the region.
#' @param ld named correlation matrix for the region's markers.
#' @param n number of individuals.
#' @param dir_run scratch directory for the tool's input/output files.
#' @param n_causal_max FINEMAP's --n-causal-max.
#' @return tibble(Gene, PIP, log10bf) or NULL if the run failed.
finemap_region <- function(tab, ld, n, dir_run, n_causal_max = 3) {
  out <- try(
    run_finemap(tab, ld, n, dir_run = dir_run,
                args = paste("--n-causal-max", n_causal_max)),
    silent = TRUE
  )
  if (inherits(out, "try-error") || !isTRUE(out$status)) {
    warning("FINEMAP failed in ", dir_run, call. = FALSE)
    return(NULL)
  }
  # region.snp columns: index, snp, snp_prob, snp_log10bf
  tibble(
    Gene    = as.character(out$snp$snp),
    PIP     = as.numeric(out$snp$snp_prob),
    log10bf = as.numeric(out$snp$snp_log10bf)
  )
}

#' Run CAVIAR on one region.
#'
#' @param n_causal CAVIAR's -c (max number of causal variants).
#' @return tibble(Gene, PIP, prob_in_set, in_tool_set) or NULL if it failed.
#'
#' CAVIAR's `*_post` columns are SNP_ID / Prob_in_pCausalSet /
#' Causal_Post._Prob, which finemapr renames to snp / snp_prob_set / snp_prob.
#' `snp_prob` (Causal_Post._Prob) is the posterior causal probability, i.e. the
#' PIP analogue; `out$set` is CAVIAR's own rho-level causal set.
caviar_region <- function(tab, ld, dir_run, n_causal = 3) {
  out <- try(
    run_caviar(tab, ld, dir_run = dir_run, args = paste("-c", n_causal)),
    silent = TRUE
  )
  if (inherits(out, "try-error") || !isTRUE(out$status)) {
    warning("CAVIAR failed in ", dir_run, call. = FALSE)
    return(NULL)
  }
  tibble(
    Gene        = as.character(out$snp$snp),
    PIP         = as.numeric(out$snp$snp_prob),
    prob_in_set = as.numeric(out$snp$snp_prob_set),
    in_tool_set = as.character(out$snp$snp) %in% out$set
  )
}

#' Greedy credible set from PIPs: highest PIP first until cumulative >= coverage.
#'
#' finemapr ships `extract_credible_set()`, but it only supports the newer
#' `finemapr()` S3 pipeline (it expects `x$snp` to be a *list* of tables with a
#' `snp_prob_cumsum` column) and errors on the `run_finemap()`/`run_caviar()`
#' return value. Computing it here keeps FINEMAP, CAVIAR and SuSiE credible
#' sets defined the same way.
credible_set_from_pip <- function(df, coverage = 0.95) {
  df <- arrange(df, desc(PIP))
  cum <- cumsum(df$PIP) / max(sum(df$PIP), .Machine$double.eps)
  # keep everything up to and including the first marker crossing `coverage`
  k <- which(cum >= coverage)[1]
  if (is.na(k)) k <- nrow(df)
  df$in_cs <- seq_len(nrow(df)) <= k
  df$cs_cumsum <- cum
  df
}

# Per-condition driver

#' Fine-map a single condition with FINEMAP and CAVIAR, region by region.
#'
#' @param df_cond rows of the training frame for ONE condition.
#' @param condition condition name (for labelling).
#' @param dataset dataset label (for scratch paths).
#' @param methods which tools to run.
#' @param region_by "signal" (default; genome-wide-significant association peaks
#'   via call_signal_loci), or "window"/"chromosome" (genome-tiling via
#'   make_regions, kept for comparison).
#' @param window_size markers per region when region_by = "window".
#' @param alpha,flank passed to call_signal_loci() when region_by = "signal":
#'   Bonferroni significance level and per-locus half-width (markers).
#' @param prune_r2 LD threshold for tag collapsing; NULL disables it.
#' @param max_markers per-region marker cap (top |z| kept).
#' @param n_causal max causal variants (FINEMAP --n-causal-max / CAVIAR -c).
#' @param coverage credible-set coverage.
#' @param work_dir scratch root for the tools' files.
#' @return list(pip_df, cs_df, meta)
finemap_caviar_condition <- function(df_cond, condition, dataset,
                                     pheno_col = "Phenotype",
                                     methods = c("finemap", "caviar"),
                                     region_by = "signal",
                                     window_size = 100,
                                     alpha = 0.05,
                                     flank = 50,
                                     prune_r2 = 0.99,
                                     max_markers = 400,
                                     n_causal = 1,
                                     coverage = 0.95,
                                     work_dir = FM_WORK_ROOT) {
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

  # Genome-wide Bonferroni z threshold (reported per locus in the meta table).
  z_thr <- stats::qnorm(1 - alpha / (2 * length(z)))

  regions <- if (region_by == "signal") {
    call_signal_loci(feats, z, alpha = alpha, flank = flank)
  } else {
    make_regions(feats, by = region_by, window_size = window_size)
  }
  region_ids <- sort(unique(regions))      # NA (markers outside any locus) dropped

  cond_tag <- gsub("[^A-Za-z0-9._-]+", "_", condition)
  pip_list <- list(); cs_list <- list(); meta_list <- list()

  for (rg in region_ids) {
    idx <- which(regions == rg)
    if (length(idx) < 2) next            # nothing to fine-map

    Xr <- X[, idx, drop = FALSE]
    zr <- z[idx]
    fr <- feats[idx]

    # Per-locus signal summary, captured on the RAW region (pre-prune/-cap).
    lead_i      <- which.max(abs(zr))
    lead_gene   <- fr[lead_i]
    lead_z      <- abs(zr[lead_i])
    lead_lod    <- lead_z^2 / (2 * log(10))
    n_sig       <- sum(abs(zr) >= z_thr)
    locus_width <- length(idx)

    Rr <- Rfast::cora(Xr, large = TRUE)
    dimnames(Rr) <- list(fr, fr)

    # --- collapse near-perfect LD to tag markers -----------------------------
    n_before <- length(fr)
    tag_map <- tibble(tag = fr, member = fr)
    if (!is.null(prune_r2)) {
      tag_map <- ld_tag_groups(Rr, zr, thr = prune_r2)
      tags <- unique(tag_map$tag)
      sel <- match(tags, fr)
      Rr <- Rr[sel, sel, drop = FALSE]
      zr <- zr[sel]
      fr <- fr[sel]
    }
    n_after_prune <- length(fr)

    # --- cap region size ----------------------------------------------------
    n_dropped_cap <- 0L
    if (length(fr) > max_markers) {
      sel <- order(abs(zr), decreasing = TRUE)[seq_len(max_markers)]
      sel <- sort(sel)
      n_dropped_cap <- length(fr) - length(sel)
      Rr <- Rr[sel, sel, drop = FALSE]
      zr <- zr[sel]
      fr <- fr[sel]
    }

    if (length(fr) < 2) next
    diag(Rr) <- 1
    tab <- tibble(snp = fr, zscore = zr)

    for (mth in methods) {
      dir_run <- file.path(work_dir, mth, dataset, cond_tag, rg)
      dir.create(dir_run, recursive = TRUE, showWarnings = FALSE)

      res <- switch(
        mth,
        finemap = finemap_region(tab, Rr, n, dir_run, n_causal_max = n_causal),
        caviar  = caviar_region(tab, Rr, dir_run, n_causal = n_causal)
      )
      if (is.null(res)) next

      res <- credible_set_from_pip(res, coverage = coverage) %>%
        mutate(Condition = condition, Method = mth, region = rg)

      # Carry the collapsed LD partners so a credible set can be read as the
      # full set of implicated genes, not just the tag markers.
      res <- res %>%
        left_join(
          tag_map %>% group_by(tag) %>%
            summarise(tagged_markers = paste(member, collapse = "|"),
                      n_tagged = n(), .groups = "drop"),
          by = c("Gene" = "tag")
        ) %>%
        mutate(tagged_markers = coalesce(tagged_markers, Gene),
               n_tagged = coalesce(n_tagged, 1L))

      pip_list[[length(pip_list) + 1]] <- res
      cs_list[[length(cs_list) + 1]] <- filter(res, in_cs)
      meta_list[[length(meta_list) + 1]] <- tibble(
        Condition = condition, Method = mth, region = rg, n = n,
        lead_gene = lead_gene, lead_z = lead_z, lead_lod = lead_lod,
        n_sig = n_sig, locus_width = locus_width,
        n_markers_raw = n_before, n_markers_after_prune = n_after_prune,
        n_dropped_by_cap = n_dropped_cap, n_markers_used = length(fr),
        n_cs = sum(res$in_cs)
      )
    }
    rm(Xr, Rr); gc(verbose = FALSE)
  }

  list(
    pip_df = bind_rows(pip_list),
    cs_df  = bind_rows(cs_list),
    meta   = bind_rows(meta_list),
    # Per-condition scan summary (method-independent): how many loci the signal
    # threshold called and the strongest marginal signal seen.
    scan   = tibble(
      Condition = condition, region_by = region_by,
      max_abs_z = if (length(z)) max(abs(z)) else NA_real_,
      z_thr = z_thr, n_loci = length(region_ids)
    )
  )
}

# Driver: fine-map every condition for a dataset

#' @param dataname feather stem, e.g. "bloom2013".
#' @param shap_name capitalised dataset name used in results/SHAP paths.
#' @param conditions optional subset of conditions (NULL = all).
#' @param ... passed through to finemap_caviar_condition().
run_conditionwise_finemapping <- function(dataname, shap_name,
                                          conditions = NULL,
                                          methods = c("finemap", "caviar"),
                                          ...) {
  df <- read_feather(sprintf(
    "data/training_data/0.5_sigma/%s_clf.feather", dataname
  )) %>%
    select(Strain, Condition, Phenotype, starts_with("Y"))   # drop latents

  if (is.null(conditions)) conditions <- sort(unique(df$Condition))

  pip_all <- list(); cs_all <- list(); meta_all <- list(); scan_all <- list()

  for (i in seq_along(conditions)) {
    cond <- conditions[i]
    message(sprintf("[%s] %2d/%d  %s", dataname, i, length(conditions), cond))
    fit <- finemap_caviar_condition(filter(df, Condition == cond), cond,
                                    dataset = shap_name, methods = methods, ...)
    scan_all[[i]] <- fit$scan            # method-independent, collect before skip
    if (nrow(fit$pip_df) == 0) next

    for (mth in unique(fit$cs_df$Method)) {
      out_dir <- file.path(FM_RESULT_ROOT, mth, shap_name)
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      write_csv(filter(fit$cs_df, Method == mth),
                file.path(out_dir, paste0(cond, "_credible_sets.csv")))
    }
    pip_all[[i]] <- fit$pip_df
    cs_all[[i]]  <- fit$cs_df
    meta_all[[i]] <- fit$meta
  }
  rm(df); gc(verbose = FALSE)

  combined_pip  <- bind_rows(pip_all)
  combined_cs   <- bind_rows(cs_all)
  combined_meta <- bind_rows(meta_all)
  combined_scan <- bind_rows(scan_all)

  for (mth in unique(combined_pip$Method)) {
    dir.create(file.path(FM_RESULT_ROOT, mth), recursive = TRUE,
               showWarnings = FALSE)
    write_csv(filter(combined_pip, Method == mth),
              sprintf("%s/%s/%s_conditionwise_pip.csv",
                      FM_RESULT_ROOT, mth, shap_name))
    write_csv(filter(combined_meta, Method == mth),
              sprintf("%s/%s/%s_conditionwise_meta.csv",
                      FM_RESULT_ROOT, mth, shap_name))
  }

  # Locus scan is method-independent (loci are chosen before either tool runs),
  # so write it once at the dataset level. Record conditions with no locus so a
  # missing condition is visible rather than silently absent from the outputs.
  dir.create(FM_RESULT_ROOT, recursive = TRUE, showWarnings = FALSE)
  if (nrow(combined_scan) > 0) {
    write_csv(combined_scan,
              sprintf("%s/%s_locus_scan.csv", FM_RESULT_ROOT, shap_name))
    no_signal <- filter(combined_scan, n_loci == 0)
    write_csv(no_signal,
              sprintf("%s/%s_no_signal_conditions.csv", FM_RESULT_ROOT, shap_name))
    if (nrow(no_signal) > 0) {
      message(sprintf("[%s] %d/%d conditions had no genome-wide-significant locus: %s",
                      shap_name, nrow(no_signal), nrow(combined_scan),
                      paste(no_signal$Condition, collapse = ", ")))
    }
  }

  list(pip = combined_pip, cs = combined_cs, meta = combined_meta,
       scan = combined_scan)
}

# Comparison with conditionwise SHAP

#' Load conditionwise SHAP, keeping only positive-SHAP gene (Y*) features and
#' harmonising the lowercase `condition` column to `Condition`.
#' (Same loader as the SuSiE script.)
load_conditionwise_shap <- function(shap_name) {
  read_parquet(sprintf(
    "data/shap/sigmas/shap_classification_0.5/%s_conditionwise_shap.parquet",
    shap_name
  )) %>%
    filter(Value > 0, grepl("^Y", Feature)) %>%
    rename(Gene = Feature, SHAP = Value, Condition = condition)
}

#' Compare one method's credible sets / PIPs against conditionwise SHAP.
#'
#' Mirrors compare_susie_shap(): (i) Spearman(PIP, SHAP) table + gt/plot,
#' (ii) per-condition PIP-vs-SHAP scatter coloured by credible-set membership,
#' (iii) per-condition Venn of credible-set genes vs top-SHAP genes, and
#' (iv) an overlap summary.
#'
#' @param fm_res result list from run_conditionwise_finemapping().
#' @param shap_name capitalised dataset name (e.g. "Bloom2013").
#' @param method "finemap" or "caviar".
#' @param top_shap_n number of top-SHAP genes per condition for the Venn.
compare_finemapping_shap <- function(fm_res, shap_name, method,
                                     top_shap_n = 50) {
  shap_file <- sprintf(
    "data/shap/sigmas/shap_classification_0.5/%s_conditionwise_shap.parquet",
    shap_name
  )
  if (!file.exists(shap_file)) {
    message(sprintf("No conditionwise SHAP for %s - skipping comparison.",
                    shap_name))
    return(invisible(NULL))
  }

  pip <- filter(fm_res$pip, Method == method)
  cs  <- filter(fm_res$cs,  Method == method)
  if (nrow(pip) == 0) {
    message("No ", method, " results for ", shap_name, " - skipping.")
    return(invisible(NULL))
  }

  shap <- load_conditionwise_shap(shap_name)
  cmp_dir <- file.path(FM_RESULT_ROOT, method, shap_name, "comparison")
  dir.create(file.path(cmp_dir, "scatter"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(cmp_dir, "venn"),    recursive = TRUE, showWarnings = FALSE)

  joined <- pip %>% inner_join(shap, by = c("Gene", "Condition"))

  conditions <- sort(unique(joined$Condition))
  cor_rows <- list(); overlap_rows <- list()

  for (cond in conditions) {
    jc <- filter(joined, Condition == cond)
    cs_genes <- cs %>% filter(Condition == cond) %>% pull(Gene) %>% unique()

    # (i) Spearman correlation between PIP and SHAP.
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
      labs(title = paste(shap_name, "-", toupper(method), "-", cond),
           x = paste(toupper(method), "PIP"), y = "SHAP value", color = NULL)
    save_plot(file.path(cmp_dir, "scatter", paste0(cond, "_pip_vs_shap.svg")),
              plot = p_sc, width = 15, height = 10, units = "cm")

    # (iii) Venn: credible-set genes vs top-SHAP genes.
    top_shap <- shap %>% filter(Condition == cond) %>%
      slice_max(SHAP, n = top_shap_n) %>% pull(Gene)
    if (length(cs_genes) > 0) {
      venn_sets <- list(`Credible sets` = cs_genes, `Top SHAP` = top_shap)
      p_v <- ggvenn(venn_sets, auto_scale = TRUE, show_percentage = FALSE)
      save_plot(file.path(cmp_dir, "venn", paste0(cond, "_venn.svg")),
                plot = p_v, width = 12, height = 8, units = "cm")
    }

    # (iv) Overlap metrics: where do credible-set genes sit in the SHAP ranking?
    cs_with_shap <- jc %>% filter(in_cs)
    shap_pct <- if (nrow(jc) >= 3 && nrow(cs_with_shap) > 0) {
      pr <- rank(jc$SHAP) / nrow(jc)
      median(pr[jc$in_cs])
    } else NA_real_
    overlap_rows[[cond]] <- tibble(
      Condition          = cond,
      n_cs_genes         = length(cs_genes),
      n_cs_genes_scored  = nrow(cs_with_shap),
      n_top_shap_hit     = length(intersect(cs_genes, top_shap)),
      median_shap_pctile = shap_pct
    )
  }

  cor_df <- bind_rows(cor_rows)
  overlap_df <- bind_rows(overlap_rows)
  write_csv(cor_df, file.path(cmp_dir, "pip_shap_spearman.csv"))
  write_csv(overlap_df, file.path(cmp_dir, "pip_shap_overlap_summary.csv"))

  tryCatch({
    gt_tbl <- cor_df %>% arrange(desc(spearman)) %>%
      gt() %>% fmt_number(columns = spearman, decimals = 3) %>%
      tab_header(title = paste(shap_name, "-", toupper(method),
                               "PIP vs SHAP (Spearman)"))
    gtsave(gt_tbl, file.path(cmp_dir, "pip_shap_spearman.html"))
  }, error = function(e) message("gt table skipped: ", conditionMessage(e)))

  p_bar <- cor_df %>% filter(!is.na(spearman)) %>%
    ggplot(aes(x = reorder(Condition, spearman), y = spearman)) +
    geom_col(fill = "#0072B2") + coord_flip() +
    labs(x = NULL, y = "Spearman(PIP, SHAP)",
         title = paste(shap_name, "-", toupper(method),
                       "per-condition PIP/SHAP rank correlation"))
  save_plot(file.path(cmp_dir, "pip_shap_spearman_bar.svg"),
            plot = p_bar, width = 16, height = 20, units = "cm")

  list(spearman = cor_df, overlap = overlap_df, joined = joined)
}

#' Cross-reference credible-set genes with previously detected QTL genes.
#' Only available for Bloom2013 (data/qtl/detected_qtl_bloom2013.csv).
compare_finemapping_qtl_bloom2013 <- function(fm_res, method,
                                              shap_name = "Bloom2013") {
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

  cs_by_cond <- fm_res$cs %>% filter(Method == method) %>%
    group_by(Condition) %>% summarise(cs_genes = list(unique(Gene)), .groups = "drop")

  summ <- cs_by_cond %>%
    left_join(qtl_genes, by = c("Condition" = "Trait")) %>%
    mutate(
      n_cs_genes  = map_int(cs_genes, length),
      n_qtl_genes = map_int(qtl_genes, ~ if (is.null(.x)) 0L else length(.x)),
      n_overlap   = map2_int(cs_genes, qtl_genes,
                             ~ if (is.null(.y)) 0L else length(intersect(.x, .y)))
    ) %>%
    select(Condition, n_cs_genes, n_qtl_genes, n_overlap)

  out_dir <- file.path(FM_RESULT_ROOT, method, shap_name, "comparison")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write_csv(summ, file.path(out_dir, "qtl_overlap_summary.csv"))
  summ
}

# Per-condition triangulation (credible-set members across evidence lines)

# Same idea as the SuSiE script: lay each credible set's member genes side by
# side across independent lines of evidence - PIP, marginal z/beta, MAF,
# conditionwise SHAP (value + rank) and QTL membership - to argue which gene is
# the likely causal one.

#' Marginal logistic-regression stats (beta, se, z, maf) for a few genes.
#'
#' Markers are collapsed to 0/1 via `x >= 1` (some Bloom markers are dosage
#' coded 0/1/2, e.g. YNL087W) to match the genotype coding used in the
#' QTL/genotype work (see `src/genotype/preocessing_genotype.R`). NOTE: the
#' fine-mapping itself is dosage-based, so a dosage marker's z here differs
#' slightly from the pipeline's value. Identical to the SuSiE script's helper.
marginal_stats_for_genes <- function(df_cond, genes, pheno_col = "Phenotype") {
  genes <- intersect(genes, colnames(df_cond))
  if (length(genes) == 0) {
    return(tibble(Gene = character(), beta = numeric(), se = numeric(),
                  z = numeric(), converged = logical(), maf = numeric()))
  }
  pheno <- df_cond[[pheno_col]]
  int <- rep(1, nrow(df_cond))
  map_dfr(genes, function(g) {
    x <- as.double(df_cond[[g]] >= 1)              # collapse 0/1/2 -> 0/1
    fit <- fastglm(cbind(int, x), pheno, family = binomial(), method = 2)
    p <- mean(x)
    tibble(Gene = g, beta = fit$coefficients[2], se = fit$se[2],
           z = fit$coefficients[2] / fit$se[2],
           converged = fit$converged, maf = min(p, 1 - p))
  })
}

#' Build the triangulation table for one condition's credible-set members.
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
    group_by(region) %>% mutate(lead = PIP == max(PIP)) %>% ungroup() %>%
    arrange(region, desc(PIP))
}

#' Genomic-position vs PIP plot: point size = SHAP, colour = region, top-PIP
#' gene per region labelled. (Analog of plot_position_pip in the SuSiE script.)
plot_position_pip <- function(tri_df, title) {
  tri_df <- tri_df %>% mutate(SHAP_size = coalesce(SHAP, 0))
  ggplot(tri_df, aes(x = genomic_pos, y = PIP)) +
    geom_point(aes(size = SHAP_size, color = region), alpha = 0.75) +
    ggrepel::geom_text_repel(data = filter(tri_df, lead),
                             aes(label = Gene), size = 2.8, max.overlaps = 20) +
    scale_size(range = c(1, 6)) +
    labs(x = "Genomic position (relative)", y = "PIP",
         color = "Region", size = "SHAP", title = title)
}

#' Triangulate every condition for one dataset/method.
run_triangulation <- function(dataname, shap_name, method) {
  out_dir <- file.path(FM_RESULT_ROOT, method, shap_name)
  cs_files <- list.files(out_dir, pattern = "_credible_sets\\.csv$",
                         full.names = TRUE)
  if (length(cs_files) == 0) {
    message("No credible-set CSVs for ", shap_name, "/", method,
            " - run run_conditionwise_finemapping() first.")
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
              plot = plot_position_pip(
                tri, paste(shap_name, "-", toupper(method), "-", cond)),
              width = 18, height = 10, units = "cm")
    tri_all[[i]] <- tri
  }
  rm(df); gc(verbose = FALSE)

  combined <- bind_rows(tri_all)
  write_csv(combined, file.path(out_dir, "triangulation_all.csv"))
  invisible(combined)
}

# Cross-method comparison (FINEMAP vs CAVIAR vs SuSiE)

#' Compare per-gene PIPs across FINEMAP, CAVIAR and (if present on disk) SuSiE.
#'
#' Reads the SuSiE PIPs written by susie_conditionwise_analysis.R
#' (results/result_susie/<shap_name>_conditionwise_pip.csv) so the three
#' methods can be put on the same axes without re-running SuSiE.
#'
#' @return list(long = per-gene PIPs by method, cor = pairwise Spearman by
#'   condition, jaccard = credible-set Jaccard by condition)
compare_methods <- function(fm_res, shap_name) {
  susie_file <- sprintf("results/result_susie/%s_conditionwise_pip.csv",
                        shap_name)

  long <- fm_res$pip %>% select(Condition, Gene, Method, PIP, in_cs)

  if (file.exists(susie_file)) {
    susie <- read_csv(susie_file, show_col_types = FALSE) %>%
      transmute(Condition, Gene, Method = "susie", PIP,
                in_cs = !is.na(cs_id))
    long <- bind_rows(long, susie)
  } else {
    message("No SuSiE PIP file at ", susie_file,
            " - comparing FINEMAP vs CAVIAR only.")
  }

  out_dir <- file.path(FM_RESULT_ROOT, "method_comparison", shap_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write_csv(long, file.path(out_dir, "pip_by_method_long.csv"))

  methods <- sort(unique(long$Method))
  pairs <- if (length(methods) >= 2) combn(methods, 2, simplify = FALSE) else list()

  cor_rows <- list(); jac_rows <- list()
  for (cond in sort(unique(long$Condition))) {
    lc <- filter(long, Condition == cond)
    for (pr in pairs) {
      a <- filter(lc, Method == pr[1]); b <- filter(lc, Method == pr[2])
      j <- inner_join(select(a, Gene, PIP_a = PIP),
                      select(b, Gene, PIP_b = PIP), by = "Gene")
      rho <- if (nrow(j) >= 3) {
        suppressWarnings(cor(j$PIP_a, j$PIP_b, method = "spearman"))
      } else NA_real_
      cor_rows[[length(cor_rows) + 1]] <- tibble(
        Condition = cond, method_a = pr[1], method_b = pr[2],
        n_common = nrow(j), spearman = rho
      )

      ga <- a %>% filter(in_cs) %>% pull(Gene) %>% unique()
      gb <- b %>% filter(in_cs) %>% pull(Gene) %>% unique()
      uni <- length(union(ga, gb))
      jac_rows[[length(jac_rows) + 1]] <- tibble(
        Condition = cond, method_a = pr[1], method_b = pr[2],
        n_cs_a = length(ga), n_cs_b = length(gb),
        n_overlap = length(intersect(ga, gb)),
        jaccard = if (uni > 0) length(intersect(ga, gb)) / uni else NA_real_
      )
    }
  }

  cor_df <- bind_rows(cor_rows); jac_df <- bind_rows(jac_rows)
  write_csv(cor_df, file.path(out_dir, "pip_spearman_between_methods.csv"))
  write_csv(jac_df, file.path(out_dir, "credible_set_jaccard.csv"))

  if (nrow(cor_df) > 0) {
    p <- cor_df %>% filter(!is.na(spearman)) %>%
      mutate(pair = paste(method_a, "vs", method_b)) %>%
      ggplot(aes(x = pair, y = spearman)) +
      geom_boxplot(outlier.alpha = 0.4, fill = "#56B4E9") +
      labs(x = NULL, y = "Spearman(PIP, PIP)",
           title = paste(shap_name, "- agreement between fine-mappers"))
    save_plot(file.path(out_dir, "pip_spearman_between_methods.svg"),
              plot = p, width = 14, height = 10, units = "cm")
  }

  list(long = long, cor = cor_df, jaccard = jac_df)
}

# Full pipeline runner

# Same three datasets as the SuSiE script. Bloom2019 = the base BY x RM cross,
# NOT the _BYxM22 / _RMxYPS163 crosses. Bloom2015 has no conditionwise SHAP
# file, so it is fine-mapped but its SHAP comparison is skipped.
DATASETS <- tibble::tribble(
  ~dataname,     ~shap_name,
  "bloom2013",   "Bloom2013",
  "bloom2015",   "Bloom2015",
  "bloom2019",   "Bloom2019"
)

#' Run the full condition-wise FINEMAP/CAVIAR + SHAP comparison over datasets.
#'
#' @param conditions optional condition subset (for quick tests); NULL = all.
#' @param methods which fine-mappers to run.
#' @param ... passed to finemap_caviar_condition() (region_by, prune_r2,
#'   max_markers, n_causal, coverage, ...).
main_conditionwise_finemapping <- function(datasets = DATASETS,
                                           conditions = NULL,
                                           methods = c("finemap", "caviar"),
                                           ...) {
  setup_finemapr()

  results <- list()
  for (k in seq_len(nrow(datasets))) {
    dn <- datasets$dataname[k]; sn <- datasets$shap_name[k]
    fm_res <- run_conditionwise_finemapping(dn, sn, conditions = conditions,
                                            methods = methods, ...)
    per_method <- list()
    for (mth in methods) {
      shap_cmp <- compare_finemapping_shap(fm_res, sn, mth)
      qtl_cmp  <- if (sn == "Bloom2013")
        compare_finemapping_qtl_bloom2013(fm_res, mth, sn) else NULL
      tri_cmp  <- run_triangulation(dn, sn, mth)
      per_method[[mth]] <- list(shap = shap_cmp, qtl = qtl_cmp,
                                triangulation = tri_cmp)
    }
    method_cmp <- compare_methods(fm_res, sn)
    results[[sn]] <- list(finemapping = fm_res, per_method = per_method,
                          method_comparison = method_cmp)
  }
  invisible(results)
}

# Usage

# One-off setup (in WSL, from the project root):
#   wsl -d Ubuntu -- bash src/finemapping/install_finemap_caviar.sh
#
# Quick subset on one dataset (recommended first run - CAVIAR is the slow one):
#   setup_finemapr()
#   res <- run_conditionwise_finemapping("bloom2013", "Bloom2013",
#                                        conditions = c("4NQO", "maltose"))
#   cmp <- compare_finemapping_shap(res, "Bloom2013", "finemap")
#   qtl <- compare_finemapping_qtl_bloom2013(res, "finemap")
#   mth <- compare_methods(res, "Bloom2013")
#
# FINEMAP only (much faster than CAVIAR):
#   res <- run_conditionwise_finemapping("bloom2013", "Bloom2013",
#                                        methods = "finemap")
#
# Full pipeline (all datasets, all conditions, both tools - slow):
#   results <- main_conditionwise_finemapping()
