# =============================================================================
# Figure 3A recreated for Bloom2015 and Bloom2019.
#
# The original panel (Bloom2013, in shap_qtl_analysis.R) plots log(global SHAP)
# against -log10(pooled odds ratio) for genes that are QTLs in more than one
# condition, sized by the number of conditions and coloured by whether the
# contingency test is significant after FDR correction.
#
# Bloom2013 takes its pleiotropy count from data/qtl/detected_qtl_bloom2013.csv,
# which does not exist for the other two panels. Each has its own published
# equivalent instead, so the count comes from a source independent of both plot
# axes in all three cases:
#   Bloom2015 - data/papers/bloom2015_detected_qtls.csv, one gene per QTL peak
#   Bloom2019 - data/papers/causal_genes_bloom2019.csv, fine-mapped causal genes
#
# Run from the project root (renv needs it):
#   "C:/Program Files/R/R-4.5.1/bin/Rscript.exe" src/Contingency/pleiotropic_odds_vs_shap_2015_2019.R
# =============================================================================

suppressPackageStartupMessages({
  # annotation packages first: AnnotationDbi and S4Vectors mask select, rename
  # and filter, so dplyr is loaded after them to win the conflicts
  library(org.Sc.sgd.db)
  library(AnnotationDbi)
  library(readr)
  library(dplyr)
  library(tibble)
  library(purrr)
  library(ggplot2)
})

source("./src/utils/theme_set.R")

OUT_DIR <- "results/result_contingency"
SHAP_FLOOR <- 1e-5            # as used throughout the contingency analyses
MIN_VARIANCE_EXPLAINED <- 0.05 # as used for Bloom2013

# =============================================================================
# Pleiotropy counts
# =============================================================================

# Bloom2015 names only the single peak gene per QTL, so counting those directly
# yields almost no recurrence (747 genes across 797 QTLs, 3 usable points).
# Bloom2013's counts instead come from every gene inside the QTL confidence
# interval, so the interval is reconstructed here from the refined QTL bounds
# and SGD gene coordinates to make the two panels comparable.
gene_coords <- local({
  loc <- AnnotationDbi::as.list(org.Sc.sgdCHRLOC)
  end <- AnnotationDbi::as.list(org.Sc.sgdCHRLOCEND)
  keep <- names(loc)[grepl("^Y", names(loc))]
  map(keep, function(g) {
    a <- loc[[g]]
    b <- end[[g]]
    if (is.null(a) || is.null(b) || all(is.na(a))) return(NULL)
    # CHRLOC is signed to encode strand; position is the magnitude
    tibble(Gene = g, chr_num = names(a)[1],
           start = min(abs(a[1]), abs(b[1])), end = max(abs(a[1]), abs(b[1])))
  }) %>% compact() %>% list_rbind() %>% filter(!is.na(chr_num))
})

ROMAN <- setNames(
  as.character(1:16),
  paste0("chr", c("I", "II", "III", "IV", "V", "VI", "VII", "VIII",
                  "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI"))
)

qtl_2015 <- read_csv("data/papers/bloom2015_detected_qtls.csv", show_col_types = FALSE) %>%
  filter(variance_explained > MIN_VARIANCE_EXPLAINED, !is.na(pos.left_refineqtl)) %>%
  mutate(chr_num = unname(ROMAN[chr])) %>%
  filter(!is.na(chr_num))

counts_2015 <- qtl_2015 %>%
  # AnnotationDbi masks dplyr::select
  dplyr::select(trait, chr_num, left = pos.left_refineqtl, right = pos.right_refineqtl) %>%
  inner_join(gene_coords, "chr_num", relationship = "many-to-many") %>%
  filter(start <= right, end >= left) %>%
  group_by(Gene) %>%
  summarise(Count = n_distinct(trait), .groups = "drop")

# Bloom2019: trait strings carry concentration and replicate, e.g.
# "Cadmium_Chloride;75uM;2". Counting them raw would score one condition tested
# at three concentrations as a pleiotropy of three, so the condition is the part
# before the first semicolon.
counts_2019 <- read_csv("data/papers/causal_genes_bloom2019.csv", show_col_types = FALSE) %>%
  rename(Gene = ORF) %>%
  filter(!is.na(Gene), grepl("^Y", Gene), !grepl("^Y", trait)) %>%
  mutate(condition = sub(";.*$", "", trait)) %>%
  group_by(Gene) %>%
  summarise(Count = n_distinct(condition), .groups = "drop")

# =============================================================================
# Panel builder
# =============================================================================

build_panel <- function(dataset, counts, shap_path, contingency_path, out_file) {
  shap <- read_csv(shap_path, show_col_types = FALSE) %>%
    filter(Value > SHAP_FLOOR) %>%
    rename(Gene = Feature)

  contingency <- read_csv(contingency_path, show_col_types = FALSE) %>%
    mutate(adj_pval = p.adjust(pval, "fdr"))

  pleiotropic <- counts %>%
    filter(Count > 1) %>%
    left_join(shap, "Gene") %>%
    left_join(contingency, "Gene") %>%
    mutate(Significance = ifelse(adj_pval > 0.05, "Not Significant", "Significant"))

  plotted <- pleiotropic %>% filter(!is.na(Value), !is.na(odds_ratio))

  message(sprintf(
    "%s: %d genes with a QTL, %d pleiotropic, %d plotted (Count %d-%d), %d significant",
    dataset, nrow(counts), sum(counts$Count > 1), nrow(plotted),
    min(plotted$Count), max(plotted$Count),
    sum(plotted$Significance == "Significant")
  ))

  p <- ggplot(plotted, aes(x = log(Value), y = -log10(odds_ratio))) +
    geom_point(aes(color = Significance, size = Count)) +
    ggrepel::geom_text_repel(mapping = aes(label = Gene), max.overlaps = 20, size = 3) +
    labs(y = "-log10(Odds Ratio)", x = "log(SHAP value)")

  save_plot(file.path(OUT_DIR, out_file), plot = p, width = 10, height = 6)
  plotted
}

panel_2015 <- build_panel(
  "Bloom2015", counts_2015,
  "data/shap/sigmas/shap_classification_0.5/shap_Bloom2015_Boosting.csv",
  file.path(OUT_DIR, "contingency_bloom2015_0.5sigma_full.csv"),
  "bloom2015_pleiotropic_odds_vs_shap.svg"
)

panel_2019 <- build_panel(
  "Bloom2019", counts_2019,
  "data/shap/sigmas/shap_classification_0.5/shap_Bloom2019_BYxRM_Boosting.csv",
  file.path(OUT_DIR, "contingency_bloom2019_0.5sigma_full.csv"),
  "bloom2019_pleiotropic_odds_vs_shap.svg"
)

write_csv(panel_2015, file.path(OUT_DIR, "bloom2015_pleiotropic_odds_vs_shap.csv"))
write_csv(panel_2019, file.path(OUT_DIR, "bloom2019_pleiotropic_odds_vs_shap.csv"))

message("Done. Figures and source data in ", OUT_DIR)
