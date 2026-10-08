# =============================================================================
# Is the prominence of pleiotropic genes in Figure 3A a consequence of joint
# training, or is it present in the data independently of the model?
#
# Figure 3A plots the global SHAP value of a model trained across all conditions
# simultaneously against an odds ratio pooled over all conditions. Both axes
# aggregate across conditions, so a gene acting in many conditions is favoured
# by both. The condition-wise contingency results give effect sizes estimated
# inside a single condition, where breadth of action cannot contribute, and so
# provide the control the figure itself cannot.
#
# Produces the supplementary table reported in the response to reviewers.
#
# Run from the project root (renv needs it):
#   "C:/Program Files/R/R-4.5.1/bin/Rscript.exe" src/Contingency/pleiotropy_effect_size_bias.R
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(gt)
})

OUT_DIR <- "results/result_contingency"
MIN_VARIANCE_EXPLAINED <- 0.05
SHAP_FLOOR <- 1e-5 # same threshold used throughout the contingency analyses

# The 39 conditions Bloom2013 actually contains. The condition-wise folder also
# holds six enrichment gene-lists (carbons/oxidative/genotoxic, _top/_bottom)
# and several Bloom2019-only conditions; including them changes the numbers.
B2013_CONDITIONS <- c(
  "4NQO", "6-azauracil", "berbamine", "CaCl2", "caffeine", "CdCl2", "cisplatin",
  "CoCl2", "congo_red", "copper", "cycloheximide", "diamide", "ethanol",
  "fluorocytosine", "fluorouracil", "formamide", "galactose", "H2O2",
  "hydroquinone", "hydroxybenzaldehyde", "hydroxyurea", "indoleacetic_acid",
  "lactate", "lactose", "LiCl", "maltose", "mannose", "menadione", "MgCl2",
  "MgSO4", "neomycin", "paraquat", "raffinose", "sds", "sorbitol", "trehalose",
  "tunicamycin", "xylose", "zeocin"
)

# =============================================================================
# Degree of pleiotropy, from the published Bloom2013 QTL mapping
# =============================================================================

# Genes are the systematic names falling inside each QTL confidence interval,
# stored pipe-separated. A gene's degree of pleiotropy is the number of distinct
# traits for which it falls inside such an interval.
qtl <- read_csv("data/qtl/detected_qtl_bloom2013.csv", show_col_types = FALSE) %>%
  select(Trait, Pheno_fraction_explained, Genes) %>%
  mutate(Genes = strsplit(as.character(Genes), split = "|", fixed = TRUE)) %>%
  unnest_longer(Genes) %>%
  filter(grepl("^Y", Genes), Pheno_fraction_explained > MIN_VARIANCE_EXPLAINED)

counts <- qtl %>%
  group_by(Genes) %>%
  summarise(n_conditions = n_distinct(Trait), .groups = "drop") %>%
  mutate(group = ifelse(n_conditions > 1, "Pleiotropic", "Single condition"))

# =============================================================================
# The three effect-size measures
# =============================================================================

# 1. Global SHAP from the model trained on all conditions jointly - the measure
#    Figure 3A places on its x axis.
shap <- read_csv(
  "data/shap/sigmas/shap_classification_0.5/shap_Bloom2013_Boosting.csv",
  show_col_types = FALSE
) %>%
  filter(Value > SHAP_FLOOR) %>%
  rename(Genes = Feature, shap = Value)

# 2. Odds ratio pooled across conditions - Figure 3A's y axis. Built by putting
#    every strain x condition row into one 2x2 table, so it is aggregated in the
#    same way the SHAP is and is not an independent control.
pooled <- read_csv(
  file.path(OUT_DIR, "contingency_bloom2013_0.5sigma_full.csv"),
  show_col_types = FALSE
) %>%
  rename(Genes = Gene) %>%
  transmute(Genes, pooled_abs_logOR = abs(log10(odds_ratio)))

# 3. Odds ratios estimated separately inside each condition.
files <- list.files(file.path(OUT_DIR, "0.5_condition-wise/Bloom2013"),
                    pattern = "csv$", full.names = TRUE)
files <- files[sub("\\.csv$", "", basename(files)) %in% B2013_CONDITIONS]
stopifnot(length(files) == length(B2013_CONDITIONS))

conditionwise <- map(files, function(f) {
  read_csv(f, show_col_types = FALSE) %>%
    mutate(Condition = sub("\\.csv$", "", basename(f)))
}) %>%
  list_rbind() %>%
  rename(Genes = Gene) %>%
  mutate(abs_logOR = abs(log10(odds_ratio)))

per_gene <- conditionwise %>%
  group_by(Genes) %>%
  summarise(
    mean_abs_logOR = mean(abs_logOR, na.rm = TRUE),
    max_abs_logOR = max(abs_logOR, na.rm = TRUE),
    n_sig_conditions = sum(adj.pval < 0.05, na.rm = TRUE),
    .groups = "drop"
  )

genes <- counts %>%
  left_join(shap, "Genes") %>%
  left_join(pooled, "Genes") %>%
  left_join(per_gene, "Genes")

# =============================================================================
# Tests
# =============================================================================

MEASURES <- tribble(
  ~column,             ~label,                                          ~block,
  "shap",              "Global SHAP value",                             "Derived from the jointly trained model",
  "pooled_abs_logOR",  "Absolute log10 odds ratio, pooled",             "Association, pooled across all conditions",
  "mean_abs_logOR",    "Mean absolute log10 odds ratio",                "Association, estimated within single conditions",
  "max_abs_logOR",     "Maximum absolute log10 odds ratio",             "Association, estimated within single conditions",
  "n_sig_conditions",  "Conditions with FDR-significant association",    "Association, estimated within single conditions"
)

summarise_measure <- function(column, label, block) {
  d <- genes %>% filter(is.finite(.data[[column]]))
  a <- d[[column]][d$group == "Pleiotropic"]
  b <- d[[column]][d$group == "Single condition"]

  mw <- suppressWarnings(wilcox.test(a, b, alternative = "two.sided"))
  sp <- suppressWarnings(cor.test(d$n_conditions, d[[column]], method = "spearman"))

  tibble(
    block = block, label = label,
    n_pleio = length(a), n_single = length(b),
    med_pleio = median(a), iqr_pleio = IQR(a),
    med_single = median(b), iqr_single = IQR(b),
    mw_p = mw$p.value,
    rho = unname(sp$estimate), rho_p = sp$p.value
  )
}

results <- pmap(MEASURES, summarise_measure) %>% list_rbind()

write_csv(results, file.path(OUT_DIR, "pleiotropy_effect_size_bias.csv"))

cat("\nGene sets\n")
cat("  QTLs explaining more than 5 percent of variance -> distinct genes:", nrow(counts), "\n")
cat("  of which pleiotropic:", sum(counts$group == "Pleiotropic"),
    " single condition:", sum(counts$group == "Single condition"), "\n")
cat("  with a SHAP value above the threshold:", sum(is.finite(genes$shap)), "\n")
cat("  with a condition-wise odds ratio:", sum(is.finite(genes$max_abs_logOR)), "\n\n")
print(as.data.frame(results %>% select(label, n_pleio, n_single, med_pleio, med_single, mw_p, rho, rho_p)),
      digits = 3)

# =============================================================================
# Supplementary table
# =============================================================================

supp_gt <- results %>%
  select(block, label, n_pleio, med_pleio, iqr_pleio, n_single, med_single, iqr_single,
         mw_p, rho, rho_p) %>%
  gt(rowname_col = "label", groupname_col = "block") %>%
  fmt_number(columns = c(med_pleio, iqr_pleio, med_single, iqr_single, rho), decimals = 3) %>%
  # the last measure is a count of conditions, not a continuous effect size
  fmt_number(
    columns = c(med_pleio, iqr_pleio, med_single, iqr_single),
    rows = label == "Conditions with FDR-significant association",
    decimals = 0
  ) %>%
  fmt_scientific(columns = c(mw_p, rho_p), decimals = 1) %>%
  cols_merge(columns = c(med_pleio, iqr_pleio), pattern = "{1} ({2})") %>%
  cols_merge(columns = c(med_single, iqr_single), pattern = "{1} ({2})") %>%
  cols_label(
    n_pleio = "n", med_pleio = "Median (IQR)",
    n_single = "n", med_single = "Median (IQR)",
    mw_p = "p", rho = "RHOSYMBOL", rho_p = "p"
  ) %>%
  tab_spanner("Pleiotropic genes", columns = c(n_pleio, med_pleio)) %>%
  tab_spanner("Single-condition genes", columns = c(n_single, med_single)) %>%
  tab_spanner("Mann-Whitney", columns = mw_p) %>%
  tab_spanner("Spearman, versus number of conditions", columns = c(rho, rho_p)) %>%
  tab_stubhead("Effect-size measure") %>%
  cols_align(align = "center", columns = -label) %>%
  tab_footnote(
    "Two-sided Mann-Whitney U test comparing pleiotropic against single-condition genes.",
    locations = cells_column_spanners(spanners = "Mann-Whitney")
  ) %>%
  tab_footnote(
    "Spearman rank correlation between the measure and the number of conditions in which the gene is a QTL, across both groups combined.",
    locations = cells_column_spanners(spanners = "Spearman, versus number of conditions")
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

CAPTION <- paste(
  "Effect size versus degree of pleiotropy for Bloom2013 QTL genes.",
  "Genes were taken from the published Bloom2013 QTL mapping",
  "(591 QTLs across 46 traits), restricted to QTLs explaining more than 5 percent",
  "of phenotypic variance and to systematic gene names lying within the QTL",
  "confidence interval, giving 381 distinct genes. A gene was called pleiotropic",
  "if it fell within such an interval for more than one condition (139 genes) and",
  "single-condition otherwise (242 genes). Of these, 155 genes carried a SHAP",
  "value above the threshold used throughout (53 pleiotropic, 102 single-condition)",
  "and 154 also had a condition-wise odds ratio (53 pleiotropic, 101",
  "single-condition); group sizes therefore differ by one between the SHAP row and",
  "the odds-ratio rows. The first two measures aggregate information across",
  "conditions and so favour genes acting in many conditions by construction; the",
  "last three are estimated within individual conditions, across the 39 conditions",
  "of the Bloom2013 panel, and cannot be inflated by breadth of action. Both",
  "groups take their maximum over the same 39 conditions, so that comparison is",
  "symmetric."
)

# gt emits no caption of its own; insert one directly after the opening
# environment, and swap the unicode gt writes for math macros so the table also
# compiles under pdfLaTeX rather than only XeLaTeX or LuaLaTeX.
as_tex <- function(gt_tbl, caption) {
  lines <- strsplit(as.character(as_latex(gt_tbl)), "\n", fixed = TRUE)[[1]]
  at <- grep("\\begin{table}", lines, fixed = TRUE)[1]
  stopifnot(!is.na(at))
  lines <- append(lines, paste0("\\caption{", caption, "}"), after = at)
  x <- paste(lines, collapse = "\n")
  x <- gsub("±", "$\\pm$", x, fixed = TRUE)
  x <- gsub("Δ", "$\\Delta$", x, fixed = TRUE)
  # placeholder so the Greek letter reaches LaTeX as a math macro rather than
  # as a unicode character gt would escape
  gsub("RHOSYMBOL", "$\\rho$", x, fixed = TRUE)
}

writeLines(as_tex(supp_gt, CAPTION), file.path(OUT_DIR, "pleiotropy_effect_size_bias.tex"))

cat("\nWrote", file.path(OUT_DIR, "pleiotropy_effect_size_bias.csv"), "and .tex\n")
