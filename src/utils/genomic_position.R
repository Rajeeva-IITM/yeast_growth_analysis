# Shared helper: map a systematic yeast ORF name to a relative genomic position.
# Used by SuSiE triangulation and the QTL/SHAP ranking analysis.

library(purrr)

#' Relative genomic position from a systematic ORF name (e.g. "YNL085W").
#'
#' Encodes chromosome (letter 2 -> number), arm (letter 3: R = +0.3, else -0.3)
#' and position on the arm (digits 4-6) into a single ordered coordinate.
get_genomic_position <- function(genename) {
  genename_vec <- strsplit(genename, "", TRUE)[[1]]           # every character
  chrom_num <- match(genename_vec[2], LETTERS)                # chromosome number
  pos_on_chrom_arm <- purrr::reduce(genename_vec[4:6], paste0) %>% as.integer
  arm_factor <- ifelse(genename_vec[3] == "R", 0.3, -0.3)
  chrom_num + pos_on_chrom_arm * 1e-4 + arm_factor
}
