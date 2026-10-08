# ---------------------------------------------------------------------------
# finemapr / FINEMAP / CAVIAR environment setup
# ---------------------------------------------------------------------------
# Sourced by finemap_caviar_conditionwise_analysis.R. Responsibilities:
#   1. install the (patched) finemapr package from tools/finemapr-src if absent
#   2. point finemapr's options() at the locally installed binaries
#   3. fail loudly and early if the binaries are missing
#
# Run src/finemapping/install_finemap_caviar.sh (in WSL) first; it fetches the
# binaries and the finemapr source into tools/ (gitignored).
#
# Platform note: on Windows the options point at .cmd shims that forward into
# WSL, because FINEMAP/CAVIAR are Linux binaries. On Linux/macOS they point at
# the binaries directly.

FINEMAPR_TOOLS_DIR <- file.path(normalizePath(".", winslash = "/"), "tools")

#' Ensure the finemapr package is installed (from the local patched clone).
ensure_finemapr <- function(src_dir = file.path(FINEMAPR_TOOLS_DIR, "finemapr-src")) {
  if (requireNamespace("finemapr", quietly = TRUE)) return(invisible(TRUE))

  if (!dir.exists(src_dir)) {
    stop("finemapr is not installed and its source is missing at ", src_dir,
         "\nRun src/finemapping/install_finemap_caviar.sh first.", call. = FALSE)
  }
  message("Installing finemapr from ", src_dir, " ...")
  install.packages(src_dir, repos = NULL, type = "source")
  invisible(TRUE)
}

#' Resolve the path to a tool, preferring the Windows .cmd shim when on Windows.
finemapr_tool_path <- function(tool = c("finemap", "CAVIAR"),
                               tools_dir = FINEMAPR_TOOLS_DIR) {
  tool <- match.arg(tool)
  bin <- file.path(tools_dir, "bin", tool)
  if (.Platform$OS.type == "windows") {
    shim <- paste0(bin, ".cmd")
    if (file.exists(shim)) return(shim)
  }
  bin
}

#' Set finemapr's options() and verify the tools are actually present.
#'
#' @param check if TRUE (default) stop when a tool is missing.
setup_finemapr <- function(check = TRUE) {
  ensure_finemapr()

  finemap_path <- finemapr_tool_path("finemap")
  caviar_path  <- finemapr_tool_path("CAVIAR")

  options(
    finemapr_finemap = finemap_path,
    finemapr_caviar  = caviar_path
  )

  if (check) {
    missing <- c(finemap = finemap_path, caviar = caviar_path)
    missing <- missing[!file.exists(missing)]
    if (length(missing) > 0) {
      stop("Fine-mapping tool(s) not found:\n",
           paste(sprintf("  %s -> %s", names(missing), missing), collapse = "\n"),
           "\nRun src/finemapping/install_finemap_caviar.sh (in WSL) first.",
           call. = FALSE)
    }
  }

  message("finemapr configured:")
  message("  FINEMAP : ", finemap_path)
  message("  CAVIAR  : ", caviar_path)
  invisible(list(finemap = finemap_path, caviar = caviar_path))
}
