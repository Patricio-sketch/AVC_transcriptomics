# =============================================================================
# Package Setup — Ischemic Stroke Transcriptomics Pipeline
# =============================================================================
#
# Run this script once on any new machine before executing script.R.
# It installs all R packages required by the pipeline from Bioconductor
# and CRAN, then verifies that every package loaded successfully.
#
# Usage:
#   Rscript configurar_ambiente.R
#   or in RStudio: source("configurar_ambiente.R")
# =============================================================================

options(repos = c(CRAN = "https://cloud.r-project.org"))

# BiocManager is required to install Bioconductor packages (GEOquery, limma,
# STRINGdb). It also handles CRAN packages, so a single call covers everything.
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

packages <- c(
  "GEOquery",   # GEO database access and expression matrix download
  "limma",      # Linear models for differential expression (microarray)
  "STRINGdb",   # STRING protein interaction database and enrichment
  "ggplot2",    # Grammar-of-graphics plotting
  "ggrepel",    # Non-overlapping text labels for ggplot2
  "dplyr",      # Data manipulation and pipelines
  "igraph",     # Graph construction and network analysis
  "ggraph",     # Graph visualization using ggplot2 grammar
  "rio",        # Unified data import/export (TSV, CSV, RDS)
  "stringr"     # String manipulation utilities
)

BiocManager::install(packages, ask = FALSE, update = FALSE)

# rstudioapi is only needed when running interactively in RStudio to detect
# the source file location. It is not required for Rscript execution.
if (!requireNamespace("rstudioapi", quietly = TRUE)) {
  install.packages("rstudioapi")
}

# Verify that all packages are available after installation
failed <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]

if (length(failed) > 0) {
  stop(
    "The following packages could not be installed: ",
    paste(failed, collapse = ", "), "\n",
    "Check your internet connection or Bioconductor version compatibility."
  )
}

message("All packages installed successfully: ", paste(packages, collapse = ", "))
message("You can now run script.R.")

# ── Optional: snapshot environment with renv ──────────────────────────────────
# Uncomment the lines below to record exact package versions for reproducibility.
# Run once after all packages are confirmed working.
#
# if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
# renv::init(bare = TRUE)
# renv::snapshot()
# message("renv.lock created. Commit renv.lock and renv/ to version control.")
