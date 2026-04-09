# =============================================================================
# Transcriptomics Pipeline — Ischemic Stroke (GSE16561)
# =============================================================================
#
# This script processes publicly available peripheral blood gene expression
# data from ischemic stroke patients (GEO dataset GSE16561) and produces:
#
#   1. Differential expression analysis (limma)           → DEG table
#   2. Volcano plot                                        → outputs/figures/
#   3. Protein-protein interaction network (STRINGdb)      → outputs/figures/
#   4. Functional enrichment (Up and Down gene sets)       → string_enrichment/
#   5. Enrichment dotplots                                 → outputs/figures/
#
# The enrichment output feeds into clinicalbert_vias_AVC.ipynb, which maps
# the biological pathways to NANDA-I nursing diagnoses via ClinicalBERT.
#
# Run: Rscript script.R
#   or in RStudio: Session > Set Working Directory > To Source File Location,
#   then source("script.R")
# =============================================================================


# ── Working directory detection ───────────────────────────────────────────────
# Resolves the project root in three environments in order:
#   1. AVC_TRANSCRIPTOMICS_DIR environment variable (CI / Docker)
#   2. --file= argument (Rscript from command line)
#   3. RStudio active document path
set_project_wd <- function() {
  env <- Sys.getenv("AVC_TRANSCRIPTOMICS_DIR", unset = "")
  if (nzchar(env) && dir.exists(env)) {
    setwd(env)
    return(invisible(TRUE))
  }

  ca <- commandArgs(trailingOnly = FALSE)
  f  <- grep("^--file=", ca, value = TRUE)
  if (length(f)) {
    setwd(dirname(normalizePath(sub("^--file=", "", f[1]), winslash = "/")))
    return(invisible(TRUE))
  }

  if (file.exists("script.R")) {
    return(invisible(TRUE))
  }

  if (interactive() &&
      requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    p <- rstudioapi::getActiveDocumentContext()$path
    if (nzchar(p) && tolower(basename(p)) == "script.r") {
      setwd(dirname(normalizePath(p, winslash = "/")))
      return(invisible(TRUE))
    }
  }

  stop(
    "Could not determine the project directory. ",
    "In RStudio: Session > Set Working Directory > To Source File Location. ",
    "Or set the environment variable AVC_TRANSCRIPTOMICS_DIR to the folder ",
    "containing script.R."
  )
}

set_project_wd()


# =============================================================================
# SECTION 1 — Data Acquisition and Preprocessing
# =============================================================================
# We download the GSE16561 expression matrix from GEO (Illumina HumanHT-12
# platform, peripheral blood, ischemic stroke vs. healthy controls).
# The matrix is cached locally so re-runs do not require internet access.
# Expression values are log2-transformed if not already, then quantile-
# normalized to make samples comparable.
# =============================================================================

library(GEOquery)
library(limma)

dir.create("geo_cache", showWarnings = FALSE)
gset <- getGEO("GSE16561", GSEMatrix = TRUE, destdir = "geo_cache")
gset <- gset[[1]]

# Expression matrix — detect scale and normalize
expr <- exprs(gset)

if (max(expr, na.rm = TRUE) > 30) {
  message("Expression values appear to be on a linear scale. Applying log2(x + 1).")
  expr <- log2(expr + 1)
}

expr <- normalizeBetweenArrays(expr, method = "quantile")
message("Quantile normalization applied.")

cat(sprintf(
  "Matrix: %d probes × %d samples  |  min=%.2f  median=%.2f  max=%.2f\n",
  nrow(expr), ncol(expr),
  min(expr, na.rm = TRUE),
  median(expr, na.rm = TRUE),
  max(expr, na.rm = TRUE)
))

# Clinical metadata — classify samples as Stroke or Control
# GSE16561 stores group labels in the 'description' column of the phenotype table
meta <- pData(gset)

if (!"description" %in% colnames(meta)) {
  stop("Column 'description' not found in sample metadata. ",
       "Check pData(gset) for the correct group label column.")
}

grupo <- ifelse(
  grepl("Stroke",  meta$description, ignore.case = TRUE), "AVC",
  ifelse(
    grepl("Control", meta$description, ignore.case = TRUE), "Controle",
    NA
  )
)
grupo <- factor(grupo)

# Remove any samples that could not be assigned to a group
unclassified <- is.na(grupo)
if (any(unclassified)) {
  cat("Removing", sum(unclassified), "sample(s) with no group assignment.\n")
}
expr  <- expr[, !unclassified]
grupo <- droplevels(grupo[!unclassified])

cat("Sample distribution:\n")
print(table(grupo))


# =============================================================================
# SECTION 2 — Differential Expression Analysis (limma)
# =============================================================================
# We fit a linear model for each probe with empirical Bayes moderation
# (eBayes), which shrinks per-gene variance estimates toward a global prior —
# particularly useful for microarray data with moderate sample sizes.
# Contrast: AVC − Controle (positive logFC = higher in stroke).
# Significance criteria: FDR (Benjamini-Hochberg) < 0.05, |logFC| > 1.
# =============================================================================

design <- model.matrix(~0 + grupo)
colnames(design) <- c("AVC", "Controle")

fit  <- lmFit(expr, design)
contr <- makeContrasts(AVCvsControle = AVC - Controle, levels = design)
fit2 <- contrasts.fit(fit, contr)
fit2 <- eBayes(fit2)

# Extract all significant DEGs (FDR < 0.05, no logFC filter yet — applied below)
deg <- topTable(fit2, adjust = "fdr", p.value = 0.05, number = Inf)

# Annotate probes with gene symbols using the platform's feature data.
# When multiple probes map to the same gene, we keep the one with the
# largest absolute logFC (most informative probe per gene).
anotar_genes <- function(deg, gset, col_symbol = NULL, col_title = NULL,
                         colapsar = TRUE) {
  anot    <- fData(gset)
  deg$ProbeID <- rownames(deg)

  gene_col <- if (!is.null(col_symbol)) {
    col_symbol
  } else {
    grep("symbol", colnames(anot), ignore.case = TRUE, value = TRUE)[1]
  }

  title_col <- if (!is.null(col_title) && col_title %in% colnames(anot)) {
    col_title
  } else {
    NULL
  }

  anot2 <- if (!is.null(title_col)) {
    data.frame(ProbeID   = rownames(anot),
               Gene      = anot[[gene_col]],
               GeneTitle = anot[[title_col]],
               stringsAsFactors = FALSE)
  } else {
    data.frame(ProbeID = rownames(anot),
               Gene    = anot[[gene_col]],
               stringsAsFactors = FALSE)
  }

  deg_anotado <- merge(deg, anot2, by = "ProbeID", all.x = TRUE)

  if (colapsar) {
    deg_anotado <- deg_anotado[order(-abs(deg_anotado$logFC)), ]
    deg_anotado <- deg_anotado[!duplicated(deg_anotado$Gene), ]
  }

  deg_anotado
}

deg_anotado <- anotar_genes(deg, gset)

cat(sprintf(
  "DEGs: %d total  |  Up: %d  |  Down: %d\n",
  nrow(deg_anotado),
  sum(deg_anotado$logFC >  1, na.rm = TRUE),
  sum(deg_anotado$logFC < -1, na.rm = TRUE)
))


# =============================================================================
# SECTION 3 — Volcano Plot
# =============================================================================
# Visual summary of the differential expression results.
# Yellow = up-regulated in stroke; purple = down-regulated; grey = not significant.
# Dashed lines mark the |logFC| > 1 and FDR < 0.05 thresholds.
# =============================================================================

library(ggplot2)
library(ggrepel)

lfc_cutoff <- 1
fdr_cutoff <- 0.05

deg_plot <- deg_anotado
deg_plot$gene_symbol  <- deg_plot$Gene
deg_plot$volcano_class <- "NS"
deg_plot$volcano_class[deg_plot$adj.P.Val < fdr_cutoff & deg_plot$logFC >  lfc_cutoff] <- "Up"
deg_plot$volcano_class[deg_plot$adj.P.Val < fdr_cutoff & deg_plot$logFC < -lfc_cutoff] <- "Down"

p_volcano <- ggplot(
    deg_plot,
    aes(x = logFC, y = -log10(adj.P.Val), color = volcano_class)
  ) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values = c(Up = "#E3B505", Down = "purple", NS = "grey70")) +
  geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff),
             linetype = "dashed", linewidth = 0.4) +
  geom_hline(yintercept = -log10(fdr_cutoff),
             linetype = "dashed", linewidth = 0.4) +
  geom_text_repel(
    data        = subset(deg_plot, volcano_class %in% c("Up", "Down")),
    aes(label   = gene_symbol),
    size        = 3.8,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "grey50"
  ) +
  labs(
    title = "Volcano Plot — Ischemic Stroke vs. Controls (GSE16561)",
    x     = "log2 Fold Change",
    y     = "-log10(FDR)",
    color = "Regulation"
  ) +
  theme_classic(base_size = 14)

dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
print(p_volcano)
ggsave("outputs/figures/Volcano.png", plot = p_volcano,
       width = 7, height = 6, dpi = 300)


# =============================================================================
# SECTION 4 — Protein-Protein Interaction Network (STRINGdb)
# =============================================================================
# We map all significant DEGs (both Up and Down) to STRING protein identifiers
# and retrieve high-confidence interactions (combined_score >= 400).
# The largest connected component of the resulting graph is plotted using a
# Fruchterman-Reingold layout. Edge thickness reflects interaction confidence.
# =============================================================================

library(dplyr)
library(STRINGdb)
library(igraph)
library(ggraph)

lfc_cutoff <- 1
fdr_cutoff <- 0.05

deg_plot <- deg_anotado
deg_plot$gene_symbol  <- deg_plot$Gene
deg_plot$volcano_class <- "NS"
deg_plot$volcano_class[deg_plot$adj.P.Val < fdr_cutoff & deg_plot$logFC >  lfc_cutoff] <- "Up"
deg_plot$volcano_class[deg_plot$adj.P.Val < fdr_cutoff & deg_plot$logFC < -lfc_cutoff] <- "Down"

up_genes   <- deg_plot %>% filter(volcano_class == "Up")   %>% pull(gene_symbol) %>% unique()
down_genes <- deg_plot %>% filter(volcano_class == "Down") %>% pull(gene_symbol) %>% unique()
deg_genes  <- sort(unique(c(up_genes, down_genes)))

cat("Up genes:", length(up_genes), "\n")
cat("Down genes:", length(down_genes), "\n")
cat("Total DEGs for PPI:", length(deg_genes), "\n")

if (length(deg_genes) < 2) stop("Too few DEGs to build a PPI network.")

# score_threshold = 400 is STRING's recommended cutoff for medium-confidence
# interactions. Increasing it (e.g., 700) gives a sparser but higher-quality graph.
options(timeout = 600)
string_db <- STRINGdb$new(
  version         = "11.5",
  species         = 9606,
  score_threshold = 400,
  input_directory = ""
)

# Normalise gene symbols before mapping: convert to upper case, strip
# surrounding whitespace, and drop blank or missing entries.
deg_genes_norm <- toupper(trimws(deg_genes))
deg_genes_norm <- deg_genes_norm[!is.na(deg_genes_norm) & deg_genes_norm != ""]
deg_genes_norm <- unique(deg_genes_norm)

node_annot <- data.frame(
  gene_symbol = deg_genes_norm,
  regulation  = ifelse(deg_genes_norm %in% toupper(trimws(up_genes)), "Up", "Down"),
  stringsAsFactors = FALSE
)

mapped <- string_db$map(node_annot, "gene_symbol", removeUnmappedRows = TRUE)

# ── Mapping diagnostics (PPI) ─────────────────────────────────────────────────
n_submitted_ppi  <- nrow(node_annot)
n_mapped_ppi     <- nrow(mapped)
mapping_rate_ppi <- if (n_submitted_ppi > 0) 100 * n_mapped_ppi / n_submitted_ppi else 0
cat("--- STRING mapping diagnostics (PPI) ---\n")
cat("Submitted:", n_submitted_ppi, "\n")
cat("Mapped:   ", n_mapped_ppi, "\n")
cat(sprintf("Rate:      %.1f%%\n", mapping_rate_ppi))
cat("Unmapped (first 20):",
    head(setdiff(toupper(node_annot$gene_symbol), toupper(mapped$gene_symbol)), 20), "\n")

if (mapping_rate_ppi < 5) {
  warning(sprintf(
    "PPI mapping rate is very low (%.1f%%). Possible causes:\n", mapping_rate_ppi),
    "  * Gene symbols may be outdated aliases not recognised by STRING v11.5\n",
    "  * There may be a mismatch between the STRING database version and the\n",
    "    gene nomenclature used in the microarray annotation file\n",
    "Suggested action: paste the unmapped symbols listed above into\n",
    "  https://string-db.org and check which identifiers STRING accepts."
  )
}

cat("Genes mapped to STRING:", n_mapped_ppi, "of", n_submitted_ppi, "\n")

if (nrow(mapped) < 2) stop("Fewer than 2 genes mapped to STRING; cannot build PPI.")

# Retrieve interactions with retry logic in case of transient network errors
MAX_ATTEMPTS <- 3
ppi_raw <- NULL
for (attempt in seq_len(MAX_ATTEMPTS)) {
  ppi_raw <- tryCatch(
    string_db$get_interactions(mapped$STRING_id),
    error = function(e) {
      message("Attempt ", attempt, " failed: ", conditionMessage(e))
      Sys.sleep(10 * attempt)
      NULL
    }
  )
  if (!is.null(ppi_raw)) break
}
if (is.null(ppi_raw)) {
  stop("PPI download failed after ", MAX_ATTEMPTS, " attempts. ",
       "Check your internet connection or try again later.")
}

ppi <- ppi_raw %>%
  filter(from %in% mapped$STRING_id, to %in% mapped$STRING_id)

cat("Edges retained (score >= 400):", nrow(ppi), "\n")
if (nrow(ppi) == 0) stop("No edges after filter. Consider lowering score_threshold.")

# Build undirected weighted graph
g <- graph_from_data_frame(
  d        = ppi %>% transmute(from, to, weight = combined_score),
  directed = FALSE,
  vertices = mapped %>% transmute(
    name        = STRING_id,
    gene_symbol = gene_symbol,
    regulation  = regulation
  )
)

# Keep only the largest connected component for clarity
components_g <- components(g)
giant_idx    <- which.max(components_g$csize)
g_cc         <- induced_subgraph(g, vids = V(g)[components_g$membership == giant_idx])

if (vcount(g_cc) < 2) {
  stop("The largest connected component has fewer than 2 nodes. ",
       "Try lowering score_threshold or broadening the DEG criteria.")
}

cat("Largest connected component — nodes:", vcount(g_cc),
    "  edges:", ecount(g_cc), "\n")

set.seed(1)
p_ppi <- ggraph(g_cc, layout = "fr") +
  geom_edge_link(aes(width = weight), alpha = 0.25) +
  scale_edge_width(range = c(0.2, 2.2), guide = "none") +
  geom_node_point(aes(color = regulation), size = 4) +
  scale_color_manual(values = c(Up = "#E3B505", Down = "purple")) +
  geom_node_text(aes(label = gene_symbol), repel = TRUE, size = 3.5) +
  theme_void(base_size = 14) +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  ggtitle("STRING PPI Network — Ischemic Stroke DEGs (score ≥ 400)")

print(p_ppi)
ggsave("outputs/figures/PPI.png", plot = p_ppi,
       width = 9, height = 7, dpi = 300, bg = "white")


# =============================================================================
# SECTION 5 — STRING Functional Enrichment
# =============================================================================
# We run enrichment separately for Up-regulated and Down-regulated gene sets.
# This matters biologically: the Up signature reflects innate immune activation
# (neutrophil degranulation, cytokine signaling), while the Down signature
# reflects post-stroke immunodepression (lymphocyte suppression, adaptive
# immunity downregulation) — two distinct and clinically relevant processes.
#
# The full DEG set (not a top-N subset) is sent to STRING to avoid
# artificially truncating the input signal.
# =============================================================================

suppressPackageStartupMessages({
  library(STRINGdb)
  library(dplyr)
  library(rio)
})

lfc_cutoff <- 1
fdr_cutoff <- 0.05

deg_filtered <- deg_anotado %>%
  filter(!is.na(Gene)) %>%
  filter(abs(logFC) > lfc_cutoff, adj.P.Val < fdr_cutoff) %>%
  mutate(regulation = ifelse(logFC > 0, "Up", "Down"))

up_genes   <- deg_filtered %>% filter(regulation == "Up")   %>% pull(Gene)
down_genes <- deg_filtered %>% filter(regulation == "Down") %>% pull(Gene)

cat("Genes sent to STRING enrichment:\n")
cat("  Up:", length(up_genes), "\n")
cat("  Down:", length(down_genes), "\n")

# Reuse the STRING connection from the PPI section if it exists; otherwise
# create a new one. For enrichment, the score threshold is handled server-side,
# so we initialize with score_threshold = 0.
if (!exists("string_db")) {
  string_db <- STRINGdb$new(
    version         = "11.5",
    species         = 9606,
    score_threshold = 0,
    input_directory = ""
  )
  message("STRINGdb connection initialized for enrichment.")
}

# Maps a vector of gene symbols to STRING protein identifiers.
# Returns only rows that could be mapped (unmapped genes are dropped).
symbols_to_string_map <- function(symbols, regulation_label, string_db) {
  # Normalise gene symbols: convert to upper case, strip surrounding
  # whitespace, and drop blank or missing entries before mapping.
  symbols_norm <- toupper(trimws(as.character(symbols)))
  symbols_norm <- symbols_norm[!is.na(symbols_norm) & symbols_norm != ""]

  df <- data.frame(
    gene_symbol = unique(symbols_norm),
    regulation  = regulation_label,
    stringsAsFactors = FALSE
  )
  if (nrow(df) == 0) return(data.frame())

  mapped <- string_db$map(df, "gene_symbol", removeUnmappedRows = TRUE)

  if (!"STRING_id" %in% colnames(mapped)) {
    stop("STRINGdb$map() did not return a STRING_id column. Columns found: ",
         paste(colnames(mapped), collapse = ", "))
  }

  # ── Mapping diagnostics ────────────────────────────────────────────────────
  n_submitted  <- nrow(df)
  n_mapped     <- nrow(mapped)
  mapping_rate <- if (n_submitted > 0) 100 * n_mapped / n_submitted else 0
  cat(sprintf("--- STRING mapping diagnostics (%s) ---\n", regulation_label))
  cat("Submitted:", n_submitted, "\n")
  cat("Mapped:   ", n_mapped, "\n")
  cat(sprintf("Rate:      %.1f%%\n", mapping_rate))
  cat("Unmapped (first 20):",
      head(setdiff(toupper(df$gene_symbol), toupper(mapped$gene_symbol)), 20),
      "\n")

  if (mapping_rate < 5) {
    warning(sprintf(
      "%s mapping rate is very low (%.1f%%). Possible causes:\n",
      regulation_label, mapping_rate),
      "  * Gene symbols may be outdated aliases not recognised by STRING v11.5\n",
      "  * There may be a mismatch between the STRING database version and the\n",
      "    gene nomenclature used in the microarray annotation file\n",
      "Suggested action: paste the unmapped symbols listed above into\n",
      "  https://string-db.org and check which identifiers STRING accepts."
    )
  }

  mapped
}

# Runs STRING enrichment for a gene set and annotates the result with
# metadata about how many genes were submitted and successfully mapped.
run_string_enrichment <- function(symbols, label, string_db) {
  mapped <- symbols_to_string_map(symbols, label, string_db)

  cat(label, "— submitted:", length(unique(symbols)),
      "  mapped:", nrow(mapped), "\n")

  if (nrow(mapped) < 2) {
    warning(label, ": fewer than 2 proteins mapped; enrichment is not meaningful.")
    return(list(mapped = mapped, enrich = data.frame()))
  }

  enr <- as.data.frame(string_db$get_enrichment(mapped$STRING_id))
  enr$set_label         <- label
  enr$n_input_symbols   <- length(unique(symbols))
  enr$n_mapped_proteins <- nrow(mapped)

  list(mapped = mapped, enrich = enr)
}

res_up   <- run_string_enrichment(up_genes,   "Up",   string_db)
res_down <- run_string_enrichment(down_genes, "Down", string_db)

enr_up   <- res_up$enrich
enr_down <- res_down$enrich

cat("\n--- Top 20 enriched terms (Up) ---\n")
print(head(enr_up,   20))
cat("\n--- Top 20 enriched terms (Down) ---\n")
print(head(enr_down, 20))

# Export mapping tables and full enrichment results
dir.create("string_enrichment", showWarnings = FALSE)
rio::export(res_up$mapped,   "string_enrichment/STRING_mapping_Up.tsv")
rio::export(res_down$mapped, "string_enrichment/STRING_mapping_Down.tsv")
rio::export(enr_up,          "string_enrichment/STRING_enrichment_Up_all.tsv")
rio::export(enr_down,        "string_enrichment/STRING_enrichment_Down_all.tsv")

# Persist key R objects for downstream use or resuming the session
dir.create("rds", showWarnings = FALSE)
saveRDS(deg_anotado, "rds/deg_anotado.rds")
saveRDS(enr_up,      "rds/enr_up.rds")
saveRDS(enr_down,    "rds/enr_down.rds")
rio::export(deg_anotado, "string_enrichment/DEGs_anotados.tsv")
message("Full annotated DEG table saved to string_enrichment/DEGs_anotados.tsv")

# The PPI graph object is saved conditionally — it may not exist if the
# PPI section encountered a network error
if (exists("g_cc")) {
  saveRDS(g_cc, "rds/ppi_graph_cc.rds")
  message("PPI graph saved to rds/ppi_graph_cc.rds")
} else {
  warning("PPI graph (g_cc) not found — the PPI section may not have completed. ",
          "rds/ppi_graph_cc.rds was not generated.")
}

gc()


# =============================================================================
# SECTION 6 — Enrichment Dotplots
# =============================================================================
# Dotplots visualize the top enriched biological terms ranked by significance.
# Only functional annotation categories are shown — PMID entries (literature
# co-citations returned by STRING) are explicitly excluded, as their
# 'description' field contains paper titles, not biological process names.
# =============================================================================

library(ggplot2)
library(dplyr)
library(stringr)
library(grid)

# Categories that contain genuine biological annotations.
# "PMID", "NetworkNeighborAL", and "TISSUES" are excluded.
BIOLOGICAL_CATEGORIES <- c(
  "Process", "Component", "Function",
  "Keyword", "COMPARTMENTS", "RCTM",
  "InterPro", "SMART", "HPO", "Pfam"
)

# Builds a dotplot of the top enriched terms.
# Point size = number of genes in the term; color = -log10(FDR).
plot_string_enrichment_dotplot <- function(enr_df, top_n = 50,
                                           title = "Enrichment Dotplot") {
  if (nrow(enr_df) == 0) {
    warning("Enrichment data frame is empty.")
    return(NULL)
  }

  enr_df <- enr_df %>% filter(category %in% BIOLOGICAL_CATEGORIES)

  if (nrow(enr_df) == 0) {
    warning("No terms from biological annotation categories found after filtering.")
    return(NULL)
  }

  top_terms <- enr_df %>%
    arrange(fdr) %>%
    slice_head(n = min(top_n, nrow(enr_df))) %>%
    mutate(
      description_wrapped = str_wrap(description, width = 65),
      description_factor  = factor(
        description_wrapped,
        levels = rev(unique(description_wrapped))
      )
    )

  ggplot(top_terms, aes(
    x     = -log10(fdr),
    y     = description_factor,
    size  = number_of_genes,
    color = -log10(fdr)
  )) +
    geom_point(alpha = 0.9) +
    scale_color_gradient(low = "#E3B505", high = "#d73027") +
    guides(
      size  = guide_legend(order = 1, override.aes = list(alpha = 1)),
      color = guide_colorbar(order = 2)
    ) +
    labs(
      x     = "-log10(FDR)",
      y     = "Biological Process / Term",
      title = title,
      size  = "Genes in term",
      color = "-log10(FDR)"
    ) +
    theme_bw(base_size = 13) +
    theme(
      axis.text.y      = element_text(size = 8, lineheight = 0.9),
      legend.position  = "right",
      legend.box       = "vertical",
      legend.spacing.y = unit(0.35, "cm"),
      plot.margin      = margin(t = 10, r = 18, b = 10, l = 14)
    )
}

# Computes a sensible figure height (inches) based on how many terms will
# be plotted and their wrapped label lengths.
compute_dotplot_height <- function(enr_df, top_n = 50) {
  enr_df  <- enr_df %>% filter(category %in% BIOLOGICAL_CATEGORIES)
  n_terms <- min(top_n, nrow(enr_df))
  if (n_terms == 0) return(10)

  n_lines <- enr_df %>%
    arrange(fdr) %>%
    slice_head(n = n_terms) %>%
    mutate(description_wrapped = str_wrap(description, width = 65)) %>%
    pull(description_wrapped) %>%
    str_count("\n") + 1

  max(10, min(28, 4 + (0.18 * sum(n_lines))))
}

dotplot_up <- plot_string_enrichment_dotplot(
  enr_up,
  top_n = 50,
  title = "STRING Enrichment — Up-regulated Genes (Ischemic Stroke)"
)
print(dotplot_up)
ggsave(
  "outputs/figures/dotplot_up_STRING.png",
  plot   = dotplot_up,
  width  = 15,
  height = compute_dotplot_height(enr_up, top_n = 50),
  dpi    = 300,
  bg     = "white"
)

dotplot_down <- plot_string_enrichment_dotplot(
  enr_down,
  top_n = 50,
  title = "STRING Enrichment — Down-regulated Genes (Ischemic Stroke)"
)
if (!is.null(dotplot_down)) {
  print(dotplot_down)
  ggsave(
    "outputs/figures/dotplot_down_STRING.png",
    plot   = dotplot_down,
    width  = 15,
    height = compute_dotplot_height(enr_down, top_n = 50),
    dpi    = 300,
    bg     = "white"
  )
}


# =============================================================================
# SECTION 7 — Import NANDA-I Mapping Results (post-notebook)
# =============================================================================
# This section is run AFTER clinicalbert_vias_AVC.ipynb has been executed.
# It reads the NANDA-I mapping CSVs generated by the notebook for any
# further R-side analysis or custom visualization.
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)

nanda_mapping <- read.csv(
  "outputs/nanda/nanda_mapping_completo_threshold65.csv",
  fileEncoding = "UTF-8"
)

cat("NANDA-I valid pairs (similarity >= 0.65):", nrow(nanda_mapping), "\n")
print(nanda_mapping)

# Full cosine similarity matrix (Up genes vs NANDA-I)
matriz_up <- read.csv(
  "outputs/nanda/matriz_similaridade_up_nanda.csv",
  row.names    = 1,
  fileEncoding = "UTF-8"
)

# Full cosine similarity matrix (Down genes vs NANDA-I)
matriz_down <- read.csv(
  "outputs/nanda/matriz_similaridade_down_nanda.csv",
  row.names    = 1,
  fileEncoding = "UTF-8"
)

# ── Summary statistics ────────────────────────────────────────────────────────
n_up   <- sum(nanda_mapping$regulacao == "Up",   na.rm = TRUE)
n_down <- sum(nanda_mapping$regulacao == "Down", na.rm = TRUE)

cat("\n--- NANDA-I Mapping Summary ---\n")
cat(sprintf(
  "Valid pairs (similarity >= 0.65)  Up: %d  |  Down: %d  |  Total: %d\n",
  n_up, n_down, nrow(nanda_mapping)
))

cat("\nTop 5 pairs by cosine similarity:\n")
top5 <- nanda_mapping %>%
  arrange(desc(similaridade)) %>%
  slice_head(n = 5)
print(top5)

cat("\nPairs per NANDA-I diagnosis (both directions):\n")
dist_diag <- nanda_mapping %>%
  count(diagnostico_nanda, name = "n_pares") %>%
  arrange(desc(n_pares))
print(dist_diag)

# ── Normalização de pares por número de vias ─────────────────────────────
# Forma correta: conta vias únicas por direção (evita confundir 9 Up vs 26 Down).
n_vias_up   <- nanda_mapping %>%
  filter(regulacao == "Up") %>%
  pull(via) %>%
  n_distinct()
n_vias_down <- nanda_mapping %>%
  filter(regulacao == "Down") %>%
  pull(via) %>%
  n_distinct()

cat(sprintf(
  "\nVias únicas — Up: %d  |  Down: %d\n",
  n_vias_up, n_vias_down
))

# Pares normalizados por via (= densidade de cobertura)
dist_norm <- nanda_mapping %>%
  group_by(diagnostico_nanda, regulacao) %>%
  summarise(n_pares_raw = n(), .groups = "drop") %>%
  mutate(
    n_vias_direcao = if_else(regulacao == "Up", n_vias_up, n_vias_down),
    pares_por_via  = n_pares_raw / n_vias_direcao,
    label_curto    = stringr::str_extract(diagnostico_nanda, "^[^:]+")
  ) %>%
  arrange(desc(pares_por_via))

cat("\nDensidade de cobertura (pares / n_vias) por diagnóstico:\n")
print(dist_norm %>% select(label_curto, regulacao,
                            n_pares_raw, n_vias_direcao,
                            pares_por_via), n = Inf)

dir.create("outputs/nanda", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
write.csv(
  dist_norm,
  "outputs/nanda/cobertura_normalizada.csv",
  row.names = FALSE, fileEncoding = "UTF-8"
)
cat("✓ Exportado → outputs/nanda/cobertura_normalizada.csv\n")

# Plot comparativo normalizado
p_norm <- ggplot(
    dist_norm,
    aes(
      x = pares_por_via,
      y = reorder(label_curto, pares_por_via),
      fill = regulacao
    )
  ) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = c(Up = "#E3B505", Down = "purple")) +
  geom_vline(xintercept = 1.0, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  annotate("text", x = 1.02, y = 0.5,
           label = "1 par / via", size = 3,
           color = "grey40", hjust = 0) +
  labs(
    title   = "Cobertura Normalizada — Pares NANDA-I por Via",
    subtitle = "Valores > 1.0 = diagnóstico aparece em média em mais de uma via por direção",
    x       = "Pares válidos / número de vias na direção",
    y       = NULL,
    fill    = "Regulação"
  ) +
  theme_bw(base_size = 12) +
  theme(axis.text.y = element_text(size = 9),
        legend.position = "right")

ggsave("outputs/figures/cobertura_normalizada.png",
  plot = p_norm, width = 12, height = 6,
  dpi = 300, bg = "white")
cat("✓ Plot exportado → outputs/figures/cobertura_normalizada.png\n")

# ── Análise de especificidade direcional ─────────────────────────────────────
# Para cada diagnóstico NANDA-I, calcula:
#   - score médio nos pares Up
#   - score médio nos pares Down
#   - índice de especificidade: (score_up - score_down) / (score_up + score_down)
#   Valores próximos de +1 = diagnóstico primariamente Up-driven
#   Valores próximos de -1 = diagnóstico primariamente Down-driven
#   Valores próximos de 0  = diagnóstico não-específico (sinal de colapso)

especificidade <- nanda_mapping %>%
  group_by(diagnostico_nanda, regulacao) %>%
  summarise(score_medio = mean(similaridade), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from  = regulacao,
    values_from = score_medio,
    names_prefix = "score_"
  ) %>%
  {
    d <- .
    if (!"score_Up" %in% names(d)) d$score_Up <- 0
    if (!"score_Down" %in% names(d)) d$score_Down <- 0
    d
  } %>%
  mutate(
    score_Up   = replace_na(score_Up, 0),
    score_Down = replace_na(score_Down, 0),
    indice_especificidade = (score_Up - score_Down) /
      (score_Up + score_Down + 1e-9)
  ) %>%
  arrange(desc(abs(indice_especificidade)))

cat("\n--- Especificidade Direcional dos Diagnósticos NANDA-I ---\n")
cat("(+1 = exclusivo Up/Inato | -1 = exclusivo Down/Adaptativo | 0 = não-específico)\n\n")
print(especificidade, n = Inf)

dir.create("outputs/nanda", recursive = TRUE, showWarnings = FALSE)
write.csv(
  especificidade,
  "outputs/nanda/especificidade_direcional.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# Plot de especificidade
p_espec <- ggplot(
    especificidade,
    aes(
      x = indice_especificidade,
      y = reorder(
        stringr::str_extract(diagnostico_nanda, "^[^:]+"),
        indice_especificidade
      ),
      fill = indice_especificidade > 0
    )
  ) +
  geom_col(show.legend = FALSE) +
  scale_fill_manual(values = c("TRUE" = "#E3B505", "FALSE" = "purple")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    title = "Especificidade Direcional dos Diagnósticos NANDA-I",
    subtitle = "Amarelo = primariamente Up/Inato  |  Roxo = primariamente Down/Adaptativo",
    x = "Índice de Especificidade  [(Up - Down) / (Up + Down)]",
    y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(axis.text.y = element_text(size = 9))

ggsave(
  "outputs/figures/especificidade_direcional.png",
  plot   = p_espec,
  width  = 12,
  height = 7,
  dpi    = 300,
  bg     = "white"
)

message("✓ Análise de especificidade direcional salva.")

# ── Barplot: valid pairs per NANDA-I diagnosis, coloured by regulation ────────
# Counts how many pathway–diagnosis pairs each NANDA-I code accumulated, split
# by regulatory direction, to highlight which diagnoses are most supported by
# the transcriptomic evidence.
plot_data <- nanda_mapping %>%
  count(diagnostico_nanda, regulacao, name = "n_pares")

# Order diagnoses on the y-axis by their total count across both directions
order_diag <- nanda_mapping %>%
  count(diagnostico_nanda, name = "total") %>%
  arrange(total) %>%
  pull(diagnostico_nanda)

plot_data <- plot_data %>%
  mutate(
    diagnostico_nanda = factor(diagnostico_nanda, levels = order_diag),
    regulacao         = factor(regulacao, levels = c("Up", "Down"))
  )

p_bar <- ggplot(
    plot_data,
    aes(x = n_pares, y = diagnostico_nanda, fill = regulacao)
  ) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = c(Up = "#E3B505", Down = "purple")) +
  labs(
    title = "Valid Pathway\u2013NANDA-I Pairs per Diagnosis (similarity \u2265 0.65)",
    x     = "Number of valid pairs",
    y     = "NANDA-I Diagnosis",
    fill  = "Regulation"
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.y     = element_text(size = 9),
    legend.position = "right"
  )

print(p_bar)
ggsave(
  "outputs/figures/barplot_nanda_summary.png",
  plot   = p_bar,
  width  = 10,
  height = 6,
  dpi    = 300,
  bg     = "white"
)

message("\n=== Pipeline complete ===")
message("All outputs saved to outputs/figures/ and outputs/nanda/")
message("Next step: review outputs/nanda/nanda_mapping_completo_threshold65.csv")
