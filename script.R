# Diretório do projeto: variável de ambiente (opcional), depois pasta do script.
set_project_wd <- function() {
  env <- Sys.getenv("AVC_TRANSCRIPTOMICS_DIR", unset = "")
  if (nzchar(env) && dir.exists(env)) {
    setwd(env)
    return(invisible(TRUE))
  }
  ca <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", ca, value = TRUE)
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
  warning(
    "Defina o diretório de trabalho para a pasta que contém script.R ",
    "(RStudio: Session > Set Working Directory > To Source File Location), ",
    "ou defina a variável de ambiente AVC_TRANSCRIPTOMICS_DIR."
  )
  invisible(FALSE)
}
set_project_wd()

library(GEOquery)
library(limma)

# baixar dataset
gset <- getGEO("GSE16561", GSEMatrix = TRUE)
gset <- gset[[1]]

# matriz de express�o
expr <- exprs(gset)

# metadados cl�nicos
meta <- pData(gset)

# ver colunas dispon�veis
colnames(meta)

# visualizar conte�do de algumas colunas
head(meta)

# no GSE16561 a condi��o geralmente est� em "characteristics_ch1"
table(meta$characteristics_ch1)

# criar vetor de grupos a partir da coluna description
grupo <- ifelse(grepl("Stroke", meta$description, ignore.case = TRUE),
                "AVC",
                ifelse(grepl("Control", meta$description, ignore.case = TRUE),
                       "Controle",
                       NA))

# transformar em fator
grupo <- factor(grupo)

# visualizar distribui��o
table(grupo)

#remover amostras n�o classificadas
validos <- !is.na(grupo)

expr <- expr[, validos]
grupo <- droplevels(grupo[validos])

#limma
design <- model.matrix(~0 + grupo)
colnames(design) <- c("AVC","Controle")

fit <- lmFit(expr, design)

contr <- makeContrasts(AVCvsControle = AVC - Controle, levels = design)

fit2 <- contrasts.fit(fit, contr)
fit2 <- eBayes(fit2)

results <- topTable(fit2, number = 20)

results

deg <- topTable(fit2, adjust="fdr", p.value=0.05, number=Inf)

anotar_genes <- function(deg, gset, col_symbol = NULL, col_title = NULL, colapsar = TRUE) {
  
  anot <- fData(gset)
  
  # garantir coluna de probe
  deg$ProbeID <- rownames(deg)
  
  # identificar automaticamente coluna de s�mbolo do gene se n�o informada
  if(is.null(col_symbol)) {
    gene_col <- grep("symbol", colnames(anot), ignore.case = TRUE, value = TRUE)[1]
  } else {
    gene_col <- col_symbol
  }
  
  # identificar coluna com t�tulo do gene (opcional)
  if(!is.null(col_title)) {
    if(col_title %in% colnames(anot)) {
      title_col <- col_title
    } else {
      title_col <- NULL
    }
  } else {
    title_col <- NULL
  }
  
  # construir tabela de anota��o
  if(!is.null(title_col)) {
    anot2 <- data.frame(
      ProbeID = rownames(anot),
      Gene = anot[[gene_col]],
      GeneTitle = anot[[title_col]],
      stringsAsFactors = FALSE
    )
  } else {
    anot2 <- data.frame(
      ProbeID = rownames(anot),
      Gene = anot[[gene_col]],
      stringsAsFactors = FALSE
    )
  }
  
  # merge com resultados
  deg_anotado <- merge(deg, anot2, by="ProbeID", all.x = TRUE)
  
  # colapsar probes duplicadas mantendo maior |logFC|
  if(colapsar) {
    deg_anotado <- deg_anotado[order(-abs(deg_anotado$logFC)), ]  # ordenar por logFC absoluto
    deg_anotado <- deg_anotado[!duplicated(deg_anotado$Gene), ]   # manter apenas a primeira ocorr�ncia de cada gene
  }
  
  return(deg_anotado)
}

deg_anotado <- anotar_genes(deg, gset)
head(deg_anotado)

View(deg_anotado)

library(ggplot2)
library(ggrepel)

# usar deg_anotado e criar coluna gene_symbol
deg_plot <- deg_anotado
deg_plot$gene_symbol <- deg_plot$Gene  # renomear para compatibilidade

# Limiares
lfc_cutoff <- 1
fdr_cutoff <- 0.05

# Classificar genes
deg_plot$volcano_class <- "NS"
deg_plot$volcano_class[deg_plot$adj.P.Val < fdr_cutoff & deg_plot$logFC >  lfc_cutoff] <- "Up"
deg_plot$volcano_class[deg_plot$adj.P.Val < fdr_cutoff & deg_plot$logFC < -lfc_cutoff] <- "Down"

# Gr�fico base
p <- ggplot(
  deg_plot,
  aes(x = logFC, y = -log10(adj.P.Val), color = volcano_class)
) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(
    values = c("Up"="blue", "Down"="purple", "NS"="grey70")
  ) +
  geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype="dashed", linewidth=0.4) +
  geom_hline(yintercept = -log10(fdr_cutoff), linetype="dashed", linewidth=0.4) +
  labs(title="VolcanoPlot - Stroke", x="log2 Fold Change", y="-log10(FDR)", color="Regulation") +
  theme_classic(base_size=14)

# Rotular apenas DEGs
p_labeled <- p +
  geom_text_repel(
    data = subset(deg_plot, volcano_class %in% c("Up","Down")),
    aes(label=gene_symbol),
    size=3.8,
    max.overlaps=Inf,
    box.padding=0.4,
    point.padding=0.3,
    segment.color="grey50"
  )

print(p_labeled)
ggsave("Volcano.png", plot=p_labeled, width=7, height=6, dpi=300)


library(dplyr)
library(STRINGdb)
library(igraph)
library(ggraph)
library(ggplot2)

# --- Preparar conjuntos de genes Up e Down ---
lfc_cutoff <- 1
fdr_cutoff <- 0.05

deg_plot <- deg_anotado
deg_plot$gene_symbol <- deg_plot$Gene  # renomeia para compatibilidade

deg_plot$volcano_class <- "NS"
deg_plot$volcano_class[deg_plot$adj.P.Val < fdr_cutoff & deg_plot$logFC >  lfc_cutoff] <- "Up"
deg_plot$volcano_class[deg_plot$adj.P.Val < fdr_cutoff & deg_plot$logFC < -lfc_cutoff] <- "Down"

up_genes   <- deg_plot %>% filter(volcano_class=="Up")   %>% pull(gene_symbol) %>% unique()
down_genes <- deg_plot %>% filter(volcano_class=="Down") %>% pull(gene_symbol) %>% unique()
deg_genes  <- sort(unique(c(up_genes, down_genes)))

cat("Up genes:", length(up_genes), "\n")
cat("Down genes:", length(down_genes), "\n")
cat("Total DEGs (Up+Down):", length(deg_genes), "\n")

if(length(deg_genes) < 2) stop("Poucos DEGs para gerar PPI")

# --- Inicializar STRINGdb ---
options(timeout=300000)
string_db <- STRINGdb$new(version="11.5", species=9606, score_threshold=0, input_directory="")

# --- Criar tabela de n�s com regula��o ---
node_annot <- data.frame(
  gene_symbol = deg_genes,
  regulation  = ifelse(deg_genes %in% up_genes, "Up", "Down"),
  stringsAsFactors = FALSE
)

# --- Mapear genes para STRING IDs ---
mapped <- string_db$map(node_annot, "gene_symbol", removeUnmappedRows = TRUE)
cat("Genes mapeados para STRING:", nrow(mapped), "de", nrow(node_annot), "\n")

if(nrow(mapped) < 2) stop("Menos de 2 genes mapeados no STRING; PPI imposs�vel")

# --- Recuperar intera��es PPI ---
ppi_raw <- string_db$get_interactions(mapped$STRING_id)

# Filtrar arestas apenas entre nossos genes e score alto
score_cutoff <- 400
ppi <- ppi_raw %>%
  filter(from %in% mapped$STRING_id, to %in% mapped$STRING_id) %>%
  filter(combined_score >= score_cutoff)

cat("Arestas ap�s filtro (score >=", score_cutoff, "):", nrow(ppi), "\n")
if(nrow(ppi)==0) stop("Nenhuma aresta ap�s filtro; diminua score_cutoff")

# --- Construir grafo igraph ---
g <- graph_from_data_frame(
  d = ppi %>% transmute(from, to, weight=combined_score),
  directed = FALSE,
  vertices = mapped %>% transmute(name=STRING_id, gene_symbol=gene_symbol, regulation=regulation)
)

# Manter apenas maior componente conexo
components <- components(g)
giant <- which.max(components$csize)
g_cc <- induced_subgraph(g, vids=V(g)[components$membership==giant])

cat("Nodes (maior CC):", vcount(g_cc), "Edges:", ecount(g_cc), "\n")

# --- Plotar rede PPI com ggraph ---
set.seed(1)
p_ppi <- ggraph(g_cc, layout="fr") +
  geom_edge_link(aes(width=weight), alpha=0.25) +
  scale_edge_width(range=c(0.2,2.2), guide="none") +
  geom_node_point(aes(color=regulation), size=4) +
  scale_color_manual(values=c(Up="blue", Down="purple")) +
  geom_node_text(aes(label=gene_symbol), repel=TRUE, size=3.5) +
  theme_void(base_size=14) +
  theme(
    plot.background = element_rect(fill="white", color=NA),
    panel.background = element_rect(fill="white", color=NA)
  ) +
  ggtitle(paste0("Rede PPI STRING (score >=", score_cutoff, ")"))

print(p_ppi)

# --- Salvar PNG com fundo branco ---
ggsave("PPI.png", plot=p_ppi, width=9, height=7, dpi=300, bg="white")


##############################
# Bibliotecas
##############################
suppressPackageStartupMessages({
  library(STRINGdb)
  library(dplyr)
  library(rio)
})

##############################
# 1) Selecionar top 50 genes Up e Down por logFC
##############################
lfc_cutoff <- 1
top_n <- 50

deg_filtered <- deg_anotado %>%
  filter(!is.na(Gene)) %>%
  filter(abs(logFC) > lfc_cutoff) %>%
  mutate(regulation = ifelse(logFC > 0, "Up", "Down"))

up_genes   <- deg_filtered %>% filter(regulation == "Up") %>% arrange(-logFC) %>% slice(1:top_n) %>% pull(Gene)
down_genes <- deg_filtered %>% filter(regulation == "Down") %>% arrange(logFC) %>% slice(1:top_n) %>% pull(Gene)

cat("Up genes:", length(up_genes), "\n")
cat("Down genes:", length(down_genes), "\n")

##############################
# 2) Inicializar STRINGdb
##############################
string_db <- STRINGdb$new(
  version = "11.5",
  species = 9606,
  score_threshold = 0,
  input_directory = ""
)

##############################
# 3) Função de mapeamento robusto para STRING
##############################
symbols_to_string_map <- function(symbols, regulation_label, string_db) {
  df <- data.frame(
    gene_symbol = unique(as.character(symbols)),
    regulation  = regulation_label,
    stringsAsFactors = FALSE
  )
  df <- df[!is.na(df$gene_symbol) & df$gene_symbol != "", , drop = FALSE]
  
  if (nrow(df) == 0) return(data.frame())
  
  mapped <- string_db$map(df, "gene_symbol", removeUnmappedRows = TRUE)
  
  if (!("STRING_id" %in% colnames(mapped))) {
    stop("STRINGdb$map() did not return a STRING_id column. Columns: ",
         paste(colnames(mapped), collapse = ", "))
  }
  
  mapped
}

##############################
# 4) Função de enriquecimento STRING
##############################
run_string_enrichment <- function(symbols, label, string_db) {
  
  mapped <- symbols_to_string_map(symbols, label, string_db)
  
  cat(label, "- input symbols:", length(unique(symbols)), "\n")
  cat(label, "- mapped proteins:", nrow(mapped), "\n")
  
  if (nrow(mapped) < 2) {
    warning(label, ": <2 mapped proteins; enrichment not meaningful.")
    return(list(mapped = mapped, enrich = data.frame()))
  }
  
  enr <- string_db$get_enrichment(mapped$STRING_id)
  enr <- as.data.frame(enr)
  
  enr$set_label        <- label
  enr$n_input_symbols  <- length(unique(symbols))
  enr$n_mapped_proteins <- nrow(mapped)
  
  list(mapped = mapped, enrich = enr)
}

##############################
# 5) Executar para Up e Down
##############################
res_up   <- run_string_enrichment(up_genes,   "Up",   string_db)
res_down <- run_string_enrichment(down_genes, "Down", string_db)

enr_up   <- res_up$enrich
enr_down <- res_down$enrich

##############################
# 6) Inspeçãoo rápida
##############################
View(enr_up, 20)
View(enr_down, 1) # Down tem apenas 1 termo enriquecido, provavelmente por poucos genes mapeados, não aparece nada porque não fez busca de enriquecimento

#down só tem um
mapped_down <- symbols_to_string_map(down_genes, "Down", string_db)
cat("Down genes mapeados no STRING:", nrow(mapped_down), "\n")
print(mapped_down)

##############################
# 7) Exportar resultados
##############################
dir.create("string_enrichment", showWarnings = FALSE)
rio::export(res_up$mapped,   "string_enrichment/STRING_mapping_Up.tsv")
rio::export(res_down$mapped, "string_enrichment/STRING_mapping_Down.tsv")
rio::export(enr_up,   "string_enrichment/STRING_enrichment_Up_all.tsv")

##############################
# 8) Salvar workspace
##############################
save.image(file = "workspace_DEG_STRING.RData")
gc()


######dotplot
library(ggplot2)
library(dplyr)
library(stringr)

plot_string_enrichment_dotplot <- function(enr_df, top_n = 50, title = "Enrichment Dotplot") {
  
  if(nrow(enr_df) == 0){
    warning("Dataframe de enriquecimento vazio.")
    return(NULL)
  }
  
  # Seleciona top N termos pelo FDR (menor primeiro)
  top_terms <- enr_df %>%
    arrange(fdr) %>%
    slice(1:top_n) %>%
    mutate(
      description_wrapped = str_wrap(description, width = 40),
      # criar fator �nico usando ordem pelo FDR
      description_factor = factor(description_wrapped, levels = rev(unique(description_wrapped)))
    )
  
  # Dotplot
  p <- ggplot(top_terms, aes(
    x = -log10(fdr),
    y = description_factor,
    size = number_of_genes,
    color = -log10(fdr)
  )) +
    geom_point(alpha = 0.9) +
    scale_color_gradient(low = "#74add1", high = "#d73027") +
    labs(
      x = "-log10(FDR)",
      y = "Biological Process / Term",
      title = title,
      size = "Genes associados",
      color = "-log10(FDR)"
    ) +
    theme_bw(base_size = 13)
  
  return(p)
}

# Uso para seu enr_up
dotplot_up <- plot_string_enrichment_dotplot(enr_up, top_n = 50, title = "STRING Enrichment - Genes Up")
print(dotplot_up)
ggsave("dotplot_up_STRING.png", plot = dotplot_up, width = 10, height = 8, dpi = 300)

###########################
#segunda etapa

library(dplyr)
library(ggplot2)

# Ler resultados
top5 <- read.csv('top5_condicoes_por_via.csv', fileEncoding = 'UTF-8')

# Ver resultado
print(top5)

# Heatmap no R (opcional — o PNG já foi gerado aqui)
matriz <- read.csv('matriz_similaridade_completa.csv',
                   row.names = 1,
                   fileEncoding = 'UTF-8')
