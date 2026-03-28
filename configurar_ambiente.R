# Instala dependências R do projeto (rode uma vez por máquina).
# Uso: Rscript configurar_ambiente.R
#   ou no RStudio: source("configurar_ambiente.R")

options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Bioconductor: GEOquery, limma, STRINGdb.
# Demais no CRAN; BiocManager::install cobre ambos.
pacotes <- c(
  "GEOquery",
  "limma",
  "STRINGdb",
  "ggplot2",
  "ggrepel",
  "dplyr",
  "igraph",
  "ggraph",
  "rio",
  "stringr",
  "rstudioapi"
)

BiocManager::install(pacotes, ask = FALSE, update = FALSE)

message("Pacotes instalados: ", paste(pacotes, collapse = ", "))
message("Concluído. Pode executar script.R.")
