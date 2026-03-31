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
  "stringr"
)

BiocManager::install(pacotes, ask = FALSE, update = FALSE)
if (!requireNamespace("rstudioapi", quietly = TRUE)) install.packages("rstudioapi") # FIX [P3]: instala rstudioapi separadamente apenas quando necessario para set_project_wd()

message("Pacotes instalados: ", paste(pacotes, collapse = ", "))
falhos <- pacotes[!sapply(pacotes, requireNamespace, quietly = TRUE)] # FIX [P2]: verifica se todos os pacotes de analise ficaram disponiveis apos a instalacao
if (length(falhos) > 0) { # FIX [P2]: interrompe a configuracao quando houver pacotes faltantes
  stop("Pacotes não instalados corretamente: ", paste(falhos, collapse = ", "))
}
message("Concluído. Pode executar script.R.")
