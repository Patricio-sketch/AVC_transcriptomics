# AVC Transcriptomics — Ischemic Stroke Pathway × NANDA-I Mapping

**Connecting peripheral-blood transcriptomics of ischemic stroke to NANDA-I nursing diagnoses via bioinformatics and clinical NLP.**

![Python 3.11](https://img.shields.io/badge/Python-3.11-blue) ![R 4.x](https://img.shields.io/badge/R-4.x-276DC3)

---

## Overview

Ischemic stroke triggers a rapid and multilayered transcriptional response in peripheral blood, encompassing innate immune activation, lymphocyte suppression, and systemic inflammatory signaling. This project bridges translational bioinformatics and clinical nursing science by systematically mapping the biological pathways enriched in peripheral blood gene expression data (GEO dataset GSE16561) to NANDA-I nursing diagnoses — the standardized clinical taxonomy used in nursing care planning worldwide. The analytical workflow combines differential expression analysis (limma), protein-protein interaction network construction (STRINGdb), and a hybrid semantic embedding strategy that fuses clinical domain knowledge from ClinicalBERT with biomedical sentence-level representations from a fine-tuned SBERT model. The end result is a ranked, threshold-filtered table of pathway–diagnosis pairs that can inform evidence-based nursing assessment in the acute stroke setting, grounded in molecular evidence rather than clinical intuition alone.

---

## Architecture

The pipeline operates in two sequential stages, designed so that each stage's outputs become the inputs to the next.

**Stage 1 — Differential expression and functional enrichment (R, `script.R`).** Raw microarray expression data from GSE16561 (Illumina HumanHT-12, peripheral blood, ischemic stroke vs. healthy controls) is downloaded via GEOquery, log2-transformed if necessary, and quantile-normalized. Differential expression is computed with limma using empirical Bayes moderation, with significant genes defined by FDR < 0.05 and |log2FC| > 1. Up- and down-regulated gene sets are analyzed separately because they represent distinct biological phenomena: the up-regulated signature reflects innate immune activation and neutrophil degranulation, while the down-regulated signature captures post-stroke adaptive immunodepression. Both sets are submitted to STRINGdb v11.5 (combined score ≥ 400) to retrieve protein-protein interaction networks and perform functional enrichment across Gene Ontology, KEGG, Reactome, and related annotation categories. All enrichment results are exported to `string_enrichment/` and serve as the structured biological vocabulary for Stage 2.

**Stage 2 — Semantic mapping via hybrid embeddings (Python, `clinicalbert_vias_AVC.ipynb`).** Each enriched biological pathway term and each NANDA-I nursing diagnosis (encoded as a full, structured natural-language description) is embedded independently using two complementary models: (1) the `[CLS]` token from `medicalai/ClinicalBERT`, a BERT model pre-trained on clinical notes and discharge summaries; and (2) sentence-level embeddings from a fine-tuned S-PubMedBert model (`pritamdeka/S-PubMedBert-MS-MARCO-SCIFACT`, domain-adapted in `finetune_sbert_nanda.ipynb`). Rather than concatenating the heterogeneous latent vectors and computing a single cosine distance — which would distort the metric geometry of each embedding space — the pipeline computes two separate cosine similarity matrices and combines them via weighted average in similarity space (late fusion, α = 0.5). This preserves the integrity of each model's learned metric while allowing their complementary signals — one capturing clinical linguistic context, the other capturing sentential semantic structure — to be combined in a geometrically stable way. Pathway–diagnosis pairs that exceed a cosine similarity threshold of 0.65 are retained as valid mappings and exported to `outputs/nanda/`.

---

## Repository Structure

```
AVC_transcriptomics/
│
├── script.R                          # R pipeline: DEA → PPI → enrichment → NANDA import
├── clinicalbert_vias_AVC.ipynb       # Main Python notebook: embeddings, fusion, NANDA mapping
├── finetune_sbert_nanda.ipynb        # Optional: fine-tunes S-PubMedBert on synthetic pairs
├── configurar_ambiente.R             # R package installation helper
├── environment.yml                   # Conda environment (Python 3.11, torch, sentence-transformers)
├── requirements.txt                  # pip requirements for Python stage
│
├── geo_cache/                        # Auto-downloaded GEO files (gitignored, ~24 MB)
│   ├── GPL6883.soft.gz
│   └── GSE16561_series_matrix.txt.gz
│
├── string_enrichment/                # STRINGdb outputs (versioned — small, essential for traceability)
│   ├── DEGs_anotados.tsv             # Full annotated DEG table (FDR < 0.05)
│   ├── STRING_enrichment_Up_all.tsv  # All enriched terms, up-regulated gene set
│   ├── STRING_enrichment_Down_all.tsv# All enriched terms, down-regulated gene set
│   ├── STRING_mapping_Up.tsv         # Gene → STRING ID mapping (up)
│   └── STRING_mapping_Down.tsv       # Gene → STRING ID mapping (down)
│
├── rds/                              # Serialized R objects for session resumption
│   ├── deg_anotado.rds               # Annotated DEG table
│   ├── enr_up.rds                    # Enrichment result (up)
│   ├── enr_down.rds                  # Enrichment result (down)
│   └── ppi_graph_cc.rds              # Largest connected component of PPI graph
│
├── outputs/
│   ├── figures/                      # Generated plots (gitignored, reproducible)
│   │   ├── Volcano.png               # Volcano plot: stroke vs. controls
│   │   ├── PPI.png                   # STRING PPI network, largest connected component
│   │   ├── dotplot_up_STRING.png     # Enrichment dotplot, up-regulated genes
│   │   ├── dotplot_down_STRING.png   # Enrichment dotplot, down-regulated genes
│   │   ├── heatmap_nanda_up.png      # Cosine similarity heatmap: up pathways × NANDA-I
│   │   ├── heatmap_nanda_down.png    # Cosine similarity heatmap: down pathways × NANDA-I
│   │   └── barplot_nanda_summary.png # Valid pairs per NANDA-I diagnosis, by regulation
│   │
│   ├── nanda/                        # NANDA-I mapping tables (gitignored, reproducible)
│   │   ├── nanda_mapping_completo_threshold65.csv  # All valid pairs (≥ 0.65), both directions
│   │   ├── nanda_mapping_up_threshold65.csv        # Valid pairs, up-regulated only
│   │   ├── nanda_mapping_down_threshold65.csv      # Valid pairs, down-regulated only
│   │   ├── nanda_ranking_top5_completo.csv         # Top-5 per pathway, no threshold (audit)
│   │   ├── top5_condicoes_por_via.csv              # Top-5 diagnoses per pathway (wide format)
│   │   ├── matriz_similaridade_up_nanda.csv        # Full cosine similarity matrix (up)
│   │   ├── matriz_similaridade_down_nanda.csv      # Full cosine similarity matrix (down)
│   │   └── matriz_similaridade_completa.csv        # Combined matrix (both directions)
│   │
│   └── sbert_nanda_finetuned/        # Fine-tuned SBERT model (versioned — ~418 MB, GPU artifact)
│       ├── model.safetensors         # Model weights
│       ├── tokenizer.json
│       ├── config.json
│       └── ...
│
└── .github/
    └── workflows/
        └── check-environment.yml     # CI: validates conda environment and package availability
```

---

## Quickstart

### R — Differential Expression and Enrichment

**Prerequisites:** R ≥ 4.0 with the following packages available from Bioconductor and CRAN.

```r
# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GEOquery", "limma"))

# Install CRAN packages
install.packages(c("STRINGdb", "ggplot2", "ggrepel", "ggraph",
                   "igraph", "dplyr", "stringr", "rio"))
```

**Run:**

```bash
Rscript script.R
```

Or, in RStudio: `Session > Set Working Directory > To Source File Location`, then `source("script.R")`.

The script downloads GSE16561 into `geo_cache/` on first run (internet access required), produces enrichment tables in `string_enrichment/`, figures in `outputs/figures/`, and R objects in `rds/`. Subsequent runs use the cached GEO files.

---

### Python — Embedding, Fusion, and NANDA-I Mapping

**Prerequisites:** Python 3.11. Install via conda (recommended) or pip.

```bash
# Option A — conda
conda env create -f environment.yml
conda activate avc_transcriptomics

# Option B — pip
pip install -r requirements.txt
```

**Execution order:**

1. **(Optional) Fine-tune the SBERT model** — run `finetune_sbert_nanda.ipynb` cell by cell. This adapts `pritamdeka/S-PubMedBert-MS-MARCO-SCIFACT` to the stroke × NANDA-I domain using synthetic training pairs and saves the result to `outputs/sbert_nanda_finetuned/`. Requires a GPU for reasonable runtime; the fine-tuned model is already present in this repository, so this step can be skipped.

2. **Run the main mapping notebook** — open and execute `clinicalbert_vias_AVC.ipynb` from top to bottom. This loads both models, generates embeddings, applies late fusion, and exports all NANDA-I mapping tables and heatmaps.

> **Security note — HuggingFace token:** The notebook contains a placeholder `HF_TOKEN` variable. The `medicalai/ClinicalBERT` model is publicly accessible and does not require authentication; a token is only needed if you access gated models. **Never hardcode a token in the notebook.** Set it as an environment variable instead:
>
> ```bash
> export HF_TOKEN=hf_your_token_here
> ```
>
> Then read it in the notebook with `os.environ.get("HF_TOKEN", "")`. Before committing, verify the token field is empty or cleared.

---

## Outputs

### NANDA-I Mapping Tables (`outputs/nanda/`)

| File | Description | Format |
|---|---|---|
| `nanda_mapping_completo_threshold65.csv` | All valid pathway–diagnosis pairs (cosine similarity ≥ 0.65), both regulatory directions | CSV |
| `nanda_mapping_up_threshold65.csv` | Valid pairs restricted to up-regulated gene pathways | CSV |
| `nanda_mapping_down_threshold65.csv` | Valid pairs restricted to down-regulated gene pathways | CSV |
| `nanda_ranking_top5_completo.csv` | Top-5 NANDA-I candidates per pathway, without threshold filtering (for calibration and audit) | CSV |
| `top5_condicoes_por_via.csv` | Top-5 diagnoses per pathway in wide format, both directions | CSV |
| `matriz_similaridade_up_nanda.csv` | Full pairwise cosine similarity matrix: up-regulated pathways × all NANDA-I diagnoses | CSV |
| `matriz_similaridade_down_nanda.csv` | Full pairwise cosine similarity matrix: down-regulated pathways × all NANDA-I diagnoses | CSV |
| `matriz_similaridade_completa.csv` | Combined similarity matrix (both regulatory directions) | CSV |

### Figures (`outputs/figures/`)

| File | Description | Format |
|---|---|---|
| `Volcano.png` | Volcano plot of all DEGs: blue = up-regulated, purple = down-regulated, grey = non-significant | PNG (300 dpi) |
| `PPI.png` | STRING PPI network (largest connected component, Fruchterman-Reingold layout, edge width = combined score) | PNG (300 dpi) |
| `dotplot_up_STRING.png` | Enrichment dotplot: top-50 biological terms enriched in up-regulated genes, ranked by FDR | PNG (300 dpi) |
| `dotplot_down_STRING.png` | Enrichment dotplot: top-50 biological terms enriched in down-regulated genes, ranked by FDR | PNG (300 dpi) |
| `heatmap_nanda_up.png` | Cosine similarity heatmap: up-regulated pathways (rows) × NANDA-I diagnoses (columns) | PNG (300 dpi) |
| `heatmap_nanda_down.png` | Cosine similarity heatmap: down-regulated pathways (rows) × NANDA-I diagnoses (columns) | PNG (300 dpi) |
| `barplot_nanda_summary.png` | Bar chart of valid pathway–diagnosis pairs per NANDA-I code, split by regulatory direction | PNG (300 dpi) |

---

## Methods Summary

**Dataset and preprocessing.** GSE16561 (NCBI GEO) contains peripheral blood gene expression profiles measured on the Illumina HumanHT-12 platform from ischemic stroke patients and healthy controls. Raw expression values are checked for scale and log2-transformed when necessary, then quantile-normalized across samples using limma's `normalizeBetweenArrays`. Probe annotations are drawn from the platform's feature data; when multiple probes map to the same gene symbol, the probe with the largest absolute log2 fold change is retained as the most informative representative. Samples are classified into groups using the phenotype description field of the GEO metadata, and any unclassifiable samples are excluded before modeling.

**Differential expression and functional enrichment.** A linear model is fitted per probe using `lmFit` with a contrast of AVC − Controle, and empirical Bayes moderation is applied via `eBayes` to stabilize variance estimates across genes — a critical adjustment for microarray datasets with moderate sample sizes. Significant DEGs are defined by Benjamini-Hochberg FDR < 0.05 and |log2FC| > 1. The up- and down-regulated gene sets are submitted separately to STRINGdb v11.5 (species 9606, combined score ≥ 400) for protein-protein interaction network retrieval and functional enrichment analysis. Biological annotation categories used for enrichment include Gene Ontology Biological Process and Molecular Function, Reactome (RCTM), InterPro, SMART, Pfam, HPO, COMPARTMENTS, and UniProt Keywords; PMID co-citation entries are explicitly excluded to avoid mixing literature-based and annotation-based evidence.

**Hybrid embedding and late fusion.** NANDA-I nursing diagnoses are encoded as full natural-language descriptions that include the diagnosis label, its defining characteristics, and related factors — a deliberate choice to maximize semantic surface for the embedding models rather than relying on short labels alone. Each pathway term and diagnosis description is embedded in two independent semantic spaces: the `[CLS]` token representation from `medicalai/ClinicalBERT` (L2-normalized), and the sentence-level embedding from a domain-adapted S-PubMedBert model. The `[CLS]` token is preferred over mean pooling for ClinicalBERT because ClinicalBERT is pre-trained with a next-sentence prediction objective that concentrates global sentence-level information in the `[CLS]` position; mean pooling, by contrast, dilutes this signal across token positions and is more appropriate for models explicitly trained with a pooling objective. The two cosine similarity matrices are averaged with equal weight (α = 0.5) in similarity space rather than in vector space, avoiding metric distortion that would arise from concatenating vectors of different dimensionality and training origin.

**Fine-tuning with synthetic supervision.** The base S-PubMedBert model is fine-tuned on 65 synthetic training pairs constructed from actual STRING enrichment pathway terms (drawn from the up- and down-regulated gene sets of GSE16561) paired with NANDA-I diagnosis descriptions. Each pair carries a float label encoding relevance: high-relevance pairs (label ≥ 0.85) capture direct immune-nursing alignment (e.g., "Neutrophil degranulation" × "Risk for infection"), moderate pairs (0.50–0.84) encode secondary clinical links, and hard-negative pairs (label ≤ 0.20) contrast immune pathways against functionally unrelated diagnoses such as mobility or communication deficits. Training uses `CosineSimilarityLoss` for 5 epochs with linear warmup (10 steps), directly optimizing the embedding space for cosine similarity prediction in the stroke–NANDA-I domain.

---

## Limitations

- **Conservative similarity threshold.** The cutoff of 0.65 was selected heuristically to exclude weak lexical co-occurrence, but it has not been calibrated against a held-out clinical gold standard. Pairs with scores between 0.60 and 0.65 may represent clinically valid mappings that are excluded by the current threshold.

- **Synthetic training pairs.** The fine-tuning supervision signal consists entirely of pairs generated by the research team, not by certified nursing clinicians. Label values represent expert biological judgment about pathway–nursing diagnosis relevance but have not been validated through structured clinical review or inter-annotator agreement protocols.

- **No external clinical validation.** The pathway–NANDA-I mappings have not been evaluated against real patient care plans, clinical chart data, or nursing expert consensus. The output should be treated as a hypothesis-generating tool rather than a clinical decision support system.

- **Single cohort, single platform.** The entire transcriptomic evidence base is derived from one GEO dataset (GSE16561) measured on one microarray platform (Illumina HumanHT-12). Findings may not generalize across stroke subtypes, disease stages, measurement platforms, or populations with different demographic profiles.

- **Equal model weighting.** The late fusion weight α = 0.5 assigns equal contribution to ClinicalBERT and S-PubMedBert without empirical optimization. Optimal weighting for this specific domain has not been determined, and different α values may yield substantially different rankings for borderline pairs.

---

## Citation / License

If you use this code or pipeline in your research, please cite:

```bibtex
@misc{avc_transcriptomics,
  author       = {Patricio-sketch},
  title        = {{AVC Transcriptomics}: Ischemic Stroke Pathway to NANDA-I Nursing Diagnosis Mapping},
  year         = {2025},
  howpublished = {\url{https://github.com/Patricio-sketch/AVC_transcriptomics}},
  note         = {Preprint / work in progress}
}
```

This project is released under the **MIT License**. See `LICENSE` for details.

> The GEO dataset GSE16561 is publicly available and subject to its own data use terms. The NANDA-I taxonomy is a proprietary classification system owned by NANDA International, Inc.; the diagnosis descriptions used in this project are derived from published literature for research purposes only.
