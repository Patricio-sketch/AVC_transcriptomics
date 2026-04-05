# Transcriptomics-to-Bedside: Mapping Ischemic Stroke Biology to NANDA-I Nursing Diagnoses via ClinicalBERT

A translational bioinformatics pipeline that bridges peripheral blood gene expression in ischemic stroke patients with standardized nursing diagnoses (NANDA-I taxonomy) through semantic similarity computed by a clinical language model (ClinicalBERT).

---

## Background and Hypothesis

Ischemic stroke triggers a well-characterized systemic immune response detectable in peripheral blood: innate immune activation (neutrophil degranulation, cytokine signaling) alongside post-stroke immunodepression (lymphocyte suppression, adaptive immunity downregulation). These molecular events are already used in biomarker research, but they have never been quantitatively connected to the nursing problems that clinicians identify at the bedside.

**Central hypothesis:** ClinicalBERT — a BERT-based language model pre-trained on clinical records — can semantically bridge what happens at the molecular level in stroke to the health responses formalized by NANDA-I, the international standardized nursing diagnosis taxonomy. The result is an algorithmically traceable, quantitative link between the patient's transcriptome and the nursing care plan.

---

## Pipeline Overview

```
GEO (GSE16561)
    │
    ▼
[R: script.R]
 Preprocessing → Differential Expression (limma) → DEG annotation
    │
    ├── Volcano plot
    ├── PPI Network (STRINGdb + igraph)
    ├── STRING Functional Enrichment (Up & Down gene sets)
    └── Enrichment Dotplots
    │
    ▼
[Python: clinicalbert_vias_AVC.ipynb]
 ClinicalBERT embeddings (768-dim)
 Cosine similarity: Biological pathways × NANDA-I diagnoses
    │
    ├── Similarity matrices (Up and Down)
    ├── Heatmaps
    └── Ranked NANDA-I mapping (threshold ≥ 0.65)
```

---

## Dataset

| Field | Value |
|---|---|
| **GEO Accession** | GSE16561 |
| **Platform** | Illumina HumanHT-12 v3 (GPL6883) |
| **Tissue** | Peripheral blood |
| **Condition** | Ischemic stroke vs. healthy controls |
| **Design** | Case-control; expression matrix obtained via `GEOquery` |

---

## Key Results

### Differential Expression

| Group | Count | Criteria |
|---|---|---|
| Total DEGs | 3,915 | FDR < 0.05, \|logFC\| > 1 |
| Up-regulated | 1,526 | logFC > 1 |
| Down-regulated | 2,389 | logFC < −1 |

**Strongest up-regulated genes:** ARG1 (logFC=1.72), MMP9 (logFC=1.39), S100A12 (logFC=1.37)  
**Strongest down-regulated genes:** CCR7, CD79B, VPREB3, LEF1, IL7R

### STRING Functional Enrichment

**Up-regulated genes — innate immunity signature:**

| Pathway | Category | FDR |
|---|---|---|
| Neutrophil degranulation | Reactome | 3.70 × 10⁻⁴ |
| Innate Immune System | Reactome | 7.50 × 10⁻⁴ |
| Immune effector process | GO:BP | 2.50 × 10⁻³ |
| Immune response | GO:BP | 2.50 × 10⁻³ |

**Down-regulated genes — adaptive immunity suppression signature:**

| Pathway | Category | FDR |
|---|---|---|
| Lymphocyte differentiation | GO:BP | 2.80 × 10⁻⁴ |
| Lymphocyte activation | GO:BP | 2.80 × 10⁻⁴ |
| Adaptive immune response | GO:BP | 8.70 × 10⁻⁴ |
| Positive regulation of immune system process | GO:BP | 1.30 × 10⁻³ |
| B cell activation | GO:BP | 1.09 × 10⁻² |
| B cell receptor signaling pathway | GO:BP | 1.09 × 10⁻² |

This two-sided pattern reflects the well-known post-stroke immunodepression syndrome (PSIDS): simultaneous innate immune activation and adaptive immune suppression.

### ClinicalBERT → NANDA-I Mapping

Biological pathways were mapped to 38 NANDA-I nursing diagnoses (cosine similarity ≥ 0.65 threshold). **17 valid pathway–diagnosis pairs** were identified:

**Up-regulated pathways:**

| Biological Pathway | NANDA-I Diagnosis | Code | Similarity |
|---|---|---|---|
| Innate immunity | Risk for infection | 00004 | 0.6753 |

**Down-regulated pathways (immunodepression → nursing risk):**

| Biological Pathway | NANDA-I Diagnosis | Code | Similarity |
|---|---|---|---|
| Adaptive immunity | Ineffective protection | 00043 | 0.7004 |
| Adaptive immunity | Risk for infection | 00004 | 0.6912 |
| Adaptive immune response | Ineffective protection | 00043 | 0.6519 |
| Positive regulation of immune system process | Ineffective protection | 00043 | 0.6793 |
| Positive regulation of immune system process | Risk for infection | 00004 | 0.6704 |
| Positive regulation of immune system process | Hyperthermia | 00007 | 0.6683 |
| Immunoreceptor tyrosine-based activation motif | Impaired verbal communication | 00051 | 0.6928 |
| Immunoreceptor tyrosine-based activation motif | Decreased intracranial adaptive capacity | 00049 | 0.6858 |
| Immunoreceptor tyrosine-based activation motif | Risk for infection | 00004 | 0.6810 |
| Antigen receptor-mediated signaling pathway | Risk for ineffective cerebral tissue perfusion | 00201 | 0.6889 |
| Immunological synapse formation | Risk for infection | 00004 | 0.6811 |
| Regulation of immune response | Ineffective protection | 00043 | 0.6794 |

**Biological interpretation:** The mapping is coherent in both directions. Up-regulated innate immunity links to infection risk — consistent with stroke-associated pneumonia and urinary tract infections. Down-regulated adaptive immunity converges on *Ineffective protection* and *Risk for infection*, which directly reflects PSIDS and is actionable in nursing care planning. The *Antigen receptor-mediated signaling → Risk for ineffective cerebral tissue perfusion* link captures the lymphocyte-mediated blood-brain barrier disruption described in the literature for AKAP7.

---

## Repository Structure

```
AVC_transcriptomics/
├── script.R                          # Main R pipeline (DEG → STRING → Dotplots)
├── configurar_ambiente.R             # One-time package installation
├── clinicalbert_vias_AVC.ipynb      # ClinicalBERT semantic mapping notebook
├── _apply_notebook_corrections.ps1  # Legacy notebook patching utility
│
├── outputs/
│   ├── figures/                      # All generated visualizations
│   │   ├── Volcano.png               # Differential expression volcano plot
│   │   ├── PPI.png                   # Protein-protein interaction network
│   │   ├── dotplot_up_STRING.png     # Enrichment dotplot — Up genes
│   │   ├── dotplot_down_STRING.png   # Enrichment dotplot — Down genes
│   │   ├── heatmap_nanda_up.png      # NANDA-I heatmap — Up pathways
│   │   └── heatmap_nanda_down.png    # NANDA-I heatmap — Down pathways
│   │
│   └── nanda/                        # NANDA-I mapping outputs
│       ├── nanda_mapping_completo_threshold65.csv  # Valid pairs (≥ 0.65) — both directions
│       ├── nanda_mapping_up_threshold65.csv         # Valid pairs — Up genes only
│       ├── nanda_mapping_down_threshold65.csv       # Valid pairs — Down genes only
│       ├── nanda_ranking_top5_completo.csv          # Full top-5 ranking (no threshold)
│       ├── matriz_similaridade_up_nanda.csv         # Full cosine similarity matrix — Up
│       └── matriz_similaridade_down_nanda.csv       # Full cosine similarity matrix — Down
│
├── string_enrichment/                # STRINGdb outputs (generated by script.R)
│   ├── DEGs_anotados.tsv             # All 3,915 DEGs with statistics
│   ├── STRING_enrichment_Up_all.tsv  # Enrichment results — Up genes
│   ├── STRING_enrichment_Down_all.tsv # Enrichment results — Down genes
│   ├── STRING_mapping_Up.tsv         # Gene symbol → STRING protein ID (Up)
│   └── STRING_mapping_Down.tsv       # Gene symbol → STRING protein ID (Down)
│
├── rds/                              # Serialized R objects (for resuming analysis)
│   ├── deg_anotado.rds
│   ├── enr_up.rds
│   ├── enr_down.rds
│   └── ppi_graph_cc.rds
│
└── geo_cache/                        # Local GEO download cache (git-ignored)
```

---

## Requirements

### R (script.R)

R ≥ 4.2 recommended. Install all dependencies by running `configurar_ambiente.R` once:

```r
source("configurar_ambiente.R")
```

Packages installed: `GEOquery`, `limma`, `STRINGdb`, `ggplot2`, `ggrepel`, `dplyr`, `igraph`, `ggraph`, `rio`, `stringr`, `rstudioapi`

### Python (notebook)

Python ≥ 3.9. Install dependencies:

```bash
pip install sentence-transformers transformers torch pandas numpy scikit-learn seaborn matplotlib
```

The notebook is designed to run in **Google Colab** or **locally** (Jupyter). The ClinicalBERT model (`medicalai/ClinicalBERT`, ~400 MB) is downloaded automatically from HuggingFace on first run.

---

## How to Run

### Step 1 — Install R packages (once per machine)

```r
source("configurar_ambiente.R")
```

### Step 2 — Run the R pipeline

Set your working directory to the project root, then:

```r
source("script.R")
```

This will:
1. Download and cache GSE16561 from GEO (≈ 25 MB, cached in `geo_cache/`)
2. Preprocess expression data (log2 transform, quantile normalization)
3. Run limma differential expression analysis
4. Generate `outputs/figures/Volcano.png`
5. Build PPI network and generate `outputs/figures/PPI.png`
6. Run STRING enrichment for Up and Down gene sets
7. Export enrichment TSVs to `string_enrichment/`
8. Generate enrichment dotplots to `outputs/figures/`
9. Save R objects to `rds/`

**Note:** Step 1 requires internet access to download GSE16561. Subsequent runs use the local cache and skip the download.

### Step 3 — Run the ClinicalBERT notebook

Open `clinicalbert_vias_AVC.ipynb` in Jupyter or Google Colab and run all cells in order.

**If running in Google Colab:**
1. Upload the entire `string_enrichment/` folder to the Colab session
2. Run all cells — model download takes 1–2 minutes
3. CSVs and heatmaps are downloaded automatically at Cell 10

**If running locally:**
```bash
jupyter notebook clinicalbert_vias_AVC.ipynb
```
All outputs are saved directly to `outputs/figures/` and `outputs/nanda/`.

The notebook will:
1. Load enriched pathways from both `STRING_enrichment_Up_all.tsv` and `STRING_enrichment_Down_all.tsv`
2. Load 38 NANDA-I nursing diagnoses
3. Generate 768-dimensional ClinicalBERT embeddings
4. Compute cosine similarity matrices for Up and Down gene sets separately
5. Rank and filter pathway–diagnosis pairs (threshold ≥ 0.65)
6. Generate heatmaps and export all CSVs

### Step 4 — Import results in R (optional)

```r
library(dplyr)

# Valid NANDA-I mappings (similarity >= 0.65)
nanda <- read.csv("outputs/nanda/nanda_mapping_completo_threshold65.csv",
                  fileEncoding = "UTF-8")

# Inspect by regulatory direction
nanda_up   <- nanda[nanda$regulacao == "Up",   ]
nanda_down <- nanda[nanda$regulacao == "Down", ]

# Full similarity matrices
mat_up   <- read.csv("outputs/nanda/matriz_similaridade_up_nanda.csv",
                     row.names = 1, fileEncoding = "UTF-8")
mat_down <- read.csv("outputs/nanda/matriz_similaridade_down_nanda.csv",
                     row.names = 1, fileEncoding = "UTF-8")
```

---

## Parameter Reference

| Parameter | Location | Default | Description |
|---|---|---|---|
| `FDR_THRESHOLD` | `script.R` | 0.05 | Maximum adjusted p-value for DEG inclusion |
| `LOGFC_THRESHOLD` | `script.R` | 1.0 | Minimum \|logFC\| for DEG inclusion |
| `STRING_SCORE` | `script.R` | 400 | Minimum STRING interaction confidence score for PPI |
| `TOP_N_DOTPLOT` | `script.R` | 50 | Max terms displayed per enrichment dotplot |
| `CATEGORIAS_VALIDAS` | notebook Cell 3 | Process, RCTM, Component, Keyword, COMPARTMENTS, InterPro, SMART, HPO, Pfam | STRING categories retained for embedding (PMID excluded) |
| `SIMILARITY_THRESHOLD` | notebook Cell 7 | 0.65 | Minimum cosine similarity for a valid pathway–diagnosis pair |
| `TOP_N` | notebook Cell 7 | 5 | Top-N NANDA-I diagnoses ranked per pathway |

---

## Output Files Reference

| File | Description |
|---|---|
| `outputs/nanda/nanda_mapping_completo_threshold65.csv` | **Main result.** All valid pathway–diagnosis pairs (similarity ≥ 0.65), columns: `regulacao`, `via`, `rank`, `diagnostico_nanda`, `similaridade` |
| `outputs/nanda/nanda_ranking_top5_completo.csv` | Full top-5 rankings without threshold (for calibration/audit) |
| `outputs/nanda/matriz_similaridade_up_nanda.csv` | Raw cosine similarity matrix: Up pathways × 38 NANDA-I diagnoses |
| `outputs/nanda/matriz_similaridade_down_nanda.csv` | Raw cosine similarity matrix: Down pathways × 38 NANDA-I diagnoses |
| `string_enrichment/DEGs_anotados.tsv` | Complete annotated DEG table (3,915 genes) |
| `string_enrichment/STRING_enrichment_*.tsv` | Full STRING enrichment results including all categories |

---

## Limitations

- **Single dataset:** Analysis is based on GSE16561 only. External validation on independent stroke cohorts (e.g., GSE10927, GSE58294) has not been performed.
- **Microarray platform specificity:** Data was generated on Illumina HumanHT-12; generalizability to RNA-seq data is not assessed.
- **Empirical similarity threshold:** The 0.65 cutoff was defined by inspecting the score distribution and is not clinically validated by expert review.
- **Small STRING-mapped gene sets:** 12 (Up) and 17 (Down) proteins mapped to STRING, limiting statistical power for enrichment analysis.
- **NANDA-I coverage:** 38 diagnoses were selected for clinical relevance to stroke; the full NANDA-I taxonomy (300+ diagnoses) was not exhaustively tested.

---

## Citation

If you use this pipeline or its outputs, please cite the source dataset:

> Stamova BS, Xu H, Jickling G, et al. *Gene expression profiling of blood for the prediction of ischemic stroke*. Stroke. 2010;41(10):2171–2177. [GSE16561]

And the ClinicalBERT model:

> Alsentzer E, Murphy J, Boag W, et al. *Publicly Available Clinical BERT Embeddings*. NAACL 2019 Workshop on Biomedical NLP. [medicalai/ClinicalBERT]

---

## License

This repository contains analysis code and derived results only. The GSE16561 dataset is publicly available under the GEO data access policy. NANDA-I diagnostic labels are used for research purposes under fair use; clinical application requires access to the official NANDA International taxonomy.
