# AVC Transcriptomics — Ischemic Stroke Pathway × NANDA-I Mapping

![Python 3.11](https://img.shields.io/badge/Python-3.11-blue) ![R 4.x](https://img.shields.io/badge/R-4.x-276DC3) ![Status](https://img.shields.io/badge/status-work%20in%20progress-yellow)

**Connecting peripheral-blood transcriptomics of ischemic stroke to NANDA-I nursing diagnoses via bioinformatics and clinical NLP.**

---

## Overview

Ischemic stroke triggers a rapid and multilayered transcriptional response in peripheral blood, encompassing innate immune activation, lymphocyte suppression, and systemic inflammatory signaling. This project bridges translational bioinformatics and clinical nursing science by systematically mapping the biological pathways enriched in peripheral blood gene expression data (GEO dataset GSE16561) to NANDA-I nursing diagnoses — the standardized clinical taxonomy used in nursing care planning worldwide.

The analytical workflow operates across **three sequential stages**. Stage 1 performs differential expression analysis (limma) and protein-protein interaction network construction (STRINGdb) in R. Stage 2 applies a hybrid semantic embedding strategy in Python, fusing clinical domain knowledge from ClinicalBERT with biomedical sentence-level representations from a fine-tuned SBERT model, and exports a ranked, threshold-filtered table of pathway–diagnosis pairs grounded in molecular evidence. Stage 3 (currently under active development) implements an LLM-as-judge cross-validation layer using **Claude Sonnet 4.6** (temperature = 0, seed = 42) to provide an independent mechanistic signal for each pathway–diagnosis pair that survives the Stage 2 threshold.

The motivation for Stage 3 is a documented **entropy collapse** in Stage 2: Shannon entropy of the per-pathway similarity score distribution over NANDA-I diagnoses reaches H/H_max = 0.960 across all 35 pathways, indicating that the embedding pipeline alone does not discriminate biologically distinct pathways with sufficient resolution for unqualified clinical claims. The stratified output of Stage 3 (high concordance / moderate / divergent) constitutes the primary validation object for the second paper in this line of research — a structured Delphi study with clinical nursing experts.

---

## Architecture

The pipeline operates in three sequential stages, designed so that each stage's outputs become the inputs to the next.

**Stage 1 — Differential expression and functional enrichment (R, `script.R`).** Raw microarray expression data from GSE16561 (Illumina HumanHT-12, peripheral blood, ischemic stroke vs. healthy controls) is downloaded via GEOquery, log2-transformed if necessary, and quantile-normalized. Differential expression is computed with limma using empirical Bayes moderation, with significant genes defined by FDR < 0.05 and |log2FC| > 1. Up- and down-regulated gene sets are analyzed separately because they represent distinct biological phenomena: the up-regulated signature reflects innate immune activation and neutrophil degranulation, while the down-regulated signature captures post-stroke adaptive immunodepression. Both sets are submitted to STRINGdb v11.5 (combined score ≥ 400) to retrieve protein-protein interaction networks and perform functional enrichment across Gene Ontology, KEGG, Reactome, and related annotation categories. All enrichment results are exported to `string_enrichment/` and serve as the structured biological vocabulary for Stage 2. A post-processing section (Section 7 of `script.R`) re-imports the Stage 2 NANDA-I mapping tables for additional R-side analyses including directional specificity indices and normalized coverage plots.

**Stage 2 — Semantic mapping via hybrid embeddings (Python, `clinicalbert_vias_AVC.ipynb`).** Each enriched biological pathway is semantically expanded into a structured natural-language description incorporating gene names and mechanistic context (188–255 characters), providing richer embedding surface than the raw pathway label alone. NANDA-I nursing diagnoses are encoded in two representations: a short label form (`nanda_labels`) for structural separation during training, and a full, structured description including defining characteristics and related factors (`nanda_full`) for embedding inference — a deliberate separation that prevents label leakage into the similarity space. Each expanded pathway and each NANDA-I description is embedded independently using two complementary models: the `[CLS]` token from `medicalai/ClinicalBERT`, and sentence-level embeddings from a fine-tuned S-PubMedBert model (`pritamdeka/S-PubMedBert-MS-MARCO-SCIFACT`, adapted in `finetune_sbert_nanda.ipynb`). The two cosine similarity matrices are averaged with equal weight (α = 0.5) in similarity space — late fusion — to preserve the metric geometry of each embedding space. The similarity threshold is calibrated by pairwise permutation (z = 18.4), corresponding to a cosine cutoff of 0.65. A Shannon entropy audit is performed per pathway on the resulting similarity distribution; the results, including `ratio_H_Hmax` and `alerta_colapso`, are exported to `discriminacao_por_via.csv`. The column `via_expandida` in the main output CSV records the full expanded pathway description used for embedding. Calibration metadata is stored in `threshold_calibration.csv` with fields `colapso_confirmado`, `entropia_ratio_medio`, `n_vias_em_colapso`, `spearman_finetuning`, and `pipeline_versao`.

**Stage 3 — LLM-based mechanistic cross-validation (Python, `llm_judge_nanda.ipynb` — to be created).** Each pathway–diagnosis pair surviving the Stage 2 threshold is submitted to **Claude Sonnet 4.6** (temperature = 0, seed = 42) for independent mechanistic evaluation. To prevent anchoring bias, the prompt is structured to elicit justification before score, in the order: mechanistic context → reasoning → score. Each response returns four structured fields: `justificativa` (free-text mechanistic rationale), `score_concordancia` (numeric, 0–10), `nivel_confianca` (categorical: high / moderate / low), and `flag_lexical` (boolean: flags pairs where similarity is driven primarily by shared terminology rather than mechanistic correspondence). Pairs are stratified into three tiers: high concordance, moderate concordance, and divergent. Estimated cost per full run: ~$5–6 with Sonnet 4.6 standard API; ~$2.50–3 with the Batch API. The stratified output (`nanda_judge_concordancia.csv`) is the primary validation object for the Delphi expert study.

---

## Repository Structure

```
AVC_transcriptomics/
│
├── script.R                          # Stage 1: DEA → PPI → enrichment → NANDA import (R)
├── clinicalbert_vias_AVC.ipynb       # Stage 2: embeddings, late fusion, NANDA mapping (Python)
├── finetune_sbert_nanda.ipynb        # Stage 2 prerequisite: fine-tunes S-PubMedBert
├── llm_judge_nanda.ipynb             # Stage 3: LLM-as-judge cross-validation (to be created)
├── configurar_ambiente.R             # R package installation helper
├── environment.yml                   # Conda environment (Python 3.11, torch, sentence-transformers)
├── requirements.txt                  # pip requirements for Python stages
├── LICENSE
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
├── rds/                              # Serialized R objects (gitignored — regenerated by script.R)
│   ├── deg_anotado.rds
│   ├── enr_up.rds
│   ├── enr_down.rds
│   └── ppi_graph_cc.rds
│
├── checkpoints/                      # Intermediate SBERT training checkpoints (gitignored)
│
├── outputs/
│   ├── figures/                      # Generated plots (gitignored — reproducible)
│   │   ├── Volcano.png
│   │   ├── PPI.png
│   │   ├── dotplot_up_STRING.png
│   │   ├── dotplot_down_STRING.png
│   │   ├── heatmap_nanda_up.png
│   │   ├── heatmap_nanda_down.png
│   │   ├── barplot_nanda_summary.png
│   │   ├── cobertura_normalizada.png
│   │   └── especificidade_direcional.png
│   │
│   ├── nanda/                        # NANDA-I mapping tables (gitignored — reproducible)
│   │   ├── nanda_mapping_completo_threshold65.csv  # All valid pairs (≥ 0.65); includes via_expandida
│   │   ├── nanda_mapping_up_threshold65.csv
│   │   ├── nanda_mapping_down_threshold65.csv
│   │   ├── nanda_ranking_top5_completo.csv
│   │   ├── top5_condicoes_por_via.csv
│   │   ├── matriz_similaridade_up_nanda.csv
│   │   ├── matriz_similaridade_down_nanda.csv
│   │   ├── matriz_similaridade_completa.csv
│   │   ├── threshold_calibration.csv               # Fields: colapso_confirmado,
│   │   │                                           #   entropia_ratio_medio, n_vias_em_colapso,
│   │   │                                           #   spearman_finetuning, pipeline_versao
│   │   ├── discriminacao_por_via.csv               # Fields include ratio_H_Hmax, alerta_colapso
│   │   ├── especificidade_direcional.csv
│   │   ├── cobertura_normalizada.csv
│   │   └── nanda_judge_concordancia.csv            # Stage 3 output (to be created)
│   │   [discriminacao_por_via_metodo_row_norm.csv — removed: mathematically incorrect method]
│   │
│   └── sbert_nanda_finetuned/        # Fine-tuned SBERT model (gitignored, ~418 MB)
│       ├── model.safetensors
│       ├── tokenizer.json
│       ├── config.json
│       ├── checkpoint_meta.json      # Pipeline version and training metadata
│       └── ...
│
└── .github/
    └── workflows/
        └── check-environment.yml     # CI: validates conda environment and package availability
```

---

## How to Run

The three stages must be executed in order. Each stage's outputs are required inputs for the next.

### Stage 1 — Differential Expression and Enrichment (R)

**Prerequisites:** R ≥ 4.0. Run `configurar_ambiente.R` once on any new machine to install all required packages.

```bash
Rscript configurar_ambiente.R   # one-time setup
Rscript script.R                # full Stage 1 + Stage 3 post-processing (Section 7)
```

Or in RStudio: `Session > Set Working Directory > To Source File Location`, then `source("script.R")`.

The script downloads GSE16561 into `geo_cache/` on first run (internet access required), produces enrichment tables in `string_enrichment/`, figures in `outputs/figures/`, and serialized objects in `rds/`. Subsequent runs use the cached GEO files. Section 7 of `script.R` imports Stage 2 outputs; it will fail gracefully if Stage 2 has not yet been run.

### Stage 2 — Embedding, Fusion, and NANDA-I Mapping (Python)

**Prerequisites:** Python 3.11. Install via conda (recommended) or pip.

```bash
# Option A — conda
conda env create -f environment.yml
conda activate avc_transcriptomics

# Option B — pip
pip install -r requirements.txt
```

**Execution order:**

1. **Fine-tune the SBERT model** — run `finetune_sbert_nanda.ipynb` cells 1–6. This adapts `pritamdeka/S-PubMedBert-MS-MARCO-SCIFACT` to the stroke × NANDA-I domain using `MultipleNegativesRankingLoss` and saves the result to `outputs/sbert_nanda_finetuned/`. Requires a GPU for reasonable runtime.

2. **Verify the checkpoint** — before running the main notebook, inspect `outputs/sbert_nanda_finetuned/checkpoint_meta.json` to confirm the correct pipeline version (`v4-expansao-semantica-MNRL`) and that training completed without errors.

3. **Run the main mapping notebook** — open and execute `clinicalbert_vias_AVC.ipynb` cells 1–13. This loads both models, generates semantically expanded pathway descriptions, computes embeddings, applies late fusion, runs the permutation-calibrated threshold (z = 18.4), performs the Shannon entropy audit, and exports all NANDA-I mapping tables and figures.

> **Security note — HuggingFace token:** Cell 1b reads `HF_TOKEN` from the environment and applies `.strip()`, so whitespace-only values are treated as missing. The `medicalai/ClinicalBERT` model is public; a token primarily improves Hub rate limits. **Never commit a real token.**
>
> **Option A — environment variable (shell):**
> ```bash
> export HF_TOKEN=hf_your_token_here
> ```
> Start Jupyter from the same terminal so the kernel inherits `HF_TOKEN`.
>
> **Option B — `python-dotenv` and a local `.env` file:** Create `.env` at the repository root:
> ```
> HF_TOKEN=hf_your_token_here
> ```
> The `.gitignore` already excludes `.env`. Before running Cell 1b, add these lines:
> ```python
> from dotenv import load_dotenv
> load_dotenv()
> ```

### Stage 3 — LLM-as-Judge Cross-Validation (Python) — *Next implementation*

Run `llm_judge_nanda.ipynb` (to be created) after Stage 2 is complete and `outputs/nanda/nanda_mapping_completo_threshold65.csv` is present. Requires an Anthropic API key. Estimated cost: ~$5–6 per full run (Sonnet 4.6 standard), ~$2.50–3 with Batch API.

---

## Outputs

### NANDA-I Mapping Tables (`outputs/nanda/`)

| File | Description | Format |
|---|---|---|
| `nanda_mapping_completo_threshold65.csv` | All valid pathway–diagnosis pairs (cosine similarity ≥ 0.65), both regulatory directions; includes `via_expandida` column | CSV |
| `nanda_mapping_up_threshold65.csv` | Valid pairs restricted to up-regulated gene pathways | CSV |
| `nanda_mapping_down_threshold65.csv` | Valid pairs restricted to down-regulated gene pathways | CSV |
| `nanda_ranking_top5_completo.csv` | Top-5 NANDA-I candidates per pathway, without threshold filtering (for calibration and audit) | CSV |
| `top5_condicoes_por_via.csv` | Top-5 diagnoses per pathway in wide format, both directions | CSV |
| `matriz_similaridade_up_nanda.csv` | Full pairwise cosine similarity matrix: up-regulated pathways × all NANDA-I diagnoses | CSV |
| `matriz_similaridade_down_nanda.csv` | Full pairwise cosine similarity matrix: down-regulated pathways × all NANDA-I diagnoses | CSV |
| `matriz_similaridade_completa.csv` | Combined similarity matrix (both regulatory directions) | CSV |
| `threshold_calibration.csv` | Permutation calibration metadata; fields include `colapso_confirmado`, `entropia_ratio_medio`, `n_vias_em_colapso`, `spearman_finetuning`, `pipeline_versao` | CSV |
| `discriminacao_por_via.csv` | Per-pathway Shannon entropy audit; fields include `ratio_H_Hmax` and `alerta_colapso` | CSV |
| `especificidade_direcional.csv` | Directional specificity index per NANDA-I diagnosis | CSV |
| `cobertura_normalizada.csv` | Normalized coverage: valid pairs per pathway per NANDA-I diagnosis | CSV |
| `nanda_judge_concordancia.csv` | Stage 3 LLM-as-judge output: fields `justificativa`, `score_concordancia`, `nivel_confianca`, `flag_lexical`, stratification tier *(to be created)* | CSV |

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
| `cobertura_normalizada.png` | Normalized coverage plot: valid pairs per pathway per NANDA-I diagnosis (directional) | PNG (300 dpi) |
| `especificidade_direcional.png` | Directional specificity index per NANDA-I diagnosis (divergence bar chart) | PNG (300 dpi) |

---

## Methods Summary

**Stage 1 — Dataset, preprocessing, and functional enrichment.** GSE16561 (NCBI GEO) contains peripheral blood gene expression profiles measured on the Illumina HumanHT-12 platform from ischemic stroke patients and healthy controls. Raw expression values are checked for scale and log2-transformed when necessary, then quantile-normalized across samples using limma's `normalizeBetweenArrays`. Probe annotations are drawn from the platform's feature data; when multiple probes map to the same gene symbol, the probe with the largest absolute log2 fold change is retained. Samples are classified using the phenotype description field of the GEO metadata, and unclassifiable samples are excluded before modeling. A linear model is fitted per probe using `lmFit` with a contrast of AVC − Controle, and empirical Bayes moderation is applied via `eBayes` to stabilize variance estimates — a critical adjustment for microarray datasets with moderate sample sizes. Significant DEGs are defined by Benjamini-Hochberg FDR < 0.05 and |log2FC| > 1. The up- and down-regulated gene sets are submitted separately to STRINGdb v11.5 (species 9606, combined score ≥ 400) for protein-protein interaction network retrieval and functional enrichment analysis. Annotation categories include Gene Ontology Biological Process and Molecular Function, Reactome (RCTM), InterPro, SMART, Pfam, HPO, COMPARTMENTS, and UniProt Keywords; PMID co-citation entries are excluded to avoid conflating literature-based and annotation-based evidence.

**Stage 2 — Semantic expansion, hybrid embedding, and threshold calibration.** Each enriched biological pathway is expanded into a structured natural-language description incorporating gene names and mechanistic context (188–255 characters), providing substantially richer embedding surface than the raw pathway label. NANDA-I nursing diagnoses are encoded as full structured descriptions including defining characteristics and related factors, separated from short label forms (`nanda_labels` vs. `nanda_full`) to prevent label leakage during fine-tuning. The S-PubMedBert base model is fine-tuned using `MultipleNegativesRankingLoss` — rather than `CosineSimilarityLoss` — because MNRL operates on relative ranking objectives and avoids the embedding compression artifact that CosineSimilarityLoss induces in models trained on small datasets. Training and validation sets are split without leakage by stratifying on coarse similarity bins (low / mid / high); the honest Spearman correlation on the held-out validation set is 0.477. Pathway and diagnosis descriptions are embedded in two independent semantic spaces: the `[CLS]` token from `medicalai/ClinicalBERT` (L2-normalized), preferred over mean pooling because ClinicalBERT's next-sentence prediction pre-training concentrates global sentence information in the `[CLS]` position; and sentence-level embeddings from the fine-tuned S-PubMedBert model. The two cosine similarity matrices are averaged with equal weight (α = 0.5) in similarity space (late fusion), preserving each model's metric geometry. The similarity threshold is calibrated by pairwise permutation (z = 18.4), yielding a cosine cutoff of 0.65 and 175 valid pairs across 15 distinct NANDA-I diagnoses. A Shannon entropy audit is performed per pathway on the resulting similarity score distribution; the ratio H/H_max = 0.960 across all 35 pathways indicates near-uniform distribution of similarity mass over NANDA-I diagnoses, motivating Stage 3 cross-validation.

**Stage 3 — LLM-as-judge mechanistic cross-validation (pending implementation).** Each pathway–diagnosis pair surviving the Stage 2 threshold is independently evaluated by Claude Sonnet 4.6 (temperature = 0, seed = 42) operating as a mechanistic judge. The prompt is structured with anti-anchoring design: the model generates a mechanistic justification before assigning a score, preventing score priming from influencing the reasoning. Each response returns four structured fields: `justificativa` (free-text mechanistic rationale linking the biological pathway to the clinical nursing diagnosis), `score_concordancia` (numeric, 0–10), `nivel_confianca` (categorical), and `flag_lexical` (boolean: detects pairs where similarity is driven primarily by shared terminology rather than mechanistic correspondence). Pairs are stratified into three tiers — high concordance, moderate concordance, and divergent — generating a filtered set that constitutes the primary evidence base for a structured Delphi expert validation study.

---

## Limitations

- **Conservative similarity threshold.** The cosine cutoff of 0.65 was calibrated by pairwise permutation (z = 18.4) but has not been validated against a held-out clinical gold standard. Pairs with scores between 0.60 and 0.65 may represent clinically valid mappings excluded by the current threshold.

- **Synthetic training pairs.** The fine-tuning supervision signal consists entirely of pairs generated by the research team, not by certified nursing clinicians. Label values represent expert biological judgment about pathway–nursing diagnosis relevance but have not been validated through structured clinical review or inter-annotator agreement protocols.

- **No external clinical validation.** The pathway–NANDA-I mappings have not been evaluated against real patient care plans, clinical chart data, or nursing expert consensus. The output should be treated as a hypothesis-generating tool rather than a clinical decision support system.

- **Single cohort, single platform.** The entire transcriptomic evidence base is derived from one GEO dataset (GSE16561) measured on one microarray platform (Illumina HumanHT-12). Findings may not generalize across stroke subtypes, disease stages, measurement platforms, or populations with different demographic profiles.

- **Entropy collapse in Stage 2.** Shannon entropy of the similarity score distribution over NANDA-I diagnoses reaches H/H_max = 0.960 across all 35 pathways, indicating that the embedding pipeline does not discriminate biologically distinct pathways with sufficient resolution for unqualified clinical claims. This limitation motivates the Stage 3 LLM-as-judge cross-validation, which provides an independent mechanistic signal for each pathway–diagnosis pair.

- **Equal fusion weights.** The late-fusion weight α = 0.5 assigns equal contribution to ClinicalBERT and S-PubMedBert without empirical optimisation for this domain. Optimal weighting for this specific biological and clinical domain has not been determined; different α values may yield substantially different rankings for borderline pairs.

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

---

> **Current pipeline version:** v4-expansao-semantica-MNRL
> Stage 3 (LLM-as-judge) is under active development and will be implemented in `llm_judge_nanda.ipynb`. The current output `nanda_mapping_completo_threshold65.csv` (175 pairs, 15 distinct NANDA-I diagnoses, threshold z = 18.4) represents the Stage 2 baseline pending Stage 3 cross-validation.
