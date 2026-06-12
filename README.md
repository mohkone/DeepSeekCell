# DeepSeekCell

**An Ontology-Guided Large Language Model Framework and Interactive Shiny Platform for Explainable Cell Type Annotation in Single-Cell RNA Sequencing**

DeepSeekCell annotates single-cell RNA-seq clusters from marker genes using an LLM-assisted, ontology-aware workflow. It returns standardized cell type labels, confidence scores, marker-based reasoning, Cell Ontology mappings, validation summaries, and exportable reports through both R functions and an interactive Shiny interface.

## Core Features

- LLM annotation from per-cluster marker genes using DeepSeek or local Ollama.
- Cell Ontology mapping with exact, synonym, and conservative fuzzy matching provenance.
- Explainability fields for marker-based reasoning, tissue consistency, possible doublets, and possible contamination.
- Offline validation metrics, quality score, ontology coverage, and missing-cluster checks.
- Shiny application for interactive annotation and report download.
- Benchmark scripts for replicated comparisons against SingleR and scType.

## Installation

```r
install.packages(c(
  "shiny", "shinythemes", "shinycssloaders", "DT", "plotly",
  "ggplot2", "dplyr", "httr2", "jsonlite", "openxlsx",
  "ontologyIndex", "cachem", "stringdist", "logger"
))
```

Optional benchmarking dependencies:

```r
install.packages(c("Seurat", "mclust", "testthat", "yaml"))
BiocManager::install(c("SingleR", "celldex", "SingleCellExperiment", "scRNAseq"))
```

## API Configuration

Set keys as environment variables instead of storing them in scripts:

```r
Sys.setenv(DEEPSEEK_API_KEY = "...")
```

Supported endpoint overrides:

- `DEEPSEEK_API_URL`, `DEEPSEEK_MODEL_ID`
- `OLLAMA_API_URL`, `OLLAMA_MODEL_ID`

Ollama can be used without an API key when a local server is running.

## R Usage

```r
source("R/utils.R")
invisible(lapply(setdiff(list.files("R", "\\.R$", full.names = TRUE), "R/utils.R"), source))

markers <- list(
  Cluster1 = c("CD3D", "CD3E", "CD8A", "NKG7"),
  Cluster2 = c("MS4A1", "CD79A", "CD74"),
  Cluster3 = c("LYZ", "S100A8", "S100A9", "FCN1")
)

result <- annotate_cell_types(
  markers = markers,
  tissue = "PBMC",
  species = "Human",
  model_name = "deepseek",
  use_ontology = TRUE,
  validate = TRUE
)

result$annotations
generate_html_report(result, "annotation_report.html")
```

For local development without network calls:

```r
source_files <- list.files("R", "\\.R$", full.names = TRUE)
source_files <- c("R/utils.R", setdiff(source_files, "R/utils.R"))
invisible(lapply(source_files, source))
testthat::test_dir("tests/testthat")
```

## Shiny App

```r
shiny::runApp("inst/shiny")
```

The app accepts marker genes for up to five clusters, calls the selected model, maps annotations to the Cell Ontology, displays confidence and validation summaries, and exports CSV, XLSX, or HTML reports.

## Benchmarking

Benchmark scripts live in `benchmarks/`. They run a closed-label,
marker-guided cluster annotation benchmark comparing DeepSeekCell with
SingleR and scType on curated PBMC, pancreas, brain, and lung datasets.
Set `DEEPSEEK_API_KEY` to include DeepSeekCell; otherwise only non-LLM
baselines that do not require an API key will run.

```r
source("benchmarks/run_benchmark.R")
main(n_replicates = 3)
```

Generated benchmark outputs are ignored by default through `.gitignore` and `.Rbuildignore`.
Key outputs include `results/benchmark_results_summary.csv`,
`results/benchmark_results_full.csv`, `results/final_benchmark_table.csv`, `results/dataset_summary.csv`,
`results/cluster_summary.csv`, and `results/benchmark_manifest.txt`.

## Reproducibility Notes

- Do not commit API keys, `.Renviron`, `.RData`, benchmark caches, or generated figures.
- Prefer the full Cell Ontology OBO file at `data/cl.obo`; a small fallback ontology is provided only for offline smoke tests and graceful failure.
- Result metadata includes model name, model ID, token usage, runtime, ontology fallback status, and schema version.

## Citation

If you use this software, cite:

> An Ontology-Guided Large Language Model Framework and Interactive Shiny Platform for Explainable Cell Type Annotation in Single-Cell RNA Sequencing.
