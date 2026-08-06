# DeepSeekCell

**An evidence-guided reliability framework for LLM-based single-cell annotation under limited refinement budgets**

DeepSeekCell formulates marker-based LLM cell type annotation as a reliability
and resource-allocation problem. A first-pass LLM generates structured
candidate labels from cluster marker genes. Deterministic marker, ontology, and
tissue evidence are then used to score reliability, detect biological conflict,
adjust confidence, and selectively allocate a limited second-pass refinement
budget. The learned Risk-k extension trains a reliability model to estimate
first-pass error risk or expected refinement benefit and ranks clusters under
the same fixed-budget objective.

The software still includes DeepSeek and local Ollama backends plus a Shiny
interface, but the main contribution is the model-agnostic formulation of
selective LLM refinement as a reliability and resource-allocation problem
rather than any single LLM or interface.

Current software version: `0.1.0`.

## Central Computational Problem

LLM-based annotation workflows often produce one label and one self-reported
confidence score. When predictions are uncertain or biologically inconsistent,
users must decide whether to trust the first pass, refine everything, or inspect
clusters manually. DeepSeekCell addresses the fixed-budget question:

> Which first-pass LLM annotations deserve additional reasoning?

## Algorithmic Components

- Candidate generation from per-cluster marker genes using DeepSeek or local Ollama.
- Deterministic marker-profile evidence scoring.
- Cell Ontology validation with exact, synonym, context-aware, and conservative fuzzy matching provenance.
- Tissue-consistency and marker-prediction consensus scoring.
- Evidence-adjusted confidence while preserving the original `LLMConfidence`.
- Frozen fixed-budget selectors: Evidence-k, Confidence-k, Random-k, NoOntology-k, Full, and None.
- Learned Risk-k selector using logistic error-risk or refinement-benefit models.
- Refinement provenance columns for first-pass labels, flagging, refinement, label changes, and reasons.
- Benchmark scripts for paired cached ablations, calibration metrics, selector efficiency, cost, and comparisons against SingleR, scType, scmap, and native CellTypist.

## Reliability Pipeline

```text
Marker genes
  -> candidate generation
  -> marker / ontology / tissue evidence scoring
  -> conflict detection
  -> evidence-adjusted confidence
  -> fixed-budget refinement selection (Evidence-k or learned Risk-k)
  -> optional second-pass reasoning
  -> final annotation with provenance
```

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
BiocManager::install(c("SingleR", "celldex", "SingleCellExperiment", "scRNAseq", "scmap"))
install.packages("reticulate")
```

CellTypist is a Python dependency. Install it into the Python environment used by
`reticulate`:

```bash
pip install celltypist anndata numpy
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
  validate = TRUE,
  calibrate_confidence = TRUE,
  self_refine = TRUE,
  refinement_strategy = "evidence",
  refinement_budget = NULL
)

result$annotations
generate_html_report(result, "annotation_report.html")
```

Use `select_refinement_candidates()` directly when you want to audit or compare
selectors without making a second LLM call:

```r
scored <- calibrate_annotation_confidence(
  annotations = result$annotations,
  markers = markers,
  tissue = "PBMC"
)

evidence_k <- select_refinement_candidates(
  scored,
  strategy = "evidence",
  budget = 3
)

confidence_k <- select_refinement_candidates(
  scored,
  strategy = "confidence",
  budget = 3
)
```

Train a learned Risk-k selector from paired benchmark debug files:

```bash
Rscript benchmarks/train_reliability_model.R results/benchmark_debug results/reliability_model_v1.1_error.rds first_pass_error
```

Training writes explainability outputs beside the model, including feature
importance, odds ratios, linear SHAP-style cluster contributions, calibration
bins, calibration metrics, feature ablations, and PDF plots. Recompute those
outputs for an existing model with:

```bash
Rscript benchmarks/explain_reliability_model.R results/reliability_model_v1.1_error.rds results/reliability_model_v1.1_error_training_features.csv
```

Then use it for risk-aware selection:

```r
model <- readRDS("results/reliability_model_v1.1_error.rds")

risk_k <- select_refinement_candidates(
  scored,
  strategy = "risk",
  budget = 3,
  reliability_model = model
)
```

For benchmark runs, set `DEEPSEEKCELL_RELIABILITY_MODEL` to the saved model RDS
to add the compute-matched `DeepSeekCell-RiskK` arm.

Assess whether the learned reliability model transfers across biological
domains with the cross-tissue generalization workflow:

```bash
Rscript benchmarks/analyse_reliability_generalization.R results/benchmark_debug results/reliability_generalization first_pass_error 10
```

This analysis does not introduce another selector. It reuses the existing
Risk-k model family to generate leave-one-tissue-out transfer results,
single-tissue transfer matrices, learning curves, domain-shift summaries,
oracle-gap analyses, calibration comparisons, and a failure taxonomy. The
protocol is documented in `benchmarks/generalization_protocol.md`.

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

The app accepts marker genes for up to five clusters, calls the selected model,
maps annotations to the Cell Ontology, displays confidence and validation
summaries, and exports CSV, XLSX, or HTML reports. The app is an implementation
layer for the reliability framework; the algorithm can also be used entirely
from R scripts.

## Benchmarking

Benchmark scripts live in `benchmarks/`. They run a closed-label,
marker-guided cluster annotation benchmark comparing DeepSeekCell with
SingleR, scType, scmap, and CellTypist on curated PBMC, pancreas, brain, and
lung datasets. DeepSeekCell is intended for cluster marker lists rather than raw
expression matrices, so the benchmark evaluates marker-driven annotation.
Optional baselines are included only when their dependencies are available;
otherwise they are skipped with an explanatory message.

CellTypist uses tissue-aware pretrained models by default: `Immune_All_Low.pkl`
for PBMC, `Adult_Human_PancreaticIslet.pkl` for pancreas, and
`Human_Lung_Atlas.pkl` for lung. Set `CELLTYPIST_MODEL` to override this
automatic choice. CellTypist is evaluated in its native cell-level workflow
where compatible human models are available; cluster-level labels are derived
by majority vote only as a secondary harmonized comparison.
Set `DEEPSEEK_API_KEY` to include DeepSeekCell; otherwise only non-LLM
baselines that do not require an API key will run.

The benchmark also writes paired ablation, calibration, selector-efficiency,
stability, runtime, cost, token, and statistical summaries including
`benchmark_pairwise_wilcoxon.csv`, `benchmark_friedman_tests.csv`,
`ablation_confidence_quality.csv`, `ablation_refinement_behavior.csv`,
`refinement_efficiency_summary.csv`, and `benchmark_llm_stability.csv`.

### Locked external validation

The reliability layer is frozen as
`DeepSeekCell reliability specification v1.0`. The machine-readable lock file
and marker audit tables are distributed with the package:

- `inst/extdata/reliability_spec_v1.0.json`
- `inst/extdata/marker_profiles_v1.0.csv`
- `inst/extdata/marker_aliases_v1.0.csv`

The frozen specification records the evidence-score weights, marker profiles,
aliases, conflict thresholds, selector rules, Cell Ontology MD5 hash, prompt
version, software version, model defaults, and primary external-validation
endpoint. It can also be inspected from R:

```r
get_reliability_spec(include_marker_profiles = FALSE)
```

The learned Risk-k extension is specified in
`inst/extdata/reliability_model_spec_v1.1.json` and described in
`benchmarks/risk_model_protocol.md`. Train it on development benchmark outputs,
freeze the saved RDS, and set `DEEPSEEKCELL_RELIABILITY_MODEL` before held-out
validation.

Use `benchmarks/external_validation_protocol.md` as the locked held-out
validation plan. Copy `benchmarks/external_validation_manifest_template.csv` to
`benchmarks/external_validation_manifest.csv`, fill in genuinely held-out
datasets, and run a preflight check:

```bash
Rscript benchmarks/run_external_validation.R benchmarks/external_validation_manifest.csv 1
```

The preflight writes `results/external_validation_lock.json` and
`results/external_validation_plan.csv` without calling an LLM. To execute the
locked validation after the prepared RDS paths are available:

```bash
set DEEPSEEKCELL_RUN_EXTERNAL_VALIDATION=true
set DEEPSEEKCELL_RELIABILITY_MODEL=results/reliability_model_v1.1_error.rds
Rscript benchmarks/run_external_validation.R benchmarks/external_validation_manifest.csv 3 deepseek
```

Each prepared external RDS should contain a named marker list `markers`, a named
truth vector `truth`, and optional `tissue`, `species`, and `purity` fields. The
runner preserves the paired design: one hashed first-pass response per
model-dataset-replicate block is reused by all selectors.

Secondary robustness checks can be run with:

```bash
Rscript benchmarks/sensitivity_analysis.R cluster_level_input.csv results/sensitivity
```

When no input is supplied, the script writes
`results/sensitivity_input_schema.csv` describing the expected cluster-level
columns.

### Fresh-clone benchmark workflow

From a fresh clone, install the R and optional Python dependencies, then obtain
the external resources used by the benchmark:

- Cell Ontology OBO file at `data/cl.obo`.
- ScType database at `scType/ScTypeDB_full.xlsx`.
- Optional CellTypist Python environment configured through `RETICULATE_PYTHON`.

Run the reproducibility smoke tests first:

```r
devtools::test()
```

For a full benchmark with new LLM calls:

```r
Sys.setenv(
  DEEPSEEK_API_KEY = "...",
  DEEPSEEKCELL_USE_LLM_CACHE = "true"
)
source("benchmarks/run_benchmark.R")
main(n_replicates = 3)
```

`DEEPSEEKCELL_USE_LLM_CACHE=true` means that existing cached first-pass and
refinement responses in `results/benchmark_debug/` are reused, and missing cache
entries are generated and saved. To reproduce a benchmark without additional
API calls, restore `results/benchmark_debug/` from a previous run or submission
package before executing the command above. The paired ablation design requires
that all DeepSeekCell arms reuse the same cached first-pass response hash within
each dataset and replicate.

For local LLM pilot checks with Ollama, start the Ollama server and run:

```bash
Rscript benchmarks/run_ollama_multimodel_pilot.R ZilionisLung 1
Rscript benchmarks/summarise_ollama_multimodel_pilots.R
```

## Archiving

The GitHub-Zenodo software archive is available at
https://doi.org/10.5281/zenodo.20680434. The all-version Zenodo concept DOI is
https://doi.org/10.5281/zenodo.20680433. The repository includes `.zenodo.json`
metadata so future Zenodo records are populated consistently.

If `reticulate::py_config()` reports a missing Python from `RETICULATE_PYTHON`,
restart R and point reticulate to a valid environment before sourcing the
benchmark:

```r
Sys.setenv(
  RETICULATE_PYTHON =
    "C:/Users/Mohamed KONE/Documents/.virtualenvs/deepseek_env/Scripts/python.exe"
)
reticulate::py_config()
reticulate::py_module_available("celltypist")
```

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
- Result metadata includes model name, model ID, token usage, runtime, ontology fallback status, refinement strategy, refinement budget, evidence flags, and schema version.
- The local `paper/` folder is ignored and should not be committed to GitHub.

## Citation

If you use this software, cite:

> Evidence-guided selective refinement improves the reliability of LLM-based single-cell annotation.
