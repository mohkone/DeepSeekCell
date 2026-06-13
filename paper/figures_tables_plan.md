# Figures and Tables Plan

## Figure 1. DeepSeekCell Workflow

Suggested caption:

> Overview of the DeepSeekCell workflow. Cluster marker genes, tissue, and species metadata are submitted through the R API or Shiny interface. A structured prompt is sent to the selected model backend. Returned annotations are parsed into a standardized schema, mapped to the Cell Ontology, validated for confidence and completeness, and exported as tabular or HTML reports with marker-based reasoning.

Panel suggestions:

- A: User input: marker genes, tissue, species, model.
- B: Prompt construction and LLM annotation.
- C: Structured parser and validation.
- D: Cell Ontology mapping and match provenance.
- E: Shiny results, explainability, and export.

## Figure 2. Interactive Shiny Platform

Suggested caption:

> DeepSeekCell Shiny interface for interactive cell type annotation. The application accepts cluster marker genes, displays predicted cell types with confidence scores and Cell Ontology identifiers, and provides explainability, validation, cost, metadata, and report export views.

Panel suggestions:

- A: Input panel and annotation setup.
- B: Results table with ontology IDs.
- C: Explainability tab showing reasoning and marker evidence.
- D: Confidence and cost/performance tabs.

## Figure 3. Macro-F1 Benchmark

Source:

`results/benchmark_macroF1.pdf`

Suggested caption:

> Macro-F1 performance of DeepSeekCell, ScType, and SingleR across six benchmark datasets. Bars represent the mean across three benchmark replicates. Error bars represent 95% confidence intervals when replicate variability is present.

## Figure 4. Exact Accuracy Benchmark

Source:

`results/benchmark_accuracy.pdf`

Suggested caption:

> Exact cluster-level annotation accuracy across benchmark datasets. DeepSeekCell achieved the highest average exact accuracy across datasets, with strongest performance in PBMC and TasicBrain.

## Figure 5. Ontology-Aware Clade Accuracy

Source:

`results/benchmark_clade_accuracy.pdf`

Suggested caption:

> Cell Ontology-aware clade accuracy across benchmark datasets. This metric credits predictions that are ontologically related to the reference label, capturing biologically meaningful agreement beyond exact string matching.

## Figure 6. Runtime

Source:

`results/benchmark_runtime.pdf`

Suggested caption:

> Runtime comparison across evaluated methods. DeepSeekCell produced rapid marker-guided annotations, while SingleR showed longer runtimes in larger reference-based comparisons.

## Table 1. Dataset Characteristics

Source:

`results/dataset_summary.csv`

Suggested columns:

- Dataset.
- Tissue.
- Species.
- Cells.
- Genes.
- Clusters.
- Mean cluster purity.
- Minimum cluster purity.

## Table 2. Main Benchmark Results

Source:

`results/final_benchmark_table.csv`

Suggested columns:

- Dataset.
- Tissue.
- Method.
- Macro-F1 mean.
- Macro-F1 95% CI.
- Accuracy mean.
- Clade accuracy mean.
- Unknown rate.
- Evaluated clusters.
- Runtime.

## Table 3. Software Features

Suggested rows:

| Feature                      | DeepSeekCell support                 |
|------------------------------|--------------------------------------|
| Marker-guided annotation     | Yes                                  |
| Interactive Shiny interface  | Yes                                  |
| Cell Ontology ID mapping     | Yes                                  |
| Ontology match provenance    | Exact, synonym, context exact, fuzzy |
| Explainable marker reasoning | Yes                                  |
| Confidence score             | Yes                                  |
| Tissue consistency flag      | Yes                                  |
| Mixed-cluster flag           | Yes                                  |
| CSV/XLSX/HTML export         | Yes                                  |
| Benchmark scripts            | Yes                                  |
| Offline fallback ontology    | Yes                                  |

## Supplementary Table S1. Full Benchmark Results

Source:

`results/benchmark_results_full.csv`

Suggested use:

Report all replicate-level results, including ARI, macro-F1, exact accuracy, balanced accuracy, clade accuracy, unknown rate, runtime, token usage, estimated cost, and cluster purity.

## Supplementary Table S2. Cluster-Level Benchmark Labels

Source:

`results/cluster_summary.csv`

Suggested use:

Report per-cluster truth labels, cluster purity, marker counts, and dataset metadata for reproducibility.
