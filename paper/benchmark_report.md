# DeepSeekCell Benchmark Report

## Purpose

This report summarizes the final benchmark used to support the DeepSeekCell manuscript. The benchmark strengthens the Shiny application paper by showing that the annotation framework is not only interactive and explainable, but also quantitatively competitive with established cell type annotation methods.

## Benchmark Configuration

| Field | Value |
|---|---|
| Benchmark mode | closed-label-marker-guided |
| Replicates | 3 |
| Seeds | 100, 200, 300 |
| Methods | DeepSeekCell, ScType, SingleR |
| Top markers per cluster | 25 |
| Seurat resolution | 0.5 |
| PCA dimensions | 50 |
| Cell Ontology terms loaded | 3,437 |
| Cache version | 2026-06-12-publication-v5 |

## Datasets

| Dataset | Tissue | Species | Cells | Genes | Clusters | Mean cluster purity |
|---|---|---|---:|---:|---:|---:|
| PBMC | PBMC | Human | 2,700 | 13,714 | 8 | 0.941 |
| BaronPancreas | Pancreas | Human | 8,569 | 20,125 | 16 | 0.952 |
| MuraroPancreas | Pancreas | Human | 2,122 | 19,059 | 11 | 0.954 |
| TasicBrain | Brain | Mouse | 1,809 | 24,058 | 15 | 0.929 |
| ZeiselBrain | Brain | Mouse | 3,005 | 20,006 | 16 | 0.773 |
| ZilionisLung | Lung | Human | 5,000 per replicate | 41,861 | 16-18 | 0.910 |

ZeiselBrain has the lowest mean cluster purity and should be described as a more difficult cluster-level benchmark. ZilionisLung is the only dataset with seed-dependent subsampling, so its confidence intervals are nonzero.

## Overall Results

| Method | Mean macro-F1 | Mean accuracy | Mean clade accuracy | Mean runtime (s) |
|---|---:|---:|---:|---:|
| DeepSeekCell | 0.818 | 0.810 | 0.868 | 1.86 |
| ScType | 0.749 | 0.773 | 0.778 | 1.38 |
| SingleR | 0.466 | 0.479 | 0.757 | 51.37 |

DeepSeekCell achieved the best overall mean macro-F1 and exact accuracy. ScType remained very competitive in pancreas datasets. SingleR performed best in ZilionisLung and showed strong ontology-aware clade accuracy in several datasets despite lower exact label matching in pancreas.

## Dataset-Level Results

| Dataset | Best macro-F1 method | Best macro-F1 | DeepSeekCell macro-F1 | Main interpretation |
|---|---|---:|---:|---|
| PBMC | DeepSeekCell | 1.000 | 1.000 | Canonical immune markers were annotated perfectly. |
| BaronPancreas | ScType | 0.947 | 0.926 | DeepSeekCell was close to ScType; pancreas marker database coverage favors ScType. |
| MuraroPancreas | ScType | 0.958 | 0.917 | Similar pattern to BaronPancreas; DeepSeekCell remains competitive. |
| TasicBrain | DeepSeekCell | 0.984 | 0.984 | Strong brain annotation with perfect clade accuracy. |
| ZeiselBrain | ScType | 0.580 | 0.465 | All methods were affected by lower cluster purity and label granularity. |
| ZilionisLung | SingleR | 0.659 | 0.619 | SingleR led by exact macro-F1; DeepSeekCell had comparable clade accuracy. |

## Publication-Ready Interpretation

The benchmark supports three main claims:

1. DeepSeekCell is quantitatively competitive with established tools across multiple tissues.
2. Ontology-aware evaluation provides additional biological interpretability beyond exact string matching.
3. The Shiny application is supported by a reproducible benchmark rather than being only a user interface demonstration.

Suggested wording:

> Across six curated scRNA-seq datasets spanning immune, pancreatic, neural, and lung tissues, DeepSeekCell achieved the highest average macro-F1 and exact accuracy among the evaluated methods. The framework performed particularly well on PBMC and Tasic brain datasets and remained competitive in pancreas and lung settings. Ontology-aware clade accuracy revealed biologically related predictions even when exact label matches differed, supporting the value of Cell Ontology integration for transparent annotation evaluation.

## Caveats to Report

- Most datasets are deterministic across replicates because the cached preprocessing and clustering are fixed. This explains zero confidence intervals for PBMC, pancreas, TasicBrain, and ZeiselBrain.
- ZilionisLung uses seed-dependent subsampling and therefore has nonzero confidence intervals.
- The benchmark is cluster-level and marker-guided; it does not replace cell-level benchmarking on raw count matrices.
- The ZeiselBrain dataset has lower mean cluster purity and should be interpreted cautiously.
- LLM reasoning is useful for review but should not be treated as independent experimental validation.

## Recommended Tables and Figures

| Item | Source | Purpose |
|---|---|---|
| Table 1 | `results/dataset_summary.csv` | Dataset characteristics and cluster purity. |
| Table 2 | `results/final_benchmark_table.csv` | Main benchmark performance. |
| Figure 1 | Workflow schematic | DeepSeekCell pipeline from markers to ontology-linked report. |
| Figure 2 | Shiny screenshots | Interactive annotation and explainability tabs. |
| Figure 3 | `results/benchmark_macroF1.pdf` | Macro-F1 comparison across datasets. |
| Figure 4 | `results/benchmark_accuracy.pdf` | Exact accuracy comparison. |
| Figure 5 | `results/benchmark_clade_accuracy.pdf` | Ontology-aware clade accuracy. |
| Figure 6 | `results/benchmark_runtime.pdf` | Runtime comparison. |

## Files Used

- `results/final_benchmark_table.csv`
- `results/benchmark_results_full.csv`
- `results/benchmark_results_summary.csv`
- `results/dataset_summary.csv`
- `results/cluster_summary.csv`
- `results/benchmark_manifest.txt`
