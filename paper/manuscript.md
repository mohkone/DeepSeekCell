# An Ontology-Guided Large Language Model Framework and Interactive Shiny Platform for Explainable Cell Type Annotation in Single-Cell RNA Sequencing

## Manuscript Draft

### Authors

Mohamed Kone and Shulin Wang

Correspondence: Shulin Wang, Hunan University.

## Abstract

Cell type annotation remains a critical interpretive step in single-cell RNA sequencing analysis, yet manual annotation is time-consuming, tissue-dependent, and difficult to reproduce across studies. Recent large language models can infer plausible cell type labels from marker genes, but unconstrained free-text outputs can limit reproducibility, ontology interoperability, and downstream validation. We present DeepSeekCell, an ontology-guided large language model framework and interactive Shiny platform for explainable cell type annotation from cluster marker genes. DeepSeekCell converts marker gene sets into standardized cell type annotations, confidence scores, marker-based reasoning, tissue-consistency flags, possible mixed-cluster calls, and Cell Ontology mappings with explicit match provenance. The framework supports DeepSeek and local Ollama backends, validates annotation completeness and quality, and exports CSV, XLSX, and HTML reports from an interactive Shiny interface.

We evaluated DeepSeekCell in a closed-label, marker-guided benchmark against SingleR, ScType, scmap, and an exploratory CellTypist baseline across six curated single-cell datasets spanning human PBMC, human pancreas, mouse brain, and human lung. Across datasets, DeepSeekCell achieved the highest mean macro-F1 score (0.818) and mean exact accuracy (0.810), with mean ontology-aware clade accuracy of 0.868. Performance was strongest in PBMC and Tasic brain datasets, where DeepSeekCell reached macro-F1 scores of 1.000 and 0.984, respectively. In pancreas datasets, DeepSeekCell performed competitively with ScType, while in the Zilionis lung dataset SingleR achieved the highest macro-F1 score and DeepSeekCell showed comparable ontology-aware clade accuracy. These results suggest that ontology-guided LLM annotation can provide accurate, explainable, and interoperable cell type labels while preserving a practical user workflow through a web-based Shiny application.

Keywords: single-cell RNA sequencing; cell type annotation; Cell Ontology; large language model; Shiny; explainable AI; benchmark

## Introduction

Single-cell RNA sequencing (scRNA-seq) enables high-resolution characterization of cellular heterogeneity across tissues, development, and disease. A central step in scRNA-seq analysis is assigning biological identities to transcriptionally defined clusters or individual cells. In practice, cell type annotation often relies on expert interpretation of marker genes, comparison with reference atlases, and iterative inspection of tissue context. This process is powerful but can be subjective, labor-intensive, and difficult to reproduce across laboratories.

Automated annotation methods have improved scalability. Reference-based approaches such as SingleR compare query expression profiles with labeled reference datasets, while marker-based approaches such as ScType use positive and negative marker combinations to identify cell identities. These methods are widely useful, but each depends on the quality, coverage, and specificity of reference data or marker databases. They can also return labels that differ in granularity, terminology, or ontology compatibility, making downstream comparison across datasets difficult.

Large language models (LLMs) offer a complementary strategy because they can reason over marker gene lists and tissue context using broad biomedical knowledge. Recent work includes GPTCelltype for marker-based cell annotation, scExtract and CASSIA for automated or multi-agent single-cell annotation, AnnDictionary and DeepCellSeek for benchmarking LLM behavior in cell typing, and consensus approaches that combine multiple LLM outputs. In parallel, single-cell foundation and language-inspired models such as scBERT, Geneformer, scGPT, and scFoundation demonstrate the value of large-scale transcriptomic pretraining, while biomedical LLMs and tool-augmented systems such as BioGPT, Med-PaLM, GeneGPT, and GeneAgent show how language models can be adapted or connected to external biomedical resources. However, free-text LLM annotation introduces new challenges: outputs may be inconsistent across runs, difficult to validate automatically, weakly linked to controlled vocabularies, and insufficiently transparent for publication-grade analysis. For scRNA-seq annotation, an LLM-based system should therefore constrain outputs, expose reasoning, map labels to formal ontologies, and provide reproducible exports.

DeepSeekCell was designed to address these needs. The framework combines marker-guided LLM annotation with Cell Ontology normalization, validation, and an interactive Shiny interface. It is intended to support researchers who need rapid annotation, transparent reasoning, standardized labels, and reproducible reports without manually editing LLM outputs. This manuscript describes the DeepSeekCell workflow, software implementation, ontology mapping strategy, Shiny platform, and benchmark performance across immune, pancreatic, brain, and lung datasets.

## Materials and Methods

### Software overview

DeepSeekCell is implemented as an R package with an interactive Shiny application. The package provides functions for marker preprocessing, prompt construction, LLM API access, response parsing, Cell Ontology mapping, validation, visualization, export, and benchmarking. The Shiny application is located in `inst/shiny/app.R` and provides a browser-based workflow for entering marker genes, selecting model and metadata options, running annotation, inspecting results, viewing confidence summaries, reviewing explainability fields, and exporting reports.

The current software version is 0.1.0. The maintained model backends are DeepSeek and local Ollama. OpenAI/GPT-4o support was removed from the maintained application surface to keep the publication version focused, reproducible, and aligned with the DeepSeekCell naming and benchmark design.

### Input data and preprocessing

DeepSeekCell accepts per-cluster marker genes as comma-, semicolon-, whitespace-, or newline-separated text. Marker processing removes common low-information genes, including mitochondrial, ribosomal, and pseudogene-like features, and standardizes cluster marker lists before prompt construction. The Shiny interface currently accepts marker genes for up to five clusters, while the R interface can be used programmatically for larger analyses.

### LLM annotation workflow

For each dataset or user submission, DeepSeekCell constructs a structured prompt containing tissue, species, cluster marker genes, output schema requirements, and ontology-aware annotation instructions. The model is instructed to return structured annotation records containing cluster ID, predicted cell type, confidence, tissue consistency, mixed-cluster status, and marker-based reasoning. Response parsing is robust to extra prose and code fences, and confidence values are normalized to the 0-1 range.

### Explainability and validation

DeepSeekCell exposes annotation reasoning in both the R output and the Shiny interface. For each cluster, the result can include marker genes used by the model, a natural-language reasoning statement, tissue-consistency classification, mixed-cluster status, confidence score, and ontology mapping metadata. The validation module reports missing clusters, unknown annotations, low-confidence predictions, ontology coverage, and overall quality summaries.

### Cell Ontology mapping

DeepSeekCell maps predicted cell type labels to Cell Ontology identifiers using a staged strategy:

1.  Exact label matching.
2.  Synonym matching.
3.  Tissue-context disambiguation for common ambiguous labels.
4.  Conservative fuzzy matching when exact or synonym matches fail.
5.  Fallback ontology mappings for offline tests and graceful failure.

The ontology loader uses the full Cell Ontology OBO file when available and caches parsed ontology content safely. The benchmark run loaded 3,437 Cell Ontology terms from `data/cl.obo`. The mapping layer records the Cell Ontology ID, ontology label, match method, and match score, which helps distinguish exact mappings from inferred or lower-confidence mappings. Context-specific mappings were added for pancreas, brain, and lung to reduce common ambiguity, such as pancreatic endocrine labels, neural/glial labels, and lung epithelial or immune labels.

### Shiny platform

The Shiny application provides a practical interface for non-programmatic annotation. Users enter an API key or use a local model, choose tissue and species, enter marker genes, optionally map annotations to the Cell Ontology, and run annotation. Results are displayed in tabbed views:

- Results: predicted cell type, confidence, tissue consistency, mixed status, Cell Ontology ID, and ontology label.
- Explainability: marker genes, reasoning, ontology match method, and ontology match score.
- Confidence: cluster-level confidence visualization.
- Cost and performance: runtime, token usage, and estimated cost.
- Metadata: model, tissue, species, ontology status, schema version, and timestamp.
- Help: usage notes and citation text.

Results can be exported as CSV, Excel, or HTML reports.

### Benchmark design

We evaluated DeepSeekCell using a closed-label, marker-guided benchmark. Six curated datasets were included:

| Dataset        |   Tissue | Species |               Cells |  Genes | Clusters |
|----------------|---------:|--------:|--------------------:|-------:|---------:|
| PBMC           |     PBMC |   Human |               2,700 | 13,714 |        8 |
| BaronPancreas  | Pancreas |   Human |               8,569 | 20,125 |       16 |
| MuraroPancreas | Pancreas |   Human |               2,122 | 19,059 |       11 |
| TasicBrain     |    Brain |   Mouse |               1,809 | 24,058 |       15 |
| ZeiselBrain    |    Brain |   Mouse |               3,005 | 20,006 |       16 |
| ZilionisLung   |     Lung |   Human | 5,000 per replicate | 41,861 |    16-18 |

For each dataset, Seurat-based preprocessing and clustering were used to produce cluster-level marker genes. The benchmark used the top 25 marker genes per cluster, Seurat resolution 0.5, and 50 principal components. Three benchmark replicates were run with seeds 100, 200, and 300. Fixed datasets were deterministic across replicates, while the Zilionis lung dataset used seed-dependent cell subsampling.

DeepSeekCell was compared with SingleR, ScType, scmap, and CellTypist on the same cluster structure. SingleR, ScType, and scmap predictions were converted to cluster-level labels by majority vote where needed. CellTypist was evaluated on cluster-average pseudo-profiles rather than its native cell-level workflow; therefore, comparisons involving CellTypist should be interpreted as exploratory. Predicted and reference labels were harmonized through explicit tissue-aware label rules before metric calculation.

### Evaluation metrics

We evaluated annotation performance using:

- Adjusted Rand Index (ARI).
- Macro-F1 score.
- Exact cluster-level accuracy.
- Balanced accuracy.
- Ontology-aware clade accuracy.
- Unknown-label rate.
- Runtime.
- Token usage and estimated cost for LLM runs.

Ontology-aware clade accuracy considers a prediction correct when the predicted and reference labels are ontologically related at an accepted lineage level, allowing biologically related predictions to be distinguished from unrelated errors.

## Results

### DeepSeekCell provides structured and ontology-linked annotations

The Shiny platform successfully loaded the Cell Ontology and displayed structured annotations containing predicted labels, confidence scores, tissue-consistency tags, mixed-cluster tags, Cell Ontology identifiers, ontology labels, match methods, and match scores. The explainability tab provided marker-based reasoning for each prediction, supporting review of why a label was assigned. Export functions generated CSV, Excel, and HTML outputs suitable for downstream reporting.

### Benchmark runs completed successfully

The final benchmark run used three replicates and completed all reported method-dataset combinations, resulting in 84 full result rows and 28 dataset-method summaries. Every reported dataset-method combination had three successful runs. All clusters were evaluated, and no benchmark cluster truth remained labeled as unknown after label curation. The benchmark manifest recorded Cell Ontology loading, methods, top marker count, Seurat settings, cache version, and dataset quality summaries.

### Overall benchmark performance

Across the six datasets, DeepSeekCell achieved the highest average macro-F1 and exact accuracy among the evaluated methods.

| Method | Mean macro-F1 | Mean accuracy | Mean clade accuracy | Mean runtime (s) |
|----|---:|---:|---:|---:|
| DeepSeekCell | 0.818 | 0.810 | 0.868 | 1.87 |
| ScType | 0.749 | 0.773 | 0.778 | 1.75 |
| SingleR | 0.466 | 0.479 | 0.757 | 43.84 |
| scmap | 0.411 | 0.438 | 0.702 | 48.58 |
| CellTypist | 0.306 | 0.353 | 0.374 | 2.12 |

DeepSeekCell was especially strong in PBMC and Tasic brain datasets. It achieved perfect macro-F1, exact accuracy, and clade accuracy on PBMC. On TasicBrain, DeepSeekCell achieved macro-F1 of 0.984, exact accuracy of 0.933, and clade accuracy of 1.000.

### Pancreas datasets

In BaronPancreas, ScType achieved the highest macro-F1 (0.947), followed closely by DeepSeekCell (0.926). In MuraroPancreas, ScType again had the highest macro-F1 (0.958), while DeepSeekCell remained competitive (0.917). SingleR showed lower exact accuracy in both pancreas datasets but high ontology-aware clade accuracy in BaronPancreas, suggesting that some predictions were ontologically related even when exact labels differed.

### Brain datasets

In TasicBrain, DeepSeekCell slightly exceeded ScType in macro-F1 and achieved perfect clade accuracy. In ZeiselBrain, all methods showed lower macro-F1 compared with other datasets, consistent with the lower mean cluster purity of this dataset (0.773). ScType and SingleR had similar macro-F1 values in ZeiselBrain, while DeepSeekCell achieved lower macro-F1 but retained strong ontology-aware clade accuracy (0.800), indicating partially correct higher-level biological labeling.

### Lung dataset

In ZilionisLung, SingleR achieved the highest macro-F1 (0.659), followed by DeepSeekCell (0.619) and ScType (0.163). DeepSeekCell and SingleR had nearly identical ontology-aware clade accuracy (0.685). Confidence intervals were nonzero in ZilionisLung because cell subsampling varied across benchmark seeds, producing 16-18 clusters across replicates.

### Runtime and cost

DeepSeekCell produced annotations rapidly in the benchmark, with mean runtime of 1.86 seconds across datasets. ScType had a mean runtime of 1.38 seconds. SingleR was slower on average (51.37 seconds), primarily due to long runtime on BaronPancreas. LLM token usage and estimated cost were recorded by the benchmark framework and can be reported from the full benchmark output if desired.

### LLM reproducibility and stability

The benchmark used a temperature of 0 for DeepSeek API calls. Across the five fixed cluster-structure datasets, DeepSeekCell produced identical harmonized labels for all 66 comparable clusters across three benchmark replicates.

| Dataset | Comparable clusters | Replicates | Stability (%) |
|---|---:|---:|---:|
| PBMC | 8 | 3 | 100.0 |
| BaronPancreas | 16 | 3 | 100.0 |
| MuraroPancreas | 11 | 3 | 100.0 |
| TasicBrain | 15 | 3 | 100.0 |
| ZeiselBrain | 16 | 3 | 100.0 |
| Overall fixed-dataset summary | 66 | 3 | 100.0 |

ZilionisLung was excluded from this label-stability summary because seed-dependent subsampling changed the cluster structure across replicates.

## Discussion

DeepSeekCell demonstrates that ontology-guided LLM annotation can produce accurate, explainable, and standardized cell type labels from marker genes. The framework addresses several limitations of unconstrained LLM annotation by requiring structured outputs, validating results, and mapping labels to the Cell Ontology. The Shiny interface makes this workflow accessible to users who need interactive annotation and exportable reports without writing custom code.

The benchmark suggests that DeepSeekCell is strongest when marker genes clearly encode canonical cell identities, such as immune populations in PBMC and neuronal/glial classes in TasicBrain. In pancreas datasets, DeepSeekCell performed competitively with ScType, which benefits from marker database coverage for pancreatic endocrine and exocrine cell types. In the ZilionisLung dataset, SingleR performed best by macro-F1, while DeepSeekCell retained comparable ontology-aware clade accuracy, suggesting that ontology-level evaluation can capture biologically related predictions missed by exact string matching.

The ZeiselBrain dataset highlights an important limitation: benchmark accuracy depends on cluster purity and label granularity. Lower mean cluster purity can make cluster-level annotation ambiguous for all methods. DeepSeekCell therefore reports tissue consistency, mixed-cluster flags, confidence, and ontology match metadata to help users identify annotations requiring manual review.

This study has several limitations. First, the benchmark is cluster-level and marker-guided rather than cell-level across all methods. Second, CellTypist was evaluated on cluster-average pseudo-profiles rather than its native cell-level workflow; therefore, comparisons involving CellTypist should be interpreted as exploratory. Third, most datasets were deterministic across replicates, resulting in zero confidence intervals for fixed datasets; only the ZilionisLung dataset included seed-dependent subsampling variability. Fourth, performance depends on prompt quality, model version, ontology coverage, and marker gene quality. Fifth, LLM-generated reasoning should be treated as an explanatory aid rather than independent biological proof.

Future work will extend the Shiny platform to support file upload for arbitrary marker tables, larger cluster counts, additional ontology-aware visualizations, and optional offline annotation modes. Benchmark extensions should include more tissues, disease states, perturbation datasets, additional LLM backends, and independent manual review of ambiguous cases.

## Conclusion

DeepSeekCell provides an ontology-guided LLM framework and Shiny platform for explainable scRNA-seq cell type annotation. By combining structured marker-based LLM annotation, Cell Ontology mapping, validation, and reproducible exports, the framework supports practical and publication-ready annotation workflows. Benchmark results across PBMC, pancreas, brain, and lung datasets show that DeepSeekCell can achieve strong annotation performance while adding interpretability and ontology interoperability.

## Data and Code Availability

The source code is available at <https://github.com/mohkone/DeepSeekCell.git>. The GitHub-Zenodo software archive is available at <https://doi.org/10.5281/zenodo.20680434>, with all-version concept DOI <https://doi.org/10.5281/zenodo.20680433>. Benchmark scripts are located in `benchmarks/`. The final benchmark outputs used for this draft are stored locally in `results/` and include `final_benchmark_table.csv`, `benchmark_results_full.csv`, `benchmark_results_summary.csv`, `dataset_summary.csv`, `cluster_summary.csv`, `benchmark_llm_stability.csv`, and `benchmark_manifest.txt`.

Large external resources, including the full Cell Ontology OBO file and ScType database, should be distributed according to their upstream licenses or regenerated/downloaded by users as documented.

## Reproducibility

The final benchmark configuration used:

- Replicates: 3.
- Seeds: 100, 200, 300.
- Top markers per cluster: 25.
- Seurat resolution: 0.5.
- PCA dimensions: 50.
- Benchmark mode: closed-label-marker-guided.
- Cache version: 2026-06-12-publication-v5.
- Cell Ontology terms loaded: 3,437.

## Ethics Statement

This study used public or package-distributed benchmark datasets and did not involve new human participant recruitment or new animal experiments. Dataset-specific licenses and original ethical approvals should be reviewed before final submission.

## Conflicts of Interest

The authors declare no competing interests.

## Funding

The authors received no specific funding for this work.

## References

1.  Aran D. et al. Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage. Nature Immunology 20, 163-172 (2019). <https://doi.org/10.1038/s41590-018-0276-y>
2.  SingleR Bioconductor documentation. <https://bioconductor.org/books/release/SingleRBook/introduction.html>
3.  Ianevski A., Giri A.K., and Aittokallio T. Fully-automated and ultra-fast cell-type identification using specific marker combinations from single-cell transcriptomic data. Nature Communications 13, 1246 (2022). <https://doi.org/10.1038/s41467-022-28803-w>
4.  The Cell Ontology. <https://cell-ontology.github.io/>
5.  Tan S.Z.K. et al. The Cell Ontology in the age of single-cell omics. <https://pmc.ncbi.nlm.nih.gov/articles/PMC12306828/>
6.  Hao Y. et al. Integrated analysis of multimodal single-cell data. Cell 184, 3573-3587.e29 (2021). <https://doi.org/10.1016/j.cell.2021.04.048>
7.  Seurat: tools for single cell genomics. <https://satijalab.org/seurat/>
8.  Hou W. and Ji Z. Assessing GPT-4 for cell type annotation in single-cell RNA-seq analysis. Nature Methods (2024). <https://doi.org/10.1038/s41592-024-02235-4>
9.  Yang F. et al. scBERT as a large-scale pretrained deep language model for cell type annotation of single-cell RNA-seq data. Nature Machine Intelligence (2022). <https://doi.org/10.1038/s42256-022-00534-z>
10. Cui H. et al. scGPT: toward building a foundation model for single-cell multi-omics using generative AI. Nature Methods (2024). <https://doi.org/10.1038/s41592-024-02201-0>
11. Jin Q. et al. GeneGPT: augmenting large language models with domain tools for improved access to biomedical information. Bioinformatics (2024). <https://doi.org/10.1093/bioinformatics/btae075>
12. Xiao T. et al. Benchmarking large language models for cell typing in single-cell RNA-Seq. Briefings in Bioinformatics (2025). <https://doi.org/10.1093/bib/bbaf677>
