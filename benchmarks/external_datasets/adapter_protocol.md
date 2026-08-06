# External Dataset Adapter Protocol

This protocol standardizes how independent studies are converted into locked
DeepSeekCell external-validation inputs. It addresses the main remaining
scientific risk: prepared RDS files must be created reproducibly and without
leakage from the benchmark results.

## Locked RDS Contract

Every prepared external dataset must be an R object with:

```r
list(
  markers = named_marker_list,
  truth = named_cluster_truth,
  tissue = "Lung",
  species = "Human",
  purity = named_cluster_purity,
  metadata = list(...)
)
```

The names of `markers`, `truth`, and `purity` must match exactly. Each marker
cluster must contain at least one marker gene.

## Registry

External studies are declared in:

```text
benchmarks/external_datasets/registry.csv
```

The registry records accession, source repository, center, laboratory,
sequencing platform, disease status, local input paths, truth-label column,
label mapping path, leakage flags, annotation source, inclusion decision, and
selection rationale.

Rows are not confirmatory until `InclusionDecision = include` and the leakage
flags pass validation.

Registry `DataPath` values may point to local files/directories or to locked
virtual sources of the form `scRNAseq::FunctionName(arguments)`. Virtual
sources are evaluated by the adapter engine and then converted through the same
common preparation workflow as local Seurat, SingleCellExperiment, h5ad, 10x,
and matrix inputs.

## Adapter Entry Points

Current study entrypoints are:

```text
prepare_wilk2020.R
prepare_segerstolpe2016.R
prepare_motor_cortex2021.R
prepare_all_external_datasets.R
```

Each study adapter calls the same common preparation engine in
`prepare_utils.R`.

The initial prepared pilot panel is generated from `prepare_all_external_datasets.R`
and currently includes:

```text
BunisHSPC
SegerstolpePancreas
LawlorPancreas
XinPancreas
DarmanisBrain
PollenGlia
WuKidneyHealthy
```

Large or unresolved candidates remain in the registry with `pending` or
`planned` status until their retrieval, gene-symbol handling, label mappings,
and leakage checks are frozen.

## Common Locked Preparation Steps

Each adapter:

1. Loads a local Seurat, SingleCellExperiment, list-RDS, h5ad, 10x directory,
   10x h5, or matrix file declared in the registry.
2. Selects the prespecified `GroundTruthColumn`.
3. Applies an optional fixed label mapping with columns
   `original_label`, `harmonized_label`, `exclude`, and `reason`.
4. Removes unannotated, doublet, low-quality, or explicitly excluded cells.
5. Applies the common Seurat workflow:
   normalization, variable features, scaling, PCA, neighbors, clustering, and
   marker extraction.
6. Selects top marker genes per cluster using the same marker filtering rules
   as the internal benchmark.
7. Computes majority cluster truth and cluster purity.
8. Writes the prepared RDS, cluster audit CSV, dataset audit CSV, hash CSV, and
   HTML audit report.
9. Validates the prepared object against the locked contract and leakage rules.

When source row names are Ensembl identifiers and no gene-symbol column is
available, the preparation engine may apply a locked species-specific symbol
mapping. The audit output records mapping success and unmapped identifiers.

## Strict Eligibility Checks

`validate_prepared_dataset.R` checks:

- object structure and named marker/truth/purity vectors;
- no duplicated cluster identifiers;
- non-empty marker lists;
- matching marker, truth, and purity names;
- purity values in `[0, 1]`;
- tissue and species metadata;
- exclusion from development datasets;
- exclusion from marker-profile development;
- exclusion from Risk-k training data;
- no marker-list inspection before freezing;
- no known donor overlap for confirmatory datasets;
- inclusion decision and confirmatory status consistency.

Datasets with uncertain independence should be marked exploratory or pending,
not silently included in the confirmatory panel.

## Blinded Two-Stage Workflow

Stage A prepares and locks external datasets:

```bash
Rscript benchmarks/external_datasets/prepare_all_external_datasets.R \
  benchmarks/external_datasets/registry.csv \
  data/external_prepared \
  benchmarks/external_validation_manifest.csv
```

To prepare only a named subset, set `DEEPSEEKCELL_EXTERNAL_DATASETS` to a
comma-separated dataset list before running the adapter.

Stage B runs the locked evaluation:

```bash
Rscript benchmarks/run_external_validation.R benchmarks/external_validation_manifest.csv 1
```

The evaluation stage reuses the prepared marker lists and uses the truth labels
only for metric computation.

## Audit Outputs

Each prepared study writes:

```text
results/external_dataset_audits/<Dataset>_cluster_audit.csv
results/external_dataset_audits/<Dataset>_dataset_audit.csv
results/external_dataset_audits/<Dataset>_hashes.csv
results/external_dataset_audits/<Dataset>_validation.csv
results/external_dataset_audits/external_dataset_audit_<Dataset>.html
```

These reports include provenance, cell counts before and after filtering, label
mapping, cluster purity, marker counts, prepared RDS MD5, and eligibility
status.

## Important Rule

Do not use benchmark results, observed difficulty, or LLM outputs to decide
which external datasets remain in the confirmatory panel. Select datasets using
biological, technical, provenance, and licensing criteria before inference, and
report all eligible locked studies.
