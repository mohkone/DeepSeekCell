# Locked External Validation Protocol

This protocol defines the confirmatory validation of DeepSeekCell reliability
specification v1.0. The goal is to test whether the frozen Evidence-k selector
generalises to datasets that were not used to choose marker profiles, aliases,
weights, thresholds, prompts, or refinement rules.

## Frozen Method

The locked method is identified as:

`DeepSeekCell reliability specification v1.0`

The machine-readable specification is stored in:

`inst/extdata/reliability_spec_v1.0.json`

The supplementary audit tables are:

`inst/extdata/marker_profiles_v1.0.csv`

`inst/extdata/marker_aliases_v1.0.csv`

No marker profile, alias, weight, threshold, prompt rule, tie-breaking rule, or
selector definition should be changed after external validation begins. Any
later methodological change should be labelled as a new specification, for
example v1.1, and evaluated separately.

The learned extension is:

`DeepSeekCell reliability model v1.1`

Its machine-readable specification is stored in:

`inst/extdata/reliability_model_spec_v1.1.json`

Risk-k models must be trained only on development benchmark outputs and frozen
as an RDS file before held-out validation. The external-validation lock file
records the model path and MD5 hash when `DEEPSEEKCELL_RELIABILITY_MODEL` is set.

## Development Benchmark

The following datasets are treated as development datasets because they were
used during method construction or internal benchmarking:

- PBMC
- BaronPancreas
- MuraroPancreas
- TasicBrain
- ZeiselBrain
- ZilionisLung

Results on these datasets remain useful for comparison, ablation, debugging,
and manuscript context, but they are not the primary evidence of locked-method
generalisability.

## Held-Out Dataset Requirements

The confirmatory external set should include at least:

- One independent immune dataset.
- One independent pancreas dataset.
- One independent neural dataset.
- One challenging disease or tumour dataset.
- Preferably one dataset from a previously unseen tissue.
- Multiple data sources or consortia, for example Human Cell Atlas,
  CELLxGENE, Tabula Sapiens, or laboratory-specific studies.
- Multiple sequencing platforms where possible, including 10x Chromium,
  Smart-seq2, Drop-seq, Seq-Well, microwell, or inDrop.
- Healthy and disease or tumour conditions when compatible labels are
  available.

Datasets must be selected before running the locked analysis. The selection
should be recorded in `benchmarks/external_validation_manifest.csv` using the
template `benchmarks/external_validation_manifest_template.csv`.

The manifest records study accession, source repository, center, laboratory,
country, sequencing platform, chemistry, disease status, condition, donor
count, cell count, expected cluster count, prospective status, and the selection
rationale. These fields are copied into the result tables and the validation
lock so stratified robustness summaries can be generated without retuning the
method.

Prepared RDS files should be generated with the registry-driven adapter system
in `benchmarks/external_datasets/`. The adapter protocol is documented in:

```text
benchmarks/external_datasets/adapter_protocol.md
```

The initial planned registry includes independent immune, pancreas, neural, and
unseen-tissue entries. Rows should remain `pending` or `planned` until local
source files, label mappings, leakage flags, and eligibility decisions are
fixed before inference.

The first prepared pilot panel currently includes seven locked studies:
BunisHSPC, SegerstolpePancreas, LawlorPancreas, XinPancreas, DarmanisBrain,
PollenGlia, and WuKidneyHealthy. The panel is intentionally labelled a pilot
until the larger 10-20 study confirmatory panel is prepared and frozen.

## Primary Endpoint

The primary endpoint is correction efficiency under a matched refinement budget:

```text
CorrectionEfficiency = (WrongToCorrect - CorrectToWrong) / NRefined
```

The primary comparison for the learned method is Risk-k versus compute-matched
Evidence-k, Random-k, and Confidence-k within each model-dataset-replicate
block. If no trained Risk-k model is supplied, the deterministic v1.0 fallback
comparison is Evidence-k versus Random-k and Confidence-k.

Primary hypothesis:

```text
Risk-k achieves higher correction efficiency than compute-matched Evidence-k,
Confidence-k and Random-k on held-out datasets.
```

## Secondary Endpoints

Secondary outcomes include:

- Wrong-to-correct revisions.
- Correct-to-wrong revisions.
- Error-recovery rate relative to FullRefined.
- Selection precision, recall, specificity, negative predictive value, and MCC.
- Final accuracy, macro-F1, balanced accuracy, ARI, and Cell Ontology clade accuracy.
- Raw versus evidence-adjusted confidence quality: Brier score, binary correctness NLL, ECE, AUROC, and AURC.
- Runtime, prompt tokens, completion tokens, total tokens, API cost, and cost per corrected error.
- Stratified correction efficiency, calibration, runtime, cost, and oracle gap
  by center, sequencing platform, disease status, prospective status, cluster
  count, cell count, and marker-list size.

## Paired Design

For every model-dataset-replicate block:

1. Generate exactly one first-pass LLM response.
2. Cache and hash the raw response.
3. Parse the cached response into first-pass annotations.
4. Apply Plain, Evidence, Calibrated, Random-k, Confidence-k, NoOntology-k, Risk-k, Evidence-k/SelfRefined, and FullRefined to the same parsed response.
5. Use the same refinement budget k for Random-k, Confidence-k, NoOntology-k, Risk-k, and Evidence-k.
6. Record selected clusters, reviewed clusters, label changes, confidence changes, tokens, latency, cost, and response hashes.

Do not independently call the first-pass LLM for different ablation arms.

By default, Evidence-k uses the locked deterministic conflict rule and may
therefore refine zero clusters on an external dataset. For a fixed-budget audit
of selector ranking quality, set either `DEEPSEEKCELL_EXTERNAL_REFINEMENT_BUDGET_K`
or `DEEPSEEKCELL_EXTERNAL_REFINEMENT_BUDGET_FRACTION` before execution. Under
that explicit audit mode, Evidence-k ranks all clusters by the frozen
evidence-conflict priority rule and selects the top-k clusters so that all
selectors receive the same non-zero budget.

## Statistical Analysis

Use dataset-replicate blocks as the primary paired unit. Report:

- Bootstrap confidence intervals over dataset-replicate blocks.
- Paired Wilcoxon signed-rank tests for Risk-k versus Evidence-k, Random-k and Confidence-k.
- Effect sizes for paired differences.
- False-discovery-rate correction for families of secondary tests.

Replicate-level results should be interpreted as stability evidence, not as
fully independent biological datasets.

## Sensitivity Analyses

Sensitivity analyses are secondary robustness checks and must not be used to
select a new best configuration after inspecting held-out performance.

Prespecified confidence-weight variants:

- Default v1.0 weights.
- Equal weights.
- No ontology evidence.
- No tissue evidence.
- Marker-dominant weights.
- LLM-confidence-dominant weights.
- Random weight vectors sampled from the simplex.

Prespecified threshold grid:

```text
tau = 0.25, 0.30, ..., 0.70
```

Report selector rank correlation, selected-cluster Jaccard overlap, correction
efficiency, recovery rate, and correct-to-wrong revisions across the grid.

## Multi-Center and Platform Robustness

After the locked validation run completes, generate reviewer-facing robustness
tables with:

```bash
Rscript benchmarks/analyse_external_validation_robustness.R \
  results/external_validation_results_full.csv \
  results/external_validation_refinement_behavior.csv \
  results/external_validation_confidence_quality.csv \
  results/external_validation_robustness
```

The script writes:

- `*_dataset_scorecard.csv`
- `*_refinement_by_domain.csv`
- `*_runtime_cost_by_domain.csv`
- `*_confidence_by_domain.csv`
- `*_selector_contrasts.csv`
- platform-efficiency, scale-efficiency, and runtime-scaling PDF plots when
  `ggplot2` is available.

These summaries answer whether the frozen reliability framework remains stable
across laboratories, repositories, tissues, technologies, disease conditions,
prospective datasets, dataset scales, and marker-list sizes. They are evidence
analyses, not additional model-selection steps.

## Benchmark Release

After the final benchmark and external validation are complete, generate a
release manifest:

```bash
Rscript benchmarks/build_benchmark_release_manifest.R results/benchmark_release_manifest.csv
```

The manifest records paths, artifact categories, file sizes, MD5 hashes, and
archive recommendations for result tables, figures, debug decisions, cached LLM
responses, validation locks, frozen specifications, and benchmark scripts. Large
prepared dataset caches are excluded by default and should be released only
when the original data license permits redistribution.

## Reporting

The manuscript should keep two questions separate:

- Annotation-method comparison: SingleR, scType, native CellTypist, scmap, and LLM-based annotation methods.
- Refinement-selector comparison: NoRefinement, Random-k, Confidence-k, NoOntology-k, Evidence-k, Risk-k, and FullRefined.

SingleR and CellTypist are annotation methods, not refinement selectors.
