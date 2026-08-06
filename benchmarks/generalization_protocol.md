# Risk-k Generalization Validation Protocol

This protocol evaluates whether the learned Risk-k reliability model transfers
across biological domains. It does not define a new selector. Instead, it tests
the existing Risk-k selector under held-out tissue, limited-data, calibration,
oracle-gap, and failure-analysis settings.

## Scientific Question

The central question is:

```text
Can a reliability model trained on one set of tissues identify high-value
second-pass refinement targets in unseen tissues?
```

This addresses the main generalization risk for the learned reliability layer:
the model may appear effective only because it was trained and evaluated on
similar biological domains.

## Input Requirements

Run the paired benchmark first so that `results/benchmark_debug` contains the
paired DeepSeekCell ablation debug files, especially:

- `*-Calibrated_debug.csv` or `*-Evidence_debug.csv`
- `*-FullRefined_debug.csv`

The debug rows should include dataset, replicate, cluster, first-pass label,
truth label, reliability features, and FullRefined outcomes. Current benchmark
outputs also record tissue and species metadata for each debug row.

## Command

```bash
Rscript benchmarks/analyse_reliability_generalization.R results/benchmark_debug results/reliability_generalization first_pass_error 10
```

Arguments:

1. Debug directory, usually `results/benchmark_debug`.
2. Output prefix, usually `results/reliability_generalization`.
3. Risk target: `first_pass_error` or `refinement_benefit`.
4. Learning-curve repeats.

## Analyses

### Leave-One-Tissue-Out Transfer

For each held-out tissue, train Risk-k on all other tissues and evaluate on the
held-out tissue. Compare Risk-k with compute-matched Evidence-k, Confidence-k,
Random-k, and FullRefined.

Primary endpoint:

```text
CorrectionEfficiency = (WrongToCorrect - CorrectToWrong) / NRefined
```

Secondary endpoints include recovery fraction, refinement budget, wrong-to-
correct revisions, correct-to-wrong revisions, selection precision, and
selection recall.

### Pairwise Tissue Transfer

Train Risk-k on one tissue and test on each other tissue. This creates a tissue
transfer matrix that can reveal asymmetric domain transfer, for example whether
immune-trained reliability features transfer to lung better than to brain.

### Reliability Learning Curves

For each held-out tissue, sample increasing fractions of the available training
blocks from the remaining tissues and evaluate held-out correction efficiency.
This estimates how much paired benchmark data Risk-k needs before performance
stabilizes.

### Domain-Shift Analysis

For every held-out tissue, compare training and test distributions for:

- raw LLM confidence;
- ontology evidence;
- marker evidence;
- best competing marker evidence;
- evidence conflict score;
- candidate count;
- conflict rate;
- prediction-label entropy;
- prediction-label Jensen-Shannon divergence.

The script also reports correlations between these shift measures and Risk-k
correction efficiency.

### Calibration Comparison

Compare three estimates of first-pass error risk:

- raw LLM risk, defined as `1 - LLMConfidence`;
- evidence-adjusted risk, defined as `1 - EvidenceAdjustedConfidence`;
- learned Risk-k probability.

The script writes equal-frequency calibration bins and Brier scores for each
score type.

### Failure Taxonomy

First-pass errors are assigned to interpretable categories:

- mixed-cell cluster;
- candidate/evidence disagreement;
- evidence conflict;
- weak marker evidence;
- ontology ambiguity;
- low LLM confidence;
- biologically related or other.

This table is intended for biological interpretation of residual errors, not
for selector training.

### Oracle-Gap Analysis

FullRefined is treated as an empirical oracle for the corrections attainable by
second-pass refinement. The analysis reports how much of that oracle correction
set is recovered by Risk-k and the compute-matched controls.

## Output Files

The output prefix generates:

```text
*_features.csv
*_leave_one_tissue_out.csv
*_pairwise_tissue_transfer.csv
*_learning_curve.csv
*_domain_shift.csv
*_domain_shift_correlations.csv
*_calibration_comparison_bins.csv
*_calibration_comparison_metrics.csv
*_failure_taxonomy_clusters.csv
*_failure_taxonomy_summary.csv
*_oracle_gap.csv
*_cross_tissue_transfer_matrix.pdf
*_learning_curve.pdf
*_oracle_gap.pdf
*_calibration_comparison.pdf
*_failure_taxonomy.pdf
```

## Interpretation

The strongest validation result would show that Risk-k maintains higher
correction efficiency than Evidence-k, Confidence-k, and Random-k when the test
tissue was absent from training. A useful negative result would identify the
types of tissue shift or failure categories where learned reliability does not
transfer, which can guide future external validation.
