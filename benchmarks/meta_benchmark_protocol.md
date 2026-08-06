# Meta-Benchmark Protocol

This protocol defines secondary analyses that ask when selective LLM refinement
is useful. It does not define a new selector, tune thresholds, or retrain the
method. It turns completed benchmark outputs into a dataset-difficulty taxonomy
and relates difficulty to accuracy, calibration, correction efficiency, runtime,
cost, and biological error modes.

## Scientific Question

The meta-benchmark asks:

```text
Under which biological and technical conditions does selective LLM refinement
provide the greatest benefit?
```

This complements the primary benchmark, which asks which method performs best
under a fixed protocol.

## Command

After running the paired benchmark, execute:

```bash
Rscript benchmarks/analyse_meta_benchmark.R results results/meta_benchmark
```

If a frozen Risk-k model is available and risk-score decision curves should be
included:

```bash
Rscript benchmarks/analyse_meta_benchmark.R \
  results \
  results/meta_benchmark \
  results/reliability_model_v1.1_error.rds
```

## Dataset Difficulty Profile

For each dataset, the script computes a difficulty profile using available
benchmark metadata:

- number of cells, genes, and clusters;
- mean, median, and minimum cluster purity;
- marker-list size;
- low-purity cluster rate;
- label granularity;
- first-pass error rate;
- unknown-label rate;
- ontology ambiguity;
- marker ambiguity;
- tissue ambiguity;
- evidence-conflict rate;
- candidate/evidence disagreement rate;
- FullRefined correction opportunity.

Two scores are reported:

- `IntrinsicDifficultyScore`: based on dataset, marker, ontology, tissue, and
  label-structure properties.
- `ObservedDifficultyScore`: augments intrinsic difficulty with first-pass
  error, unknown rate, and calibration loss.

Datasets are assigned a low, medium, or high difficulty class using tertiles of
the observed score.

## Benefit and Difficulty Relationships

The script joins difficulty profiles to:

- annotation metrics for every method;
- selector-level refinement behaviour;
- confidence-quality summaries.

This supports analyses such as:

```text
Difficulty -> Evidence-k correction efficiency
Difficulty -> Risk-k correction efficiency
Difficulty -> Plain annotation accuracy
Difficulty -> calibration error
```

## Descriptive Benefit Predictor

When FullRefined outcomes are available, the script fits a descriptive logistic
model for:

```text
P(refinement is beneficial | reliability and difficulty features)
```

This model is used only for interpretation. It is not a deployed selector and
must not be used to retune Risk-k or Evidence-k after inspecting benchmark
performance.

## Pareto Analysis

The script identifies non-dominated operating points for:

- Macro-F1 versus API cost;
- Macro-F1 versus runtime;
- Macro-F1 versus token use;
- accuracy versus second-pass calls;
- correction efficiency versus refinement runtime;
- correction efficiency versus refinement tokens;
- oracle recovery versus second-pass calls.

These tables let users choose an operating point under practical compute or API
constraints.

## Biological Error Taxonomy

First-pass errors are categorized into biologically interpretable modes:

- unknown or abstained predictions;
- insufficient marker evidence;
- evidence conflict;
- tissue mismatch;
- ontology ambiguity;
- developmental-stage confusion;
- immune subtype confusion;
- sibling or related Cell Ontology terms;
- tumour/normal confusion;
- ontology synonym or label-mapping ambiguity;
- biologically related or other.

This analysis is intended to explain residual errors and refinement successes,
not to construct a new rule system.

## Reliability and Decision Curves

Reliability curves compare observed error rates across bins of:

- raw LLM confidence;
- evidence-adjusted confidence;
- ontology evidence;
- marker evidence;
- tissue evidence;
- evidence-conflict score;
- marker margin;
- Risk-k probability, when a frozen Risk-k model is supplied.

Decision curves evaluate net correction and oracle recovery as the refinement
budget increases from 0% to 100% of clusters. This estimates when additional
refinement begins to yield diminishing returns.

## Output Files

The output prefix generates:

```text
*_dataset_difficulty.csv
*_method_gains_vs_difficulty.csv
*_selector_benefit_vs_difficulty.csv
*_pareto_annotation.csv
*_pareto_refinement.csv
*_biological_error_taxonomy_clusters.csv
*_biological_error_taxonomy_summary.csv
*_reliability_curves.csv
*_decision_curves.csv
*_benefit_predictor_coefficients.csv
*_benefit_predictor_predictions.csv
*_difficulty_vs_plain_macroF1.pdf
*_difficulty_vs_selector_efficiency.pdf
*_pareto_macroF1_cost.pdf
*_reliability_curves.pdf
*_decision_curves.pdf
*_biological_error_taxonomy.pdf
```

## Interpretation

The strongest manuscript use is not another claim that one method is best.
Instead, the meta-benchmark should support statements such as:

```text
Selective refinement provided the largest benefit on datasets with low marker
specificity, high ontology ambiguity, and high first-pass unknown rates, whereas
well-separated immune datasets required little or no second-pass reasoning.
```
