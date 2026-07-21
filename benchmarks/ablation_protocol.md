# DeepSeekCell Ablation Benchmark Protocol

This protocol is intended for the revised methodological benchmark. It keeps
the algorithmic effect separate from stochastic variation in API responses.

## Core rule

For each dataset and replicate, generate the first-pass LLM response once,
save the raw response text, record its MD5 hash, and reuse that exact response
for every DeepSeekCell ablation arm.

Do not call the first-pass LLM independently for each arm.

## Ablation arms

| Arm | First-pass LLM | Evidence score | Confidence replacement | Second LLM call |
| --- | ---: | ---: | ---: | ---: |
| DeepSeek-Plain | Yes | No | No | No |
| DeepSeekCell-Evidence | Yes | Yes | No | No |
| DeepSeekCell-Calibrated | Yes | Yes | Yes | No |
| DeepSeekCell-RandomK | Yes | Yes | Yes | Random-k |
| DeepSeekCell-ConfidenceK | Yes | Yes | Yes | Lowest-confidence-k |
| DeepSeekCell-NoOntologyK | Yes | Marker/tissue only | Yes | Ontology-disabled evidence-k |
| DeepSeekCell-SelfRefined | Yes | Yes | Yes | Selective |
| DeepSeekCell-FullRefined | Yes | Yes | Yes | All clusters |

Expected behavior:

- Plain, Evidence, and Calibrated must have identical labels.
- Evidence may add conflict and marker-support fields but must not change labels
  or raw confidence.
- Calibrated may change confidence but must not change labels.
- RandomK, ConfidenceK, NoOntologyK, SelfRefined, and FullRefined are allowed
  to change labels through second-pass refinement. RandomK, ConfidenceK, and
  NoOntologyK use the same per-dataset refinement budget as SelfRefined.

## Primary annotation metrics

- Accuracy
- Macro-F1
- Balanced accuracy
- Adjusted Rand index
- Cell Ontology clade accuracy
- Unknown-label rate

The primary gain is:

```text
Delta Macro-F1 = MacroF1(SelfRefined) - MacroF1(Plain)
```

The selector-specific gains are:

```text
Delta_selector    = Metric(SelfRefined) - Metric(RandomK)
Delta_confidence  = Metric(SelfRefined) - Metric(ConfidenceK)
Delta_ontology    = Metric(SelfRefined) - Metric(NoOntologyK)
Delta_full_bound  = Metric(FullRefined) - Metric(SelfRefined)
```

## Refinement behavior metrics

- Number and fraction of clusters flagged
- Selection precision
- Selection recall
- Selection specificity
- Selection negative predictive value
- Selection Matthews correlation coefficient
- Conflict precision and recall for the evidence selector
- Number of refinement calls
- Number of labels changed
- Wrong-to-correct changes
- Correct-to-wrong changes
- Net correction rate

Definitions:

```text
ConflictPrecision = flagged initially incorrect / all flagged
ConflictRecall    = flagged initially incorrect / all initially incorrect
NetCorrectionRate = (wrong_to_correct - correct_to_wrong) / refined_clusters
```

## Confidence quality metrics

Compare raw `LLMConfidence` with final evidence-adjusted `Confidence`.

- Brier score
- Expected calibration error
- Binary correctness negative log-likelihood with clipped probabilities
- AUROC for distinguishing correct from incorrect predictions
- Area under the risk-coverage curve, including the origin at zero coverage
- Coverage and accuracy above confidence thresholds 0.7, 0.8, and 0.9

Use the term "evidence-adjusted confidence" unless probability calibration is
empirically demonstrated with reliability curves or calibration metrics.

## Compute and cost metrics

- First-pass runtime
- Refinement runtime
- Total runtime
- Prompt tokens
- Completion tokens
- Total cost
- Number of second-pass calls
- Cost per corrected error

## Compute-matched controls

To show that selective evidence-guided refinement is responsible for any gain,
the executable benchmark includes compute-matched controls:

- Random-k refinement
- Confidence-only-k refinement
- Ontology-disabled evidence-k refinement
- Evidence-conflict-k refinement
- Full refinement

Here `k` is the number of clusters selected by the evidence-conflict detector.
Random-k, confidence-only-k, and ontology-disabled evidence-k refine exactly
`k` clusters within the same dataset and replicate. Full refinement is retained
as an upper-bound control and therefore uses a larger second-pass budget.
