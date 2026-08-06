# Risk-Aware Reliability Model Protocol

This protocol defines the v1.1 learned extension to the frozen deterministic
Evidence-k reliability layer.

## Scientific Motivation

Evidence-k is a fixed deterministic selector. It ranks clusters by evidence
conflict and asks which first-pass annotations deserve a second LLM pass. The
learned Risk-k extension turns this into a supervised decision problem:

```text
P(first-pass annotation is incorrect | reliability features)
```

or, when FullRefined training outcomes are available:

```text
P(annotation improves after refinement | reliability features)
```

The first target learns error risk. The second target learns expected refinement
benefit, which is more directly aligned with correction efficiency.

## Optimization View

For a fixed refinement budget `k`, the ideal selector solves:

```text
S* = argmax_{|S| <= k} sum_{i in S} E[Delta_i]
```

where `Delta_i` is the expected improvement obtained by refining cluster `i`.

Risk-k approximates this objective by ranking clusters according to a learned
probability:

```text
score_i = P(error_i = 1 | x_i)
```

or:

```text
score_i = P(refinement_benefit_i = 1 | x_i)
```

and selecting the top `k` clusters.

## Features

The model uses only first-pass and deterministic evidence features:

- Raw LLM confidence.
- Ontology evidence score.
- Marker evidence score for the predicted label.
- Best marker-profile score.
- Tissue-consistency evidence.
- Consensus evidence.
- Evidence conflict score.
- Marker margin.
- Mixed-profile flag.
- Requires-refinement flag.
- Candidate-label count.
- Whether the first-pass label matches the strongest evidence-supported label.
- Whether the predicted label appears in the candidate set.
- Whether the evidence-supported label appears in the candidate set.
- Candidate/evidence disagreement.
- Candidate agreement score.

The machine-readable feature list is stored in:

`inst/extdata/reliability_model_spec_v1.1.json`

## Training

Train on development benchmark outputs only:

```bash
Rscript benchmarks/train_reliability_model.R results/benchmark_debug results/reliability_model_v1.1_error.rds first_pass_error
```

If FullRefined outcomes are available and the goal is expected correction rather
than error detection:

```bash
Rscript benchmarks/train_reliability_model.R results/benchmark_debug results/reliability_model_v1.1_benefit.rds refinement_benefit
```

The default implementation uses logistic regression for interpretability and
minimal dependency burden. More flexible learners can be evaluated later, but
they should be labelled as separate model variants.

## Locked Evaluation

Before held-out validation:

1. Train the model on development benchmark outputs only.
2. Save the model RDS.
3. Record the model MD5 hash.
4. Set `DEEPSEEKCELL_RELIABILITY_MODEL` to the model path.
5. Run external validation without retraining or retuning.

```bash
set DEEPSEEKCELL_RELIABILITY_MODEL=results/reliability_model_v1.1_error.rds
set DEEPSEEKCELL_RUN_EXTERNAL_VALIDATION=true
Rscript benchmarks/run_external_validation.R benchmarks/external_validation_manifest.csv 3 deepseek
```

## Primary Comparisons

Risk-k should be evaluated under the same per-dataset-replicate budget as
Evidence-k, Confidence-k, and Random-k.

Primary endpoint:

```text
CorrectionEfficiency = (WrongToCorrect - CorrectToWrong) / NRefined
```

Primary comparisons:

- Risk-k versus Evidence-k.
- Risk-k versus Confidence-k.
- Risk-k versus Random-k.

The manuscript should report the deterministic Evidence-k results as the frozen
v1.0 baseline and Risk-k as the learned v1.1 method.
