# R/reliability_model.R
# Learned reliability models for risk-aware selective refinement.

RELIABILITY_MODEL_FEATURES <- c(
  "LLMConfidence",
  "OntologyEvidenceScore",
  "MarkerEvidenceScore",
  "BestMarkerEvidenceScore",
  "TissueEvidenceScore",
  "ConsensusEvidenceScore",
  "EvidenceConflictScore",
  "MarkerMargin",
  "IsMixedNumeric",
  "RequiresRefinementNumeric",
  "CandidateCount",
  "PredictedMatchesEvidenceBest",
  "PredictedInCandidateSet",
  "EvidenceBestInCandidateSet",
  "CandidateEvidenceDisagreement",
  "CandidateAgreementScore"
)

#' Extract reliability-model features from first-pass annotations
#'
#' Converts first-pass annotations and deterministic evidence scores into a
#' numeric feature matrix suitable for learning first-pass error risk or
#' expected refinement benefit.
#'
#' @param annotations Annotation data frame.
#' @param markers Optional named marker list. When supplied and evidence columns
#' are absent, deterministic evidence scores are computed first.
#' @param tissue Optional tissue context.
#' @param truth Optional named truth vector keyed by cluster id. When supplied,
#' `InitiallyCorrect` and `FirstPassIncorrect` are added.
#' @param full_refined_annotations Optional FullRefined annotations used to
#' derive `FullRefinedCorrect` and `RefinementBeneficial`.
#' @param use_ontology_evidence Whether ontology evidence contributes if
#' evidence scores must be computed.
#' @return Data frame with numeric reliability features and optional targets.
#' @export
extract_reliability_features <- function(annotations,
                                         markers = NULL,
                                         tissue = NULL,
                                         truth = NULL,
                                         full_refined_annotations = NULL,
                                         use_ontology_evidence = TRUE) {
  if (!is.data.frame(annotations) || nrow(annotations) == 0) {
    return(data.frame())
  }

  needs_evidence <- !all(c(
    "OntologyEvidenceScore",
    "MarkerEvidenceScore",
    "BestMarkerEvidenceScore",
    "TissueEvidenceScore",
    "ConsensusEvidenceScore",
    "EvidenceConflictScore",
    "RequiresRefinement"
  ) %in% names(annotations))

  if (isTRUE(needs_evidence) && !is.null(markers)) {
    annotations <- calibrate_annotation_confidence(
      annotations,
      markers,
      tissue = tissue,
      use_ontology_evidence = use_ontology_evidence
    )
  }

  n <- nrow(annotations)
  clusters <- as.character(annotations$Cluster %||% paste0("Cluster", seq_len(n)))
  llm_confidence <- as_confidence(annotations$LLMConfidence %||% annotations$Confidence)
  ontology_score <- as_confidence(annotations$OntologyEvidenceScore %||% rep(0.5, n))
  marker_score <- as_confidence(annotations$MarkerEvidenceScore %||% rep(0, n))
  best_marker_score <- as_confidence(annotations$BestMarkerEvidenceScore %||% rep(0, n))
  tissue_score <- as_confidence(annotations$TissueEvidenceScore %||% rep(0.6, n))
  consensus_score <- as_confidence(annotations$ConsensusEvidenceScore %||% rep(1, n))
  conflict_score <- suppressWarnings(as.numeric(annotations$EvidenceConflictScore %||% rep(0, n)))
  conflict_score[is.na(conflict_score)] <- 0
  marker_margin <- pmax(best_marker_score - marker_score, 0)

  candidate_features <- .candidate_reliability_features(annotations)

  out <- data.frame(
    Cluster = clusters,
    LLMConfidence = llm_confidence,
    OntologyEvidenceScore = ontology_score,
    MarkerEvidenceScore = marker_score,
    BestMarkerEvidenceScore = best_marker_score,
    TissueEvidenceScore = tissue_score,
    ConsensusEvidenceScore = consensus_score,
    EvidenceConflictScore = conflict_score,
    MarkerMargin = marker_margin,
    IsMixedNumeric = as.numeric(as_flag(annotations$IsMixed %||% rep(FALSE, n))),
    RequiresRefinementNumeric = as.numeric(as_flag(annotations$RequiresRefinement %||% rep(FALSE, n))),
    candidate_features,
    stringsAsFactors = FALSE
  )

  out$IsMixedNumeric[is.na(out$IsMixedNumeric)] <- 0
  out$RequiresRefinementNumeric[is.na(out$RequiresRefinementNumeric)] <- 0

  if (!is.null(truth)) {
    truth_aligned <- truth[clusters]
    initially_correct <- normalize_cell_type(annotations$CellType) ==
      normalize_cell_type(truth_aligned)
    out$InitiallyCorrect <- initially_correct
    out$FirstPassIncorrect <- !initially_correct
  }

  if (!is.null(full_refined_annotations) && !is.null(truth)) {
    refined_idx <- match(clusters, full_refined_annotations$Cluster)
    refined_labels <- rep(NA_character_, length(clusters))
    valid <- !is.na(refined_idx)
    refined_labels[valid] <- as.character(full_refined_annotations$CellType[refined_idx[valid]])
    full_correct <- normalize_cell_type(refined_labels) == normalize_cell_type(truth[clusters])
    out$FullRefinedCorrect <- full_correct
    out$RefinementBeneficial <- !out$InitiallyCorrect & full_correct
    out$RefinementHarmful <- out$InitiallyCorrect & !full_correct
  }

  out
}

#' Train a learned reliability model
#'
#' Fits a transparent logistic model for either first-pass error risk or
#' expected refinement benefit. The returned object can be passed to
#' `select_refinement_candidates(strategy = "risk")`.
#'
#' @param features Data frame from [extract_reliability_features()] or another
#' table containing the required numeric feature columns.
#' @param target One of `"first_pass_error"` or `"refinement_benefit"`.
#' @param target_column Optional explicit binary target column.
#' @return Object of class `deepseekcell_reliability_model`.
#' @export
train_reliability_model <- function(features,
                                    target = c("first_pass_error", "refinement_benefit"),
                                    target_column = NULL,
                                    feature_columns = RELIABILITY_MODEL_FEATURES) {
  target <- match.arg(target)
  if (!is.data.frame(features) || nrow(features) == 0) {
    stop("features must be a non-empty data frame.", call. = FALSE)
  }

  y <- .resolve_reliability_target(features, target, target_column)
  feature_columns <- unique(as.character(feature_columns))
  feature_columns <- feature_columns[nzchar(feature_columns)]
  if (length(feature_columns) == 0) {
    stop("feature_columns must contain at least one feature.", call. = FALSE)
  }

  x <- .prepare_reliability_feature_frame(
    features,
    feature_columns = feature_columns
  )
  feature_medians <- attr(x, "feature_medians")
  keep <- !is.na(y)
  y <- as.integer(y[keep])
  x <- x[keep, , drop = FALSE]
  attr(x, "feature_medians") <- feature_medians

  if (length(y) == 0) {
    stop("No non-missing target values are available.", call. = FALSE)
  }

  feature_means <- vapply(x, mean, numeric(1), na.rm = TRUE)
  feature_sds <- vapply(x, stats::sd, numeric(1), na.rm = TRUE)
  feature_means[is.na(feature_means)] <- 0
  feature_sds[is.na(feature_sds)] <- 0

  fallback_rate <- mean(y)
  fit <- NULL
  formula <- stats::as.formula(paste("Target ~", paste(feature_columns, collapse = " + ")))

  if (length(unique(y)) > 1 && nrow(x) > length(feature_columns)) {
    training_df <- cbind(Target = y, x)
    fit <- tryCatch(
      stats::glm(formula, data = training_df, family = stats::binomial()),
      error = function(e) NULL
    )
  }

  model <- list(
    model_id = paste0("DeepSeekCell reliability model v1.1-", target),
    target = target,
    target_column = target_column,
    model_type = "logistic_regression",
    feature_columns = feature_columns,
    feature_medians = feature_medians,
    feature_means = feature_means,
    feature_sds = feature_sds,
    fallback_rate = fallback_rate,
    formula = deparse(formula),
    fit = fit,
    n_training_rows = length(y),
    positive_rate = fallback_rate,
    trained_at = Sys.time()
  )
  class(model) <- "deepseekcell_reliability_model"
  model
}

#' Predict reliability risk or expected refinement benefit
#'
#' @param model A model returned by [train_reliability_model()].
#' @param features Data frame containing reliability feature columns.
#' @return Numeric vector of predicted risk/benefit probabilities in `[0, 1]`.
#' @export
predict_reliability_risk <- function(model, features) {
  if (!inherits(model, "deepseekcell_reliability_model")) {
    stop("model must be a deepseekcell_reliability_model object.", call. = FALSE)
  }
  if (!is.data.frame(features) || nrow(features) == 0) {
    return(numeric())
  }

  feature_columns <- model$feature_columns %||% RELIABILITY_MODEL_FEATURES
  x <- .prepare_reliability_feature_frame(
    features,
    feature_columns = feature_columns,
    feature_medians = model$feature_medians
  )

  if (is.null(model$fit)) {
    return(rep(pmin(pmax(model$fallback_rate, 0), 1), nrow(x)))
  }

  pred <- suppressWarnings(
    stats::predict(model$fit, newdata = x, type = "response")
  )
  pred <- as.numeric(pred)
  pred[is.na(pred)] <- model$fallback_rate
  pmin(pmax(pred, 0), 1)
}

#' Explain a learned reliability model
#'
#' Reports logistic coefficients, standardized coefficients, odds ratios and
#' Wald confidence intervals for a trained Risk-k model.
#'
#' @param model A model returned by [train_reliability_model()].
#' @param conf_level Confidence level for Wald intervals.
#' @return Data frame of feature-level model explanations.
#' @export
explain_reliability_model <- function(model, conf_level = 0.95) {
  if (!inherits(model, "deepseekcell_reliability_model")) {
    stop("model must be a deepseekcell_reliability_model object.", call. = FALSE)
  }
  if (is.null(model$fit)) {
    return(data.frame())
  }

  coef_table <- tryCatch(
    summary(model$fit)$coefficients,
    error = function(e) NULL
  )
  if (is.null(coef_table) || nrow(coef_table) == 0) {
    return(data.frame())
  }

  features <- setdiff(rownames(coef_table), "(Intercept)")
  feature_sds <- model$feature_sds %||% stats::setNames(rep(1, length(features)), features)
  feature_means <- model$feature_means %||% stats::setNames(rep(0, length(features)), features)
  z_value <- stats::qnorm(1 - (1 - conf_level) / 2)

  estimate <- coef_table[features, "Estimate"]
  std_error <- coef_table[features, "Std. Error"]
  sds <- feature_sds[features]
  sds[is.na(sds)] <- 0
  means <- feature_means[features]
  means[is.na(means)] <- 0
  standardized <- estimate * sds

  out <- data.frame(
    ModelID = model$model_id %||% NA_character_,
    Target = model$target %||% NA_character_,
    Feature = features,
    Coefficient = as.numeric(estimate),
    StdError = as.numeric(std_error),
    ZValue = as.numeric(coef_table[features, "z value"]),
    PValue = as.numeric(coef_table[features, "Pr(>|z|)"]),
    OddsRatio = .safe_exp(estimate),
    OddsRatioLower = .safe_exp(estimate - z_value * std_error),
    OddsRatioUpper = .safe_exp(estimate + z_value * std_error),
    FeatureMean = as.numeric(means),
    FeatureSD = as.numeric(sds),
    StandardizedCoefficient = as.numeric(standardized),
    StandardizedOddsRatio = .safe_exp(standardized),
    Direction = ifelse(estimate >= 0, "increases_predicted_risk", "decreases_predicted_risk"),
    stringsAsFactors = FALSE
  )

  out <- out[order(abs(out$StandardizedCoefficient), decreasing = TRUE), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Compute linear SHAP-style reliability contributions
#'
#' For logistic regression, the per-feature log-odds contribution is computed as
#' `beta_j * (x_ij - mean_j)`. This is a transparent linear analogue of SHAP for
#' the fitted Risk-k model.
#'
#' @param model A model returned by [train_reliability_model()].
#' @param features Reliability feature data frame.
#' @return Long data frame of cluster-feature log-odds contributions.
#' @export
compute_reliability_contributions <- function(model, features) {
  if (!inherits(model, "deepseekcell_reliability_model")) {
    stop("model must be a deepseekcell_reliability_model object.", call. = FALSE)
  }
  if (is.null(model$fit) || !is.data.frame(features) || nrow(features) == 0) {
    return(data.frame())
  }

  feature_columns <- model$feature_columns %||% RELIABILITY_MODEL_FEATURES
  x <- .prepare_reliability_feature_frame(
    features,
    feature_columns = feature_columns,
    feature_medians = model$feature_medians
  )
  beta <- stats::coef(model$fit)
  intercept <- beta[["(Intercept)"]] %||% 0
  beta <- beta[feature_columns]
  beta[is.na(beta)] <- 0

  means <- model$feature_means[feature_columns] %||%
    vapply(x, mean, numeric(1), na.rm = TRUE)
  means[is.na(means)] <- 0

  centered <- sweep(as.matrix(x), 2, means, FUN = "-")
  contribution_matrix <- sweep(centered, 2, beta, FUN = "*")
  baseline_log_odds <- as.numeric(intercept + sum(beta * means))
  predicted_log_odds <- baseline_log_odds + rowSums(contribution_matrix)
  predicted_risk <- stats::plogis(predicted_log_odds)
  baseline_risk <- stats::plogis(baseline_log_odds)

  metadata <- features[intersect(
    c("Dataset", "Replicate", "LLMBackend", "LLMModelID", "Cluster"),
    names(features)
  )]
  if (!"Cluster" %in% names(metadata)) {
    metadata$Cluster <- paste0("Row", seq_len(nrow(features)))
  }

  rows <- lapply(feature_columns, function(feature) {
    data.frame(
      metadata,
      Feature = feature,
      FeatureValue = x[[feature]],
      BaselineValue = means[[feature]],
      ContributionLogOdds = as.numeric(contribution_matrix[, feature]),
      AbsContributionLogOdds = abs(as.numeric(contribution_matrix[, feature])),
      PredictedRisk = predicted_risk,
      BaselineRisk = baseline_risk,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

#' Summarise global reliability contribution importance
#'
#' @param contributions Output from [compute_reliability_contributions()].
#' @return Feature-level mean absolute and signed contributions.
#' @export
summarise_reliability_contributions <- function(contributions) {
  if (!is.data.frame(contributions) || nrow(contributions) == 0) {
    return(data.frame())
  }

  out <- stats::aggregate(
    cbind(
      MeanAbsContribution = contributions$AbsContributionLogOdds,
      MeanContribution = contributions$ContributionLogOdds
    ),
    by = list(Feature = contributions$Feature),
    FUN = mean,
    na.rm = TRUE
  )
  out <- out[order(out$MeanAbsContribution, decreasing = TRUE), , drop = FALSE]
  out$Rank <- seq_len(nrow(out))
  rownames(out) <- NULL
  out
}

#' Evaluate reliability-model calibration
#'
#' @param model A model returned by [train_reliability_model()].
#' @param features Reliability feature data frame with target columns.
#' @param target One of `"first_pass_error"` or `"refinement_benefit"`.
#' @param target_column Optional explicit binary target column.
#' @param n_bins Number of equal-frequency calibration bins.
#' @return List with `metrics`, `bins`, and per-row `predictions`.
#' @export
evaluate_reliability_calibration <- function(model,
                                             features,
                                             target = c("first_pass_error", "refinement_benefit"),
                                             target_column = NULL,
                                             n_bins = 10) {
  target <- match.arg(target)
  if (!inherits(model, "deepseekcell_reliability_model")) {
    stop("model must be a deepseekcell_reliability_model object.", call. = FALSE)
  }
  if (!is.data.frame(features) || nrow(features) == 0) {
    return(list(metrics = data.frame(), bins = data.frame(), predictions = data.frame()))
  }

  observed <- .resolve_reliability_target(features, target, target_column)
  predicted <- predict_reliability_risk(model, features)
  ok <- !is.na(observed) & !is.na(predicted)
  observed <- as.numeric(observed[ok])
  predicted <- predicted[ok]

  if (length(observed) == 0) {
    return(list(metrics = data.frame(), bins = data.frame(), predictions = data.frame()))
  }

  bins <- .reliability_calibration_bins(observed, predicted, n_bins)
  ece <- sum((bins$N / sum(bins$N)) * abs(bins$ObservedEventRate - bins$MeanPredictedRisk))
  slope <- .calibration_slope_intercept(observed, predicted)

  metrics <- data.frame(
    ModelID = model$model_id %||% NA_character_,
    Target = target,
    N = length(observed),
    EventRate = mean(observed),
    MeanPredictedRisk = mean(predicted),
    Brier = mean((predicted - observed)^2),
    BinaryNLL = .binary_nll(observed, predicted),
    ECE = ece,
    CalibrationIntercept = slope$intercept,
    CalibrationSlope = slope$slope,
    stringsAsFactors = FALSE
  )

  predictions <- data.frame(
    Cluster = (features$Cluster %||% paste0("Row", seq_len(nrow(features))))[ok],
    Observed = observed,
    PredictedRisk = predicted,
    stringsAsFactors = FALSE
  )

  list(metrics = metrics, bins = bins, predictions = predictions)
}

#' Run reliability feature-ablation analysis
#'
#' Trains full and feature-group-ablated Risk-k models, then evaluates
#' calibration and fixed-budget selection behavior on a validation feature table.
#'
#' @param training_features Development feature table.
#' @param evaluation_features Optional evaluation feature table. Defaults to
#' `training_features`.
#' @param target One of `"first_pass_error"` or `"refinement_benefit"`.
#' @param feature_groups Named list of feature groups to remove.
#' @param budget Optional fixed k per block. If `NULL`, each block uses the
#' number of deterministic v1.0 `RequiresRefinement` clusters.
#' @return Data frame comparing full and ablated models.
#' @export
run_reliability_feature_ablation <- function(training_features,
                                             evaluation_features = NULL,
                                             target = c("first_pass_error", "refinement_benefit"),
                                             feature_groups = NULL,
                                             budget = NULL) {
  target <- match.arg(target)
  if (is.null(evaluation_features)) {
    evaluation_features <- training_features
  }
  if (is.null(feature_groups)) {
    feature_groups <- .default_reliability_feature_groups()
  }

  variants <- c(list(full = character()), feature_groups)
  rows <- lapply(names(variants), function(name) {
    removed <- intersect(variants[[name]], RELIABILITY_MODEL_FEATURES)
    feature_columns <- setdiff(RELIABILITY_MODEL_FEATURES, removed)
    model <- train_reliability_model(
      training_features,
      target = target,
      feature_columns = feature_columns
    )
    calibration <- evaluate_reliability_calibration(
      model,
      evaluation_features,
      target = target
    )$metrics
    selection <- .evaluate_risk_selection(
      model,
      evaluation_features,
      budget = budget
    )

    data.frame(
      Variant = name,
      RemovedFeatures = paste(removed, collapse = ";"),
      NFeatures = length(feature_columns),
      calibration,
      selection,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.resolve_reliability_target <- function(features, target, target_column = NULL) {
  if (!is.null(target_column)) {
    if (!target_column %in% names(features)) {
      stop("target_column is not present in features: ", target_column, call. = FALSE)
    }
    return(.as_binary_target(features[[target_column]]))
  }

  if (identical(target, "first_pass_error")) {
    if ("FirstPassIncorrect" %in% names(features)) {
      return(.as_binary_target(features$FirstPassIncorrect))
    }
    if ("InitiallyCorrect" %in% names(features)) {
      return(!.as_binary_target(features$InitiallyCorrect))
    }
    stop(
      "First-pass error training requires FirstPassIncorrect or InitiallyCorrect.",
      call. = FALSE
    )
  }

  if ("RefinementBeneficial" %in% names(features)) {
    return(.as_binary_target(features$RefinementBeneficial))
  }
  if (all(c("InitiallyCorrect", "FullRefinedCorrect") %in% names(features))) {
    initial <- .as_binary_target(features$InitiallyCorrect)
    full <- .as_binary_target(features$FullRefinedCorrect)
    return(!initial & full)
  }

  stop(
    "Refinement-benefit training requires RefinementBeneficial or both ",
    "InitiallyCorrect and FullRefinedCorrect.",
    call. = FALSE
  )
}

.prepare_reliability_feature_frame <- function(features,
                                               feature_columns = RELIABILITY_MODEL_FEATURES,
                                               feature_medians = NULL) {
  missing <- setdiff(feature_columns, names(features))
  for (column in missing) {
    features[[column]] <- NA_real_
  }

  x <- features[feature_columns]
  for (column in names(x)) {
    x[[column]] <- suppressWarnings(as.numeric(x[[column]]))
  }

  if (is.null(feature_medians)) {
    feature_medians <- vapply(x, function(column) {
      value <- stats::median(column, na.rm = TRUE)
      if (is.na(value)) 0 else value
    }, numeric(1))
  }

  for (column in names(x)) {
    missing_values <- is.na(x[[column]])
    if (any(missing_values)) {
      x[[column]][missing_values] <- feature_medians[[column]] %||% 0
    }
  }

  attr(x, "feature_medians") <- feature_medians
  x
}

.candidate_reliability_features <- function(annotations) {
  n <- nrow(annotations)
  candidates <- annotations$CandidateCellTypes %||% rep("", n)
  candidate_lists <- lapply(candidates, .split_candidate_labels)
  predicted <- as.character(annotations$CellType %||% rep("", n))
  best <- as.character(annotations$EvidenceBestCellType %||% rep("", n))

  predicted_matches_best <- vapply(seq_len(n), function(i) {
    length(.match_profile_name(predicted[[i]], best[[i]])) > 0
  }, logical(1))

  predicted_in_candidates <- vapply(seq_len(n), function(i) {
    .candidate_contains_label(candidate_lists[[i]], predicted[[i]])
  }, logical(1))

  evidence_best_in_candidates <- vapply(seq_len(n), function(i) {
    .candidate_contains_label(candidate_lists[[i]], best[[i]])
  }, logical(1))

  candidate_count <- vapply(candidate_lists, length, integer(1))
  candidate_agreement <- ifelse(
    predicted_matches_best,
    1,
    ifelse(evidence_best_in_candidates, 0.5, 0)
  )
  candidate_agreement[candidate_count == 0] <- 0.5

  data.frame(
    CandidateCount = candidate_count,
    PredictedMatchesEvidenceBest = as.numeric(predicted_matches_best),
    PredictedInCandidateSet = as.numeric(predicted_in_candidates),
    EvidenceBestInCandidateSet = as.numeric(evidence_best_in_candidates),
    CandidateEvidenceDisagreement = as.numeric(evidence_best_in_candidates & !predicted_matches_best),
    CandidateAgreementScore = candidate_agreement,
    stringsAsFactors = FALSE
  )
}

.split_candidate_labels <- function(x) {
  x <- as.character(x %||% "")
  if (!nzchar(trimws(x))) {
    return(character())
  }
  out <- unlist(strsplit(x, "\\s*;\\s*|\\s*,\\s*|\\s*\\|\\s*", perl = TRUE))
  out <- trimws(out)
  unique(out[nzchar(out) & !is.na(out)])
}

.candidate_contains_label <- function(candidates, label) {
  if (length(candidates) == 0 || is.null(label) || is.na(label) || !nzchar(label)) {
    return(FALSE)
  }
  candidate_norm <- normalize_cell_type(candidates)
  label_norm <- normalize_cell_type(label)
  any(candidate_norm == label_norm) ||
    any(vapply(candidate_norm, function(candidate) {
      grepl(candidate, label_norm, fixed = TRUE) || grepl(label_norm, candidate, fixed = TRUE)
    }, logical(1)))
}

.as_binary_target <- function(x) {
  if (is.logical(x)) {
    return(x)
  }
  if (is.numeric(x)) {
    return(x > 0)
  }
  tolower(as.character(x)) %in% c("true", "1", "yes", "y")
}

.safe_exp <- function(x) {
  exp(pmin(pmax(as.numeric(x), -50), 50))
}

.binary_nll <- function(observed, predicted, eps = 1e-6) {
  predicted <- pmin(pmax(predicted, eps), 1 - eps)
  -mean(observed * log(predicted) + (1 - observed) * log(1 - predicted))
}

.reliability_calibration_bins <- function(observed, predicted, n_bins = 10) {
  n_bins <- max(1, suppressWarnings(as.integer(n_bins[1] %||% 10)))
  n <- length(observed)
  ord <- order(predicted)
  bin_id <- integer(n)
  bin_id[ord] <- ceiling(seq_len(n) / n * n_bins)
  bin_id <- pmax(pmin(bin_id, n_bins), 1)

  rows <- lapply(sort(unique(bin_id)), function(bin) {
    idx <- bin_id == bin
    data.frame(
      Bin = bin,
      N = sum(idx),
      MeanPredictedRisk = mean(predicted[idx]),
      ObservedEventRate = mean(observed[idx]),
      MinPredictedRisk = min(predicted[idx]),
      MaxPredictedRisk = max(predicted[idx]),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

.calibration_slope_intercept <- function(observed, predicted) {
  out <- list(intercept = NA_real_, slope = NA_real_)
  if (length(unique(observed)) < 2 || length(unique(predicted)) < 2) {
    return(out)
  }

  logit_pred <- stats::qlogis(pmin(pmax(predicted, 1e-6), 1 - 1e-6))
  fit <- tryCatch(
    stats::glm(observed ~ logit_pred, family = stats::binomial()),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(out)
  }

  coefs <- stats::coef(fit)
  out$intercept <- as.numeric(coefs[[1]] %||% NA_real_)
  out$slope <- as.numeric(coefs[[2]] %||% NA_real_)
  out
}

.default_reliability_feature_groups <- function() {
  list(
    no_llm_confidence = "LLMConfidence",
    no_ontology = "OntologyEvidenceScore",
    no_marker = c(
      "MarkerEvidenceScore",
      "BestMarkerEvidenceScore",
      "EvidenceConflictScore",
      "MarkerMargin",
      "PredictedMatchesEvidenceBest",
      "EvidenceBestInCandidateSet",
      "CandidateEvidenceDisagreement",
      "CandidateAgreementScore"
    ),
    no_tissue = "TissueEvidenceScore",
    no_consensus = c("ConsensusEvidenceScore", "RequiresRefinementNumeric"),
    no_candidates = c(
      "CandidateCount",
      "PredictedInCandidateSet",
      "EvidenceBestInCandidateSet",
      "CandidateEvidenceDisagreement",
      "CandidateAgreementScore"
    )
  )
}

.evaluate_risk_selection <- function(model, features, budget = NULL) {
  if (!all(c("InitiallyCorrect", "FullRefinedCorrect") %in% names(features))) {
    return(data.frame(
      NBlocks = NA_integer_,
      NRefined = NA_integer_,
      SelectionPrecision = NA_real_,
      SelectionRecall = NA_real_,
      WrongToCorrect = NA_integer_,
      CorrectToWrong = NA_integer_,
      CorrectionEfficiency = NA_real_
    ))
  }

  risk <- predict_reliability_risk(model, features)
  initially_correct <- .as_binary_target(features$InitiallyCorrect)
  full_correct <- .as_binary_target(features$FullRefinedCorrect)
  block <- .reliability_block_id(features)
  selected <- rep(FALSE, nrow(features))

  for (block_id in unique(block)) {
    idx <- which(block == block_id)
    k <- if (is.null(budget)) {
      requires <- suppressWarnings(as.numeric(features$RequiresRefinementNumeric[idx] %||% 0))
      sum(requires > 0, na.rm = TRUE)
    } else {
      suppressWarnings(as.integer(budget[1]))
    }
    if (is.na(k) || k <= 0) {
      next
    }
    k <- min(k, length(idx))
    ord <- idx[order(risk[idx], decreasing = TRUE, na.last = TRUE)]
    selected[head(ord, k)] <- TRUE
  }

  n_refined <- sum(selected)
  initially_wrong <- !initially_correct
  wrong_to_correct <- selected & initially_wrong & full_correct
  correct_to_wrong <- selected & initially_correct & !full_correct

  data.frame(
    NBlocks = length(unique(block)),
    NRefined = n_refined,
    SelectionPrecision = if (n_refined > 0) {
      sum(selected & initially_wrong, na.rm = TRUE) / n_refined
    } else {
      NA_real_
    },
    SelectionRecall = if (sum(initially_wrong, na.rm = TRUE) > 0) {
      sum(selected & initially_wrong, na.rm = TRUE) / sum(initially_wrong, na.rm = TRUE)
    } else {
      NA_real_
    },
    WrongToCorrect = sum(wrong_to_correct, na.rm = TRUE),
    CorrectToWrong = sum(correct_to_wrong, na.rm = TRUE),
    CorrectionEfficiency = if (n_refined > 0) {
      (sum(wrong_to_correct, na.rm = TRUE) - sum(correct_to_wrong, na.rm = TRUE)) / n_refined
    } else {
      NA_real_
    }
  )
}

.reliability_block_id <- function(features) {
  fields <- intersect(c("Dataset", "Replicate", "LLMBackend", "LLMModelID"), names(features))
  if (length(fields) == 0) {
    return(rep("all", nrow(features)))
  }
  do.call(paste, c(features[fields], sep = "|"))
}
