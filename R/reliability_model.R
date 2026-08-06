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
                                    target_column = NULL) {
  target <- match.arg(target)
  if (!is.data.frame(features) || nrow(features) == 0) {
    stop("features must be a non-empty data frame.", call. = FALSE)
  }

  y <- .resolve_reliability_target(features, target, target_column)
  x <- .prepare_reliability_feature_frame(features)
  keep <- !is.na(y)
  y <- as.integer(y[keep])
  x <- x[keep, , drop = FALSE]

  if (length(y) == 0) {
    stop("No non-missing target values are available.", call. = FALSE)
  }

  fallback_rate <- mean(y)
  fit <- NULL
  formula <- stats::as.formula(paste("Target ~", paste(RELIABILITY_MODEL_FEATURES, collapse = " + ")))

  if (length(unique(y)) > 1 && nrow(x) > length(RELIABILITY_MODEL_FEATURES)) {
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
    feature_columns = RELIABILITY_MODEL_FEATURES,
    feature_medians = attr(x, "feature_medians"),
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

  x <- .prepare_reliability_feature_frame(
    features,
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

.prepare_reliability_feature_frame <- function(features, feature_medians = NULL) {
  missing <- setdiff(RELIABILITY_MODEL_FEATURES, names(features))
  for (column in missing) {
    features[[column]] <- NA_real_
  }

  x <- features[RELIABILITY_MODEL_FEATURES]
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
