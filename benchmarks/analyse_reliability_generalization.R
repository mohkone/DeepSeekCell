# benchmarks/analyse_reliability_generalization.R
#
# Scientific validation analyses for Risk-k generalization. These analyses do
# not introduce a new selector. They ask whether the existing Risk-k model
# transfers across biological tissues, how much training data it needs, how
# close it gets to the FullRefined oracle, and why it fails.
#
# Usage:
#   Rscript benchmarks/analyse_reliability_generalization.R results/benchmark_debug results/reliability_generalization

Sys.setenv(DEEPSEEKCELL_RUN_RELIABILITY_TRAINING_ON_SOURCE = "false")

script_path <- tryCatch(sys.frame(1)$ofile, error = function(e) NA_character_)
if (is.null(script_path) || is.na(script_path) || !nzchar(script_path)) {
  args0 <- commandArgs(trailingOnly = FALSE)
  file_args <- args0[grepl("^--file=", args0)]
  script_path <- if (length(file_args) > 0) sub("^--file=", "", file_args[1]) else ""
}
repo_root <- if (nzchar(script_path)) {
  normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}
setwd(repo_root)

source("R/api.R")
source("R/utils.R")
source("R/ontology.R")
source("R/refinement.R")
source("R/reliability_model.R")
source("benchmarks/train_reliability_model.R")

as_binary <- function(x) {
  if (is.logical(x)) return(x)
  if (is.numeric(x)) return(x > 0)
  tolower(as.character(x)) %in% c("true", "1", "yes", "y")
}

safe_mean <- function(x) {
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

column_or_default <- function(x, column, default = NA_real_) {
  if (!is.data.frame(x) || !column %in% names(x)) {
    return(rep(default, nrow(x)))
  }
  x[[column]]
}

scalar_or_default <- function(x, default = NA_real_) {
  if (is.null(x) || length(x) == 0 || is.na(x[1])) {
    return(default)
  }
  x[1]
}

block_key <- function(features) {
  fields <- intersect(c("Dataset", "Replicate", "LLMBackend", "LLMModelID"), names(features))
  if (length(fields) == 0) {
    return(rep("all", nrow(features)))
  }
  do.call(paste, c(features[fields], sep = "|"))
}

ensure_generalization_features <- function(features) {
  features <- add_training_domain_metadata(features)
  if (!"Tissue" %in% names(features)) {
    features$Tissue <- infer_tissue_from_dataset(features$Dataset)
  }
  if (!"EvidenceAdjustedConfidence" %in% names(features)) {
    features$EvidenceAdjustedConfidence <- column_or_default(features, "Confidence", NA_real_)
  }
  if (!"FullRefinedCorrect" %in% names(features)) {
    features$FullRefinedCorrect <- NA
  }
  features$InitiallyCorrect <- as_binary(features$InitiallyCorrect)
  features$FirstPassIncorrect <- !features$InitiallyCorrect
  features$FullRefinedCorrect <- as_binary(features$FullRefinedCorrect)
  features$BlockID <- block_key(features)
  features
}

selector_indices <- function(block_df,
                             selector,
                             risk = NULL,
                             k = NULL,
                             seed = 1) {
  n <- nrow(block_df)
  if (n == 0) return(integer())

  evidence_flags <- suppressWarnings(
    as.numeric(column_or_default(block_df, "RequiresRefinementNumeric", 0))
  ) > 0
  evidence_flags[is.na(evidence_flags)] <- FALSE
  evidence_k <- sum(evidence_flags)
  if (is.null(k)) {
    k <- if (identical(selector, "FullRefined")) n else evidence_k
  }
  k <- min(max(suppressWarnings(as.integer(k %||% 0)), 0), n)
  if (k == 0 && !identical(selector, "FullRefined")) return(integer())

  if (identical(selector, "Evidence-k")) {
    idx <- which(evidence_flags)
    score <- suppressWarnings(as.numeric(column_or_default(block_df, "EvidenceConflictScore", 0)))
    score[is.na(score)] <- 0
    return(idx[order(score[idx], decreasing = TRUE)])
  }

  if (identical(selector, "Risk-k")) {
    score <- suppressWarnings(as.numeric(risk))
    score[is.na(score)] <- -Inf
    return(head(order(score, decreasing = TRUE), k))
  }

  if (identical(selector, "Confidence-k")) {
    confidence <- as_confidence(column_or_default(
      block_df,
      "LLMConfidence",
      column_or_default(block_df, "Confidence", 0.5)
    ))
    return(head(order(confidence, decreasing = FALSE, na.last = TRUE), k))
  }

  if (identical(selector, "Random-k")) {
    set.seed(seed)
    return(sort(sample(seq_len(n), k)))
  }

  if (identical(selector, "FullRefined")) {
    return(seq_len(n))
  }

  stop("Unknown selector: ", selector, call. = FALSE)
}

evaluate_selector_on_features <- function(features,
                                          selector,
                                          risk = NULL,
                                          seed = 1001) {
  features <- ensure_generalization_features(features)
  selected <- rep(FALSE, nrow(features))

  for (block in unique(features$BlockID)) {
    idx <- which(features$BlockID == block)
    block_df <- features[idx, , drop = FALSE]
    block_risk <- if (!is.null(risk)) risk[idx] else NULL
    local_idx <- selector_indices(
      block_df,
      selector = selector,
      risk = block_risk,
      seed = seed + sum(utf8ToInt(block))
    )
    selected[idx[local_idx]] <- TRUE
  }

  initially_correct <- as_binary(features$InitiallyCorrect)
  full_correct <- as_binary(features$FullRefinedCorrect)
  initially_wrong <- !initially_correct
  wrong_to_correct <- selected & initially_wrong & full_correct
  correct_to_wrong <- selected & initially_correct & !full_correct
  full_wrong_to_correct <- initially_wrong & full_correct
  full_correct_to_wrong <- initially_correct & !full_correct

  n_refined <- sum(selected)
  full_net <- sum(full_wrong_to_correct, na.rm = TRUE) -
    sum(full_correct_to_wrong, na.rm = TRUE)
  net <- sum(wrong_to_correct, na.rm = TRUE) -
    sum(correct_to_wrong, na.rm = TRUE)

  data.frame(
    Selector = selector,
    NBlocks = length(unique(features$BlockID)),
    NClusters = nrow(features),
    NRefined = n_refined,
    NInitiallyIncorrect = sum(initially_wrong, na.rm = TRUE),
    WrongToCorrect = sum(wrong_to_correct, na.rm = TRUE),
    CorrectToWrong = sum(correct_to_wrong, na.rm = TRUE),
    NetCorrections = net,
    FullWrongToCorrect = sum(full_wrong_to_correct, na.rm = TRUE),
    FullCorrectToWrong = sum(full_correct_to_wrong, na.rm = TRUE),
    FullNetCorrections = full_net,
    CorrectionEfficiency = if (n_refined > 0) net / n_refined else NA_real_,
    RecoveryFraction = if (full_net > 0) net / full_net else NA_real_,
    RelativeRefinementBudget = if (nrow(features) > 0) n_refined / nrow(features) else NA_real_,
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
    stringsAsFactors = FALSE
  )
}

train_reliability_safely <- function(features,
                                     target = c("first_pass_error", "refinement_benefit")) {
  target <- match.arg(target)
  y <- if (identical(target, "refinement_benefit") &&
           "RefinementBeneficial" %in% names(features)) {
    as_binary(features$RefinementBeneficial)
  } else {
    as_binary(features$FirstPassIncorrect)
  }
  if (length(unique(y[!is.na(y)])) < 2) {
    return(NULL)
  }
  train_reliability_model(features, target = target)
}

cross_tissue_transfer <- function(features,
                                  target = c("first_pass_error", "refinement_benefit")) {
  target <- match.arg(target)
  features <- ensure_generalization_features(features)
  tissues <- sort(unique(features$Tissue))
  rows <- list()

  for (test_tissue in tissues) {
    train <- features[features$Tissue != test_tissue, , drop = FALSE]
    test <- features[features$Tissue == test_tissue, , drop = FALSE]
    model <- train_reliability_safely(train, target)
    if (is.null(model) || nrow(test) == 0) next
    risk <- predict_reliability_risk(model, test)

    for (selector in c("Risk-k", "Evidence-k", "Confidence-k", "Random-k", "FullRefined")) {
      selector_risk <- if (identical(selector, "Risk-k")) risk else NULL
      summary <- evaluate_selector_on_features(test, selector, risk = selector_risk)
      summary$TransferMode <- "leave_one_tissue_out"
      summary$TrainTissues <- paste(sort(unique(train$Tissue)), collapse = ";")
      summary$TestTissue <- test_tissue
      summary$Target <- target
      rows[[length(rows) + 1]] <- summary
    }
  }

  if (length(rows) == 0) data.frame() else do.call(rbind, rows)
}

pairwise_tissue_transfer <- function(features,
                                     target = c("first_pass_error", "refinement_benefit")) {
  target <- match.arg(target)
  features <- ensure_generalization_features(features)
  tissues <- sort(unique(features$Tissue))
  rows <- list()

  for (train_tissue in tissues) {
    for (test_tissue in setdiff(tissues, train_tissue)) {
      train <- features[features$Tissue == train_tissue, , drop = FALSE]
      test <- features[features$Tissue == test_tissue, , drop = FALSE]
      model <- train_reliability_safely(train, target)
      if (is.null(model) || nrow(test) == 0) next
      risk <- predict_reliability_risk(model, test)
      summary <- evaluate_selector_on_features(test, "Risk-k", risk = risk)
      summary$TransferMode <- "single_tissue_to_tissue"
      summary$TrainTissue <- train_tissue
      summary$TestTissue <- test_tissue
      summary$Target <- target
      rows[[length(rows) + 1]] <- summary
    }
  }

  if (length(rows) == 0) data.frame() else do.call(rbind, rows)
}

learning_curve_analysis <- function(features,
                                    fractions = c(0.1, 0.2, 0.4, 0.6, 0.8, 1),
                                    repeats = 10,
                                    target = c("first_pass_error", "refinement_benefit"),
                                    seed = 2026) {
  target <- match.arg(target)
  features <- ensure_generalization_features(features)
  tissues <- sort(unique(features$Tissue))
  rows <- list()

  for (test_tissue in tissues) {
    train_pool <- features[features$Tissue != test_tissue, , drop = FALSE]
    test <- features[features$Tissue == test_tissue, , drop = FALSE]
    train_blocks <- unique(train_pool$BlockID)
    if (length(train_blocks) == 0 || nrow(test) == 0) next

    for (fraction in fractions) {
      for (repeat_id in seq_len(repeats)) {
        set.seed(seed + repeat_id + round(fraction * 1000) + sum(utf8ToInt(test_tissue)))
        n_blocks <- max(1, ceiling(length(train_blocks) * fraction))
        sampled_blocks <- sample(train_blocks, n_blocks)
        train <- train_pool[train_pool$BlockID %in% sampled_blocks, , drop = FALSE]
        model <- train_reliability_safely(train, target)
        if (is.null(model)) next
        risk <- predict_reliability_risk(model, test)
        summary <- evaluate_selector_on_features(test, "Risk-k", risk = risk)
        summary$TestTissue <- test_tissue
        summary$TrainingFraction <- fraction
        summary$Repeat <- repeat_id
        summary$NTrainingBlocks <- length(sampled_blocks)
        summary$NTrainingRows <- nrow(train)
        summary$Target <- target
        rows[[length(rows) + 1]] <- summary
      }
    }
  }

  if (length(rows) == 0) data.frame() else do.call(rbind, rows)
}

domain_shift_analysis <- function(features, transfer_results) {
  features <- ensure_generalization_features(features)
  tissues <- sort(unique(features$Tissue))
  numeric_features <- intersect(
    c(
      "LLMConfidence", "OntologyEvidenceScore", "MarkerEvidenceScore",
      "BestMarkerEvidenceScore", "TissueEvidenceScore", "ConsensusEvidenceScore",
      "EvidenceConflictScore", "MarkerMargin", "CandidateCount"
    ),
    names(features)
  )

  rows <- lapply(tissues, function(test_tissue) {
    train <- features[features$Tissue != test_tissue, , drop = FALSE]
    test <- features[features$Tissue == test_tissue, , drop = FALSE]
    shifts <- vapply(numeric_features, function(feature) {
      abs(safe_mean(as.numeric(train[[feature]])) - safe_mean(as.numeric(test[[feature]])))
    }, numeric(1))

    transfer_eff <- transfer_results[
      transfer_results$Selector == "Risk-k" &
        transfer_results$TestTissue == test_tissue,
      ,
      drop = FALSE
    ]

    shift_value <- function(name) {
      if (name %in% names(shifts)) unname(shifts[[name]]) else NA_real_
    }
    risk_eff <- if (nrow(transfer_eff) > 0) {
      scalar_or_default(transfer_eff$CorrectionEfficiency)
    } else {
      NA_real_
    }
    risk_recovery <- if (nrow(transfer_eff) > 0) {
      scalar_or_default(transfer_eff$RecoveryFraction)
    } else {
      NA_real_
    }

    data.frame(
      TestTissue = test_tissue,
      NTrain = nrow(train),
      NTest = nrow(test),
      ConfidenceShift = shift_value("LLMConfidence"),
      OntologyEvidenceShift = shift_value("OntologyEvidenceScore"),
      MarkerEvidenceShift = shift_value("MarkerEvidenceScore"),
      BestMarkerEvidenceShift = shift_value("BestMarkerEvidenceScore"),
      ConflictScoreShift = shift_value("EvidenceConflictScore"),
      CandidateCountShift = shift_value("CandidateCount"),
      TrainConflictRate = safe_mean(as.numeric(column_or_default(train, "RequiresRefinementNumeric", 0) > 0)),
      TestConflictRate = safe_mean(as.numeric(column_or_default(test, "RequiresRefinementNumeric", 0) > 0)),
      PredictionEntropyTrain = prediction_entropy(column_or_default(
        train, "RawPrediction", column_or_default(train, "Cluster", "")
      )),
      PredictionEntropyTest = prediction_entropy(column_or_default(
        test, "RawPrediction", column_or_default(test, "Cluster", "")
      )),
      PredictionJSD = prediction_jsd(
        column_or_default(train, "RawPrediction", column_or_default(train, "Cluster", "")),
        column_or_default(test, "RawPrediction", column_or_default(test, "Cluster", ""))
      ),
      RiskCorrectionEfficiency = risk_eff,
      RiskRecoveryFraction = risk_recovery,
      stringsAsFactors = FALSE
    )
  })

  shift <- do.call(rbind, rows)
  numeric_shift <- setdiff(names(shift), c("TestTissue"))
  correlations <- do.call(rbind, lapply(numeric_shift, function(column) {
    if (all(is.na(shift[[column]])) || all(is.na(shift$RiskCorrectionEfficiency))) {
      return(NULL)
    }
    data.frame(
      ShiftMetric = column,
      CorrelationWithRiskEfficiency = suppressWarnings(
        stats::cor(shift[[column]], shift$RiskCorrectionEfficiency, use = "complete.obs")
      ),
      stringsAsFactors = FALSE
    )
  }))
  if (is.null(correlations)) {
    correlations <- data.frame()
  }

  list(shift = shift, correlations = correlations)
}

prediction_entropy <- function(labels) {
  labels <- labels[!is.na(labels) & nzchar(as.character(labels))]
  if (length(labels) == 0) return(NA_real_)
  p <- prop.table(table(labels))
  -sum(p * log2(p))
}

prediction_jsd <- function(a, b) {
  a <- a[!is.na(a) & nzchar(as.character(a))]
  b <- b[!is.na(b) & nzchar(as.character(b))]
  labels <- union(unique(a), unique(b))
  if (length(labels) == 0) return(NA_real_)
  pa <- tabulate(match(a, labels), nbins = length(labels))
  pb <- tabulate(match(b, labels), nbins = length(labels))
  pa <- pa / sum(pa)
  pb <- pb / sum(pb)
  m <- (pa + pb) / 2
  kl <- function(p, q) {
    idx <- p > 0
    sum(p[idx] * log2(p[idx] / q[idx]))
  }
  (kl(pa, m) + kl(pb, m)) / 2
}

calibration_comparison <- function(features,
                                   target = c("first_pass_error", "refinement_benefit")) {
  target <- match.arg(target)
  features <- ensure_generalization_features(features)
  rows <- list()
  predictions <- list()

  for (test_tissue in sort(unique(features$Tissue))) {
    train <- features[features$Tissue != test_tissue, , drop = FALSE]
    test <- features[features$Tissue == test_tissue, , drop = FALSE]
    model <- train_reliability_safely(train, target)
    if (is.null(model) || nrow(test) == 0) next

    observed <- as.numeric(as_binary(test$FirstPassIncorrect))
    scores <- list(
      `Raw LLM risk` = 1 - as_confidence(column_or_default(test, "LLMConfidence", 0.5)),
      `Evidence-adjusted risk` = 1 - as_confidence(column_or_default(
        test,
        "EvidenceAdjustedConfidence",
        column_or_default(test, "Confidence", 0.5)
      )),
      `Risk-k probability` = predict_reliability_risk(model, test)
    )

    for (score_name in names(scores)) {
      score <- pmin(pmax(scores[[score_name]], 0), 1)
      bins <- calibration_bins(observed, score, n_bins = 5)
      bins$ScoreType <- score_name
      bins$TestTissue <- test_tissue
      rows[[length(rows) + 1]] <- bins
      predictions[[length(predictions) + 1]] <- data.frame(
        TestTissue = test_tissue,
        ScoreType = score_name,
        Observed = observed,
        PredictedRisk = score,
        Brier = mean((score - observed)^2, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }

  bins <- if (length(rows) == 0) data.frame() else do.call(rbind, rows)
  pred <- if (length(predictions) == 0) data.frame() else do.call(rbind, predictions)
  metrics <- if (nrow(pred) == 0) {
    data.frame()
  } else {
    stats::aggregate(
      cbind(Brier = pred$Brier),
      by = list(ScoreType = pred$ScoreType),
      FUN = mean,
      na.rm = TRUE
    )
  }

  list(bins = bins, metrics = metrics, predictions = pred)
}

calibration_bins <- function(observed, predicted, n_bins = 5) {
  ok <- !is.na(observed) & !is.na(predicted)
  observed <- observed[ok]
  predicted <- predicted[ok]
  if (length(observed) == 0) return(data.frame())
  ord <- order(predicted)
  bin <- integer(length(predicted))
  bin[ord] <- ceiling(seq_along(predicted) / length(predicted) * n_bins)
  bin <- pmax(pmin(bin, n_bins), 1)

  do.call(rbind, lapply(sort(unique(bin)), function(b) {
    idx <- bin == b
    data.frame(
      Bin = b,
      N = sum(idx),
      MeanPredictedRisk = mean(predicted[idx]),
      ObservedEventRate = mean(observed[idx]),
      Brier = mean((predicted[idx] - observed[idx])^2),
      stringsAsFactors = FALSE
    )
  }))
}

failure_taxonomy <- function(features) {
  features <- ensure_generalization_features(features)
  incorrect <- features[as_binary(features$FirstPassIncorrect), , drop = FALSE]
  if (nrow(incorrect) == 0) {
    return(list(cluster = data.frame(), summary = data.frame()))
  }

  category <- vapply(seq_len(nrow(incorrect)), function(i) {
    row <- incorrect[i, , drop = FALSE]
    row_num <- function(column, default) {
      value <- suppressWarnings(as.numeric(column_or_default(row, column, default)))
      ifelse(is.na(value), default, value)
    }
    if (row_num("IsMixedNumeric", 0) > 0) return("mixed_cell_cluster")
    if (row_num("CandidateEvidenceDisagreement", 0) > 0) return("candidate_evidence_disagreement")
    if (row_num("RequiresRefinementNumeric", 0) > 0 ||
        row_num("EvidenceConflictScore", 0) >= 0.25) return("evidence_conflict")
    if (row_num("MarkerEvidenceScore", 0) < 0.2 &&
        row_num("BestMarkerEvidenceScore", 0) < 0.2) return("weak_marker_evidence")
    if (row_num("OntologyEvidenceScore", 1) < 0.5) return("ontology_ambiguity")
    if (row_num("LLMConfidence", 1) < 0.5) return("low_llm_confidence")
    "biologically_related_or_other"
  }, character(1))

  cluster <- data.frame(
    incorrect[intersect(
      c("Dataset", "Tissue", "Replicate", "LLMBackend", "Cluster",
        "RawPrediction", "RawTruth", "HarmonisedPrediction", "HarmonisedTruth"),
      names(incorrect)
    )],
    FailureCategory = category,
    stringsAsFactors = FALSE
  )
  summary <- as.data.frame(table(FailureCategory = category), stringsAsFactors = FALSE)
  summary$Fraction <- summary$Freq / sum(summary$Freq)
  summary <- summary[order(summary$Freq, decreasing = TRUE), , drop = FALSE]
  rownames(summary) <- NULL

  list(cluster = cluster, summary = summary)
}

plot_generalization_outputs <- function(prefix,
                                        transfer,
                                        pairwise,
                                        learning,
                                        oracle_gap,
                                        calibration_bins,
                                        taxonomy_summary) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    return(invisible(FALSE))
  }

  if (nrow(pairwise) > 0) {
    p <- ggplot2::ggplot(
      pairwise,
      ggplot2::aes(
        x = .data$TrainTissue,
        y = .data$TestTissue,
        fill = .data$CorrectionEfficiency
      )
    ) +
      ggplot2::geom_tile(color = "white") +
      ggplot2::scale_fill_gradient(low = "#F7FBFF", high = "#2C7FB8", na.value = "grey90") +
      ggplot2::labs(x = "Train tissue", y = "Test tissue", fill = "Efficiency") +
      ggplot2::theme_minimal(base_size = 11)
    ggplot2::ggsave(paste0(prefix, "_cross_tissue_transfer_matrix.pdf"), p, width = 6, height = 5)
  }

  if (nrow(learning) > 0) {
    p <- ggplot2::ggplot(
      learning,
      ggplot2::aes(
        x = .data$TrainingFraction,
        y = .data$CorrectionEfficiency,
        color = .data$TestTissue
      )
    ) +
      ggplot2::stat_summary(fun = mean, geom = "line") +
      ggplot2::stat_summary(fun = mean, geom = "point") +
      ggplot2::labs(x = "Training fraction", y = "Correction efficiency", color = "Held-out tissue") +
      ggplot2::theme_minimal(base_size = 11)
    ggplot2::ggsave(paste0(prefix, "_learning_curve.pdf"), p, width = 6.5, height = 4.5)
  }

  if (nrow(oracle_gap) > 0) {
    p <- ggplot2::ggplot(
      oracle_gap,
      ggplot2::aes(x = .data$Selector, y = .data$RecoveryFraction, fill = .data$Selector)
    ) +
      ggplot2::geom_col(show.legend = FALSE) +
      ggplot2::facet_wrap(~TestTissue) +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::labs(x = NULL, y = "Recovery fraction vs FullRefined oracle") +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
    ggplot2::ggsave(paste0(prefix, "_oracle_gap.pdf"), p, width = 8, height = 5)
  }

  if (nrow(calibration_bins) > 0) {
    p <- ggplot2::ggplot(
      calibration_bins,
      ggplot2::aes(
        x = .data$MeanPredictedRisk,
        y = .data$ObservedEventRate,
        color = .data$ScoreType
      )
    ) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
      ggplot2::geom_point(ggplot2::aes(size = .data$N)) +
      ggplot2::geom_line() +
      ggplot2::facet_wrap(~TestTissue) +
      ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
      ggplot2::labs(x = "Predicted risk", y = "Observed error rate", color = NULL) +
      ggplot2::theme_minimal(base_size = 11)
    ggplot2::ggsave(paste0(prefix, "_calibration_comparison.pdf"), p, width = 8, height = 6)
  }

  if (nrow(taxonomy_summary) > 0) {
    p <- ggplot2::ggplot(
      taxonomy_summary,
      ggplot2::aes(
        x = stats::reorder(.data$FailureCategory, .data$Fraction),
        y = .data$Fraction
      )
    ) +
      ggplot2::geom_col(fill = "#2C7FB8") +
      ggplot2::coord_flip() +
      ggplot2::labs(x = NULL, y = "Fraction of first-pass errors") +
      ggplot2::theme_minimal(base_size = 11)
    ggplot2::ggsave(paste0(prefix, "_failure_taxonomy.pdf"), p, width = 7, height = 4.5)
  }

  invisible(TRUE)
}

run_generalization_analysis <- function(debug_dir = file.path("results", "benchmark_debug"),
                                        output_prefix = file.path("results", "reliability_generalization"),
                                        target = c("first_pass_error", "refinement_benefit"),
                                        learning_repeats = 10) {
  target <- match.arg(target)
  features <- ensure_generalization_features(build_training_features(debug_dir))
  dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)

  transfer <- cross_tissue_transfer(features, target = target)
  pairwise <- pairwise_tissue_transfer(features, target = target)
  learning <- learning_curve_analysis(features, repeats = learning_repeats, target = target)
  shift <- domain_shift_analysis(features, transfer)
  calibration <- calibration_comparison(features, target = target)
  taxonomy <- failure_taxonomy(features)
  oracle_gap <- transfer[transfer$Selector %in% c(
    "Risk-k", "Evidence-k", "Confidence-k", "Random-k", "FullRefined"
  ), , drop = FALSE]

  utils::write.csv(features, paste0(output_prefix, "_features.csv"), row.names = FALSE)
  utils::write.csv(transfer, paste0(output_prefix, "_leave_one_tissue_out.csv"), row.names = FALSE)
  utils::write.csv(pairwise, paste0(output_prefix, "_pairwise_tissue_transfer.csv"), row.names = FALSE)
  utils::write.csv(learning, paste0(output_prefix, "_learning_curve.csv"), row.names = FALSE)
  utils::write.csv(shift$shift, paste0(output_prefix, "_domain_shift.csv"), row.names = FALSE)
  utils::write.csv(shift$correlations, paste0(output_prefix, "_domain_shift_correlations.csv"), row.names = FALSE)
  utils::write.csv(calibration$bins, paste0(output_prefix, "_calibration_comparison_bins.csv"), row.names = FALSE)
  utils::write.csv(calibration$metrics, paste0(output_prefix, "_calibration_comparison_metrics.csv"), row.names = FALSE)
  utils::write.csv(taxonomy$cluster, paste0(output_prefix, "_failure_taxonomy_clusters.csv"), row.names = FALSE)
  utils::write.csv(taxonomy$summary, paste0(output_prefix, "_failure_taxonomy_summary.csv"), row.names = FALSE)
  utils::write.csv(oracle_gap, paste0(output_prefix, "_oracle_gap.csv"), row.names = FALSE)

  plot_generalization_outputs(
    output_prefix,
    transfer = transfer,
    pairwise = pairwise,
    learning = learning,
    oracle_gap = oracle_gap,
    calibration_bins = calibration$bins,
    taxonomy_summary = taxonomy$summary
  )

  invisible(list(
    features = features,
    transfer = transfer,
    pairwise = pairwise,
    learning = learning,
    domain_shift = shift,
    calibration = calibration,
    taxonomy = taxonomy,
    oracle_gap = oracle_gap
  ))
}

should_auto_run_generalization_analysis <- function() {
  value <- tolower(Sys.getenv("DEEPSEEKCELL_RUN_GENERALIZATION_ON_SOURCE", unset = "true"))
  value %in% c("1", "true", "yes", "y")
}

if (should_auto_run_generalization_analysis()) {
  args <- commandArgs(trailingOnly = TRUE)
  debug_dir <- if (length(args) >= 1) args[1] else file.path("results", "benchmark_debug")
  output_prefix <- if (length(args) >= 2) args[2] else file.path("results", "reliability_generalization")
  target <- if (length(args) >= 3) args[3] else "first_pass_error"
  repeats <- if (length(args) >= 4) as.integer(args[4]) else 10
  run_generalization_analysis(debug_dir, output_prefix, target, repeats)
}
