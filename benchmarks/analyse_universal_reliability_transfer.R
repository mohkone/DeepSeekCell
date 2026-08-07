# benchmarks/analyse_universal_reliability_transfer.R
#
# Universal Reliability Benchmark for frozen Risk-k models. The analysis asks
# whether a single frozen reliability model can identify risky annotations
# across annotation backends/models without retraining or retuning.
#
# Usage:
#   Rscript benchmarks/analyse_universal_reliability_transfer.R \
#     results/benchmark_debug \
#     results/reliability_model_v1.1_error.rds \
#     results/universal_reliability_transfer \
#     deepseek-chat

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

source("R/utils.R")
source("R/ontology.R")
source("R/refinement.R")
source("R/reliability_model.R")

as_number <- function(x) {
  suppressWarnings(as.numeric(x))
}

as_flag_universal <- function(x) {
  if (is.logical(x)) return(x)
  if (is.numeric(x)) return(x != 0)
  tolower(trimws(as.character(x))) %in% c("true", "1", "yes", "y")
}

safe_divide <- function(x, y) {
  out <- rep(NA_real_, length(x))
  ok <- !is.na(x) & !is.na(y) & y != 0
  out[ok] <- x[ok] / y[ok]
  out
}

ensure_column <- function(x, column, default = NA) {
  if (!column %in% names(x)) {
    x[[column]] <- default
  }
  x
}

read_debug_arm <- function(debug_dir, pattern) {
  files <- list.files(debug_dir, pattern = pattern, full.names = TRUE, recursive = FALSE)
  if (length(files) == 0) {
    return(data.frame())
  }
  rows <- lapply(files, function(path) {
    x <- tryCatch(utils::read.csv(path, stringsAsFactors = FALSE), error = function(e) data.frame())
    if (nrow(x) == 0) return(data.frame())
    x$SourceFile <- basename(path)
    fill_model_identity(x)
  })
  rows <- rows[vapply(rows, nrow, integer(1)) > 0]
  if (length(rows) == 0) {
    return(data.frame())
  }
  all_names <- unique(unlist(lapply(rows, names)))
  rows <- lapply(rows, function(x) {
    missing <- setdiff(all_names, names(x))
    for (column in missing) x[[column]] <- NA
    x[all_names]
  })
  do.call(rbind, rows)
}

fill_model_identity <- function(x) {
  x <- ensure_column(x, "Method", "")
  x <- ensure_column(x, "LLMBackend", NA_character_)
  x <- ensure_column(x, "LLMModelID", NA_character_)

  method <- as.character(x$Method)
  source <- as.character(x$SourceFile %||% "")
  backend <- as.character(x$LLMBackend)
  model <- as.character(x$LLMModelID)

  missing_backend <- is.na(backend) | !nzchar(backend)
  backend[missing_backend & grepl("DeepSeek", method, ignore.case = TRUE)] <- "deepseek"
  backend[missing_backend & grepl("Ollama", method, ignore.case = TRUE)] <- "ollama"
  backend[missing_backend & grepl("DeepSeek", source, ignore.case = TRUE)] <- "deepseek"
  backend[missing_backend & grepl("Ollama", source, ignore.case = TRUE)] <- "ollama"
  backend[is.na(backend) | !nzchar(backend)] <- "unknown"

  missing_model <- is.na(model) | !nzchar(model)
  model[missing_model & backend == "deepseek"] <- "deepseek-chat"
  model[missing_model & grepl("llama3_2|llama3\\.2", method, ignore.case = TRUE)] <- "llama3.2:latest"
  model[missing_model & grepl("mistral", method, ignore.case = TRUE)] <- "mistral:latest"
  model[missing_model & grepl("gemma2_2b|gemma2:2b", method, ignore.case = TRUE)] <- "gemma2:2b"
  model[missing_model & backend == "ollama"] <- "ollama-default"
  model[is.na(model) | !nzchar(model)] <- method[is.na(model) | !nzchar(model)]
  model[is.na(model) | !nzchar(model)] <- "unknown"

  x$LLMBackend <- backend
  x$LLMModelID <- model
  x
}

block_key <- function(x) {
  for (column in c("Dataset", "Replicate", "LLMBackend", "LLMModelID")) {
    x <- ensure_column(x, column, "")
  }
  paste(x$Dataset, x$Replicate, x$LLMBackend, x$LLMModelID, sep = "|")
}

row_key <- function(x) {
  for (column in c("Dataset", "Replicate", "LLMBackend", "LLMModelID", "Cluster")) {
    x <- ensure_column(x, column, "")
  }
  paste(x$Dataset, x$Replicate, x$LLMBackend, x$LLMModelID, x$Cluster, sep = "|")
}

prefer_calibrated_then_evidence <- function(calibrated, evidence) {
  if (nrow(calibrated) == 0) return(evidence)
  if (nrow(evidence) == 0) return(calibrated)
  cal_key <- row_key(calibrated)
  ev_key <- row_key(evidence)
  rbind(calibrated, evidence[!ev_key %in% cal_key, names(calibrated), drop = FALSE])
}

build_universal_features <- function(debug_dir) {
  calibrated <- read_debug_arm(debug_dir, "Calibrated_debug\\.csv$")
  evidence <- read_debug_arm(debug_dir, "Evidence_debug\\.csv$")
  full <- read_debug_arm(debug_dir, "FullRefined_debug\\.csv$")
  base <- prefer_calibrated_then_evidence(calibrated, evidence)
  if (nrow(base) == 0) {
    stop("No Calibrated/Evidence debug files found in ", debug_dir, call. = FALSE)
  }

  base <- ensure_column(base, "Cluster", seq_len(nrow(base)))
  base <- ensure_column(base, "HarmonisedPrediction", NA_character_)
  base <- ensure_column(base, "HarmonisedTruth", NA_character_)
  base <- ensure_column(base, "RawPrediction", base$HarmonisedPrediction)
  base <- ensure_column(base, "RawTruth", base$HarmonisedTruth)
  base <- ensure_column(base, "Confidence", base$LLMConfidence %||% 0.5)
  base <- ensure_column(base, "LLMConfidence", base$Confidence)
  base <- ensure_column(base, "RequiresRefinement", FALSE)
  base <- ensure_column(base, "EvidenceConflict", FALSE)
  base <- ensure_column(base, "CandidateCellTypes", "")
  base <- ensure_column(base, "IsMixed", FALSE)
  base$CellType <- ifelse(
    !is.na(base$HarmonisedPrediction) & nzchar(as.character(base$HarmonisedPrediction)),
    as.character(base$HarmonisedPrediction),
    as.character(base$RawPrediction)
  )

  base$InitiallyCorrect <- normalize_cell_type(base$HarmonisedPrediction) ==
    normalize_cell_type(base$HarmonisedTruth)
  base$FirstPassIncorrect <- !base$InitiallyCorrect
  base$BlockID <- block_key(base)
  base$RowID <- row_key(base)

  if (nrow(full) > 0) {
    full <- ensure_column(full, "HarmonisedPrediction", NA_character_)
    full <- ensure_column(full, "HarmonisedTruth", NA_character_)
    full$RowID <- row_key(full)
    idx <- match(base$RowID, full$RowID)
    full_pred <- full$HarmonisedPrediction[idx]
    full_truth <- full$HarmonisedTruth[idx]
    full_correct <- normalize_cell_type(full_pred) == normalize_cell_type(full_truth)
    full_correct[is.na(idx)] <- NA
    base$FullRefinedCorrect <- full_correct
    base$HasFullRefinedOutcome <- !is.na(idx)
    base$RefinementBeneficial <- !base$InitiallyCorrect & full_correct
    base$RefinementHarmful <- base$InitiallyCorrect & !full_correct
  } else {
    base$FullRefinedCorrect <- NA
    base$HasFullRefinedOutcome <- FALSE
    base$RefinementBeneficial <- NA
    base$RefinementHarmful <- NA
  }

  feature_df <- extract_reliability_features(base)
  metadata_cols <- intersect(
    c(
      "Dataset", "Replicate", "Cluster", "Method", "LLMBackend", "LLMModelID",
      "Tissue", "Species", "ClusterPurity", "RefinementSelector", "RefinementBudgetK",
      "RawPrediction", "RawTruth", "HarmonisedPrediction", "HarmonisedTruth",
      "InitiallyCorrect", "FirstPassIncorrect", "FullRefinedCorrect",
      "HasFullRefinedOutcome", "RefinementBeneficial", "RefinementHarmful",
      "FirstPassResponseHash", "SourceFile", "BlockID", "RowID"
    ),
    names(base)
  )
  out <- cbind(base[metadata_cols], feature_df[setdiff(names(feature_df), "Cluster")])
  out$RiskTargetFirstPassError <- as.integer(as_flag_universal(out$FirstPassIncorrect))
  out
}

rank_evidence <- function(block) {
  requires <- as.numeric(as_flag_universal(block$RequiresRefinementNumeric %||% block$RequiresRefinement %||% FALSE))
  conflict <- as_number(block$EvidenceConflictScore %||% 0)
  best <- as_number(block$BestMarkerEvidenceScore %||% 0)
  margin <- as_number(block$MarkerMargin %||% 0)
  requires[is.na(requires)] <- 0
  conflict[is.na(conflict)] <- 0
  best[is.na(best)] <- 0
  margin[is.na(margin)] <- 0
  10 * requires + conflict + 0.1 * best + 0.05 * margin
}

resolve_budget <- function(n, explicit_k = NULL, fraction = 0.2, min_k = 1) {
  if (!is.null(explicit_k) && !is.na(explicit_k) && explicit_k > 0) {
    return(min(n, as.integer(explicit_k)))
  }
  min(n, max(min_k, ceiling(n * fraction)))
}

select_by_strategy <- function(block, strategy, risk = NULL, k, seed = 1001) {
  n <- nrow(block)
  if (n == 0 || k <= 0) return(integer())
  k <- min(k, n)
  if (identical(strategy, "Risk-k")) {
    score <- as_number(risk)
    score[is.na(score)] <- -Inf
    return(head(order(score, decreasing = TRUE), k))
  }
  if (identical(strategy, "Evidence-k")) {
    return(head(order(rank_evidence(block), decreasing = TRUE), k))
  }
  if (identical(strategy, "Confidence-k")) {
    confidence <- as_confidence(block$LLMConfidence %||% block$Confidence %||% 0.5)
    return(head(order(confidence, decreasing = FALSE, na.last = TRUE), k))
  }
  if (identical(strategy, "Random-k")) {
    set.seed(seed)
    return(sort(sample(seq_len(n), k)))
  }
  stop("Unknown universal selector strategy: ", strategy, call. = FALSE)
}

auc_binary <- function(observed, score) {
  observed <- as.integer(observed)
  ok <- !is.na(observed) & !is.na(score)
  observed <- observed[ok]
  score <- score[ok]
  if (length(unique(observed)) < 2) return(NA_real_)
  pos <- score[observed == 1]
  neg <- score[observed == 0]
  ranks <- rank(c(pos, neg), ties.method = "average")
  n_pos <- length(pos)
  n_neg <- length(neg)
  (sum(ranks[seq_len(n_pos)]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

binary_nll <- function(observed, predicted, eps = 1e-6) {
  observed <- as.integer(observed)
  predicted <- pmin(pmax(predicted, eps), 1 - eps)
  -mean(observed * log(predicted) + (1 - observed) * log(1 - predicted), na.rm = TRUE)
}

ece_equal_frequency <- function(observed, predicted, n_bins = 10) {
  observed <- as.integer(observed)
  ok <- !is.na(observed) & !is.na(predicted)
  observed <- observed[ok]
  predicted <- predicted[ok]
  if (length(observed) == 0) return(NA_real_)
  ord <- order(predicted)
  bin <- integer(length(predicted))
  bin[ord] <- ceiling(seq_along(predicted) / length(predicted) * n_bins)
  bin <- pmax(pmin(bin, n_bins), 1)
  sum(vapply(sort(unique(bin)), function(b) {
    idx <- bin == b
    (sum(idx) / length(observed)) * abs(mean(observed[idx]) - mean(predicted[idx]))
  }, numeric(1)))
}

calibration_bins <- function(features, n_bins = 10) {
  rows <- list()
  for (model_id in sort(unique(features$LLMModelID))) {
    x <- features[features$LLMModelID == model_id, , drop = FALSE]
    observed <- as.integer(as_flag_universal(x$FirstPassIncorrect))
    predicted <- x$FrozenRisk
    ok <- !is.na(observed) & !is.na(predicted)
    observed <- observed[ok]
    predicted <- predicted[ok]
    if (length(observed) == 0) next
    ord <- order(predicted)
    bin <- integer(length(predicted))
    bin[ord] <- ceiling(seq_along(predicted) / length(predicted) * n_bins)
    bin <- pmax(pmin(bin, n_bins), 1)
    rows[[length(rows) + 1]] <- do.call(rbind, lapply(sort(unique(bin)), function(b) {
      idx <- bin == b
      data.frame(
        LLMModelID = model_id,
        Bin = b,
        N = sum(idx),
        MeanPredictedRisk = mean(predicted[idx]),
        ObservedErrorRate = mean(observed[idx]),
        MinPredictedRisk = min(predicted[idx]),
        MaxPredictedRisk = max(predicted[idx]),
        stringsAsFactors = FALSE
      )
    }))
  }
  if (length(rows) == 0) data.frame() else do.call(rbind, rows)
}

calibration_slope_intercept <- function(observed, predicted) {
  observed <- as.integer(observed)
  ok <- !is.na(observed) & !is.na(predicted)
  observed <- observed[ok]
  predicted <- predicted[ok]
  if (length(unique(observed)) < 2 || length(unique(predicted)) < 2) {
    return(c(intercept = NA_real_, slope = NA_real_))
  }
  fit <- tryCatch(
    stats::glm(observed ~ stats::qlogis(pmin(pmax(predicted, 1e-6), 1 - 1e-6)), family = stats::binomial()),
    error = function(e) NULL
  )
  if (is.null(fit)) return(c(intercept = NA_real_, slope = NA_real_))
  coef <- stats::coef(fit)
  c(intercept = unname(coef[[1]]), slope = unname(coef[[2]]))
}

summarise_transfer_metrics <- function(features, train_model_id) {
  rows <- lapply(sort(unique(features$LLMModelID)), function(model_id) {
    x <- features[features$LLMModelID == model_id, , drop = FALSE]
    observed <- as.integer(as_flag_universal(x$FirstPassIncorrect))
    predicted <- x$FrozenRisk
    slope <- calibration_slope_intercept(observed, predicted)
    data.frame(
      TrainModelID = train_model_id,
      TestModelID = model_id,
      NBlocks = length(unique(x$BlockID)),
      NClusters = nrow(x),
      EventRate = mean(observed, na.rm = TRUE),
      MeanPredictedRisk = mean(predicted, na.rm = TRUE),
      AUROC = auc_binary(observed, predicted),
      Brier = mean((predicted - observed)^2, na.rm = TRUE),
      BinaryNLL = binary_nll(observed, predicted),
      ECE = ece_equal_frequency(observed, predicted),
      CalibrationIntercept = slope[["intercept"]],
      CalibrationSlope = slope[["slope"]],
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  source <- out[out$TestModelID == train_model_id, , drop = FALSE]
  if (nrow(source) > 0) {
    out$DeltaAUROCFromSource <- out$AUROC - source$AUROC[[1]]
    out$DeltaBrierFromSource <- out$Brier - source$Brier[[1]]
    out$DeltaECEFromSource <- out$ECE - source$ECE[[1]]
  } else {
    out$DeltaAUROCFromSource <- NA_real_
    out$DeltaBrierFromSource <- NA_real_
    out$DeltaECEFromSource <- NA_real_
  }
  out
}

evaluate_selector_transfer <- function(features,
                                       budget_k = NULL,
                                       budget_fraction = 0.2,
                                       budget_min = 1,
                                       seed = 1001) {
  strategies <- c("Risk-k", "Evidence-k", "Confidence-k", "Random-k")
  rows <- list()
  selections <- list()
  for (block in unique(features$BlockID)) {
    idx <- which(features$BlockID == block)
    block_df <- features[idx, , drop = FALSE]
    k <- resolve_budget(nrow(block_df), explicit_k = budget_k, fraction = budget_fraction, min_k = budget_min)
    for (strategy in strategies) {
      local <- select_by_strategy(
        block_df,
        strategy = strategy,
        risk = block_df$FrozenRisk,
        k = k,
        seed = seed + sum(utf8ToInt(block)) + match(strategy, strategies)
      )
      selected <- rep(FALSE, nrow(block_df))
      selected[local] <- TRUE
      initial_correct <- as_flag_universal(block_df$InitiallyCorrect)
      initial_wrong <- !initial_correct
      has_full <- as_flag_universal(block_df$HasFullRefinedOutcome)
      full_correct <- as_flag_universal(block_df$FullRefinedCorrect)
      wrong_to_correct <- selected & has_full & initial_wrong & full_correct
      correct_to_wrong <- selected & has_full & initial_correct & !full_correct
      full_wrong_to_correct <- has_full & initial_wrong & full_correct
      full_correct_to_wrong <- has_full & initial_correct & !full_correct
      n_refined <- sum(selected)
      net <- sum(wrong_to_correct, na.rm = TRUE) - sum(correct_to_wrong, na.rm = TRUE)
      full_net <- sum(full_wrong_to_correct, na.rm = TRUE) - sum(full_correct_to_wrong, na.rm = TRUE)

      rows[[length(rows) + 1]] <- data.frame(
        Dataset = block_df$Dataset[[1]],
        Replicate = block_df$Replicate[[1]],
        LLMBackend = block_df$LLMBackend[[1]],
        LLMModelID = block_df$LLMModelID[[1]],
        BlockID = block,
        Selector = strategy,
        BudgetK = k,
        NClusters = nrow(block_df),
        NRefined = n_refined,
        NInitiallyIncorrect = sum(initial_wrong, na.rm = TRUE),
        SelectionPrecision = if (n_refined > 0) sum(selected & initial_wrong, na.rm = TRUE) / n_refined else NA_real_,
        SelectionRecall = if (sum(initial_wrong, na.rm = TRUE) > 0) {
          sum(selected & initial_wrong, na.rm = TRUE) / sum(initial_wrong, na.rm = TRUE)
        } else {
          NA_real_
        },
        HasFullRefinedOutcome = any(has_full, na.rm = TRUE),
        WrongToCorrect = sum(wrong_to_correct, na.rm = TRUE),
        CorrectToWrong = sum(correct_to_wrong, na.rm = TRUE),
        NetCorrections = net,
        FullWrongToCorrect = sum(full_wrong_to_correct, na.rm = TRUE),
        FullCorrectToWrong = sum(full_correct_to_wrong, na.rm = TRUE),
        FullNetCorrections = full_net,
        CorrectionEfficiency = if (n_refined > 0 && any(has_full)) net / n_refined else NA_real_,
        WrongToCorrectEfficiency = if (n_refined > 0 && any(has_full)) sum(wrong_to_correct, na.rm = TRUE) / n_refined else NA_real_,
        RecoveryFraction = if (full_net > 0) net / full_net else NA_real_,
        RelativeRefinementBudget = if (nrow(block_df) > 0) n_refined / nrow(block_df) else NA_real_,
        stringsAsFactors = FALSE
      )

      if (length(local) > 0) {
        selections[[length(selections) + 1]] <- data.frame(
          Dataset = block_df$Dataset[local],
          Replicate = block_df$Replicate[local],
          LLMBackend = block_df$LLMBackend[local],
          LLMModelID = block_df$LLMModelID[local],
          Selector = strategy,
          Cluster = as.character(block_df$Cluster[local]),
          FrozenRisk = block_df$FrozenRisk[local],
          FirstPassIncorrect = block_df$FirstPassIncorrect[local],
          RefinementBeneficial = block_df$RefinementBeneficial[local],
          stringsAsFactors = FALSE
        )
      }
    }
  }

  behavior <- if (length(rows) == 0) data.frame() else do.call(rbind, rows)
  selected <- if (length(selections) == 0) data.frame() else do.call(rbind, selections)
  list(behavior = behavior, selected = selected)
}

summarise_selector_behavior <- function(behavior) {
  if (nrow(behavior) == 0) return(data.frame())
  rows <- lapply(split(behavior, paste(behavior$LLMModelID, behavior$Selector, sep = "\r")), function(x) {
    n_refined <- sum(x$NRefined, na.rm = TRUE)
    wrong_to_correct <- sum(x$WrongToCorrect, na.rm = TRUE)
    correct_to_wrong <- sum(x$CorrectToWrong, na.rm = TRUE)
    full_net <- sum(x$FullNetCorrections, na.rm = TRUE)
    net <- wrong_to_correct - correct_to_wrong
    data.frame(
      LLMModelID = x$LLMModelID[[1]],
      Selector = x$Selector[[1]],
      NBlocks = length(unique(x$BlockID)),
      NClusters = sum(x$NClusters, na.rm = TRUE),
      NRefined = n_refined,
      WrongToCorrect = wrong_to_correct,
      CorrectToWrong = correct_to_wrong,
      NetCorrections = net,
      CorrectionEfficiency = if (n_refined > 0) net / n_refined else NA_real_,
      WrongToCorrectEfficiency = if (n_refined > 0) wrong_to_correct / n_refined else NA_real_,
      RecoveryFraction = if (full_net > 0) net / full_net else NA_real_,
      MeanSelectionPrecision = mean(x$SelectionPrecision, na.rm = TRUE),
      MeanSelectionRecall = mean(x$SelectionRecall, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

selector_overlap <- function(selected) {
  if (nrow(selected) == 0) return(data.frame())
  groups <- unique(selected[c("Dataset", "Replicate", "Selector")])
  rows <- list()
  for (i in seq_len(nrow(groups))) {
    g <- groups[i, , drop = FALSE]
    x <- selected[
      selected$Dataset == g$Dataset & selected$Replicate == g$Replicate & selected$Selector == g$Selector,
      ,
      drop = FALSE
    ]
    models <- sort(unique(x$LLMModelID))
    if (length(models) < 2) next
    pairs <- utils::combn(models, 2, simplify = FALSE)
    for (pair in pairs) {
      a <- unique(x$Cluster[x$LLMModelID == pair[[1]]])
      b <- unique(x$Cluster[x$LLMModelID == pair[[2]]])
      rows[[length(rows) + 1]] <- data.frame(
        Dataset = g$Dataset,
        Replicate = g$Replicate,
        Selector = g$Selector,
        ModelA = pair[[1]],
        ModelB = pair[[2]],
        NSelectedA = length(a),
        NSelectedB = length(b),
        Intersection = length(intersect(a, b)),
        Union = length(union(a, b)),
        Jaccard = safe_divide(length(intersect(a, b)), length(union(a, b))),
        OverlapOfSmallerSet = safe_divide(length(intersect(a, b)), min(length(a), length(b))),
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) == 0) data.frame() else do.call(rbind, rows)
}

plot_universal_outputs <- function(prefix, transfer, selector_summary, calibration, overlap) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(FALSE))
  if (nrow(transfer) > 0) {
    p <- ggplot2::ggplot(
      transfer,
      ggplot2::aes(x = .data$TestModelID, y = .data$TrainModelID, fill = .data$AUROC)
    ) +
      ggplot2::geom_tile(color = "white") +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", .data$AUROC)), size = 3) +
      ggplot2::scale_fill_gradient(low = "#F7FBFF", high = "#2C7FB8", na.value = "grey90") +
      ggplot2::labs(x = "Test annotation model", y = "Frozen Risk-k source", fill = "AUROC") +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
    ggplot2::ggsave(paste0(prefix, "_auroc_transfer_heatmap.pdf"), p, width = 7, height = 3.8)
  }
  if (nrow(selector_summary) > 0) {
    p <- ggplot2::ggplot(
      selector_summary,
      ggplot2::aes(x = .data$LLMModelID, y = .data$CorrectionEfficiency, fill = .data$Selector)
    ) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::labs(x = "Annotation model", y = "Net correction efficiency", fill = "Selector") +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
    ggplot2::ggsave(paste0(prefix, "_selector_efficiency_by_model.pdf"), p, width = 8, height = 4.8)
  }
  if (nrow(calibration) > 0) {
    p <- ggplot2::ggplot(
      calibration,
      ggplot2::aes(x = .data$MeanPredictedRisk, y = .data$ObservedErrorRate, color = .data$LLMModelID)
    ) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
      ggplot2::geom_point(ggplot2::aes(size = .data$N)) +
      ggplot2::geom_line() +
      ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
      ggplot2::labs(x = "Predicted error risk", y = "Observed error rate", color = "Model") +
      ggplot2::theme_minimal(base_size = 11)
    ggplot2::ggsave(paste0(prefix, "_risk_calibration_by_model.pdf"), p, width = 6.5, height = 5.5)
  }
  if (nrow(overlap) > 0) {
    p <- ggplot2::ggplot(
      overlap,
      ggplot2::aes(x = .data$ModelA, y = .data$ModelB, fill = .data$Jaccard)
    ) +
      ggplot2::geom_tile(color = "white") +
      ggplot2::facet_wrap(~Selector) +
      ggplot2::scale_fill_gradient(low = "#F7FBFF", high = "#2C7FB8", na.value = "grey90") +
      ggplot2::labs(x = NULL, y = NULL, fill = "Jaccard") +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
    ggplot2::ggsave(paste0(prefix, "_selector_overlap_by_model.pdf"), p, width = 8, height = 5)
  }
  invisible(TRUE)
}

run_universal_reliability_transfer <- function(debug_dir,
                                               model_rds,
                                               output_prefix,
                                               train_model_id = "deepseek-chat",
                                               budget_k = NULL,
                                               budget_fraction = as.numeric(Sys.getenv("DEEPSEEKCELL_UNIVERSAL_BUDGET_FRACTION", unset = "0.2")),
                                               budget_min = as.integer(Sys.getenv("DEEPSEEKCELL_UNIVERSAL_BUDGET_MIN", unset = "1"))) {
  if (!file.exists(model_rds)) {
    stop("Frozen Risk-k model not found: ", model_rds, call. = FALSE)
  }
  model <- readRDS(model_rds)
  if (!inherits(model, "deepseekcell_reliability_model")) {
    stop("model_rds must contain a deepseekcell_reliability_model object.", call. = FALSE)
  }
  features <- build_universal_features(debug_dir)
  features$FrozenRisk <- predict_reliability_risk(model, features)

  transfer <- summarise_transfer_metrics(features, train_model_id)
  calibration <- calibration_bins(features)
  selectors <- evaluate_selector_transfer(
    features,
    budget_k = budget_k,
    budget_fraction = budget_fraction,
    budget_min = budget_min
  )
  selector_summary <- summarise_selector_behavior(selectors$behavior)
  overlap <- selector_overlap(selectors$selected)

  dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(features, paste0(output_prefix, "_features.csv"), row.names = FALSE)
  utils::write.csv(transfer, paste0(output_prefix, "_risk_transfer_metrics.csv"), row.names = FALSE)
  utils::write.csv(calibration, paste0(output_prefix, "_risk_calibration_bins.csv"), row.names = FALSE)
  utils::write.csv(selectors$behavior, paste0(output_prefix, "_selector_behavior_by_block.csv"), row.names = FALSE)
  utils::write.csv(selector_summary, paste0(output_prefix, "_selector_summary_by_model.csv"), row.names = FALSE)
  utils::write.csv(selectors$selected, paste0(output_prefix, "_selected_clusters.csv"), row.names = FALSE)
  utils::write.csv(overlap, paste0(output_prefix, "_selector_overlap_by_model.csv"), row.names = FALSE)

  plot_universal_outputs(output_prefix, transfer, selector_summary, calibration, overlap)

  invisible(list(
    features = features,
    transfer = transfer,
    calibration = calibration,
    selector_behavior = selectors$behavior,
    selector_summary = selector_summary,
    selected_clusters = selectors$selected,
    selector_overlap = overlap
  ))
}

args <- commandArgs(trailingOnly = TRUE)
debug_dir <- if (length(args) >= 1) args[[1]] else file.path("results", "benchmark_debug")
model_rds <- if (length(args) >= 2) args[[2]] else file.path("results", "reliability_model_v1.1_error.rds")
output_prefix <- if (length(args) >= 3) args[[3]] else file.path("results", "universal_reliability_transfer")
train_model_id <- if (length(args) >= 4) args[[4]] else "deepseek-chat"
budget_k <- if (length(args) >= 5 && nzchar(args[[5]])) as.integer(args[[5]]) else NULL

run_universal_reliability_transfer(
  debug_dir = debug_dir,
  model_rds = model_rds,
  output_prefix = output_prefix,
  train_model_id = train_model_id,
  budget_k = budget_k
)
