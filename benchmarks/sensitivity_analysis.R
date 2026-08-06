# benchmarks/sensitivity_analysis.R
#
# Secondary robustness analyses for the frozen reliability layer. This script
# does not change DeepSeekCell reliability specification v1.0. It evaluates how
# selector rankings and confidence-quality summaries behave under prespecified
# weight and threshold perturbations.
#
# Usage:
#   Rscript benchmarks/sensitivity_analysis.R cluster_level_input.csv results/sensitivity
#
# Required input columns:
#   Dataset, Replicate, Cluster, LLMConfidence, OntologyEvidenceScore,
#   MarkerEvidenceScore, BestMarkerEvidenceScore, TissueEvidenceScore,
#   ConsensusEvidenceScore, InitiallyCorrect
#
# Optional columns:
#   FullRefinedCorrect, LLMBackend, LLMModelID
#
# If FullRefinedCorrect is supplied, the script estimates wrong-to-correct and
# correct-to-wrong outcomes for each simulated selector by treating FullRefined
# as the available second-pass outcome for all clusters.

Sys.setenv(DEEPSEEKCELL_RUN_BENCHMARK_ON_SOURCE = "false")

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
source("benchmarks/ablation.R")

required_columns <- c(
  "Dataset",
  "Replicate",
  "Cluster",
  "LLMConfidence",
  "OntologyEvidenceScore",
  "MarkerEvidenceScore",
  "BestMarkerEvidenceScore",
  "TissueEvidenceScore",
  "ConsensusEvidenceScore",
  "InitiallyCorrect"
)

write_input_schema <- function(path = file.path("results", "sensitivity_input_schema.csv")) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  schema <- data.frame(
    column = c(required_columns, "FullRefinedCorrect", "LLMBackend", "LLMModelID"),
    required = c(rep(TRUE, length(required_columns)), FALSE, FALSE, FALSE),
    description = c(
      "Dataset identifier.",
      "Replicate identifier.",
      "Cluster identifier.",
      "Raw first-pass LLM confidence.",
      "Ontology evidence score from frozen v1.0 scoring.",
      "Marker evidence score for the predicted first-pass label.",
      "Best marker-profile score among supported profiles.",
      "Tissue-consistency evidence score.",
      "Consensus score used by confidence recalibration.",
      "Whether the first-pass label is correct after label harmonisation.",
      "Whether the FullRefined label is correct after label harmonisation.",
      "LLM backend identifier, optional.",
      "LLM model identifier, optional."
    ),
    stringsAsFactors = FALSE
  )
  utils::write.csv(schema, path, row.names = FALSE)
  message("Wrote expected input schema to ", path)
}

normalize_weights <- function(weights) {
  weights <- as.numeric(weights)
  weights[is.na(weights) | weights < 0] <- 0
  if (sum(weights) == 0) {
    weights[] <- 1 / length(weights)
  } else {
    weights <- weights / sum(weights)
  }
  weights
}

weight_library <- function(n_random = 500, seed = 1001) {
  base <- RELIABILITY_CONFIDENCE_WEIGHTS
  named <- list(
    default = base,
    equal = stats::setNames(rep(0.2, 5), names(base)),
    no_ontology = c(
      llm_confidence = 0.45,
      ontology_evidence = 0.00,
      marker_evidence = 0.30,
      tissue_evidence = 0.15,
      consensus_evidence = 0.10
    ),
    no_tissue = c(
      llm_confidence = 0.40,
      ontology_evidence = 0.25,
      marker_evidence = 0.25,
      tissue_evidence = 0.00,
      consensus_evidence = 0.10
    ),
    marker_dominant = c(
      llm_confidence = 0.20,
      ontology_evidence = 0.20,
      marker_evidence = 0.45,
      tissue_evidence = 0.10,
      consensus_evidence = 0.05
    ),
    llm_confidence_dominant = c(
      llm_confidence = 0.60,
      ontology_evidence = 0.15,
      marker_evidence = 0.15,
      tissue_evidence = 0.05,
      consensus_evidence = 0.05
    )
  )
  named <- lapply(named, normalize_weights)

  set.seed(seed)
  random <- replicate(
    n_random,
    normalize_weights(stats::rexp(length(base), rate = 1)),
    simplify = FALSE
  )
  names(random) <- paste0("simplex_random_", seq_along(random))
  random <- lapply(random, stats::setNames, names(base))

  c(named, random)
}

block_id <- function(x) {
  backend <- if ("LLMBackend" %in% names(x)) x$LLMBackend else ""
  paste(x$Dataset, x$Replicate, backend, sep = "|")
}

as_logical_column <- function(x) {
  if (is.logical(x)) return(x)
  tolower(as.character(x)) %in% c("true", "1", "yes", "y")
}

compute_composite_confidence <- function(df, weights) {
  weights <- normalize_weights(weights)
  components <- data.frame(
    llm_confidence = as_confidence(df$LLMConfidence),
    ontology_evidence = as_confidence(df$OntologyEvidenceScore),
    marker_evidence = as_confidence(df$MarkerEvidenceScore),
    tissue_evidence = as_confidence(df$TissueEvidenceScore),
    consensus_evidence = as_confidence(df$ConsensusEvidenceScore),
    stringsAsFactors = FALSE
  )
  as.numeric(as.matrix(components[names(weights)]) %*% weights)
}

binary_brier <- function(correct, confidence) {
  mean((as.numeric(correct) - confidence)^2, na.rm = TRUE)
}

binary_nll <- function(correct, confidence, eps = 1e-6) {
  confidence <- pmin(pmax(confidence, eps), 1 - eps)
  correct <- as.numeric(correct)
  -mean(correct * log(confidence) + (1 - correct) * log(1 - confidence), na.rm = TRUE)
}

expected_calibration_error <- function(correct, confidence, n_bins = 10) {
  ok <- !is.na(correct) & !is.na(confidence)
  correct <- as.numeric(correct[ok])
  confidence <- confidence[ok]
  if (length(correct) == 0) return(NA_real_)

  breaks <- unique(stats::quantile(confidence, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE))
  if (length(breaks) <= 2) {
    breaks <- seq(0, 1, length.out = min(n_bins, length(correct)) + 1)
  }
  bins <- cut(confidence, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  sum(vapply(stats::na.omit(unique(bins)), function(bin) {
    idx <- bins == bin
    mean(idx) * abs(mean(correct[idx]) - mean(confidence[idx]))
  }, numeric(1)))
}

select_by_threshold <- function(df, tau, k = NULL) {
  best <- suppressWarnings(as.numeric(df$BestMarkerEvidenceScore))
  predicted <- suppressWarnings(as.numeric(df$MarkerEvidenceScore))
  conflict <- best >= tau &
    predicted <= pmax(
      RELIABILITY_THRESHOLDS$strong_marker_floor,
      best - RELIABILITY_THRESHOLDS$strong_marker_margin
    )
  conflict[is.na(conflict)] <- FALSE
  score <- pmax(best - predicted, 0)
  candidates <- which(conflict)
  if (length(candidates) == 0) return(integer())
  candidates <- candidates[order(score[candidates], decreasing = TRUE)]
  if (!is.null(k)) candidates <- head(candidates, k)
  candidates
}

summarise_selected <- function(df, idx, label) {
  initially_correct <- as_logical_column(df$InitiallyCorrect)
  has_full <- "FullRefinedCorrect" %in% names(df)
  full_correct <- if (has_full) as_logical_column(df$FullRefinedCorrect) else rep(NA, nrow(df))

  n_refined <- length(idx)
  wrong_to_correct <- if (has_full) sum(!initially_correct[idx] & full_correct[idx], na.rm = TRUE) else NA_integer_
  correct_to_wrong <- if (has_full) sum(initially_correct[idx] & !full_correct[idx], na.rm = TRUE) else NA_integer_
  initially_wrong <- sum(!initially_correct, na.rm = TRUE)
  flagged_wrong <- sum(!initially_correct[idx], na.rm = TRUE)
  flagged_correct <- sum(initially_correct[idx], na.rm = TRUE)

  data.frame(
    Analysis = label,
    NClusters = nrow(df),
    NRefined = n_refined,
    InitiallyWrong = initially_wrong,
    SelectionPrecision = if (n_refined > 0) flagged_wrong / n_refined else NA_real_,
    SelectionRecall = if (initially_wrong > 0) flagged_wrong / initially_wrong else NA_real_,
    WrongToCorrect = wrong_to_correct,
    CorrectToWrong = correct_to_wrong,
    CorrectionEfficiency = if (n_refined > 0 && has_full) {
      (wrong_to_correct - correct_to_wrong) / n_refined
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )
}

run_sensitivity_analysis <- function(input_csv,
                                     output_prefix = file.path("results", "sensitivity"),
                                     n_random_weights = 500) {
  if (!file.exists(input_csv)) {
    write_input_schema()
    stop("Input CSV not found: ", input_csv, call. = FALSE)
  }

  df <- utils::read.csv(input_csv, stringsAsFactors = FALSE)
  missing <- setdiff(required_columns, names(df))
  if (length(missing) > 0) {
    write_input_schema()
    stop("Input CSV is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  df$BlockID <- block_id(df)
  dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)

  weights <- weight_library(n_random = n_random_weights)
  correct <- as_logical_column(df$InitiallyCorrect)
  confidence_summary <- do.call(rbind, lapply(names(weights), function(name) {
    confidence <- compute_composite_confidence(df, weights[[name]])
    data.frame(
      WeightSet = name,
      Brier = binary_brier(correct, confidence),
      BinaryCorrectnessNLL = binary_nll(correct, confidence),
      ECE = expected_calibration_error(correct, confidence),
      MeanConfidence = mean(confidence, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))

  default_confidence <- compute_composite_confidence(df, weights$default)
  default_order <- order(default_confidence, decreasing = FALSE)
  ranking_summary <- do.call(rbind, lapply(names(weights), function(name) {
    confidence <- compute_composite_confidence(df, weights[[name]])
    ord <- order(confidence, decreasing = FALSE)
    k <- max(1, ceiling(0.1 * nrow(df)))
    data.frame(
      WeightSet = name,
      SpearmanWithDefault = suppressWarnings(stats::cor(default_confidence, confidence, method = "spearman")),
      JaccardTop10Pct = length(intersect(head(default_order, k), head(ord, k))) /
        length(union(head(default_order, k), head(ord, k))),
      stringsAsFactors = FALSE
    )
  }))

  threshold_grid <- seq(0.25, 0.70, by = 0.05)
  threshold_summary <- do.call(rbind, lapply(threshold_grid, function(tau) {
    per_block <- do.call(rbind, lapply(split(df, df$BlockID), function(block_df) {
      idx <- select_by_threshold(block_df, tau)
      summarise_selected(block_df, idx, paste0("tau_", sprintf("%.2f", tau)))
    }))
    data.frame(
      Threshold = tau,
      NBlocks = length(unique(df$BlockID)),
      NRefined = sum(per_block$NRefined, na.rm = TRUE),
      SelectionPrecision = weighted.mean(
        per_block$SelectionPrecision,
        w = pmax(per_block$NRefined, 1),
        na.rm = TRUE
      ),
      SelectionRecall = mean(per_block$SelectionRecall, na.rm = TRUE),
      WrongToCorrect = sum(per_block$WrongToCorrect, na.rm = TRUE),
      CorrectToWrong = sum(per_block$CorrectToWrong, na.rm = TRUE),
      CorrectionEfficiency = {
        n_refined <- sum(per_block$NRefined, na.rm = TRUE)
        if (n_refined > 0) {
          (sum(per_block$WrongToCorrect, na.rm = TRUE) -
             sum(per_block$CorrectToWrong, na.rm = TRUE)) / n_refined
        } else {
          NA_real_
        }
      },
      stringsAsFactors = FALSE
    )
  }))

  utils::write.csv(confidence_summary, paste0(output_prefix, "_weight_confidence.csv"), row.names = FALSE)
  utils::write.csv(ranking_summary, paste0(output_prefix, "_weight_rankings.csv"), row.names = FALSE)
  utils::write.csv(threshold_summary, paste0(output_prefix, "_threshold_grid.csv"), row.names = FALSE)

  invisible(list(
    confidence = confidence_summary,
    rankings = ranking_summary,
    thresholds = threshold_summary
  ))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  write_input_schema()
  message("No input supplied; wrote schema only.")
} else {
  output_prefix <- if (length(args) >= 2) args[2] else file.path("results", "sensitivity")
  n_random_weights <- if (length(args) >= 3) as.integer(args[3]) else 500
  run_sensitivity_analysis(args[1], output_prefix, n_random_weights)
}
