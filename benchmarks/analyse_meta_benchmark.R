# benchmarks/analyse_meta_benchmark.R
#
# Meta-benchmark analyses for DeepSeekCell. This script does not define a new
# selector or retune the method. It asks which dataset characteristics explain
# annotation accuracy, calibration, correction efficiency, and refinement value.
#
# Usage:
#   Rscript benchmarks/analyse_meta_benchmark.R results results/meta_benchmark
#   Rscript benchmarks/analyse_meta_benchmark.R results results/meta_benchmark results/reliability_model_v1.1_error.rds

repo_root <- normalizePath(getwd(), mustWork = TRUE)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

as_number <- function(x) {
  suppressWarnings(as.numeric(x))
}

as_flag_meta <- function(x) {
  if (is.logical(x)) return(x)
  if (is.numeric(x)) return(x != 0)
  tolower(trimws(as.character(x))) %in% c("true", "1", "yes", "y")
}

safe_mean <- function(x) {
  x <- as_number(x)
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_median <- function(x) {
  x <- as_number(x)
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  stats::median(x, na.rm = TRUE)
}

safe_divide <- function(num, den) {
  ifelse(!is.na(den) & den != 0, num / den, NA_real_)
}

value_or_unknown <- function(x) {
  x <- as.character(x)
  x[is.na(x) | !nzchar(trimws(x)) | x %in% c("NA", "NaN")] <- "Unknown"
  trimws(x)
}

ensure_column <- function(x, column, default = NA) {
  if (!column %in% names(x)) {
    x[[column]] <- default
  }
  x
}

column_or_default <- function(x, column, default = NA) {
  if (!is.data.frame(x) || !column %in% names(x)) {
    return(rep(default, nrow(x)))
  }
  x[[column]]
}

read_optional_csv <- function(path) {
  if (!file.exists(path)) return(data.frame())
  utils::read.csv(path, stringsAsFactors = FALSE)
}

bind_rows_base <- function(rows) {
  rows <- rows[vapply(rows, function(x) is.data.frame(x) && nrow(x) > 0, logical(1))]
  if (length(rows) == 0) return(data.frame())
  all_names <- unique(unlist(lapply(rows, names)))
  rows <- lapply(rows, function(x) {
    missing <- setdiff(all_names, names(x))
    for (column in missing) x[[column]] <- NA
    x[all_names]
  })
  do.call(rbind, rows)
}

read_debug_features <- function(debug_dir) {
  if (!dir.exists(debug_dir)) return(data.frame())
  files <- list.files(debug_dir, pattern = "-Calibrated_debug\\.csv$", full.names = TRUE)
  if (length(files) == 0) {
    files <- list.files(debug_dir, pattern = "-Evidence_debug\\.csv$", full.names = TRUE)
  }
  rows <- lapply(files, function(path) {
    x <- utils::read.csv(path, stringsAsFactors = FALSE)
    x$SourceFile <- basename(path)
    x
  })
  bind_rows_base(rows)
}

load_feature_table <- function(results_dir) {
  candidates <- file.path(
    results_dir,
    c(
      "reliability_generalization_features.csv",
      "reliability_model_v1.1_error_training_features.csv"
    )
  )
  for (path in candidates) {
    if (file.exists(path)) return(utils::read.csv(path, stringsAsFactors = FALSE))
  }
  read_debug_features(file.path(results_dir, "benchmark_debug"))
}

normalize_01 <- function(x) {
  x <- as_number(x)
  if (length(x) == 0 || all(is.na(x))) return(rep(NA_real_, length(x)))
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || identical(rng[1], rng[2])) {
    return(rep(0.5, length(x)))
  }
  (x - rng[1]) / (rng[2] - rng[1])
}

row_mean_na <- function(df) {
  if (!is.data.frame(df) || ncol(df) == 0) return(rep(NA_real_, nrow(df)))
  apply(df, 1, function(x) {
    x <- as_number(x)
    if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  })
}

block_id <- function(x) {
  for (column in c("Dataset", "Replicate", "LLMBackend", "LLMModelID")) {
    x <- ensure_column(x, column, "")
  }
  paste(x$Dataset, x$Replicate, x$LLMBackend, x$LLMModelID, sep = "|")
}

selector_label <- function(x) {
  selector <- if ("RefinementSelector" %in% names(x)) value_or_unknown(x$RefinementSelector) else rep("Unknown", nrow(x))
  method <- if ("Method" %in% names(x)) as.character(x$Method) else rep("Unknown", nrow(x))
  out <- ifelse(selector != "Unknown" & selector != "none", selector, method)
  out <- gsub("^DeepSeekCell-", "", out)
  out <- gsub("^DeepSeek-", "", out)
  tolower(out)
}

annotation_block_key <- function(x) {
  paste(x$Dataset, x$Replicate, sep = "|")
}

summarise_by_dataset <- function(x, fun) {
  if (!is.data.frame(x) || nrow(x) == 0 || !"Dataset" %in% names(x)) return(data.frame())
  bind_rows_base(lapply(split(x, x$Dataset), fun))
}

build_dataset_difficulty <- function(dataset_summary,
                                     cluster_summary,
                                     benchmark_results,
                                     refinement,
                                     confidence,
                                     features) {
  dataset_summary <- ensure_column(dataset_summary, "Dataset")
  cluster_summary <- ensure_column(cluster_summary, "Dataset")
  benchmark_results <- ensure_column(benchmark_results, "Dataset")
  refinement <- ensure_column(refinement, "Dataset")
  confidence <- ensure_column(confidence, "Dataset")
  features <- ensure_column(features, "Dataset")

  dataset_agg <- summarise_by_dataset(dataset_summary, function(x) {
    data.frame(
      Dataset = x$Dataset[[1]],
      Tissue = value_or_unknown(x$Tissue)[[1]],
      Species = value_or_unknown(x$Species)[[1]],
      NCells = safe_mean(x$NCells),
      NGenes = safe_mean(x$NGenes),
      NClusters = safe_mean(x$NClusters),
      MeanClusterPurity = safe_mean(x$MeanClusterPurity),
      MedianClusterPurity = safe_mean(x$MedianClusterPurity),
      MinClusterPurity = safe_mean(x$MinClusterPurity),
      TopMarkers = safe_mean(x$TopMarkers),
      stringsAsFactors = FALSE
    )
  })

  cluster_agg <- summarise_by_dataset(cluster_summary, function(x) {
    truth <- value_or_unknown(x$Truth)
    n_clusters <- nrow(x)
    unique_truth <- length(unique(truth))
    purity <- as_number(x$ClusterPurity)
    data.frame(
      Dataset = x$Dataset[[1]],
      ClusterRows = n_clusters,
      MeanMarkerCount = safe_mean(x$NMarkers),
      MedianMarkerCount = safe_median(x$NMarkers),
      LowPurityClusterRate = mean(purity < 0.8, na.rm = TRUE),
      VeryLowPurityClusterRate = mean(purity < 0.6, na.rm = TRUE),
      UniqueTruthLabels = unique_truth,
      LabelGranularity = safe_divide(unique_truth, n_clusters),
      DuplicateTruthLabelRate = 1 - safe_divide(unique_truth, n_clusters),
      stringsAsFactors = FALSE
    )
  })

  feature_agg <- summarise_by_dataset(features, function(x) {
    first_wrong <- if ("FirstPassIncorrect" %in% names(x)) {
      as_flag_meta(x$FirstPassIncorrect)
    } else if (all(c("HarmonisedPrediction", "HarmonisedTruth") %in% names(x))) {
      as.character(x$HarmonisedPrediction) != as.character(x$HarmonisedTruth)
    } else {
      rep(NA, nrow(x))
    }
    data.frame(
      Dataset = x$Dataset[[1]],
      LLMUncertainty = 1 - safe_mean(x$LLMConfidence),
      OntologyAmbiguity = 1 - safe_mean(x$OntologyEvidenceScore),
      MarkerAmbiguity = 1 - safe_mean(x$MarkerEvidenceScore),
      BestMarkerSupport = safe_mean(x$BestMarkerEvidenceScore),
      MarkerMarginMean = safe_mean(x$MarkerMargin),
      TissueAmbiguity = 1 - safe_mean(x$TissueEvidenceScore),
      ConsensusAmbiguity = 1 - safe_mean(x$ConsensusEvidenceScore),
      EvidenceConflictRate = mean(as_number(x$RequiresRefinementNumeric) > 0, na.rm = TRUE),
      CandidateDisagreementRate = mean(as_number(x$CandidateEvidenceDisagreement) > 0, na.rm = TRUE),
      FirstPassErrorRate = mean(first_wrong, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })

  plain <- benchmark_results[
    column_or_default(benchmark_results, "Method", "") %in% c("DeepSeek-Plain", "DeepSeek Plain") |
      column_or_default(benchmark_results, "RefinementSelector", "") %in% c("none", "None"),
    ,
    drop = FALSE
  ]
  plain <- plain[grepl("DeepSeek", plain$Method), , drop = FALSE]
  plain_agg <- summarise_by_dataset(plain, function(x) {
    data.frame(
      Dataset = x$Dataset[[1]],
      PlainMacroF1 = safe_mean(x$MacroF1),
      PlainAccuracy = safe_mean(x$Accuracy),
      PlainBalancedAcc = safe_mean(x$BalancedAcc),
      PlainCladeAcc = safe_mean(x$CladeAcc),
      PlainUnknownRate = safe_mean(x$UnknownRate),
      stringsAsFactors = FALSE
    )
  })

  full <- refinement[refinement$RefinementSelector %in% c("full", "FullRefined"), , drop = FALSE]
  full_agg <- summarise_by_dataset(full, function(x) {
    data.frame(
      Dataset = x$Dataset[[1]],
      FullWrongToCorrect = sum(as_number(x$WrongToCorrect), na.rm = TRUE),
      FullCorrectToWrong = sum(as_number(x$CorrectToWrong), na.rm = TRUE),
      OracleCorrectionOpportunity = sum(as_number(x$WrongToCorrect), na.rm = TRUE) -
        sum(as_number(x$CorrectToWrong), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })

  conf_plain <- confidence[confidence$Method %in% c("DeepSeek-Plain", "DeepSeek Plain"), , drop = FALSE]
  conf_agg <- summarise_by_dataset(conf_plain, function(x) {
    data.frame(
      Dataset = x$Dataset[[1]],
      PlainBrier = safe_mean(x$Brier),
      PlainECE = safe_mean(x$ECE),
      PlainAURC = safe_mean(x$AURC),
      stringsAsFactors = FALSE
    )
  })

  out <- Reduce(function(a, b) merge(a, b, by = "Dataset", all = TRUE), list(
    dataset_agg, cluster_agg, feature_agg, plain_agg, full_agg, conf_agg
  ))

  intrinsic_components <- data.frame(
    ClusterComplexity = normalize_01(log1p(out$NClusters)),
    CellScale = normalize_01(log1p(out$NCells)),
    LowPurity = 1 - as_number(out$MeanClusterPurity),
    MinimumPurityPenalty = 1 - as_number(out$MinClusterPurity),
    MarkerAmbiguity = as_number(out$MarkerAmbiguity),
    OntologyAmbiguity = as_number(out$OntologyAmbiguity),
    TissueAmbiguity = as_number(out$TissueAmbiguity),
    EvidenceConflictRate = as_number(out$EvidenceConflictRate),
    CandidateDisagreementRate = as_number(out$CandidateDisagreementRate),
    LabelGranularity = as_number(out$LabelGranularity),
    stringsAsFactors = FALSE
  )
  observed_components <- data.frame(
    IntrinsicDifficulty = row_mean_na(intrinsic_components),
    PlainError = 1 - as_number(out$PlainAccuracy),
    PlainMacroF1Loss = 1 - as_number(out$PlainMacroF1),
    PlainUnknownRate = as_number(out$PlainUnknownRate),
    PlainCalibrationError = normalize_01(out$PlainBrier),
    FirstPassErrorRate = as_number(out$FirstPassErrorRate),
    stringsAsFactors = FALSE
  )

  out$IntrinsicDifficultyScore <- row_mean_na(intrinsic_components)
  out$ObservedDifficultyScore <- row_mean_na(observed_components)
  q <- stats::quantile(out$ObservedDifficultyScore, probs = c(1 / 3, 2 / 3), na.rm = TRUE)
  out$DifficultyClass <- cut(
    out$ObservedDifficultyScore,
    breaks = c(-Inf, q[[1]], q[[2]], Inf),
    labels = c("Low", "Medium", "High"),
    include.lowest = TRUE
  )

  component_out <- intrinsic_components
  names(component_out) <- paste0("Component_", names(component_out))
  cbind(out, component_out)
}

compute_method_gains <- function(benchmark_results, difficulty) {
  if (nrow(benchmark_results) == 0) return(data.frame())
  benchmark_results$BlockKey <- annotation_block_key(benchmark_results)
  benchmark_results$Selector <- selector_label(benchmark_results)
  plain <- benchmark_results[
    benchmark_results$Method %in% c("DeepSeek-Plain", "DeepSeek Plain") |
      benchmark_results$Selector %in% c("plain", "none"),
    ,
    drop = FALSE
  ]
  plain <- plain[grepl("DeepSeek", plain$Method), , drop = FALSE]
  plain_key <- plain$BlockKey
  idx <- match(benchmark_results$BlockKey, plain_key)
  out <- benchmark_results
  out$PlainMacroF1Block <- as_number(plain$MacroF1)[idx]
  out$PlainAccuracyBlock <- as_number(plain$Accuracy)[idx]
  out$PlainBrierBlock <- NA_real_
  out$DeltaMacroF1VsPlain <- as_number(out$MacroF1) - out$PlainMacroF1Block
  out$DeltaAccuracyVsPlain <- as_number(out$Accuracy) - out$PlainAccuracyBlock
  merge(out, difficulty, by = "Dataset", all.x = TRUE)
}

compute_selector_difficulty_links <- function(refinement, difficulty) {
  if (nrow(refinement) == 0) return(data.frame())
  refinement$Selector <- selector_label(refinement)
  out <- merge(refinement, difficulty, by = "Dataset", all.x = TRUE)
  out$NetCorrections <- as_number(out$WrongToCorrect) - as_number(out$CorrectToWrong)
  out$CorrectionEfficiency <- as_number(out$CorrectionEfficiency)
  out$RecoveryFraction <- as_number(out$RecoveryFraction)
  out
}

is_pareto_front <- function(df, maximize, minimize) {
  if (nrow(df) == 0) return(logical())
  keep <- rep(TRUE, nrow(df))
  max_mat <- as.matrix(data.frame(lapply(df[maximize], as_number)))
  min_mat <- as.matrix(data.frame(lapply(df[minimize], as_number)))
  for (i in seq_len(nrow(df))) {
    for (j in seq_len(nrow(df))) {
      if (i == j) next
      better_or_equal_max <- all(max_mat[j, ] >= max_mat[i, ], na.rm = TRUE)
      better_or_equal_min <- all(min_mat[j, ] <= min_mat[i, ], na.rm = TRUE)
      strictly_better <- any(max_mat[j, ] > max_mat[i, ], na.rm = TRUE) ||
        any(min_mat[j, ] < min_mat[i, ], na.rm = TRUE)
      if (isTRUE(better_or_equal_max && better_or_equal_min && strictly_better)) {
        keep[i] <- FALSE
        break
      }
    }
  }
  keep
}

compute_pareto_tables <- function(benchmark_results, refinement) {
  annotation <- benchmark_results
  if (nrow(annotation) > 0) {
    annotation$Selector <- selector_label(annotation)
    annotation$MacroF1 <- as_number(annotation$MacroF1)
    annotation$Accuracy <- as_number(annotation$Accuracy)
    annotation$CostUSD <- as_number(annotation$CostUSD)
    annotation$RuntimeSec <- as_number(annotation$RuntimeSec)
    annotation$Tokens <- as_number(annotation$Tokens)
    annotation$SecondPassCalls <- as_number(annotation$SecondPassCalls)
    annotation$ParetoMacroF1Cost <- is_pareto_front(annotation, "MacroF1", "CostUSD")
    annotation$ParetoMacroF1Runtime <- is_pareto_front(annotation, "MacroF1", "RuntimeSec")
    annotation$ParetoMacroF1Tokens <- is_pareto_front(annotation, "MacroF1", "Tokens")
    annotation$ParetoAccuracyCalls <- is_pareto_front(annotation, "Accuracy", "SecondPassCalls")
  }

  refinement_out <- refinement
  if (nrow(refinement_out) > 0) {
    refinement_out$Selector <- selector_label(refinement_out)
    for (column in c(
      "CorrectionEfficiency", "RecoveryFraction", "RefinementRuntimeSec",
      "RefinementTokens", "RefinementCostUSD", "SecondPassCalls", "NRefined"
    )) {
      refinement_out[[column]] <- as_number(refinement_out[[column]])
    }
    refinement_out$ParetoEfficiencyRuntime <- is_pareto_front(
      refinement_out,
      "CorrectionEfficiency",
      "RefinementRuntimeSec"
    )
    refinement_out$ParetoEfficiencyTokens <- is_pareto_front(
      refinement_out,
      "CorrectionEfficiency",
      "RefinementTokens"
    )
    refinement_out$ParetoRecoveryCalls <- is_pareto_front(
      refinement_out,
      "RecoveryFraction",
      "SecondPassCalls"
    )
  }

  list(annotation = annotation, refinement = refinement_out)
}

normalize_label <- function(x) {
  x <- tolower(as.character(x))
  x <- gsub("[^a-z0-9]+", " ", x)
  trimws(gsub("\\s+", " ", x))
}

cell_family <- function(label) {
  x <- normalize_label(label)
  out <- rep("other", length(x))
  out[grepl("\\bt cell\\b|\\bt lymphocyte\\b|\\bcd4\\b|\\bcd8\\b", x)] <- "t_cell"
  out[grepl("\\bb cell\\b|\\bplasma\\b", x)] <- "b_cell"
  out[grepl("\\bnk\\b|natural killer", x)] <- "nk_cell"
  out[grepl("monocyte|macrophage|myeloid", x)] <- "myeloid"
  out[grepl("neutrophil", x)] <- "neutrophil"
  out[grepl("dendritic", x)] <- "dendritic"
  out[grepl("alpha cell|beta cell|delta cell|gamma|acinar|ductal|islet", x)] <- "pancreas"
  out[grepl("neuron|astro|oligo|microglia|ependymal", x)] <- "brain"
  out[grepl("epithelial|basal|goblet|secretory|alveolar|ionocyte", x)] <- "epithelial"
  out[grepl("endothelial", x)] <- "endothelial"
  out[grepl("fibroblast|stellate", x)] <- "stromal"
  out
}

contains_developmental_stage <- function(x) {
  grepl("progenitor|precursor|immature|mature|stem|development", normalize_label(x))
}

contains_tumour <- function(x) {
  grepl("tumou?r|cancer|malignant|neoplastic|carcinoma|normal", normalize_label(x))
}

classify_error_type <- function(row) {
  field <- function(name, default = "") {
    if (!name %in% names(row)) return(default)
    value <- as.character(row[[name]][[1]])
    if (is.na(value) || !nzchar(value)) default else value
  }
  num_field <- function(name, default = NA_real_) {
    value <- if (name %in% names(row)) as_number(row[[name]][[1]]) else default
    if (is.na(value)) default else value
  }
  pred <- field("HarmonisedPrediction", field("RawPrediction", ""))
  truth <- field("HarmonisedTruth", field("RawTruth", ""))
  pred_norm <- normalize_label(pred)
  truth_norm <- normalize_label(truth)
  marker <- num_field("MarkerEvidenceScore")
  best_marker <- num_field("BestMarkerEvidenceScore")
  ontology <- num_field("OntologyEvidenceScore")
  tissue <- num_field("TissueEvidenceScore")
  conflict <- isTRUE(as_flag_meta(field("EvidenceConflict", "FALSE"))) ||
    isTRUE(num_field("EvidenceConflictScore", 0) >= 0.25)

  if (is.na(pred_norm) || !nzchar(pred_norm) || pred_norm %in% c("unknown", "unassigned", "ambiguous")) {
    return("unknown_or_abstained")
  }
  if (contains_tumour(pred) || contains_tumour(truth)) {
    if (cell_family(pred) != cell_family(truth) || pred_norm != truth_norm) {
      return("tumour_normal_confusion")
    }
  }
  if (contains_developmental_stage(pred) || contains_developmental_stage(truth)) {
    return("developmental_stage_confusion")
  }
  if (!is.na(tissue) && tissue < 0.5) {
    return("tissue_mismatch")
  }
  if (!is.na(ontology) && ontology < 0.5) {
    return("ontology_ambiguity")
  }
  if (isTRUE(conflict)) {
    return("evidence_conflict")
  }
  if (!is.na(marker) && !is.na(best_marker) && marker < 0.2 && best_marker < 0.2) {
    return("insufficient_marker_evidence")
  }
  pred_family <- cell_family(pred)
  truth_family <- cell_family(truth)
  if (pred_family == truth_family && pred_family != "other") {
    if (pred_family %in% c("t_cell", "b_cell", "nk_cell", "myeloid", "dendritic", "neutrophil")) {
      return("immune_subtype_confusion")
    }
    return("sibling_or_related_cell_type")
  }
  if (!is.na(marker) && marker < 0.2 && !is.na(ontology) && ontology >= 0.8) {
    return("ontology_synonym_or_label_mapping")
  }
  "biologically_related_or_other"
}

build_error_taxonomy <- function(features) {
  if (nrow(features) == 0) return(list(clusters = data.frame(), summary = data.frame()))
  for (column in c("HarmonisedPrediction", "HarmonisedTruth", "RawPrediction", "RawTruth")) {
    features <- ensure_column(features, column, NA_character_)
  }
  incorrect <- features[
    value_or_unknown(features$HarmonisedPrediction) != value_or_unknown(features$HarmonisedTruth),
    ,
    drop = FALSE
  ]
  if (nrow(incorrect) == 0) return(list(clusters = data.frame(), summary = data.frame()))
  category <- vapply(seq_len(nrow(incorrect)), function(i) {
    classify_error_type(incorrect[i, , drop = FALSE])
  }, character(1))
  keep <- intersect(
    c(
      "Dataset", "Tissue", "Species", "Replicate", "Cluster", "RawPrediction",
      "RawTruth", "HarmonisedPrediction", "HarmonisedTruth", "LLMConfidence",
      "Confidence", "OntologyEvidenceScore", "MarkerEvidenceScore",
      "BestMarkerEvidenceScore", "TissueEvidenceScore", "EvidenceConflictScore"
    ),
    names(incorrect)
  )
  clusters <- incorrect[keep]
  clusters$ErrorCategory <- category
  summary <- as.data.frame(table(ErrorCategory = category), stringsAsFactors = FALSE)
  summary$Fraction <- summary$Freq / sum(summary$Freq)
  summary <- summary[order(summary$Freq, decreasing = TRUE), , drop = FALSE]
  rownames(summary) <- NULL
  list(clusters = clusters, summary = summary)
}

equal_frequency_bins <- function(score, n_bins = 10) {
  score <- as_number(score)
  ok <- !is.na(score)
  bin <- rep(NA_integer_, length(score))
  if (!any(ok)) return(bin)
  ord <- order(score[ok])
  local_bin <- ceiling(seq_along(ord) / length(ord) * n_bins)
  local_bin <- pmax(pmin(local_bin, n_bins), 1)
  bin[which(ok)[ord]] <- local_bin
  bin
}

build_reliability_curves <- function(features) {
  if (nrow(features) == 0) return(data.frame())
  correct <- if ("InitiallyCorrect" %in% names(features)) {
    as_flag_meta(features$InitiallyCorrect)
  } else if ("FirstPassIncorrect" %in% names(features)) {
    !as_flag_meta(features$FirstPassIncorrect)
  } else {
    value_or_unknown(features$HarmonisedPrediction) == value_or_unknown(features$HarmonisedTruth)
  }
  axes <- intersect(
    c(
      "LLMConfidence", "Confidence", "OntologyEvidenceScore", "MarkerEvidenceScore",
      "BestMarkerEvidenceScore", "TissueEvidenceScore", "ConsensusEvidenceScore",
      "EvidenceConflictScore", "MarkerMargin", "RequiresRefinementNumeric",
      "PredictedReliabilityRisk"
    ),
    names(features)
  )
  rows <- list()
  for (axis in axes) {
    score <- as_number(features[[axis]])
    bin <- equal_frequency_bins(score, n_bins = 10)
    for (b in sort(stats::na.omit(unique(bin)))) {
      idx <- bin == b
      rows[[length(rows) + 1]] <- data.frame(
        Axis = axis,
        Bin = b,
        N = sum(idx, na.rm = TRUE),
        MeanScore = mean(score[idx], na.rm = TRUE),
        ObservedCorrectRate = mean(correct[idx], na.rm = TRUE),
        ObservedErrorRate = mean(!correct[idx], na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }
  bind_rows_base(rows)
}

add_risk_predictions_if_available <- function(features, model_path = "") {
  if (!nzchar(model_path) || !file.exists(model_path) || nrow(features) == 0) {
    return(features)
  }
  Sys.setenv(DEEPSEEKCELL_RUN_RELIABILITY_TRAINING_ON_SOURCE = "false")
  source("R/api.R")
  source("R/utils.R")
  source("R/ontology.R")
  source("R/refinement.R")
  source("R/reliability_model.R")
  model <- tryCatch(readRDS(model_path), error = function(e) NULL)
  if (is.null(model)) return(features)
  features$PredictedReliabilityRisk <- tryCatch(
    predict_reliability_risk(model, features),
    error = function(e) rep(NA_real_, nrow(features))
  )
  features
}

decision_curve_scores <- function(features) {
  scores <- list(
    ConfidenceRisk = 1 - as_number(column_or_default(features, "LLMConfidence", 0.5)),
    EvidenceConflict = as_number(column_or_default(features, "EvidenceConflictScore", 0)),
    MarkerMargin = as_number(column_or_default(features, "MarkerMargin", 0)),
    Random = stats::runif(nrow(features))
  )
  if ("PredictedReliabilityRisk" %in% names(features)) {
    scores$RiskProbability <- as_number(features$PredictedReliabilityRisk)
  }
  scores
}

build_decision_curves <- function(features,
                                  budget_grid = c(0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1),
                                  seed = 2026) {
  if (nrow(features) == 0 || !"FullRefinedCorrect" %in% names(features)) {
    return(data.frame())
  }
  features$BlockID <- block_id(features)
  initially_correct <- if ("InitiallyCorrect" %in% names(features)) {
    as_flag_meta(features$InitiallyCorrect)
  } else if ("FirstPassIncorrect" %in% names(features)) {
    !as_flag_meta(features$FirstPassIncorrect)
  } else {
    value_or_unknown(features$HarmonisedPrediction) == value_or_unknown(features$HarmonisedTruth)
  }
  full_correct <- as_flag_meta(features$FullRefinedCorrect)
  initially_wrong <- !initially_correct
  set.seed(seed)
  scores <- decision_curve_scores(features)
  rows <- list()

  for (score_name in names(scores)) {
    score <- scores[[score_name]]
    score[is.na(score)] <- -Inf
    for (fraction in budget_grid) {
      selected <- rep(FALSE, nrow(features))
      for (block in unique(features$BlockID)) {
        idx <- which(features$BlockID == block)
        k <- ceiling(length(idx) * fraction)
        if (k <= 0) next
        local_order <- order(score[idx], decreasing = TRUE)
        selected[idx[head(local_order, k)]] <- TRUE
      }
      full_wrong_to_correct <- initially_wrong & full_correct
      full_correct_to_wrong <- initially_correct & !full_correct
      full_net <- sum(full_wrong_to_correct, na.rm = TRUE) - sum(full_correct_to_wrong, na.rm = TRUE)
      wrong_to_correct <- selected & initially_wrong & full_correct
      correct_to_wrong <- selected & initially_correct & !full_correct
      net <- sum(wrong_to_correct, na.rm = TRUE) - sum(correct_to_wrong, na.rm = TRUE)
      n_refined <- sum(selected)
      rows[[length(rows) + 1]] <- data.frame(
        Score = score_name,
        BudgetFraction = fraction,
        NRefined = n_refined,
        WrongToCorrect = sum(wrong_to_correct, na.rm = TRUE),
        CorrectToWrong = sum(correct_to_wrong, na.rm = TRUE),
        NetCorrections = net,
        CorrectionEfficiency = if (n_refined > 0) net / n_refined else NA_real_,
        RecoveryFraction = if (full_net > 0) net / full_net else NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }
  bind_rows_base(rows)
}

fit_meta_benefit_model <- function(features, difficulty) {
  if (nrow(features) == 0 || !"RefinementBeneficial" %in% names(features)) {
    return(list(coefficients = data.frame(), predictions = data.frame()))
  }
  x <- merge(features, difficulty[c("Dataset", "IntrinsicDifficultyScore", "ObservedDifficultyScore")],
             by = "Dataset", all.x = TRUE)
  y <- as_flag_meta(x$RefinementBeneficial)
  feature_cols <- intersect(
    c(
      "LLMConfidence", "OntologyEvidenceScore", "MarkerEvidenceScore",
      "BestMarkerEvidenceScore", "TissueEvidenceScore", "ConsensusEvidenceScore",
      "EvidenceConflictScore", "MarkerMargin", "RequiresRefinementNumeric",
      "CandidateEvidenceDisagreement", "IntrinsicDifficultyScore",
      "ObservedDifficultyScore"
    ),
    names(x)
  )
  if (length(unique(y[!is.na(y)])) < 2 || length(feature_cols) == 0) {
    return(list(coefficients = data.frame(), predictions = data.frame()))
  }
  model_df <- data.frame(Target = as.integer(y), x[feature_cols], stringsAsFactors = FALSE)
  for (column in feature_cols) {
    model_df[[column]] <- as_number(model_df[[column]])
    model_df[[column]][is.na(model_df[[column]])] <- safe_mean(model_df[[column]])
    if (is.na(model_df[[column]][1])) model_df[[column]] <- 0
  }
  formula <- stats::as.formula(paste("Target ~", paste(feature_cols, collapse = " + ")))
  fit <- tryCatch(stats::glm(formula, data = model_df, family = stats::binomial()), error = function(e) NULL)
  if (is.null(fit)) {
    return(list(coefficients = data.frame(), predictions = data.frame()))
  }
  coef_table <- summary(fit)$coefficients
  coefficients <- data.frame(
    Feature = rownames(coef_table),
    Estimate = coef_table[, "Estimate"],
    StdError = coef_table[, "Std. Error"],
    Z = coef_table[, "z value"],
    PValue = coef_table[, "Pr(>|z|)"],
    OddsRatio = exp(coef_table[, "Estimate"]),
    stringsAsFactors = FALSE
  )
  predictions <- data.frame(
    Dataset = x$Dataset,
    Replicate = x$Replicate %||% NA,
    Cluster = x$Cluster %||% NA,
    RefinementBeneficial = y,
    PredictedBenefitProbability = as.numeric(stats::predict(fit, type = "response")),
    stringsAsFactors = FALSE
  )
  list(coefficients = coefficients, predictions = predictions)
}

plot_meta_outputs <- function(prefix,
                              difficulty,
                              selector_links,
                              pareto,
                              reliability_curves,
                              decision_curves,
                              taxonomy_summary) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    return(invisible(FALSE))
  }

  if (nrow(difficulty) > 0) {
    p <- ggplot2::ggplot(
      difficulty,
      ggplot2::aes(
        x = .data$IntrinsicDifficultyScore,
        y = .data$PlainMacroF1,
        color = .data$DifficultyClass,
        label = .data$Dataset
      )
    ) +
      ggplot2::geom_point(size = 2.5) +
      ggplot2::geom_text(vjust = -0.6, show.legend = FALSE) +
      ggplot2::labs(x = "Intrinsic difficulty score", y = "Plain Macro-F1", color = "Difficulty") +
      ggplot2::theme_minimal(base_size = 11)
    ggplot2::ggsave(paste0(prefix, "_difficulty_vs_plain_macroF1.pdf"), p, width = 7, height = 5)
  }

  if (nrow(selector_links) > 0) {
    keep <- selector_links[selector_links$Selector %in% c("evidence-k", "risk", "risk-k", "confidence-k", "random-k"), , drop = FALSE]
    if (nrow(keep) > 0) {
      p <- ggplot2::ggplot(
        keep,
        ggplot2::aes(
          x = .data$ObservedDifficultyScore,
          y = .data$CorrectionEfficiency,
          color = .data$Selector
        )
      ) +
        ggplot2::geom_point(size = 2) +
        ggplot2::geom_smooth(method = "lm", se = FALSE) +
        ggplot2::coord_cartesian(ylim = c(0, 1)) +
        ggplot2::labs(x = "Observed difficulty score", y = "Correction efficiency", color = "Selector") +
        ggplot2::theme_minimal(base_size = 11)
      ggplot2::ggsave(paste0(prefix, "_difficulty_vs_selector_efficiency.pdf"), p, width = 7, height = 5)
    }
  }

  if (nrow(pareto$annotation) > 0) {
    p <- ggplot2::ggplot(
      pareto$annotation,
      ggplot2::aes(x = .data$CostUSD, y = .data$MacroF1, color = .data$ParetoMacroF1Cost)
    ) +
      ggplot2::geom_point(size = 2) +
      ggplot2::labs(x = "API cost (USD)", y = "Macro-F1", color = "Pareto front") +
      ggplot2::theme_minimal(base_size = 11)
    ggplot2::ggsave(paste0(prefix, "_pareto_macroF1_cost.pdf"), p, width = 6.5, height = 5)
  }

  if (nrow(reliability_curves) > 0) {
    p <- ggplot2::ggplot(
      reliability_curves,
      ggplot2::aes(x = .data$MeanScore, y = .data$ObservedErrorRate)
    ) +
      ggplot2::geom_point(ggplot2::aes(size = .data$N), color = "#2C7FB8") +
      ggplot2::geom_line(color = "#2C7FB8") +
      ggplot2::facet_wrap(~Axis, scales = "free_x") +
      ggplot2::labs(x = "Mean score", y = "Observed error rate", size = "Clusters") +
      ggplot2::theme_minimal(base_size = 10)
    ggplot2::ggsave(paste0(prefix, "_reliability_curves.pdf"), p, width = 9, height = 7)
  }

  if (nrow(decision_curves) > 0) {
    p <- ggplot2::ggplot(
      decision_curves,
      ggplot2::aes(x = .data$BudgetFraction, y = .data$RecoveryFraction, color = .data$Score)
    ) +
      ggplot2::geom_line() +
      ggplot2::geom_point() +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::labs(x = "Refinement budget fraction", y = "Oracle recovery fraction", color = "Ranking score") +
      ggplot2::theme_minimal(base_size = 11)
    ggplot2::ggsave(paste0(prefix, "_decision_curves.pdf"), p, width = 7, height = 5)
  }

  if (nrow(taxonomy_summary) > 0) {
    p <- ggplot2::ggplot(
      taxonomy_summary,
      ggplot2::aes(x = stats::reorder(.data$ErrorCategory, .data$Fraction), y = .data$Fraction)
    ) +
      ggplot2::geom_col(fill = "#2C7FB8") +
      ggplot2::coord_flip() +
      ggplot2::labs(x = NULL, y = "Fraction of first-pass errors") +
      ggplot2::theme_minimal(base_size = 11)
    ggplot2::ggsave(paste0(prefix, "_biological_error_taxonomy.pdf"), p, width = 7, height = 4.5)
  }

  invisible(TRUE)
}

run_meta_benchmark <- function(results_dir = "results",
                               output_prefix = file.path("results", "meta_benchmark"),
                               reliability_model_path = Sys.getenv("DEEPSEEKCELL_RELIABILITY_MODEL", unset = "")) {
  dataset_summary <- read_optional_csv(file.path(results_dir, "dataset_summary.csv"))
  cluster_summary <- read_optional_csv(file.path(results_dir, "cluster_summary.csv"))
  benchmark_results <- read_optional_csv(file.path(results_dir, "benchmark_results_full.csv"))
  refinement <- read_optional_csv(file.path(results_dir, "ablation_refinement_behavior.csv"))
  confidence <- read_optional_csv(file.path(results_dir, "ablation_confidence_quality.csv"))
  features <- load_feature_table(results_dir)
  features <- add_risk_predictions_if_available(features, reliability_model_path)

  if (nrow(dataset_summary) == 0 || nrow(benchmark_results) == 0) {
    stop("Meta-benchmark requires at least dataset_summary.csv and benchmark_results_full.csv.", call. = FALSE)
  }

  dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)

  difficulty <- build_dataset_difficulty(
    dataset_summary = dataset_summary,
    cluster_summary = cluster_summary,
    benchmark_results = benchmark_results,
    refinement = refinement,
    confidence = confidence,
    features = features
  )
  method_gains <- compute_method_gains(benchmark_results, difficulty)
  selector_links <- compute_selector_difficulty_links(refinement, difficulty)
  pareto <- compute_pareto_tables(benchmark_results, refinement)
  taxonomy <- build_error_taxonomy(features)
  reliability_curves <- build_reliability_curves(features)
  decision_curves <- build_decision_curves(features)
  benefit_model <- fit_meta_benefit_model(features, difficulty)

  utils::write.csv(difficulty, paste0(output_prefix, "_dataset_difficulty.csv"), row.names = FALSE)
  utils::write.csv(method_gains, paste0(output_prefix, "_method_gains_vs_difficulty.csv"), row.names = FALSE)
  utils::write.csv(selector_links, paste0(output_prefix, "_selector_benefit_vs_difficulty.csv"), row.names = FALSE)
  utils::write.csv(pareto$annotation, paste0(output_prefix, "_pareto_annotation.csv"), row.names = FALSE)
  utils::write.csv(pareto$refinement, paste0(output_prefix, "_pareto_refinement.csv"), row.names = FALSE)
  utils::write.csv(taxonomy$clusters, paste0(output_prefix, "_biological_error_taxonomy_clusters.csv"), row.names = FALSE)
  utils::write.csv(taxonomy$summary, paste0(output_prefix, "_biological_error_taxonomy_summary.csv"), row.names = FALSE)
  utils::write.csv(reliability_curves, paste0(output_prefix, "_reliability_curves.csv"), row.names = FALSE)
  utils::write.csv(decision_curves, paste0(output_prefix, "_decision_curves.csv"), row.names = FALSE)
  utils::write.csv(benefit_model$coefficients, paste0(output_prefix, "_benefit_predictor_coefficients.csv"), row.names = FALSE)
  utils::write.csv(benefit_model$predictions, paste0(output_prefix, "_benefit_predictor_predictions.csv"), row.names = FALSE)

  plot_meta_outputs(
    output_prefix,
    difficulty = difficulty,
    selector_links = selector_links,
    pareto = pareto,
    reliability_curves = reliability_curves,
    decision_curves = decision_curves,
    taxonomy_summary = taxonomy$summary
  )

  invisible(list(
    difficulty = difficulty,
    method_gains = method_gains,
    selector_links = selector_links,
    pareto = pareto,
    taxonomy = taxonomy,
    reliability_curves = reliability_curves,
    decision_curves = decision_curves,
    benefit_model = benefit_model
  ))
}

args <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) >= 1) args[1] else "results"
output_prefix <- if (length(args) >= 2) args[2] else file.path(results_dir, "meta_benchmark")
reliability_model_path <- if (length(args) >= 3) args[3] else Sys.getenv("DEEPSEEKCELL_RELIABILITY_MODEL", unset = "")

run_meta_benchmark(results_dir, output_prefix, reliability_model_path)
