# benchmarks/analyse_external_validation_robustness.R
#
# Reviewer-facing robustness summaries for locked external validation. This
# script does not train a model, tune thresholds, or define a new selector. It
# stratifies completed validation outputs by study center, sequencing platform,
# disease condition, prospective status, tissue, and dataset scale.
#
# Usage:
#   Rscript benchmarks/analyse_external_validation_robustness.R \
#     results/external_validation_results_full.csv \
#     results/external_validation_refinement_behavior.csv \
#     results/external_validation_confidence_quality.csv \
#     results/external_validation_robustness

as_number <- function(x) {
  suppressWarnings(as.numeric(x))
}

read_csv_required <- function(path, description) {
  if (!file.exists(path)) {
    stop(description, " not found: ", path, call. = FALSE)
  }
  utils::read.csv(path, stringsAsFactors = FALSE)
}

value_or_unrecorded <- function(x) {
  x <- as.character(x)
  x[is.na(x) | !nzchar(trimws(x)) | x %in% c("NA", "NaN")] <- "Unrecorded"
  trimws(x)
}

ensure_column <- function(x, column, default = NA) {
  if (!column %in% names(x)) {
    x[[column]] <- default
  }
  x
}

rbind_nonempty <- function(rows) {
  rows <- rows[vapply(rows, function(x) is.data.frame(x) && nrow(x) > 0, logical(1))]
  if (length(rows) == 0) {
    return(data.frame())
  }
  do.call(rbind, rows)
}

block_id <- function(x) {
  for (column in c("Dataset", "Replicate", "LLMBackend", "LLMModelID")) {
    x <- ensure_column(x, column, "")
  }
  paste(x$Dataset, x$Replicate, x$LLMBackend, x$LLMModelID, sep = "|")
}

method_label <- function(x) {
  if ("RefinementSelector" %in% names(x)) {
    selector <- value_or_unrecorded(x$RefinementSelector)
    selector[selector == "Unrecorded"] <- NA_character_
  } else {
    selector <- rep(NA_character_, nrow(x))
  }
  method <- if ("Method" %in% names(x)) as.character(x$Method) else rep("Unknown", nrow(x))
  out <- ifelse(!is.na(selector) & nzchar(selector), selector, method)
  out <- gsub("^DeepSeekCell-", "", out)
  out <- gsub("^DeepSeek-", "", out)
  out
}

bin_cluster_count <- function(x) {
  x <- as_number(x)
  cut(
    x,
    breaks = c(-Inf, 5, 20, 100, 500, Inf),
    labels = c("1-5", "6-20", "21-100", "101-500", ">500"),
    right = TRUE
  )
}

bin_cell_count <- function(x) {
  x <- as_number(x)
  cut(
    x,
    breaks = c(-Inf, 1000, 10000, 100000, Inf),
    labels = c("<1k", "1k-10k", "10k-100k", ">100k"),
    right = TRUE
  )
}

bin_marker_count <- function(x) {
  x <- as_number(x)
  cut(
    x,
    breaks = c(-Inf, 10, 25, 50, Inf),
    labels = c("<=10", "11-25", "26-50", ">50"),
    right = TRUE
  )
}

prepare_robustness_table <- function(x) {
  for (column in c(
    "Center", "Laboratory", "Country", "SourceRepository",
    "SequencingPlatform", "Chemistry", "DiseaseStatus", "Condition",
    "Tissue", "Species", "IsUnseenTissue", "IsProspectiveDataset"
  )) {
    x <- ensure_column(x, column, "Unrecorded")
    x[[column]] <- value_or_unrecorded(x[[column]])
  }

  if (!"ClusterCount" %in% names(x)) {
    x$ClusterCount <- if ("NClusters" %in% names(x)) x$NClusters else NA_real_
  }
  if (!"MarkerGenesMean" %in% names(x)) {
    x$MarkerGenesMean <- NA_real_
  }
  if (!"CellCount" %in% names(x)) {
    x$CellCount <- NA_real_
  }

  x$Selector <- method_label(x)
  x$BlockID <- block_id(x)
  x$ClusterCountBin <- value_or_unrecorded(as.character(bin_cluster_count(x$ClusterCount)))
  x$CellCountBin <- value_or_unrecorded(as.character(bin_cell_count(x$CellCount)))
  x$MarkerGenesMeanBin <- value_or_unrecorded(as.character(bin_marker_count(x$MarkerGenesMean)))
  x
}

summarise_refinement_group <- function(df, group_column) {
  if (nrow(df) == 0 || !group_column %in% names(df)) {
    return(data.frame())
  }

  split_key <- paste(df[[group_column]], df$Selector, sep = "\r")
  rows <- lapply(split(df, split_key), function(x) {
    n_refined <- sum(as_number(x$NRefined), na.rm = TRUE)
    wrong_to_correct <- sum(as_number(x$WrongToCorrect), na.rm = TRUE)
    correct_to_wrong <- sum(as_number(x$CorrectToWrong), na.rm = TRUE)
    full_wrong_to_correct <- sum(as_number(x$FullWrongToCorrect), na.rm = TRUE)
    full_correct_to_wrong <- sum(as_number(x$FullCorrectToWrong), na.rm = TRUE)
    full_net <- full_wrong_to_correct - full_correct_to_wrong
    net <- wrong_to_correct - correct_to_wrong

    out <- data.frame(
      Grouping = group_column,
      Group = x[[group_column]][[1]],
      Selector = x$Selector[[1]],
      NBlocks = length(unique(x$BlockID)),
      NDatasets = length(unique(x$Dataset)),
      NClusters = sum(as_number(x$NClusters), na.rm = TRUE),
      NRefined = n_refined,
      WrongToCorrect = wrong_to_correct,
      CorrectToWrong = correct_to_wrong,
      NetCorrections = net,
      CorrectionEfficiency = if (n_refined > 0) net / n_refined else NA_real_,
      RecoveryFraction = if (full_net > 0) net / full_net else NA_real_,
      RelativeRefinementBudget = if (sum(as_number(x$NClusters), na.rm = TRUE) > 0) {
        n_refined / sum(as_number(x$NClusters), na.rm = TRUE)
      } else {
        NA_real_
      },
      MeanSelectionPrecision = mean(as_number(x$SelectionPrecision), na.rm = TRUE),
      MeanSelectionRecall = mean(as_number(x$SelectionRecall), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    out
  })

  do.call(rbind, rows)
}

summarise_metric_group <- function(df, group_column, metric_columns, table_name) {
  if (nrow(df) == 0 || !group_column %in% names(df)) {
    return(data.frame())
  }
  metric_columns <- intersect(metric_columns, names(df))
  if (length(metric_columns) == 0) {
    return(data.frame())
  }

  split_key <- paste(df[[group_column]], df$Selector, sep = "\r")
  rows <- lapply(split(df, split_key), function(x) {
    out <- data.frame(
      Table = table_name,
      Grouping = group_column,
      Group = x[[group_column]][[1]],
      Selector = x$Selector[[1]],
      NBlocks = length(unique(x$BlockID)),
      NDatasets = length(unique(x$Dataset)),
      stringsAsFactors = FALSE
    )
    for (metric in metric_columns) {
      out[[paste0(metric, "_mean")]] <- mean(as_number(x[[metric]]), na.rm = TRUE)
      out[[paste0(metric, "_median")]] <- stats::median(as_number(x[[metric]]), na.rm = TRUE)
    }
    out
  })
  do.call(rbind, rows)
}

selector_contrasts <- function(refinement) {
  if (nrow(refinement) == 0) {
    return(data.frame())
  }

  groups <- intersect(
    c("Tissue", "Center", "SequencingPlatform", "DiseaseStatus", "ClusterCountBin"),
    names(refinement)
  )
  rows <- list()

  for (group in groups) {
    summary <- summarise_refinement_group(refinement, group)
    if (nrow(summary) == 0) next
    for (level in unique(summary$Group)) {
      x <- summary[summary$Group == level, , drop = FALSE]
      risk <- x[x$Selector %in% c("risk", "risk-k", "Risk-k", "DeepSeekCell-RiskK"), , drop = FALSE]
      evidence <- x[x$Selector %in% c("evidence", "evidence-k", "self-refined", "SelfRefined", "Evidence-k"), , drop = FALSE]
      confidence <- x[x$Selector %in% c("confidence", "confidence-k", "Confidence-k"), , drop = FALSE]
      random <- x[x$Selector %in% c("random-k", "Random-k"), , drop = FALSE]
      if (nrow(risk) == 0) next
      rows[[length(rows) + 1]] <- data.frame(
        Grouping = group,
        Group = level,
        RiskEfficiency = risk$CorrectionEfficiency[[1]],
        EvidenceEfficiency = if (nrow(evidence) > 0) evidence$CorrectionEfficiency[[1]] else NA_real_,
        ConfidenceEfficiency = if (nrow(confidence) > 0) confidence$CorrectionEfficiency[[1]] else NA_real_,
        RandomEfficiency = if (nrow(random) > 0) random$CorrectionEfficiency[[1]] else NA_real_,
        RiskMinusEvidence = if (nrow(evidence) > 0) {
          risk$CorrectionEfficiency[[1]] - evidence$CorrectionEfficiency[[1]]
        } else {
          NA_real_
        },
        RiskMinusConfidence = if (nrow(confidence) > 0) {
          risk$CorrectionEfficiency[[1]] - confidence$CorrectionEfficiency[[1]]
        } else {
          NA_real_
        },
        RiskMinusRandom = if (nrow(random) > 0) {
          risk$CorrectionEfficiency[[1]] - random$CorrectionEfficiency[[1]]
        } else {
          NA_real_
        },
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(rows) == 0) data.frame() else do.call(rbind, rows)
}

dataset_scorecard <- function(results, refinement) {
  datasets <- unique(c(results$Dataset, refinement$Dataset))
  rows <- lapply(datasets, function(dataset) {
    r <- results[results$Dataset == dataset, , drop = FALSE]
    f <- refinement[refinement$Dataset == dataset, , drop = FALSE]
    meta_source <- if (nrow(r) > 0) r[1, , drop = FALSE] else f[1, , drop = FALSE]
    data.frame(
      Dataset = dataset,
      Tissue = meta_source$Tissue[[1]],
      Center = meta_source$Center[[1]],
      SequencingPlatform = meta_source$SequencingPlatform[[1]],
      DiseaseStatus = meta_source$DiseaseStatus[[1]],
      IsProspectiveDataset = meta_source$IsProspectiveDataset[[1]],
      ClusterCount = meta_source$ClusterCount[[1]],
      CellCount = meta_source$CellCount[[1]],
      MarkerGenesMean = meta_source$MarkerGenesMean[[1]],
      NResultRows = nrow(r),
      NRefinementRows = nrow(f),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

plot_robustness_outputs <- function(prefix, refinement_summary, runtime_summary) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    return(invisible(FALSE))
  }

  platform <- refinement_summary[
    refinement_summary$Grouping == "SequencingPlatform" &
      refinement_summary$Selector %in% c("risk", "Risk-k", "evidence", "self-refined", "random-k", "confidence-k"),
    ,
    drop = FALSE
  ]
  if (nrow(platform) > 0) {
    p <- ggplot2::ggplot(
      platform,
      ggplot2::aes(x = .data$Group, y = .data$CorrectionEfficiency, fill = .data$Selector)
    ) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::labs(x = "Sequencing platform", y = "Correction efficiency", fill = "Selector") +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
    ggplot2::ggsave(paste0(prefix, "_platform_efficiency.pdf"), p, width = 8, height = 5)
  }

  scale <- refinement_summary[
    refinement_summary$Grouping == "ClusterCountBin",
    ,
    drop = FALSE
  ]
  if (nrow(scale) > 0) {
    p <- ggplot2::ggplot(
      scale,
      ggplot2::aes(x = .data$Group, y = .data$CorrectionEfficiency, fill = .data$Selector)
    ) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::labs(x = "Cluster-count bin", y = "Correction efficiency", fill = "Selector") +
      ggplot2::theme_minimal(base_size = 11)
    ggplot2::ggsave(paste0(prefix, "_scale_efficiency.pdf"), p, width = 8, height = 5)
  }

  runtime <- runtime_summary[
    runtime_summary$Grouping == "ClusterCountBin",
    ,
    drop = FALSE
  ]
  if (nrow(runtime) > 0 && "RuntimeSec_mean" %in% names(runtime)) {
    p <- ggplot2::ggplot(
      runtime,
      ggplot2::aes(x = .data$Group, y = .data$RuntimeSec_mean, fill = .data$Selector)
    ) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::labs(x = "Cluster-count bin", y = "Mean runtime (s)", fill = "Selector") +
      ggplot2::theme_minimal(base_size = 11)
    ggplot2::ggsave(paste0(prefix, "_runtime_scaling.pdf"), p, width = 8, height = 5)
  }

  invisible(TRUE)
}

run_external_robustness_analysis <- function(results_csv,
                                             refinement_csv,
                                             confidence_csv,
                                             output_prefix) {
  results <- prepare_robustness_table(read_csv_required(results_csv, "External validation results"))
  refinement <- prepare_robustness_table(read_csv_required(refinement_csv, "External validation refinement behavior"))
  confidence <- prepare_robustness_table(read_csv_required(confidence_csv, "External validation confidence quality"))
  dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)

  groupings <- c(
    "Tissue", "Center", "Laboratory", "Country", "SourceRepository",
    "SequencingPlatform", "Chemistry", "DiseaseStatus", "Condition",
    "IsUnseenTissue", "IsProspectiveDataset", "ClusterCountBin",
    "CellCountBin", "MarkerGenesMeanBin"
  )
  groupings <- unique(groupings)

  refinement_summary <- rbind_nonempty(lapply(
    groupings,
    function(group) summarise_refinement_group(refinement, group)
  ))
  runtime_summary <- rbind_nonempty(lapply(
    groupings,
    function(group) {
      summarise_metric_group(
        results,
        group,
        metric_columns = c(
        "RuntimeSec", "CostUSD", "Tokens", "FirstPassRuntimeSec",
        "RefinementRuntimeSec", "FirstPassTokens", "RefinementTokens",
        "SecondPassCalls", "MacroF1", "Accuracy", "BalancedAcc", "CladeAcc"
        ),
        table_name = "annotation_runtime_performance"
      )
    }
  ))
  confidence_summary <- rbind_nonempty(lapply(
    groupings,
    function(group) {
      summarise_metric_group(
        confidence,
        group,
        metric_columns = c(
        "Brier", "BinaryCorrectnessNLL", "NLL", "ECE", "AUROC", "AURC",
        "MeanConfidence", "Accuracy"
        ),
        table_name = "confidence_quality"
      )
    }
  ))
  contrasts <- selector_contrasts(refinement)
  scorecard <- dataset_scorecard(results, refinement)

  utils::write.csv(scorecard, paste0(output_prefix, "_dataset_scorecard.csv"), row.names = FALSE)
  utils::write.csv(refinement_summary, paste0(output_prefix, "_refinement_by_domain.csv"), row.names = FALSE)
  utils::write.csv(runtime_summary, paste0(output_prefix, "_runtime_cost_by_domain.csv"), row.names = FALSE)
  utils::write.csv(confidence_summary, paste0(output_prefix, "_confidence_by_domain.csv"), row.names = FALSE)
  utils::write.csv(contrasts, paste0(output_prefix, "_selector_contrasts.csv"), row.names = FALSE)

  plot_robustness_outputs(output_prefix, refinement_summary, runtime_summary)

  invisible(list(
    scorecard = scorecard,
    refinement = refinement_summary,
    runtime = runtime_summary,
    confidence = confidence_summary,
    contrasts = contrasts
  ))
}

args <- commandArgs(trailingOnly = TRUE)
results_csv <- if (length(args) >= 1) args[1] else file.path("results", "external_validation_results_full.csv")
refinement_csv <- if (length(args) >= 2) args[2] else file.path("results", "external_validation_refinement_behavior.csv")
confidence_csv <- if (length(args) >= 3) args[3] else file.path("results", "external_validation_confidence_quality.csv")
output_prefix <- if (length(args) >= 4) args[4] else file.path("results", "external_validation_robustness")

run_external_robustness_analysis(results_csv, refinement_csv, confidence_csv, output_prefix)
