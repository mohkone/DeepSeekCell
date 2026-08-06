# benchmarks/train_reliability_model.R
#
# Train a v1.1 risk-aware reliability model from paired benchmark debug files.
# This script uses development benchmark outputs only. Do not train on held-out
# external-validation datasets.
#
# Typical workflow:
#   1. Run the paired benchmark with FullRefined enabled.
#   2. Rscript benchmarks/train_reliability_model.R results/benchmark_debug results/reliability_model_v1.1_error.rds first_pass_error
#   3. Set DEEPSEEKCELL_RELIABILITY_MODEL to the saved RDS before running Risk-k.

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

read_debug_files <- function(debug_dir, suffix_pattern) {
  files <- list.files(
    debug_dir,
    pattern = suffix_pattern,
    full.names = TRUE,
    recursive = FALSE
  )
  if (length(files) == 0) {
    return(data.frame())
  }
  rows <- lapply(files, function(path) {
    x <- utils::read.csv(path, stringsAsFactors = FALSE)
    x$SourceFile <- basename(path)
    x
  })
  all_names <- unique(unlist(lapply(rows, names)))
  rows <- lapply(rows, function(x) {
    missing <- setdiff(all_names, names(x))
    for (column in missing) {
      x[[column]] <- NA
    }
    x[all_names]
  })
  do.call(rbind, rows)
}

build_training_features <- function(debug_dir = file.path("results", "benchmark_debug")) {
  calibrated <- read_debug_files(debug_dir, "-Calibrated_debug\\.csv$")
  evidence <- read_debug_files(debug_dir, "-Evidence_debug\\.csv$")
  full <- read_debug_files(debug_dir, "-FullRefined_debug\\.csv$")

  base <- if (nrow(calibrated) > 0) calibrated else evidence
  base <- add_training_domain_metadata(base)
  if (nrow(base) == 0) {
    stop(
      "No calibrated/evidence debug files found in ", debug_dir,
      ". Run the paired benchmark first.",
      call. = FALSE
    )
  }

  key_cols <- c("Dataset", "Replicate", "Cluster")
  missing_keys <- setdiff(key_cols, names(base))
  if (length(missing_keys) > 0) {
    stop("Debug files are missing key columns: ", paste(missing_keys, collapse = ", "), call. = FALSE)
  }

  base$InitiallyCorrect <- base$HarmonisedPrediction == base$HarmonisedTruth
  base$FirstPassIncorrect <- !base$InitiallyCorrect

  if (nrow(full) > 0 && all(key_cols %in% names(full))) {
    full <- add_training_domain_metadata(full)
    full$FullRefinedCorrect <- full$HarmonisedPrediction == full$HarmonisedTruth
    full_key <- paste(full$Dataset, full$Replicate, full$Cluster, sep = "|")
    base_key <- paste(base$Dataset, base$Replicate, base$Cluster, sep = "|")
    idx <- match(base_key, full_key)
    base$FullRefinedCorrect <- full$FullRefinedCorrect[idx]
    base$RefinementBeneficial <- !base$InitiallyCorrect & base$FullRefinedCorrect
    base$RefinementHarmful <- base$InitiallyCorrect & !base$FullRefinedCorrect
  }

  feature_df <- extract_reliability_features(base)
  metadata_cols <- intersect(
    c(
      "Dataset", "Replicate", "Cluster", "Method", "LLMBackend", "LLMModelID",
      "Tissue", "Species", "Confidence", "RefinementSelector",
      "RefinementBudgetK", "RequiresRefinement", "EvidenceConflict",
      "RawPrediction", "RawTruth", "HarmonisedPrediction", "HarmonisedTruth",
      "InitiallyCorrect", "FirstPassIncorrect", "FullRefinedCorrect",
      "RefinementBeneficial", "RefinementHarmful", "SourceFile"
    ),
    names(base)
  )
  metadata <- base[metadata_cols]
  if ("Confidence" %in% names(metadata)) {
    metadata$EvidenceAdjustedConfidence <- metadata$Confidence
  }
  cbind(metadata, feature_df[setdiff(names(feature_df), "Cluster")])
}

add_training_domain_metadata <- function(x,
                                         results_path = file.path("results", "benchmark_results_full.csv")) {
  if (!is.data.frame(x) || nrow(x) == 0 || !"Dataset" %in% names(x)) {
    return(x)
  }

  if (!"Tissue" %in% names(x)) {
    x$Tissue <- NA_character_
  }
  if (!"Species" %in% names(x)) {
    x$Species <- NA_character_
  }

  if (file.exists(results_path)) {
    meta <- tryCatch(
      utils::read.csv(results_path, stringsAsFactors = FALSE),
      error = function(e) data.frame()
    )
    if (nrow(meta) > 0 && all(c("Dataset", "Tissue", "Species") %in% names(meta))) {
      meta <- unique(meta[c("Dataset", "Tissue", "Species")])
      idx <- match(x$Dataset, meta$Dataset)
      missing_tissue <- is.na(x$Tissue) | !nzchar(as.character(x$Tissue))
      missing_species <- is.na(x$Species) | !nzchar(as.character(x$Species))
      x$Tissue[missing_tissue] <- meta$Tissue[idx[missing_tissue]]
      x$Species[missing_species] <- meta$Species[idx[missing_species]]
    }
  }

  missing_tissue <- is.na(x$Tissue) | !nzchar(as.character(x$Tissue))
  x$Tissue[missing_tissue] <- infer_tissue_from_dataset(x$Dataset[missing_tissue])

  missing_species <- is.na(x$Species) | !nzchar(as.character(x$Species))
  x$Species[missing_species] <- ifelse(
    grepl("Tasic|Zeisel", x$Dataset[missing_species], ignore.case = TRUE),
    "Mouse",
    "Human"
  )

  x
}

infer_tissue_from_dataset <- function(dataset) {
  dataset <- as.character(dataset)
  out <- rep("Unknown", length(dataset))
  out[grepl("PBMC", dataset, ignore.case = TRUE)] <- "PBMC"
  out[grepl("Pancreas|Baron|Muraro", dataset, ignore.case = TRUE)] <- "Pancreas"
  out[grepl("Brain|Tasic|Zeisel", dataset, ignore.case = TRUE)] <- "Brain"
  out[grepl("Lung|Zilionis", dataset, ignore.case = TRUE)] <- "Lung"
  out
}

write_explainability_outputs <- function(model, features, output_rds) {
  prefix <- sub("\\.rds$", "", output_rds, ignore.case = TRUE)

  importance <- explain_reliability_model(model)
  contributions <- compute_reliability_contributions(model, features)
  global_contributions <- summarise_reliability_contributions(contributions)
  calibration <- evaluate_reliability_calibration(
    model,
    features,
    target = model$target %||% "first_pass_error"
  )
  ablation <- run_reliability_feature_ablation(
    training_features = features,
    evaluation_features = features,
    target = model$target %||% "first_pass_error"
  )

  utils::write.csv(
    importance,
    paste0(prefix, "_feature_importance.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    global_contributions,
    paste0(prefix, "_global_contributions.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    contributions,
    paste0(prefix, "_cluster_contributions.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    calibration$metrics,
    paste0(prefix, "_calibration_metrics.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    calibration$bins,
    paste0(prefix, "_calibration_bins.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    ablation,
    paste0(prefix, "_feature_ablation.csv"),
    row.names = FALSE
  )

  write_explainability_plots(
    importance,
    global_contributions,
    calibration$bins,
    prefix
  )
}

write_explainability_plots <- function(importance,
                                       global_contributions,
                                       calibration_bins,
                                       prefix) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    return(invisible(FALSE))
  }

  if (nrow(importance) > 0) {
    plot_df <- importance[order(abs(importance$StandardizedCoefficient)), , drop = FALSE]
    p <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(
        x = stats::reorder(.data$Feature, .data$StandardizedCoefficient),
        y = .data$StandardizedCoefficient,
        fill = .data$StandardizedCoefficient > 0
      )
    ) +
      ggplot2::geom_col(show.legend = FALSE) +
      ggplot2::coord_flip() +
      ggplot2::labs(
        x = NULL,
        y = "Standardized logistic coefficient",
        title = "Risk-k feature importance"
      ) +
      ggplot2::theme_minimal(base_size = 11)
    ggplot2::ggsave(paste0(prefix, "_feature_importance.pdf"), p, width = 7, height = 5)
  }

  if (nrow(global_contributions) > 0) {
    plot_df <- global_contributions[order(global_contributions$MeanAbsContribution), , drop = FALSE]
    p <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(
        x = stats::reorder(.data$Feature, .data$MeanAbsContribution),
        y = .data$MeanAbsContribution
      )
    ) +
      ggplot2::geom_col(fill = "#2C7FB8") +
      ggplot2::coord_flip() +
      ggplot2::labs(
        x = NULL,
        y = "Mean absolute log-odds contribution",
        title = "Global linear SHAP-style importance"
      ) +
      ggplot2::theme_minimal(base_size = 11)
    ggplot2::ggsave(paste0(prefix, "_contribution_importance.pdf"), p, width = 7, height = 5)
  }

  if (nrow(calibration_bins) > 0) {
    p <- ggplot2::ggplot(
      calibration_bins,
      ggplot2::aes(x = .data$MeanPredictedRisk, y = .data$ObservedEventRate)
    ) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
      ggplot2::geom_point(ggplot2::aes(size = .data$N), color = "#2C7FB8") +
      ggplot2::geom_line(color = "#2C7FB8") +
      ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
      ggplot2::labs(
        x = "Mean predicted risk",
        y = "Observed event rate",
        size = "Clusters",
        title = "Risk-k calibration"
      ) +
      ggplot2::theme_minimal(base_size = 11)
    ggplot2::ggsave(paste0(prefix, "_calibration_curve.pdf"), p, width = 5.5, height = 5)
  }

  invisible(TRUE)
}

train_and_save_reliability_model <- function(debug_dir,
                                             output_rds,
                                             target = c("first_pass_error", "refinement_benefit")) {
  target <- match.arg(target)
  features <- build_training_features(debug_dir)
  model <- train_reliability_model(features, target = target)
  dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)
  saveRDS(model, output_rds)

  feature_csv <- sub("\\.rds$", "_training_features.csv", output_rds, ignore.case = TRUE)
  summary_csv <- sub("\\.rds$", "_summary.csv", output_rds, ignore.case = TRUE)

  utils::write.csv(features, feature_csv, row.names = FALSE)
  utils::write.csv(
    data.frame(
      ModelID = model$model_id,
      Target = model$target,
      ModelType = model$model_type,
      NTrainingRows = model$n_training_rows,
      PositiveRate = model$positive_rate,
      HasFittedModel = !is.null(model$fit),
      stringsAsFactors = FALSE
    ),
    summary_csv,
    row.names = FALSE
  )
  write_explainability_outputs(model, features, output_rds)

  message("Saved reliability model to ", output_rds)
  message("Saved training features to ", feature_csv)
  invisible(model)
}

should_auto_run_reliability_training <- function() {
  value <- tolower(Sys.getenv("DEEPSEEKCELL_RUN_RELIABILITY_TRAINING_ON_SOURCE", unset = "true"))
  value %in% c("1", "true", "yes", "y")
}

if (should_auto_run_reliability_training()) {
  args <- commandArgs(trailingOnly = TRUE)
  debug_dir <- if (length(args) >= 1) args[1] else file.path("results", "benchmark_debug")
  output_rds <- if (length(args) >= 2) args[2] else file.path("results", "reliability_model_v1.1_error.rds")
  target <- if (length(args) >= 3) args[3] else "first_pass_error"

  train_and_save_reliability_model(debug_dir, output_rds, target)
}
