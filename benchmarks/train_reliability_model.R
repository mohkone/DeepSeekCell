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
      "RawPrediction", "RawTruth", "HarmonisedPrediction", "HarmonisedTruth",
      "InitiallyCorrect", "FirstPassIncorrect", "FullRefinedCorrect",
      "RefinementBeneficial", "RefinementHarmful", "SourceFile"
    ),
    names(base)
  )
  metadata <- base[metadata_cols]
  cbind(metadata, feature_df[setdiff(names(feature_df), "Cluster")])
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

  message("Saved reliability model to ", output_rds)
  message("Saved training features to ", feature_csv)
  invisible(model)
}

args <- commandArgs(trailingOnly = TRUE)
debug_dir <- if (length(args) >= 1) args[1] else file.path("results", "benchmark_debug")
output_rds <- if (length(args) >= 2) args[2] else file.path("results", "reliability_model_v1.1_error.rds")
target <- if (length(args) >= 3) args[3] else "first_pass_error"

train_and_save_reliability_model(debug_dir, output_rds, target)
