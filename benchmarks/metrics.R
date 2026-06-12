# benchmarks/metrics.R

suppressPackageStartupMessages({
  library(mclust)
})

# -----------------------------------------------------------------------------
# Package ontology bridge
# -----------------------------------------------------------------------------

source_deepseekcell_core <- function() {
  r_dir <- normalizePath("R", winslash = "/", mustWork = TRUE)
  core_files <- file.path(
    r_dir,
    c("utils.R", "ontology.R", "ontology_loader.R")
  )

  missing_files <- core_files[!file.exists(core_files)]

  if (length(missing_files) > 0) {
    stop(
      "Could not locate package ontology source files: ",
      paste(missing_files, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(lapply(core_files, source, local = FALSE))
}

validate_benchmark_ontology <- function(ontology) {
  required_fields <- c("dataframe", "ancestor_cache")
  missing_fields <- required_fields[!vapply(required_fields, function(field) {
    !is.null(ontology[[field]])
  }, logical(1))]

  if (length(missing_fields) > 0) {
    stop(
      "Benchmark ontology did not use the DeepSeekCell package loader. ",
      "Missing field(s): ",
      paste(missing_fields, collapse = ", "),
      ". Restart R or re-source benchmarks/run_benchmark.R after updating.",
      call. = FALSE
    )
  }

  if (!is.data.frame(ontology$dataframe) || nrow(ontology$dataframe) == 0) {
    stop(
      "Benchmark ontology has no term dataframe. ",
      "Check data/cl.obo and the package ontology loader.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

load_benchmark_ontology <- function(local_path = "data/cl.obo") {
  source_deepseekcell_core()

  if (!file.exists(local_path)) {
    stop("Cell Ontology file not found: ", local_path, call. = FALSE)
  }

  ontology <- load_cell_ontology(local_path)
  validate_benchmark_ontology(ontology)
  n_terms <- nrow(ontology$dataframe)
  message("Benchmark ontology loaded: ", n_terms, " CL terms")
  ontology
}

# -----------------------------------------------------------------------------
# Label normalization
# -----------------------------------------------------------------------------

normalize_label <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("\\+", " positive ", x)
  x <- gsub("[/\\-]", " ", x)
  x <- gsub("\\(.*?\\)", "", x)
  x <- gsub("[^a-z0-9\\s]", " ", x)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

# -----------------------------------------------------------------------------
# Cell Ontology mapping
# -----------------------------------------------------------------------------

name_to_cl_id <- function(cell_type, ont_data, tissue = NULL) {
  if (
    is.na(cell_type) ||
    cell_type == "Unknown" ||
    !nzchar(cell_type)
  ) {
    return(NA_character_)
  }

  mapped <- tryCatch(
    map_to_cell_ontology(cell_type, ont_data, tissue = tissue),
    error = function(e) NULL
  )

  if (is.null(mapped) || nrow(mapped) == 0 || is.na(mapped$CL_ID[1])) {
    return(NA_character_)
  }

  as.character(mapped$CL_ID[1])
}

ontology_related <- function(pred_cl, true_cl, ont_data) {
  if (is.na(pred_cl) || is.na(true_cl)) {
    return(FALSE)
  }

  are_related_terms(pred_cl, true_cl, ont_data)
}

clade_accuracy <- function(pred_std, true_std, ont_data, tissue = NULL) {
  pred_cl <- vapply(
    pred_std,
    name_to_cl_id,
    character(1),
    ont_data = ont_data,
    tissue = tissue
  )

  true_cl <- vapply(
    true_std,
    name_to_cl_id,
    character(1),
    ont_data = ont_data,
    tissue = tissue
  )

  valid <- !is.na(pred_cl) & !is.na(true_cl)

  if (sum(valid) == 0) {
    return(NA_real_)
  }

  related <- mapply(
    ontology_related,
    pred_cl[valid],
    true_cl[valid],
    MoreArgs = list(ont_data = ont_data)
  )

  mean(related)
}

# -----------------------------------------------------------------------------
# Safe ARI
# -----------------------------------------------------------------------------

safe_ari <- function(pred, truth) {
  pred <- as.character(pred)
  truth <- as.character(truth)

  valid <- !is.na(pred) &
    !is.na(truth) &
    pred != "" &
    truth != ""

  pred <- pred[valid]
  truth <- truth[valid]

  if (length(pred) < 3) {
    return(NA_real_)
  }

  if (length(unique(pred)) < 2 || length(unique(truth)) < 2) {
    return(NA_real_)
  }

  if (
    length(unique(pred)) == length(pred) &&
    length(unique(truth)) == length(truth)
  ) {
    return(mean(pred == truth))
  }

  ari <- tryCatch(
    mclust::adjustedRandIndex(pred, truth),
    error = function(e) NA_real_
  )

  if (is.nan(ari)) {
    return(NA_real_)
  }

  ari
}

# -----------------------------------------------------------------------------
# Evaluation metrics
# -----------------------------------------------------------------------------

evaluate_metrics <- function(pred, truth, ont_data, tissue = NULL) {
  pred <- trimws(as.character(pred))
  truth <- trimws(as.character(truth))

  pred[is.na(pred) | pred == ""] <- "Unknown"
  truth[is.na(truth) | truth == ""] <- "Unknown"

  valid <- truth != "Unknown"

  pred_valid <- pred[valid]
  truth_valid <- truth[valid]

  if (length(pred_valid) < 2) {
    return(c(
      ARI = NA_real_,
      MacroF1 = NA_real_,
      Accuracy = NA_real_,
      BalancedAcc = NA_real_,
      CladeAcc = NA_real_,
      UnknownRate = mean(pred == "Unknown"),
      EvaluatedClusters = length(pred_valid)
    ))
  }

  ari <- safe_ari(pred_valid, truth_valid)
  acc <- mean(pred_valid == truth_valid)
  classes <- unique(truth_valid)

  f1s <- sapply(classes, function(cls) {
    tp <- sum(pred_valid == cls & truth_valid == cls)
    fp <- sum(pred_valid == cls & truth_valid != cls)
    fn <- sum(pred_valid != cls & truth_valid == cls)

    precision <- if ((tp + fp) > 0) tp / (tp + fp) else 0
    recall <- if ((tp + fn) > 0) tp / (tp + fn) else 0

    if ((precision + recall) > 0) {
      2 * precision * recall / (precision + recall)
    } else {
      0
    }
  })

  recalls <- sapply(classes, function(cls) {
    tp <- sum(pred_valid == cls & truth_valid == cls)
    fn <- sum(pred_valid != cls & truth_valid == cls)

    if ((tp + fn) > 0) tp / (tp + fn) else 0
  })

  clade_acc <- clade_accuracy(
    pred_valid,
    truth_valid,
    ont_data,
    tissue = tissue
  )

  metrics <- c(
    ARI = ari,
    MacroF1 = mean(f1s, na.rm = TRUE),
    Accuracy = acc,
    BalancedAcc = mean(recalls, na.rm = TRUE),
    CladeAcc = clade_acc,
    UnknownRate = mean(pred_valid == "Unknown"),
    EvaluatedClusters = length(pred_valid)
  )

  metrics[is.nan(metrics)] <- NA_real_
  metrics
}
