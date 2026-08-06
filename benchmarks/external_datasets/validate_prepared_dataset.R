# benchmarks/external_datasets/validate_prepared_dataset.R
#
# Validate prepared external datasets against the locked RDS contract expected
# by benchmarks/run_external_validation.R.

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

as_external_flag <- function(x) {
  if (is.logical(x)) return(x)
  if (is.numeric(x)) return(x != 0)
  tolower(trimws(as.character(x))) %in% c("true", "1", "yes", "y")
}

is_blank_external <- function(x) {
  is.null(x) ||
    length(x) == 0 ||
    all(is.na(x)) ||
    !nzchar(trimws(as.character(x[1])))
}

external_required_registry_columns <- function() {
  c(
    "Dataset", "Domain", "Tissue", "Species", "PreparedRdsPath",
    "GroundTruthColumn", "IsDevelopmentDataset", "ConfirmatoryStatus",
    "AppearedInMarkerProfileDevelopment", "LabelsInformedHarmonizationRules",
    "MarkerListsInspectedBeforeFreeze", "SameStudyInRiskTraining",
    "DonorOverlapKnown", "AnnotationSource", "InclusionDecision",
    "SelectionRationale"
  )
}

validate_external_registry <- function(registry) {
  if (!is.data.frame(registry) || nrow(registry) == 0) {
    stop("registry must be a non-empty data frame.", call. = FALSE)
  }
  missing <- setdiff(external_required_registry_columns(), names(registry))
  if (length(missing) > 0) {
    stop("External dataset registry is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (anyDuplicated(registry$Dataset)) {
    duplicated <- unique(registry$Dataset[duplicated(registry$Dataset)])
    stop("External dataset registry has duplicate Dataset values: ", paste(duplicated, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

validate_external_dataset <- function(x, manifest_row = NULL, strict_confirmatory = TRUE) {
  errors <- character()
  warnings <- character()

  add_error <- function(message) {
    errors <<- c(errors, message)
  }
  add_warning <- function(message) {
    warnings <<- c(warnings, message)
  }

  if (!is.list(x)) {
    add_error("Prepared dataset must be a list.")
  }
  if (!is.list(x$markers)) {
    add_error("markers must be a named list.")
  } else {
    marker_names <- names(x$markers)
    if (is.null(marker_names) || any(!nzchar(trimws(marker_names)))) {
      add_error("markers must have non-empty cluster names.")
    }
    if (anyDuplicated(marker_names)) {
      add_error("markers must not contain duplicated cluster names.")
    }
    if (length(x$markers) == 0) {
      add_error("markers must contain at least one cluster.")
    }
    if (any(lengths(x$markers) == 0)) {
      add_error("every marker cluster must contain at least one gene.")
    }
  }

  if (is.null(x$truth)) {
    add_error("truth must be present.")
  } else {
    truth_names <- names(x$truth)
    if (is.null(truth_names) || any(!nzchar(trimws(truth_names)))) {
      add_error("truth must be a named vector keyed by cluster id.")
    }
    if (anyDuplicated(truth_names)) {
      add_error("truth must not contain duplicated cluster names.")
    }
    if (is.list(x$markers) && !is.null(names(x$markers)) && !setequal(names(x$markers), truth_names)) {
      add_error("names(markers) and names(truth) must match exactly.")
    }
    if (any(is.na(x$truth) | !nzchar(trimws(as.character(x$truth))))) {
      add_error("truth labels must be non-empty after preparation.")
    }
  }

  if (!is.null(x$purity)) {
    if (is.null(names(x$purity))) {
      add_error("purity must be named when supplied.")
    } else if (is.list(x$markers) && !setequal(names(x$markers), names(x$purity))) {
      add_error("names(purity) must match names(markers) when purity is supplied.")
    }
    purity <- suppressWarnings(as.numeric(x$purity))
    if (any(!is.na(purity) & (purity < 0 | purity > 1))) {
      add_error("purity values must be in [0, 1].")
    }
  } else {
    add_warning("purity is missing; benchmark results will lack cluster-purity metadata.")
  }

  if (is_blank_external(x$tissue)) {
    add_error("tissue must be supplied.")
  }
  if (is_blank_external(x$species)) {
    add_error("species must be supplied.")
  }

  metadata <- x$metadata %||% list()
  if (!is.list(metadata)) {
    add_error("metadata must be a list when supplied.")
    metadata <- list()
  }

  if (!is.null(manifest_row)) {
    row <- as.list(manifest_row[1, , drop = FALSE])
    confirmatory <- tolower(as.character(row$ConfirmatoryStatus %||% "")) == "confirmatory"
    included <- tolower(as.character(row$InclusionDecision %||% "")) %in% c("include", "included", "yes", "true")

    if (as_external_flag(row$IsDevelopmentDataset %||% FALSE)) {
      add_error("Dataset is flagged as a development dataset.")
    }
    if (as_external_flag(row$AppearedInMarkerProfileDevelopment %||% FALSE)) {
      add_error("Dataset appeared in marker-profile development.")
    }
    if (as_external_flag(row$SameStudyInRiskTraining %||% FALSE)) {
      add_error("Dataset appears in Risk-k training data.")
    }
    if (as_external_flag(row$MarkerListsInspectedBeforeFreeze %||% FALSE)) {
      add_error("Marker lists were inspected before freezing the external panel.")
    }
    if (isTRUE(confirmatory) && !isTRUE(included)) {
      add_error("Confirmatory datasets must have InclusionDecision = include.")
    }
    if (isTRUE(strict_confirmatory) && isTRUE(confirmatory)) {
      uncertain <- c(
        LabelsInformedHarmonizationRules = as_external_flag(row$LabelsInformedHarmonizationRules %||% FALSE),
        DonorOverlapKnown = as_external_flag(row$DonorOverlapKnown %||% FALSE)
      )
      if (any(uncertain)) {
        add_error(paste0(
          "Confirmatory dataset has leakage/overlap uncertainty: ",
          paste(names(uncertain)[uncertain], collapse = ", ")
        ))
      }
    }
  }

  valid <- length(errors) == 0
  result <- list(
    valid = valid,
    errors = errors,
    warnings = warnings,
    n_clusters = if (is.list(x$markers)) length(x$markers) else NA_integer_,
    n_marker_genes = if (is.list(x$markers)) sum(lengths(x$markers), na.rm = TRUE) else NA_integer_,
    mean_marker_genes = if (is.list(x$markers) && length(x$markers) > 0) mean(lengths(x$markers)) else NA_real_,
    min_purity = if (!is.null(x$purity)) min(as.numeric(x$purity), na.rm = TRUE) else NA_real_,
    mean_purity = if (!is.null(x$purity)) mean(as.numeric(x$purity), na.rm = TRUE) else NA_real_
  )
  class(result) <- "deepseekcell_external_validation_result"
  result
}

validation_result_to_data_frame <- function(result, dataset = NA_character_) {
  data.frame(
    Dataset = dataset,
    Valid = isTRUE(result$valid),
    Errors = paste(result$errors, collapse = " | "),
    Warnings = paste(result$warnings, collapse = " | "),
    NClusters = result$n_clusters %||% NA_integer_,
    NMarkerGenes = result$n_marker_genes %||% NA_integer_,
    MeanMarkerGenes = result$mean_marker_genes %||% NA_real_,
    MinPurity = result$min_purity %||% NA_real_,
    MeanPurity = result$mean_purity %||% NA_real_,
    stringsAsFactors = FALSE
  )
}

should_run_external_dataset_validator <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  any(grepl("validate_prepared_dataset\\.R$", args))
}

if (should_run_external_dataset_validator()) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1) {
    stop("Usage: Rscript validate_prepared_dataset.R prepared_dataset.rds [registry.csv] [dataset_id]", call. = FALSE)
  }
  dataset <- readRDS(args[1])
  manifest_row <- NULL
  dataset_id <- NA_character_
  if (length(args) >= 2 && file.exists(args[2])) {
    registry <- utils::read.csv(args[2], stringsAsFactors = FALSE)
    validate_external_registry(registry)
    dataset_id <- if (length(args) >= 3) args[3] else dataset$dataset_name %||% dataset$metadata$dataset %||% NA_character_
    manifest_row <- registry[registry$Dataset == dataset_id, , drop = FALSE]
    if (nrow(manifest_row) == 0) {
      stop("Dataset not found in registry: ", dataset_id, call. = FALSE)
    }
  }
  result <- validate_external_dataset(dataset, manifest_row)
  print(validation_result_to_data_frame(result, dataset_id))
  if (!isTRUE(result$valid)) quit(status = 1)
}
