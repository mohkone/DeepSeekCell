# benchmarks/external_datasets/prepare_all_external_datasets.R
#
# Prepare all eligible external datasets listed in registry.csv. Rows without a
# local DataPath, MatrixPath, or supported virtual source are skipped rather
# than silently entering the confirmatory panel.

source(file.path("benchmarks", "external_datasets", "prepare_utils.R"))

row_has_input <- function(row) {
  clean_path <- function(x) {
    x <- as.character(x %||% "")
    x[is.na(x)] <- ""
    trimws(x[[1]])
  }
  path <- clean_path(row$DataPath[[1]] %||% "")
  matrix_path <- clean_path(row$MatrixPath[[1]] %||% "")
  if (is_external_virtual_source(path) || is_external_virtual_source(matrix_path)) {
    return(TRUE)
  }
  (nzchar(path) && (file.exists(path) || dir.exists(path))) ||
    (nzchar(matrix_path) && (file.exists(matrix_path) || dir.exists(matrix_path)))
}

registry_row_is_selected <- function(row, include_pending = FALSE) {
  decision <- tolower(trimws(as.character(row$InclusionDecision[[1]] %||% "")))
  status <- tolower(trimws(as.character(row$ConfirmatoryStatus[[1]] %||% "")))
  decision %in% c("include", "included", "yes", "true") ||
    (isTRUE(include_pending) && status %in% c("planned", "exploratory", "confirmatory", "pending"))
}

external_registry_to_manifest <- function(registry_rows, output_path) {
  manifest_cols <- c(
    "Dataset", "StudyAccession", "SourceRepository", "Tissue", "Species",
    "Center", "Laboratory", "Country", "SequencingPlatform", "Chemistry",
    "DiseaseStatus", "Condition", "DonorCount", "CellCount", "ExpectedClusters",
    "DataPath", "PreparedRdsPath", "GroundTruthColumn", "ClusterColumn",
    "MarkerSource", "PreprocessingScript", "IsDevelopmentDataset",
    "IsUnseenTissue", "IsProspectiveDataset", "SelectionRationale", "Notes"
  )
  for (column in manifest_cols) {
    if (!column %in% names(registry_rows)) registry_rows[[column]] <- ""
  }
  manifest <- registry_rows[manifest_cols]
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(manifest, output_path, row.names = FALSE)
  invisible(manifest)
}

prepare_all_external_datasets <- function(registry_path = file.path("benchmarks", "external_datasets", "registry.csv"),
                                          output_dir = file.path("data", "external_prepared"),
                                          manifest_output = file.path("benchmarks", "external_validation_manifest.csv"),
                                          include_pending = FALSE,
                                          strict_confirmatory = TRUE) {
  registry <- read_external_registry(registry_path)
  dataset_filter <- trimws(strsplit(Sys.getenv("DEEPSEEKCELL_EXTERNAL_DATASETS", unset = ""), ",", fixed = TRUE)[[1]])
  dataset_filter <- dataset_filter[nzchar(dataset_filter)]
  if (length(dataset_filter) > 0) {
    missing_filter <- setdiff(dataset_filter, registry$Dataset)
    if (length(missing_filter) > 0) {
      stop(
        "DEEPSEEKCELL_EXTERNAL_DATASETS contains unknown dataset(s): ",
        paste(missing_filter, collapse = ", "),
        call. = FALSE
      )
    }
    registry <- registry[registry$Dataset %in% dataset_filter, , drop = FALSE]
  }
  audit_dir <- file.path("results", "external_dataset_audits")
  dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

  rows <- list()
  successful <- list()

  for (i in seq_len(nrow(registry))) {
    row <- registry[i, , drop = FALSE]
    dataset <- row$Dataset[[1]]

    if (!registry_row_is_selected(row, include_pending = include_pending)) {
      rows[[length(rows) + 1]] <- data.frame(
        Dataset = dataset,
        Status = "skipped_not_selected",
        PreparedRdsPath = row$PreparedRdsPath[[1]],
        Message = "InclusionDecision is not include.",
        stringsAsFactors = FALSE
      )
      next
    }

    if (!row_has_input(row)) {
      rows[[length(rows) + 1]] <- data.frame(
        Dataset = dataset,
        Status = "skipped_missing_input",
        PreparedRdsPath = row$PreparedRdsPath[[1]],
        Message = "Provide DataPath or MatrixPath before preparation.",
        stringsAsFactors = FALSE
      )
      next
    }

    result <- tryCatch(
      prepare_external_dataset_from_registry(
        dataset,
        registry_path = registry_path,
        output_dir = output_dir,
        strict_confirmatory = strict_confirmatory
      ),
      error = function(e) e
    )

    if (inherits(result, "error")) {
      rows[[length(rows) + 1]] <- data.frame(
        Dataset = dataset,
        Status = "failed",
        PreparedRdsPath = row$PreparedRdsPath[[1]],
        Message = conditionMessage(result),
        stringsAsFactors = FALSE
      )
    } else {
      rows[[length(rows) + 1]] <- data.frame(
        Dataset = dataset,
        Status = "prepared",
        PreparedRdsPath = result$outputs$prepared_rds,
        Message = "Prepared and validated.",
        stringsAsFactors = FALSE
      )
      row$PreparedRdsPath <- result$outputs$prepared_rds
      successful[[length(successful) + 1]] <- row
    }
  }

  summary <- do.call(rbind, rows)
  summary_path <- file.path(audit_dir, "prepare_all_external_datasets_summary.csv")
  utils::write.csv(summary, summary_path, row.names = FALSE)

  if (length(successful) > 0) {
    successful_rows <- do.call(rbind, successful)
    external_registry_to_manifest(successful_rows, manifest_output)
  }

  message("Wrote preparation summary to ", summary_path)
  if (length(successful) > 0) {
    message("Wrote external validation manifest to ", manifest_output)
  } else {
    message("No datasets were prepared; external validation manifest was not updated.")
  }

  invisible(list(summary = summary, manifest = if (length(successful) > 0) manifest_output else NULL))
}

if (any(grepl("prepare_all_external_datasets\\.R$", commandArgs(trailingOnly = FALSE)))) {
  args <- commandArgs(trailingOnly = TRUE)
  registry_path <- if (length(args) >= 1) args[1] else file.path("benchmarks", "external_datasets", "registry.csv")
  output_dir <- if (length(args) >= 2) args[2] else file.path("data", "external_prepared")
  manifest_output <- if (length(args) >= 3) args[3] else file.path("benchmarks", "external_validation_manifest.csv")
  include_pending <- if (length(args) >= 4) tolower(args[4]) %in% c("1", "true", "yes", "y") else FALSE
  prepare_all_external_datasets(
    registry_path = registry_path,
    output_dir = output_dir,
    manifest_output = manifest_output,
    include_pending = include_pending
  )
}
