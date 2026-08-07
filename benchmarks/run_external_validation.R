# benchmarks/run_external_validation.R
#
# Preflight and execute the locked held-out validation for DeepSeekCell
# reliability specification v1.0. The default mode is preflight: it validates
# the manifest, writes a validation lock file, and does not spend API credits.
# To execute the benchmark, set DEEPSEEKCELL_RUN_EXTERNAL_VALIDATION=true.

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

source("benchmarks/run_benchmark.R")

SPEC_PATH <- file.path("inst", "extdata", "reliability_spec_v1.0.json")
MANIFEST_TEMPLATE <- file.path("benchmarks", "external_validation_manifest_template.csv")
DEVELOPMENT_DATASETS <- c(
  "PBMC",
  "BaronPancreas",
  "MuraroPancreas",
  "TasicBrain",
  "ZeiselBrain",
  "ZilionisLung"
)

OPTIONAL_MANIFEST_COLUMNS <- c(
  "StudyAccession",
  "SourceRepository",
  "Center",
  "Laboratory",
  "Country",
  "SequencingPlatform",
  "Chemistry",
  "DiseaseStatus",
  "Condition",
  "DonorCount",
  "CellCount",
  "ExpectedClusters",
  "IsProspectiveDataset"
)

load_locked_spec <- function(spec_path = SPEC_PATH) {
  if (!file.exists(spec_path)) {
    stop("Frozen reliability spec is missing: ", spec_path, call. = FALSE)
  }
  spec <- jsonlite::read_json(spec_path, simplifyVector = TRUE)
  if (!identical(spec$spec_id, "DeepSeekCell reliability specification v1.0")) {
    stop("Unexpected reliability spec id: ", spec$spec_id, call. = FALSE)
  }
  spec
}

validate_external_manifest <- function(manifest_path) {
  if (!file.exists(manifest_path)) {
    stop(
      "External validation manifest is missing: ", manifest_path,
      ". Copy ", MANIFEST_TEMPLATE, " and complete it first.",
      call. = FALSE
    )
  }

  manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
  required <- c(
    "Dataset", "Tissue", "Species", "PreparedRdsPath",
    "IsDevelopmentDataset", "IsUnseenTissue", "SelectionRationale"
  )
  missing <- setdiff(required, names(manifest))
  if (length(missing) > 0) {
    stop("Manifest is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  manifest$Dataset <- trimws(manifest$Dataset)
  manifest$Tissue <- trimws(manifest$Tissue)
  manifest$Species <- trimws(manifest$Species)
  manifest$PreparedRdsPath <- trimws(manifest$PreparedRdsPath)
  for (column in OPTIONAL_MANIFEST_COLUMNS) {
    if (!column %in% names(manifest)) {
      manifest[[column]] <- NA_character_
    }
    manifest[[column]] <- trimws(as.character(manifest[[column]]))
  }

  manifest <- manifest[nzchar(manifest$Dataset), , drop = FALSE]
  if (nrow(manifest) == 0) {
    stop("Manifest does not contain any dataset rows.", call. = FALSE)
  }

  dev_flag <- tolower(as.character(manifest$IsDevelopmentDataset)) %in% c("true", "1", "yes", "y")
  known_dev <- manifest$Dataset %in% DEVELOPMENT_DATASETS
  if (any(dev_flag | known_dev)) {
    blocked <- manifest$Dataset[dev_flag | known_dev]
    stop(
      "External validation manifest includes development dataset(s): ",
      paste(blocked, collapse = ", "),
      call. = FALSE
    )
  }

  if (anyDuplicated(manifest$Dataset)) {
    duplicated_names <- unique(manifest$Dataset[duplicated(manifest$Dataset)])
    stop("Duplicate dataset names in manifest: ", paste(duplicated_names, collapse = ", "), call. = FALSE)
  }

  manifest$PreparedRdsExists <- file.exists(manifest$PreparedRdsPath)
  manifest
}

external_manifest_metadata <- function(manifest) {
  keep <- intersect(
    c(
      "Dataset", "StudyAccession", "SourceRepository", "Center", "Laboratory",
      "Country", "SequencingPlatform", "Chemistry", "DiseaseStatus",
      "Condition", "DonorCount", "CellCount", "ExpectedClusters",
      "IsUnseenTissue", "IsProspectiveDataset", "MarkerSource",
      "SelectionRationale", "Notes"
    ),
    names(manifest)
  )
  meta <- manifest[keep]
  meta$ExternalValidationDataset <- TRUE
  meta
}

add_external_validation_metadata <- function(x, manifest) {
  if (!is.data.frame(x) || nrow(x) == 0 || !"Dataset" %in% names(x)) {
    return(x)
  }
  meta <- external_manifest_metadata(manifest)
  idx <- match(x$Dataset, meta$Dataset)
  for (column in setdiff(names(meta), "Dataset")) {
    value <- meta[[column]][idx]
    if (column %in% names(x)) {
      missing <- is.na(x[[column]]) | !nzchar(as.character(x[[column]]))
      x[[column]][missing] <- value[missing]
    } else {
      x[[column]] <- value
    }
  }
  x
}

external_dataset_scale <- function(data) {
  markers <- data$markers %||% list()
  marker_lengths <- lengths(markers)
  data.frame(
    ClusterCount = length(markers),
    MarkerGeneCountTotal = sum(marker_lengths, na.rm = TRUE),
    MarkerGenesMean = if (length(marker_lengths) > 0) mean(marker_lengths, na.rm = TRUE) else NA_real_,
    MarkerGenesMedian = if (length(marker_lengths) > 0) stats::median(marker_lengths, na.rm = TRUE) else NA_real_,
    MarkerGenesMax = if (length(marker_lengths) > 0) max(marker_lengths, na.rm = TRUE) else NA_real_,
    stringsAsFactors = FALSE
  )
}

external_refinement_budget_k <- function(n_clusters) {
  fixed_k <- Sys.getenv("DEEPSEEKCELL_EXTERNAL_REFINEMENT_BUDGET_K", unset = "")
  if (nzchar(fixed_k)) {
    k <- suppressWarnings(as.integer(fixed_k))
    if (!is.na(k)) {
      return(min(max(k, 0), n_clusters))
    }
  }

  fraction <- Sys.getenv("DEEPSEEKCELL_EXTERNAL_REFINEMENT_BUDGET_FRACTION", unset = "")
  if (!nzchar(fraction)) {
    return(NULL)
  }

  fraction <- suppressWarnings(as.numeric(fraction))
  if (is.na(fraction) || fraction <= 0) {
    return(NULL)
  }
  fraction <- min(fraction, 1)
  min_k <- suppressWarnings(as.integer(Sys.getenv(
    "DEEPSEEKCELL_EXTERNAL_REFINEMENT_BUDGET_MIN",
    unset = "1"
  )))
  if (is.na(min_k)) {
    min_k <- 1
  }

  min(max(ceiling(fraction * n_clusters), min_k), n_clusters)
}

add_external_dataset_scale <- function(x, scale) {
  if (!is.data.frame(x) || nrow(x) == 0) {
    return(x)
  }
  for (column in names(scale)) {
    x[[column]] <- scale[[column]][[1]]
  }
  x
}

load_prepared_external_dataset <- function(row) {
  path <- row$PreparedRdsPath[[1]]
  if (!file.exists(path)) {
    stop("PreparedRdsPath does not exist for ", row$Dataset[[1]], ": ", path, call. = FALSE)
  }

  data <- readRDS(path)
  required <- c("markers", "truth")
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop(
      "Prepared RDS for ", row$Dataset[[1]], " is missing: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  data$tissue <- data$tissue %||% row$Tissue[[1]]
  data$species <- data$species %||% row$Species[[1]]
  data$purity <- data$purity %||% NULL
  data$external_metadata <- as.list(row[intersect(names(row), OPTIONAL_MANIFEST_COLUMNS)])

  if (is.null(names(data$markers)) || any(!nzchar(names(data$markers)))) {
    stop("markers must be a named list keyed by cluster id.", call. = FALSE)
  }
  if (is.null(names(data$truth))) {
    stop("truth must be a named character vector keyed by cluster id.", call. = FALSE)
  }
  if (!all(names(data$markers) %in% names(data$truth))) {
    stop("truth labels must cover every marker cluster.", call. = FALSE)
  }

  data
}

write_external_validation_lock <- function(manifest,
                                           spec,
                                           output_path = file.path("results", "external_validation_lock.json")) {
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  git_commit <- tryCatch(
    system2("git", c("rev-parse", "--short=12", "HEAD"), stdout = TRUE),
    error = function(e) NA_character_
  )
  risk_model_path <- Sys.getenv("DEEPSEEKCELL_RELIABILITY_MODEL", unset = "")
  risk_model <- if (nzchar(risk_model_path) && file.exists(risk_model_path)) {
    tryCatch(readRDS(risk_model_path), error = function(e) NULL)
  } else {
    NULL
  }
  lock <- list(
    created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    spec_id = spec$spec_id,
    spec_version = spec$spec_version,
    spec_md5 = unname(tools::md5sum(SPEC_PATH)),
    software_commit = unname(git_commit),
    risk_model_path = if (nzchar(risk_model_path)) risk_model_path else NA_character_,
    risk_model_md5 = if (nzchar(risk_model_path) && file.exists(risk_model_path)) {
      unname(tools::md5sum(risk_model_path))
    } else {
      NA_character_
    },
    risk_model_id = risk_model$model_id %||% NA_character_,
    risk_model_target = risk_model$target %||% NA_character_,
    risk_model_training_rows = risk_model$n_training_rows %||% NA_integer_,
    external_refinement_budget_k = Sys.getenv("DEEPSEEKCELL_EXTERNAL_REFINEMENT_BUDGET_K", unset = NA_character_),
    external_refinement_budget_fraction = Sys.getenv("DEEPSEEKCELL_EXTERNAL_REFINEMENT_BUDGET_FRACTION", unset = NA_character_),
    external_refinement_budget_min = Sys.getenv("DEEPSEEKCELL_EXTERNAL_REFINEMENT_BUDGET_MIN", unset = NA_character_),
    manifest_md5 = unname(tools::md5sum(attr(manifest, "manifest_path"))),
    datasets = manifest$Dataset,
    tissues = manifest$Tissue,
    centers = manifest$Center,
    laboratories = manifest$Laboratory,
    sequencing_platforms = manifest$SequencingPlatform,
    disease_status = manifest$DiseaseStatus,
    prospective_datasets = manifest$IsProspectiveDataset,
    prepared_rds_paths = manifest$PreparedRdsPath,
    prepared_rds_md5 = vapply(
      manifest$PreparedRdsPath,
      function(path) if (file.exists(path)) unname(tools::md5sum(path)) else NA_character_,
      character(1)
    ),
    primary_endpoint = spec$locked_external_validation$primary_endpoint,
    primary_hypothesis = if (!is.null(risk_model)) {
      paste(
        "Risk-k has higher correction efficiency than matched Evidence-k,",
        "Confidence-k and Random-k on held-out dataset-replicate blocks."
      )
    } else {
      spec$locked_external_validation$primary_hypothesis
    },
    primary_comparisons = if (!is.null(risk_model)) {
      c("Risk-k vs Evidence-k", "Risk-k vs Confidence-k", "Risk-k vs Random-k")
    } else {
      c("Evidence-k vs Confidence-k", "Evidence-k vs Random-k")
    }
  )
  jsonlite::write_json(lock, output_path, pretty = TRUE, auto_unbox = TRUE, null = "null")
  invisible(lock)
}

run_locked_external_validation <- function(manifest_path,
                                           n_replicates = 1,
                                           model_key = "deepseek",
                                           api_key = NULL) {
  spec <- load_locked_spec()
  manifest <- validate_external_manifest(manifest_path)
  attr(manifest, "manifest_path") <- manifest_path
  model_config <- MODELS[[model_key]]
  if (is.null(model_config)) {
    stop("Unknown model_key: ", model_key, call. = FALSE)
  }
  if (is.null(api_key) || identical(api_key, "")) {
    for (env_name in model_config$api_key_env %||% character()) {
      value <- Sys.getenv(env_name, unset = "")
      if (nzchar(value)) {
        api_key <- value
        break
      }
    }
  }

  write_external_validation_lock(manifest, spec)
  utils::write.csv(
    manifest,
    file.path("results", "external_validation_plan.csv"),
    row.names = FALSE
  )

  run_flag <- tolower(Sys.getenv("DEEPSEEKCELL_RUN_EXTERNAL_VALIDATION", unset = "false"))
  run_enabled <- run_flag %in% c("1", "true", "yes", "y")
  if (!isTRUE(run_enabled)) {
    message("Preflight complete. Set DEEPSEEKCELL_RUN_EXTERNAL_VALIDATION=true to execute.")
    message("Prepared RDS available for ", sum(manifest$PreparedRdsExists), " of ", nrow(manifest), " dataset(s).")
    return(invisible(list(spec = spec, manifest = manifest)))
  }

  if (
    isTRUE(model_config$requires_api_key) &&
      (is.null(api_key) || identical(api_key, ""))
  ) {
    env_hint <- paste(model_config$api_key_env %||% character(), collapse = ", ")
    stop(
      "API key is required for model: ", model_config$name,
      if (nzchar(env_hint)) paste0(". Set one of: ", env_hint) else "",
      call. = FALSE
    )
  }

  if (!all(manifest$PreparedRdsExists)) {
    missing <- manifest$Dataset[!manifest$PreparedRdsExists]
    stop(
      "Cannot run external validation until PreparedRdsPath exists for: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  ont_data <- load_benchmark_ontology(ONTOLOGY_FILE)
  all_results <- list()
  all_confidence <- list()
  all_reliability <- list()
  all_refinement <- list()

  for (replicate in seq_len(n_replicates)) {
    for (i in seq_len(nrow(manifest))) {
      row <- manifest[i, , drop = FALSE]
      dataset <- load_prepared_external_dataset(row)
      dataset_scale <- external_dataset_scale(dataset)
      refinement_budget_k <- external_refinement_budget_k(length(dataset$markers))
      message("External validation replicate ", replicate, " - ", row$Dataset)
      result <- run_llm_ablation_wrapper(
        dataset_name = row$Dataset,
        data = dataset,
        ont_data = ont_data,
        replicate = replicate,
        model_key = model_key,
        api_key = api_key,
        method_prefix = model_config$name %||% model_key,
        cache_slug = paste0("external_", model_key),
        include_full_refinement = TRUE,
        fixed_refinement_budget_k = refinement_budget_k
      )
      result$results$Replicate <- replicate
      result$results <- add_external_dataset_scale(result$results, dataset_scale)
      result$confidence_quality <- add_external_dataset_scale(result$confidence_quality, dataset_scale)
      result$reliability <- add_external_dataset_scale(result$reliability, dataset_scale)
      result$refinement_behavior <- add_external_dataset_scale(result$refinement_behavior, dataset_scale)

      all_results[[length(all_results) + 1]] <- result$results
      all_confidence[[length(all_confidence) + 1]] <- result$confidence_quality
      all_reliability[[length(all_reliability) + 1]] <- result$reliability
      all_refinement[[length(all_refinement) + 1]] <- result$refinement_behavior
    }
  }

  results <- add_external_validation_metadata(dplyr::bind_rows(all_results), manifest)
  confidence <- add_external_validation_metadata(dplyr::bind_rows(all_confidence), manifest)
  reliability <- add_external_validation_metadata(dplyr::bind_rows(all_reliability), manifest)
  refinement <- add_external_validation_metadata(dplyr::bind_rows(all_refinement), manifest)

  utils::write.csv(results, "results/external_validation_results_full.csv", row.names = FALSE)
  utils::write.csv(confidence, "results/external_validation_confidence_quality.csv", row.names = FALSE)
  utils::write.csv(reliability, "results/external_validation_reliability_bins.csv", row.names = FALSE)
  utils::write.csv(refinement, "results/external_validation_refinement_behavior.csv", row.names = FALSE)

  if (nrow(refinement) > 0) {
    efficiency <- summarise_refinement_efficiency(refinement)
    utils::write.csv(
      efficiency,
      "results/external_validation_refinement_efficiency_summary.csv",
      row.names = FALSE
    )
  }

  invisible(list(
    results = results,
    confidence_quality = confidence,
    reliability = reliability,
    refinement_behavior = refinement
  ))
}

args <- commandArgs(trailingOnly = TRUE)
manifest_path <- if (length(args) >= 1) args[1] else MANIFEST_TEMPLATE
n_replicates <- if (length(args) >= 2) as.integer(args[2]) else 1
model_key <- if (length(args) >= 3) args[3] else "deepseek"

run_locked_external_validation(
  manifest_path = manifest_path,
  n_replicates = n_replicates,
  model_key = model_key
)
