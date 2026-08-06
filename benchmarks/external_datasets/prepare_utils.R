# benchmarks/external_datasets/prepare_utils.R
#
# Shared utilities for converting independent studies into the locked external
# validation RDS contract:
# list(markers, truth, tissue, species, purity, metadata).

external_script_path <- tryCatch(sys.frame(1)$ofile, error = function(e) NA_character_)
if (is.null(external_script_path) || is.na(external_script_path) || !nzchar(external_script_path)) {
  args0 <- commandArgs(trailingOnly = FALSE)
  file_args <- args0[grepl("^--file=", args0)]
  external_script_path <- if (length(file_args) > 0) sub("^--file=", "", file_args[1]) else ""
}
external_repo_root <- if (nzchar(external_script_path)) {
  normalizePath(file.path(dirname(external_script_path), "..", ".."), mustWork = FALSE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}
setwd(external_repo_root)

source(file.path("benchmarks", "external_datasets", "validate_prepared_dataset.R"))

external_default_config <- function() {
  list(
    resolution = as.numeric(Sys.getenv("DEEPSEEKCELL_EXTERNAL_SEURAT_RESOLUTION", unset = "0.5")),
    npcs = as.integer(Sys.getenv("DEEPSEEKCELL_EXTERNAL_NPCS", unset = "50")),
    top_markers = as.integer(Sys.getenv("DEEPSEEKCELL_EXTERNAL_TOP_MARKERS", unset = "25")),
    min_pct = as.numeric(Sys.getenv("DEEPSEEKCELL_EXTERNAL_MIN_PCT", unset = "0.25")),
    logfc_threshold = as.numeric(Sys.getenv("DEEPSEEKCELL_EXTERNAL_LOGFC_THRESHOLD", unset = "0.25")),
    seed = as.integer(Sys.getenv("DEEPSEEKCELL_EXTERNAL_SEED", unset = "100"))
  )
}

read_external_registry <- function(registry_path = file.path("benchmarks", "external_datasets", "registry.csv")) {
  if (!file.exists(registry_path)) {
    stop("External dataset registry not found: ", registry_path, call. = FALSE)
  }
  registry <- utils::read.csv(registry_path, stringsAsFactors = FALSE)
  validate_external_registry(registry)
  registry
}

get_external_registry_row <- function(dataset_id,
                                      registry_path = file.path("benchmarks", "external_datasets", "registry.csv")) {
  registry <- read_external_registry(registry_path)
  row <- registry[registry$Dataset == dataset_id, , drop = FALSE]
  if (nrow(row) != 1) {
    stop("Expected exactly one registry row for dataset: ", dataset_id, call. = FALSE)
  }
  row
}

external_bool <- function(x) {
  as_external_flag(x)
}

external_clean_label <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x)] <- ""
  x
}

is_unannotated_external_label <- function(x) {
  x <- tolower(trimws(as.character(x)))
  is.na(x) |
    !nzchar(x) |
    x %in% c("na", "nan", "none", "unknown", "unassigned", "unannotated", "ambiguous") |
    grepl("unclassified|unresolved|doublet|multiplet|low quality|remove|discard", x)
}

read_label_mapping <- function(path) {
  if (is_blank_external(path) || !file.exists(path)) {
    return(data.frame(
      original_label = character(),
      harmonized_label = character(),
      exclude = logical(),
      reason = character(),
      stringsAsFactors = FALSE
    ))
  }
  mapping <- utils::read.csv(path, stringsAsFactors = FALSE)
  required <- c("original_label", "harmonized_label")
  missing <- setdiff(required, names(mapping))
  if (length(missing) > 0) {
    stop("Label mapping is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (!"exclude" %in% names(mapping)) mapping$exclude <- FALSE
  if (!"reason" %in% names(mapping)) mapping$reason <- ""
  mapping$original_label <- external_clean_label(mapping$original_label)
  mapping$harmonized_label <- external_clean_label(mapping$harmonized_label)
  mapping$exclude <- external_bool(mapping$exclude)
  mapping
}

apply_external_label_mapping <- function(labels, mapping) {
  raw <- external_clean_label(labels)
  out <- raw
  excluded <- is_unannotated_external_label(raw)
  exclusion_reason <- ifelse(excluded, "unannotated_or_low_quality", "")

  if (nrow(mapping) > 0) {
    idx <- match(raw, mapping$original_label)
    mapped <- !is.na(idx)
    out[mapped] <- mapping$harmonized_label[idx[mapped]]
    mapped_excluded <- rep(FALSE, length(raw))
    mapped_excluded[mapped] <- mapping$exclude[idx[mapped]]
    excluded[mapped] <- excluded[mapped] | mapped_excluded[mapped]
    exclusion_reason[mapped_excluded] <- mapping$reason[idx[mapped_excluded]]
  }

  out[is_unannotated_external_label(out)] <- "Unknown"
  excluded <- excluded | out == "Unknown"
  list(
    labels = out,
    excluded = excluded,
    exclusion_reason = exclusion_reason,
    mapping = mapping
  )
}

read_external_table <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("tsv", "txt")) {
    utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  }
}

read_external_matrix <- function(path) {
  if (is_blank_external(path)) return(NULL)
  if (!file.exists(path) && !dir.exists(path)) {
    stop("Matrix/DataPath not found: ", path, call. = FALSE)
  }

  ext <- tolower(tools::file_ext(path))
  if (dir.exists(path)) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("Seurat is required to read 10x directories.", call. = FALSE)
    }
    return(Seurat::Read10X(path))
  }
  if (ext == "rds") {
    return(readRDS(path))
  }
  if (ext %in% c("csv", "tsv", "txt")) {
    table <- read_external_table(path)
    rownames(table) <- table[[1]]
    table[[1]] <- NULL
    return(as.matrix(table))
  }
  if (ext == "h5" && requireNamespace("Seurat", quietly = TRUE)) {
    return(Seurat::Read10X_h5(path))
  }
  if (ext == "h5ad" && requireNamespace("zellkonverter", quietly = TRUE)) {
    return(zellkonverter::readH5AD(path))
  }
  stop(
    "Unsupported input format: ", path,
    ". Supported formats: Seurat/SCE/list RDS, matrix CSV/TSV, 10x directory, 10x h5, or h5ad with zellkonverter.",
    call. = FALSE
  )
}

coerce_external_to_seurat <- function(object, row) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat is required for external dataset preparation.", call. = FALSE)
  }
  if (inherits(object, "Seurat")) {
    return(object)
  }
  if (inherits(object, "SingleCellExperiment")) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      stop("SummarizedExperiment is required to convert SingleCellExperiment inputs.", call. = FALSE)
    }
    assay_name <- row$AssayName %||% "counts"
    counts <- SummarizedExperiment::assay(object, assay_name)
    meta <- as.data.frame(SummarizedExperiment::colData(object))
    return(Seurat::CreateSeuratObject(counts = counts, meta.data = meta))
  }
  if (is.list(object) && all(c("counts", "metadata") %in% names(object))) {
    return(Seurat::CreateSeuratObject(counts = object$counts, meta.data = as.data.frame(object$metadata)))
  }
  if (is.matrix(object) || inherits(object, "Matrix")) {
    metadata <- NULL
    metadata_path <- row$MetadataPath %||% ""
    if (!is_blank_external(metadata_path) && file.exists(metadata_path)) {
      metadata <- read_external_table(metadata_path)
      rownames(metadata) <- metadata[[1]]
      metadata[[1]] <- NULL
    }
    return(Seurat::CreateSeuratObject(counts = object, meta.data = metadata))
  }
  stop("Cannot coerce external input to a Seurat object.", call. = FALSE)
}

load_external_seurat <- function(row) {
  path <- row$DataPath %||% ""
  matrix_path <- row$MatrixPath %||% ""
  object <- if (!is_blank_external(path)) {
    read_external_matrix(path)
  } else if (!is_blank_external(matrix_path)) {
    read_external_matrix(matrix_path)
  } else {
    stop("Registry row must provide DataPath or MatrixPath before preparation.", call. = FALSE)
  }
  coerce_external_to_seurat(object, row)
}

standardize_external_markers <- function(markers, top_n) {
  if (!"gene" %in% names(markers)) {
    markers$gene <- rownames(markers)
  }
  effect_col <- intersect(c("avg_log2FC", "avg_logFC"), names(markers))[1]
  if (is.na(effect_col)) {
    effect_col <- "pct.1"
  }
  markers <- markers[
    !grepl(
      "^MT-|^MTRNR|^RP[SL]|^MALAT|^XIST|^HB[AB]|^RP11-|^CTD-|^AC[0-9]|^AL[0-9]|^LINC",
      markers$gene,
      ignore.case = TRUE
    ),
    ,
    drop = FALSE
  ]
  if ("p_val_adj" %in% names(markers)) {
    markers <- markers[is.na(markers$p_val_adj) | markers$p_val_adj < 0.05, , drop = FALSE]
  }
  markers <- markers[order(markers$cluster, -markers$pct.1, -markers[[effect_col]]), , drop = FALSE]
  markers_list <- split(markers$gene, markers$cluster)
  markers_list <- lapply(markers_list, function(x) head(unique(external_clean_label(x)), top_n))
  markers_list[lengths(markers_list) > 0]
}

process_external_seurat <- function(seu,
                                    truth_vector,
                                    row,
                                    config = external_default_config()) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat is required for external dataset preparation.", call. = FALSE)
  }
  set.seed(config$seed)
  truth_vector <- external_clean_label(truth_vector)
  names(truth_vector) <- names(truth_vector) %||% colnames(seu)
  truth_vector <- truth_vector[colnames(seu)]

  seu <- Seurat::NormalizeData(seu, verbose = FALSE)
  seu <- Seurat::FindVariableFeatures(seu, verbose = FALSE)
  seu <- Seurat::ScaleData(seu, verbose = FALSE)
  seu <- Seurat::RunPCA(seu, npcs = config$npcs, verbose = FALSE)
  dims <- seq_len(min(config$npcs, 30))
  seu <- Seurat::FindNeighbors(seu, dims = dims, verbose = FALSE)
  seu <- Seurat::FindClusters(seu, resolution = config$resolution, verbose = FALSE)

  markers <- Seurat::FindAllMarkers(
    seu,
    only.pos = TRUE,
    min.pct = config$min_pct,
    logfc.threshold = config$logfc_threshold,
    test.use = "wilcox",
    verbose = FALSE
  )
  markers_list <- standardize_external_markers(markers, top_n = config$top_markers)

  seu$external_truth_label <- truth_vector[colnames(seu)]

  cluster_truth <- stats::setNames(vapply(names(markers_list), function(cl) {
    cells <- Seurat::WhichCells(seu, idents = cl)
    labs <- truth_vector[cells]
    labs <- labs[!is.na(labs) & nzchar(labs)]
    if (length(labs) == 0) return("Unknown")
    names(sort(table(labs), decreasing = TRUE))[1]
  }, character(1)), names(markers_list))

  purity <- stats::setNames(vapply(names(markers_list), function(cl) {
    cells <- Seurat::WhichCells(seu, idents = cl)
    labs <- truth_vector[cells]
    labs <- labs[!is.na(labs) & nzchar(labs)]
    if (length(labs) == 0) return(NA_real_)
    max(table(labs)) / length(labs)
  }, numeric(1)), names(markers_list))

  cluster_summary <- data.frame(
    Dataset = row$Dataset[[1]],
    Tissue = row$Tissue[[1]],
    Species = row$Species[[1]],
    Cluster = names(markers_list),
    Truth = as.character(cluster_truth[names(markers_list)]),
    NMarkers = lengths(markers_list),
    ClusterPurity = as.numeric(purity[names(markers_list)]),
    stringsAsFactors = FALSE
  )

  dataset_summary <- data.frame(
    Dataset = row$Dataset[[1]],
    Tissue = row$Tissue[[1]],
    Species = row$Species[[1]],
    NCells = ncol(seu),
    NGenes = nrow(seu),
    NClusters = length(markers_list),
    MeanClusterPurity = mean(purity, na.rm = TRUE),
    MedianClusterPurity = stats::median(purity, na.rm = TRUE),
    MinClusterPurity = min(purity, na.rm = TRUE),
    TopMarkers = config$top_markers,
    SeuratResolution = config$resolution,
    NPCs = config$npcs,
    stringsAsFactors = FALSE
  )

  list(
    markers = markers_list,
    truth = cluster_truth[names(markers_list)],
    tissue = row$Tissue[[1]],
    species = row$Species[[1]],
    purity = purity[names(markers_list)],
    dataset_name = row$Dataset[[1]],
    metadata = as.list(row),
    audit = list(
      cluster_summary = cluster_summary,
      dataset_summary = dataset_summary,
      config = config
    )
  )
}

write_external_audit_html <- function(prepared, audit, output_path) {
  html_escape <- function(x) {
    x <- as.character(x)
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    x <- gsub(">", "&gt;", x, fixed = TRUE)
    x
  }
  metadata <- prepared$metadata %||% list()
  cluster_summary <- prepared$audit$cluster_summary
  dataset_summary <- prepared$audit$dataset_summary
  lines <- c(
    "<!doctype html><html><head><meta charset='utf-8'>",
    "<title>DeepSeekCell External Dataset Audit</title>",
    "<style>body{font-family:Arial,sans-serif;max-width:1100px;margin:2rem auto;line-height:1.4}table{border-collapse:collapse;width:100%;margin:1rem 0}td,th{border:1px solid #ddd;padding:6px}th{background:#f4f4f4;text-align:left}.ok{color:#0a7f2e}.bad{color:#b00020}</style>",
    "</head><body>",
    paste0("<h1>External dataset audit: ", html_escape(metadata$Dataset %||% prepared$dataset_name), "</h1>"),
    "<h2>Study provenance</h2><table>",
    paste0("<tr><th>Field</th><th>Value</th></tr>",
           paste(vapply(names(metadata), function(name) {
             paste0("<tr><td>", html_escape(name), "</td><td>", html_escape(metadata[[name]]), "</td></tr>")
           }, character(1)), collapse = "")),
    "</table>",
    "<h2>Dataset summary</h2><table>",
    paste0("<tr>", paste0("<th>", names(dataset_summary), "</th>", collapse = ""), "</tr>"),
    paste0("<tr>", paste0("<td>", html_escape(dataset_summary[1, ]), "</td>", collapse = ""), "</tr>"),
    "</table>",
    "<h2>Cluster summary</h2><table>",
    paste0("<tr>", paste0("<th>", names(cluster_summary), "</th>", collapse = ""), "</tr>"),
    paste(apply(cluster_summary, 1, function(row) {
      paste0("<tr>", paste0("<td>", html_escape(row), "</td>", collapse = ""), "</tr>")
    }), collapse = ""),
    "</table>",
    "</body></html>"
  )
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, output_path, useBytes = TRUE)
  invisible(output_path)
}

write_prepared_external_dataset <- function(prepared, row, output_dir = file.path("data", "external_prepared")) {
  dataset <- row$Dataset[[1]]
  output_path <- row$PreparedRdsPath[[1]]
  if (is_blank_external(output_path)) {
    output_path <- file.path(output_dir, paste0(dataset, "_prepared.rds"))
  }
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(prepared, output_path)

  audit_dir <- file.path("results", "external_dataset_audits")
  dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
  cluster_csv <- file.path(audit_dir, paste0(dataset, "_cluster_audit.csv"))
  dataset_csv <- file.path(audit_dir, paste0(dataset, "_dataset_audit.csv"))
  hash_csv <- file.path(audit_dir, paste0(dataset, "_hashes.csv"))
  html_path <- file.path(audit_dir, paste0("external_dataset_audit_", dataset, ".html"))

  utils::write.csv(prepared$audit$cluster_summary, cluster_csv, row.names = FALSE)
  utils::write.csv(prepared$audit$dataset_summary, dataset_csv, row.names = FALSE)
  hashes <- data.frame(
    Dataset = dataset,
    File = c(output_path, row$DataPath[[1]] %||% NA_character_, row$LabelMappingPath[[1]] %||% NA_character_),
    Role = c("prepared_rds", "source_data", "label_mapping"),
    Exists = file.exists(c(output_path, row$DataPath[[1]] %||% "", row$LabelMappingPath[[1]] %||% "")),
    MD5 = vapply(
      c(output_path, row$DataPath[[1]] %||% "", row$LabelMappingPath[[1]] %||% ""),
      function(path) if (nzchar(path) && file.exists(path)) unname(tools::md5sum(path)) else NA_character_,
      character(1)
    ),
    stringsAsFactors = FALSE
  )
  utils::write.csv(hashes, hash_csv, row.names = FALSE)
  write_external_audit_html(prepared, prepared$audit, html_path)

  list(
    prepared_rds = output_path,
    cluster_audit_csv = cluster_csv,
    dataset_audit_csv = dataset_csv,
    hash_csv = hash_csv,
    audit_html = html_path
  )
}

prepare_external_dataset_from_registry <- function(dataset_id,
                                                   registry_path = file.path("benchmarks", "external_datasets", "registry.csv"),
                                                   output_dir = file.path("data", "external_prepared"),
                                                   strict_confirmatory = TRUE) {
  row <- get_external_registry_row(dataset_id, registry_path)
  config <- external_default_config()
  seu <- load_external_seurat(row)
  label_column <- row$GroundTruthColumn[[1]]
  if (is_blank_external(label_column) || !label_column %in% names(seu@meta.data)) {
    stop("GroundTruthColumn not found in Seurat metadata for ", dataset_id, ": ", label_column, call. = FALSE)
  }
  raw_labels <- seu@meta.data[[label_column]]
  names(raw_labels) <- rownames(seu@meta.data)
  mapping <- read_label_mapping(row$LabelMappingPath[[1]] %||% "")
  mapped <- apply_external_label_mapping(raw_labels, mapping)
  keep <- !mapped$excluded
  if (sum(keep) == 0) {
    stop("No annotated cells remain after applying exclusion rules for ", dataset_id, call. = FALSE)
  }
  seu <- subset(seu, cells = names(raw_labels)[keep])
  truth <- mapped$labels[keep]
  names(truth) <- names(raw_labels)[keep]

  prepared <- process_external_seurat(seu, truth, row, config)
  prepared$metadata$source_cell_count <- length(raw_labels)
  prepared$metadata$post_filter_cell_count <- length(truth)
  prepared$metadata$excluded_cell_count <- sum(!keep)
  prepared$metadata$excluded_labels <- paste(sort(unique(raw_labels[!keep])), collapse = ";")
  prepared$metadata$label_mapping_path <- row$LabelMappingPath[[1]] %||% ""
  prepared$metadata$frozen_spec_version <- "DeepSeekCell reliability specification v1.0"

  validation <- validate_external_dataset(prepared, row, strict_confirmatory = strict_confirmatory)
  prepared$audit$validation <- validation
  outputs <- write_prepared_external_dataset(prepared, row, output_dir)
  validation_csv <- file.path("results", "external_dataset_audits", paste0(dataset_id, "_validation.csv"))
  utils::write.csv(validation_result_to_data_frame(validation, dataset_id), validation_csv, row.names = FALSE)

  if (!isTRUE(validation$valid)) {
    stop("Prepared dataset failed validation: ", paste(validation$errors, collapse = " | "), call. = FALSE)
  }

  invisible(list(dataset = prepared, outputs = outputs, validation = validation))
}
