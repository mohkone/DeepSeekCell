# benchmarks/build_benchmark_release_manifest.R
#
# Build a machine-readable manifest for benchmark-release archives. The script
# does not copy or publish data. It inventories result tables, plots, debug
# decisions, cached LLM outputs, dataset cache files, frozen specifications, and
# validation plans with file sizes and MD5 hashes.
#
# Usage:
#   Rscript benchmarks/build_benchmark_release_manifest.R results/benchmark_release_manifest.csv

repo_root <- normalizePath(getwd(), mustWork = TRUE)

relative_path <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  root <- normalizePath(repo_root, winslash = "/", mustWork = TRUE)
  sub(paste0("^", gsub("([\\^$.|?*+(){}\\[\\]])", "\\\\\\1", root), "/?"), "", path)
}

file_category <- function(path) {
  rel <- relative_path(path)
  if (grepl("^inst/extdata/", rel)) return("frozen_specification")
  if (grepl("^benchmark_cache/", rel)) return("prepared_dataset_cache")
  if (grepl("^results/benchmark_debug/", rel)) return("paired_debug_decisions")
  if (grepl("^results/.+first_pass", rel, ignore.case = TRUE)) return("cached_llm_response")
  if (grepl("^results/external_validation", rel)) return("external_validation_output")
  if (grepl("^results/reliability_", rel)) return("reliability_model_or_analysis")
  if (grepl("^results/.*\\.pdf$", rel)) return("benchmark_figure")
  if (grepl("^results/.*\\.csv$", rel)) return("benchmark_table")
  if (grepl("^results/.*\\.json$", rel)) return("validation_lock")
  if (grepl("^benchmarks/", rel)) return("benchmark_protocol_or_script")
  "other"
}

release_recommendation <- function(category) {
  category %in% c(
    "frozen_specification",
    "paired_debug_decisions",
    "cached_llm_response",
    "external_validation_output",
    "reliability_model_or_analysis",
    "benchmark_figure",
    "benchmark_table",
    "validation_lock",
    "benchmark_protocol_or_script"
  )
}

list_release_files <- function(include_dataset_cache = FALSE) {
  roots <- c("results", "inst/extdata", "benchmarks")
  files <- unlist(lapply(roots, function(root) {
    if (!dir.exists(root)) return(character())
    list.files(root, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  }))
  files <- files[file.exists(files) & !dir.exists(files)]
  files <- files[!grepl("\\.log$|\\.tmp$", files, ignore.case = TRUE)]

  if (dir.exists("benchmark_cache")) {
    cache_files <- list.files("benchmark_cache", recursive = TRUE, full.names = TRUE)
    cache_files <- cache_files[file.exists(cache_files) & !dir.exists(cache_files)]
    files <- c(files, cache_files)
  }

  categories <- vapply(files, file_category, character(1))
  if (!isTRUE(include_dataset_cache)) {
    files <- files[categories != "prepared_dataset_cache"]
  }
  unique(files)
}

build_release_manifest <- function(output_csv = file.path("results", "benchmark_release_manifest.csv"),
                                   include_dataset_cache = FALSE) {
  files <- list_release_files(include_dataset_cache = include_dataset_cache)
  dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)
  if (length(files) == 0) {
    manifest <- data.frame()
  } else {
    info <- file.info(files)
    category <- vapply(files, file_category, character(1))
    manifest <- data.frame(
      RelativePath = vapply(files, relative_path, character(1)),
      Category = category,
      FileSizeBytes = as.numeric(info$size),
      MD5 = unname(tools::md5sum(files)),
      ModifiedTime = format(info$mtime, "%Y-%m-%dT%H:%M:%S%z"),
      RecommendedForArchive = release_recommendation(category),
      Notes = ifelse(
        category == "prepared_dataset_cache",
        "Large derived RDS cache; include only if redistribution is permitted.",
        ""
      ),
      stringsAsFactors = FALSE
    )
    manifest <- manifest[order(manifest$Category, manifest$RelativePath), , drop = FALSE]
  }

  utils::write.csv(manifest, output_csv, row.names = FALSE)
  summary_csv <- sub("\\.csv$", "_summary.csv", output_csv, ignore.case = TRUE)
  if (nrow(manifest) > 0) {
    summary <- stats::aggregate(
      FileSizeBytes ~ Category + RecommendedForArchive,
      data = manifest,
      FUN = sum
    )
    names(summary)[names(summary) == "FileSizeBytes"] <- "TotalSizeBytes"
    counts <- stats::aggregate(
      RelativePath ~ Category + RecommendedForArchive,
      data = manifest,
      FUN = length
    )
    names(counts)[names(counts) == "RelativePath"] <- "NFiles"
    summary <- merge(summary, counts, by = c("Category", "RecommendedForArchive"), all = TRUE)
  } else {
    summary <- data.frame()
  }
  utils::write.csv(summary, summary_csv, row.names = FALSE)
  message("Wrote benchmark release manifest to ", output_csv)
  message("Wrote benchmark release summary to ", summary_csv)
  invisible(list(manifest = manifest, summary = summary))
}

args <- commandArgs(trailingOnly = TRUE)
output_csv <- if (length(args) >= 1) args[1] else file.path("results", "benchmark_release_manifest.csv")
include_dataset_cache <- if (length(args) >= 2) {
  tolower(args[2]) %in% c("1", "true", "yes", "y")
} else {
  FALSE
}

build_release_manifest(output_csv, include_dataset_cache = include_dataset_cache)
