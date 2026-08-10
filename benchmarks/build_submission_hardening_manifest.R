#!/usr/bin/env Rscript

# Build a reviewer-facing submission hardening manifest from frozen local
# outputs. This script performs no API calls and does not rerun benchmarks.

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

results_dir <- Sys.getenv("DEEPSEEKCELL_RESULTS_DIR", unset = "results")
paper_dir <- Sys.getenv(
  "DEEPSEEKCELL_PAPER_DIR",
  unset = file.path("paper", "sn-article-template")
)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

run_cmd <- function(command, args = character()) {
  out <- tryCatch(
    system2(command, args, stdout = TRUE, stderr = TRUE),
    error = function(e) NA_character_
  )
  paste(out, collapse = "\n")
}

git_commit <- run_cmd("git", c("rev-parse", "HEAD"))
git_short <- run_cmd("git", c("rev-parse", "--short=12", "HEAD"))
git_status <- run_cmd("git", c("status", "--short"))

md5_file <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  unname(tools::md5sum(path))
}

sha256_file <- function(path) {
  if (!file.exists(path) || !requireNamespace("openssl", quietly = TRUE)) {
    return(NA_character_)
  }
  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  as.character(openssl::sha256(con))
}

file_record <- function(path, role) {
  exists <- file.exists(path)
  info <- if (exists) file.info(path) else data.frame(size = NA_real_, mtime = as.POSIXct(NA))
  data.frame(
    Role = role,
    Path = normalizePath(path, winslash = "/", mustWork = FALSE),
    Exists = exists,
    SizeBytes = if (exists) info$size else NA_real_,
    ModifiedTime = if (exists) as.character(info$mtime) else NA_character_,
    MD5 = md5_file(path),
    SHA256 = sha256_file(path),
    stringsAsFactors = FALSE
  )
}

expected_outputs <- c(
  file.path(results_dir, "benchmark_results_full.csv"),
  file.path(results_dir, "benchmark_results_summary.csv"),
  file.path(results_dir, "ablation_refinement_behavior.csv"),
  file.path(results_dir, "ablation_confidence_quality.csv"),
  file.path(results_dir, "refinement_efficiency_summary.csv"),
  file.path(results_dir, "refinement_efficiency.pdf"),
  file.path(results_dir, "celltypist_native_cell_metrics.csv"),
  file.path(results_dir, "openai_gpt5_3rep_selector_summary.csv"),
  file.path(results_dir, "openai_gpt5_3rep_performance_summary.csv"),
  file.path(results_dir, "openai_gpt5_3rep_hash_audit.csv"),
  file.path(results_dir, "openai_gpt5_fullrefined_harmful_revisions.csv"),
  file.path(results_dir, "cross_model_selector_summary.csv"),
  file.path(results_dir, "cross_model_selector_net_table.csv"),
  file.path(results_dir, "cross_model_selector_main_table.csv"),
  file.path(results_dir, "cross_model_selector_main_table.tex"),
  file.path(results_dir, "cross_model_selector_contrasts.csv"),
  file.path(results_dir, "cross_model_selector_net_efficiency.pdf"),
  file.path(results_dir, "external_validation_results_full.csv"),
  file.path(results_dir, "external_validation_refinement_behavior.csv"),
  file.path(results_dir, "external_validation_meta_analysis_selector_forest.pdf"),
  file.path(results_dir, "biological_validation_summary_by_method.csv"),
  file.path(results_dir, "biological_validation_failure_taxonomy.csv"),
  file.path(results_dir, "submission_reproducibility_report.json"),
  file.path(results_dir, "benchmark_release_manifest.csv"),
  file.path("inst", "extdata", "reliability_spec_v1.0.json"),
  file.path("inst", "extdata", "marker_profiles_v1.0.csv"),
  file.path("inst", "extdata", "marker_aliases_v1.0.csv"),
  file.path(paper_dir, "deepseekcell_bmc_bioinformatics.tex"),
  file.path(paper_dir, "deepseekcell_bmc_refs.bib"),
  file.path(paper_dir, "cross_model_selector_net_efficiency.pdf"),
  file.path(paper_dir, "cross_model_selector_main_table.csv"),
  file.path("paper", "bmc_cover_letter_resubmission.md"),
  file.path(paper_dir, "DeepSeekCell.pdf")
)

manifest <- do.call(rbind, lapply(expected_outputs, function(path) {
  file_record(path, role = ifelse(grepl("^paper", path), "manuscript", "result_or_spec"))
}))

utils::write.csv(
  manifest,
  file.path(results_dir, "submission_hardening_manifest.csv"),
  row.names = FALSE
)

session_path <- file.path(results_dir, "submission_hardening_sessionInfo.txt")
capture.output(utils::sessionInfo(), file = session_path)

scan_claims <- function(path) {
  if (!file.exists(path)) return(data.frame())
  lines <- readLines(path, warn = FALSE)
  terms <- c(
    "superior", "consistently", "universally", "robust", "significantly",
    "generalizes", "model-agnostic", "improves", "outperforms"
  )
  rows <- list()
  for (i in seq_along(lines)) {
    for (term in terms) {
      if (grepl(term, lines[[i]], ignore.case = TRUE, fixed = FALSE)) {
        rows[[length(rows) + 1]] <- data.frame(
          File = normalizePath(path, winslash = "/", mustWork = FALSE),
          Line = i,
          Term = term,
          Text = trimws(lines[[i]]),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(rows) == 0) data.frame() else do.call(rbind, rows)
}

claim_paths <- c(
  file.path(paper_dir, "deepseekcell_bmc_bioinformatics.tex"),
  file.path("paper", "bmc_cover_letter_resubmission.md"),
  "README.md"
)
claim_audit <- do.call(rbind, lapply(claim_paths, scan_claims))
utils::write.csv(
  claim_audit,
  file.path(results_dir, "submission_claim_audit.csv"),
  row.names = FALSE
)

editor_response_audit <- data.frame(
  PriorConcern = c(
    "No clear new algorithmic contribution",
    "Insufficient comparison with LLM-based approaches",
    "CellTypist not evaluated natively",
    "Results insufficient to support benefit",
    "Reproducibility and provenance"
  ),
  SubmissionEvidence = c(
    "Algorithmic selective-refinement formulation, Evidence-k, Risk-k, compute-matched controls",
    "DeepSeek, GPT-5, Llama 3.2, Gemma2-2B and Mistral cross-model synthesis",
    "Native cell-level CellTypist output retained separately with cluster majority vote only as secondary comparison",
    "Paired cached ablation, GPT-5 three-replicate benchmark, external validation, biological validation and failure taxonomy",
    "Frozen reliability specification, first-pass hashes, pricing source, result hashes and submission hardening manifest"
  ),
  PrimaryOutput = c(
    "benchmarks/run_benchmark.R; refinement_efficiency_summary.csv",
    "cross_model_selector_main_table.csv; cross_model_selector_net_efficiency.pdf",
    "celltypist_native_cell_metrics.csv",
    "ablation_refinement_behavior.csv; openai_gpt5_3rep_selector_summary.csv; biological_validation_summary_by_method.csv",
    "inst/extdata/reliability_spec_v1.0.json; submission_hardening_manifest.csv"
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(
  editor_response_audit,
  file.path(results_dir, "bmc_editor_response_audit.csv"),
  row.names = FALSE
)

summary_lines <- c(
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste("Git commit:", git_commit),
  paste("Git short commit:", git_short),
  paste("Working tree tracked changes:", ifelse(nzchar(git_status), "yes", "no")),
  paste("Manifest rows:", nrow(manifest)),
  paste("Missing expected outputs:", sum(!manifest$Exists)),
  paste("Claim-audit rows:", nrow(claim_audit)),
  paste("OpenSSL SHA256 available:", requireNamespace("openssl", quietly = TRUE))
)
writeLines(summary_lines, file.path(results_dir, "submission_hardening_manifest_summary.txt"))

cat(paste(summary_lines, collapse = "\n"), "\n")
