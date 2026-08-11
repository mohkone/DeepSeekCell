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

normalize_optional <- function(path) {
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

read_csv_if_exists <- function(path) {
  if (!file.exists(path)) return(data.frame())
  tryCatch(
    utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) data.frame()
  )
}

safe_numeric <- function(x) {
  out <- suppressWarnings(as.numeric(x))
  if (length(out) == 0) NA_real_ else out[[1]]
}

get_table_value <- function(path, row_column, row_value, value_column) {
  x <- read_csv_if_exists(path)
  if (nrow(x) == 0 || !(row_column %in% names(x)) || !(value_column %in% names(x))) {
    return(NA_real_)
  }
  rows <- which(x[[row_column]] == row_value)
  if (length(rows) != 1) return(NA_real_)
  safe_numeric(x[[value_column]][[rows]])
}

run_cmd_status <- function(command, args = character(), wd = NULL) {
  old_wd <- getwd()
  if (!is.null(wd)) setwd(wd)
  on.exit(setwd(old_wd), add = TRUE)

  out <- tryCatch(
    system2(command, args, stdout = TRUE, stderr = TRUE),
    error = function(e) {
      structure(conditionMessage(e), status = 1)
    }
  )
  status <- attr(out, "status") %||% 0
  list(status = status, output = paste(out, collapse = "\n"))
}

build_clean_manuscript <- function(tex_path, pdf_path) {
  tex_engine <- Sys.which("pdflatex")
  bibtex_engine <- Sys.which("bibtex")
  clean_build_requested <- tolower(Sys.getenv("DEEPSEEKCELL_CLEAN_BUILD", unset = "false")) %in%
    c("1", "true", "yes", "y")

  if (!clean_build_requested) {
    return(list(
      attempted = FALSE,
      status = "not_requested",
      reason = "Set DEEPSEEKCELL_CLEAN_BUILD=true to delete existing PDF/auxiliary files and rebuild from TeX.",
      tex_engine = tex_engine,
      bibtex_engine = bibtex_engine,
      output = ""
    ))
  }

  if (!file.exists(tex_path)) {
    return(list(
      attempted = FALSE,
      status = "missing_tex",
      reason = "Manuscript .tex file is missing.",
      tex_engine = tex_engine,
      bibtex_engine = bibtex_engine,
      output = ""
    ))
  }

  if (!nzchar(tex_engine)) {
    return(list(
      attempted = FALSE,
      status = "unavailable_tex_engine",
      reason = "pdflatex was not found on PATH; clean rebuild cannot be verified on this machine.",
      tex_engine = tex_engine,
      bibtex_engine = bibtex_engine,
      output = ""
    ))
  }

  base <- tools::file_path_sans_ext(basename(tex_path))
  wd <- dirname(tex_path)
  aux_ext <- c("aux", "bbl", "bcf", "blg", "fdb_latexmk", "fls", "lof", "log", "lot", "out", "run.xml", "toc")
  aux_paths <- file.path(wd, paste(base, aux_ext, sep = "."))
  unlink(c(pdf_path, aux_paths[file.exists(aux_paths)]), force = TRUE)

  commands <- list(
    run_cmd_status(tex_engine, c("-interaction=nonstopmode", "-halt-on-error", basename(tex_path)), wd = wd)
  )
  if (file.exists(file.path(wd, paste0(base, ".aux"))) && nzchar(bibtex_engine)) {
    commands[[length(commands) + 1]] <- run_cmd_status(bibtex_engine, base, wd = wd)
  }
  commands[[length(commands) + 1]] <- run_cmd_status(tex_engine, c("-interaction=nonstopmode", "-halt-on-error", basename(tex_path)), wd = wd)
  commands[[length(commands) + 1]] <- run_cmd_status(tex_engine, c("-interaction=nonstopmode", "-halt-on-error", basename(tex_path)), wd = wd)

  output <- paste(vapply(commands, `[[`, character(1), "output"), collapse = "\n\n---\n\n")
  utils::writeLines(output, file.path(results_dir, "submission_clean_build.log"))

  failed <- any(vapply(commands, function(x) x$status != 0, logical(1)))
  log_path <- file.path(wd, paste0(base, ".log"))
  log_text <- if (file.exists(log_path)) paste(readLines(log_path, warn = FALSE), collapse = "\n") else ""
  unresolved <- grepl(
    "undefined references|undefined citations|Citation .* undefined|Reference .* undefined|Rerun to get cross-references right",
    log_text,
    ignore.case = TRUE
  )

  if (failed) {
    status <- "failed"
    reason <- "At least one LaTeX command returned a non-zero exit status."
  } else if (!file.exists(pdf_path)) {
    status <- "failed"
    reason <- "LaTeX commands completed but the expected PDF was not created."
  } else if (unresolved) {
    status <- "failed_unresolved_references"
    reason <- "LaTeX log still contains unresolved citation/reference warnings after rebuild."
  } else {
    status <- "success"
    reason <- "PDF was rebuilt from a clean LaTeX state."
  }

  list(
    attempted = TRUE,
    status = status,
    reason = reason,
    tex_engine = tex_engine,
    bibtex_engine = bibtex_engine,
    output = output
  )
}

manuscript_tex <- file.path(paper_dir, "deepseekcell_bmc_bioinformatics.tex")
manuscript_pdf <- file.path(paper_dir, "DeepSeekCell.pdf")
build_source_paths <- unique(c(
  file.path(paper_dir, c(
    "deepseekcell_bmc_bioinformatics.tex",
    "deepseekcell_bmc_refs.bib",
    "sn-jnl.cls",
    "cross_model_selector_net_efficiency.pdf",
    "cross_model_selector_main_table.csv",
    "cross_model_selector_main_table.tex",
    "refinement_efficiency.pdf",
    "refinement_efficiency_summary.csv",
    "ablation_refinement_behavior.csv",
    "ablation_confidence_quality.csv",
    "benchmark_results_full.csv",
    "benchmark_results_summary.csv",
    "benchmark_macroF1.pdf",
    "benchmark_accuracy.pdf",
    "benchmark_clade_accuracy.pdf",
    "benchmark_runtime.pdf",
    "celltypist_native_cell_metrics.csv"
  )),
  Sys.glob(file.path(paper_dir, "global_*.pdf")),
  Sys.glob(file.path(paper_dir, "Fig", "*"))
))

build_result <- build_clean_manuscript(manuscript_tex, manuscript_pdf)

existing_source_paths <- build_source_paths[file.exists(build_source_paths)]
source_info <- if (length(existing_source_paths) > 0) file.info(existing_source_paths) else data.frame()
source_max_mtime <- if (length(existing_source_paths) > 0) max(source_info$mtime, na.rm = TRUE) else as.POSIXct(NA)
pdf_exists <- file.exists(manuscript_pdf)
pdf_mtime <- if (pdf_exists) file.info(manuscript_pdf)$mtime else as.POSIXct(NA)
pdf_newer_than_sources <- isTRUE(pdf_exists && length(existing_source_paths) > 0 && pdf_mtime >= source_max_mtime)
stale_source_paths <- if (pdf_exists && length(existing_source_paths) > 0) {
  existing_source_paths[source_info$mtime > pdf_mtime]
} else {
  character()
}
clean_build_verified <- identical(build_result$status, "success") && pdf_newer_than_sources

build_audit <- data.frame(
  ManuscriptTex = normalize_optional(manuscript_tex),
  ManuscriptPdf = normalize_optional(manuscript_pdf),
  TexExists = file.exists(manuscript_tex),
  PdfExists = pdf_exists,
  CleanBuildRequested = tolower(Sys.getenv("DEEPSEEKCELL_CLEAN_BUILD", unset = "false")) %in% c("1", "true", "yes", "y"),
  BuildAttempted = build_result$attempted,
  BuildStatus = build_result$status,
  BuildReason = build_result$reason,
  TeXEngineAvailable = nzchar(build_result$tex_engine),
  TeXEngine = ifelse(nzchar(build_result$tex_engine), normalize_optional(build_result$tex_engine), NA_character_),
  BibTeXEngineAvailable = nzchar(build_result$bibtex_engine),
  BibTeXEngine = ifelse(nzchar(build_result$bibtex_engine), normalize_optional(build_result$bibtex_engine), NA_character_),
  SourceFilesChecked = length(existing_source_paths),
  MissingSourceFiles = paste(normalize_optional(setdiff(build_source_paths, existing_source_paths)), collapse = "; "),
  SourceMaxModifiedTime = ifelse(is.na(source_max_mtime), NA_character_, as.character(source_max_mtime)),
  PdfModifiedTime = ifelse(is.na(pdf_mtime), NA_character_, as.character(pdf_mtime)),
  PdfNewerThanSources = pdf_newer_than_sources,
  CleanBuildVerified = clean_build_verified,
  StaleSourceCount = length(stale_source_paths),
  StaleSourceExamples = paste(normalize_optional(head(stale_source_paths, 8)), collapse = "; "),
  PDF_MD5 = md5_file(manuscript_pdf),
  PDF_SHA256 = sha256_file(manuscript_pdf),
  stringsAsFactors = FALSE
)
utils::write.csv(build_audit, file.path(results_dir, "submission_build_audit.csv"), row.names = FALSE)

visual_audit <- data.frame(
  PdfPath = normalize_optional(manuscript_pdf),
  RenderAttempted = FALSE,
  RenderStatus = if (!pdf_exists) {
    "missing_pdf"
  } else if (!clean_build_verified) {
    "skipped_until_clean_build_verified"
  } else {
    "not_requested"
  },
  PageCount = NA_integer_,
  ManualInspectionRequired = TRUE,
  Notes = if (!clean_build_verified) {
    "Visual PDF audit intentionally deferred until a clean, source-current PDF has been rebuilt."
  } else {
    "Render and inspect every page before submission."
  },
  stringsAsFactors = FALSE
)
utils::write.csv(visual_audit, file.path(results_dir, "submission_pdf_visual_audit.csv"), row.names = FALSE)

numeric_rows <- list()
add_numeric_row <- function(claim_id, claim_text, source_artifact, source_filter,
                            source_column, source_value, expected_value,
                            tolerance = 1e-6, unit = "number",
                            severity = "FAIL") {
  passed <- !is.na(source_value) && abs(source_value - expected_value) <= tolerance
  status <- if (passed) "PASS" else severity
  numeric_rows[[length(numeric_rows) + 1]] <<- data.frame(
    ClaimID = claim_id,
    ClaimText = claim_text,
    SourceArtifact = normalize_optional(source_artifact),
    SourceFilter = source_filter,
    SourceColumn = source_column,
    SourceValue = source_value,
    ExpectedValue = expected_value,
    Tolerance = tolerance,
    Unit = unit,
    Status = status,
    Verified = passed,
    stringsAsFactors = FALSE
  )
}

paper_refinement_summary <- file.path(paper_dir, "refinement_efficiency_summary.csv")
current_refinement_summary <- file.path(results_dir, "refinement_efficiency_summary.csv")
cross_model_main_table <- file.path(results_dir, "cross_model_selector_main_table.csv")
cross_model_table <- read_csv_if_exists(cross_model_main_table)

add_numeric_row(
  "DEEPSEEK_RECOVERY_59_1",
  "Evidence-guided refinement recovered 59.1% of full-refinement corrections.",
  paper_refinement_summary,
  "Method == 'DeepSeekCell-SelfRefined'",
  "RecoveryFraction",
  get_table_value(paper_refinement_summary, "Method", "DeepSeekCell-SelfRefined", "RecoveryFraction"),
  13 / 22,
  unit = "fraction"
)
add_numeric_row(
  "DEEPSEEK_RELATIVE_BUDGET_6_4",
  "Evidence-guided refinement used 6.4% as many refinement calls as full refinement.",
  paper_refinement_summary,
  "Method == 'DeepSeekCell-SelfRefined'",
  "RelativeRefinementBudget",
  get_table_value(paper_refinement_summary, "Method", "DeepSeekCell-SelfRefined", "RelativeRefinementBudget"),
  16 / 249,
  unit = "fraction"
)
add_numeric_row(
  "DEEPSEEK_EVIDENCE_EFFICIENCY_81_3",
  "Evidence-k correction efficiency was 81.3%.",
  paper_refinement_summary,
  "Method == 'DeepSeekCell-SelfRefined'",
  "CorrectionEfficiency",
  get_table_value(paper_refinement_summary, "Method", "DeepSeekCell-SelfRefined", "CorrectionEfficiency"),
  13 / 16,
  unit = "fraction"
)
add_numeric_row(
  "DEEPSEEK_RANDOM_EFFICIENCY_56_3",
  "Random-k correction efficiency was 56.3%.",
  paper_refinement_summary,
  "Method == 'DeepSeekCell-RandomK'",
  "CorrectionEfficiency",
  get_table_value(paper_refinement_summary, "Method", "DeepSeekCell-RandomK", "CorrectionEfficiency"),
  9 / 16,
  unit = "fraction"
)
add_numeric_row(
  "DEEPSEEK_CONFIDENCE_EFFICIENCY_43_8",
  "Confidence-k correction efficiency was 43.8%.",
  paper_refinement_summary,
  "Method == 'DeepSeekCell-ConfidenceK'",
  "CorrectionEfficiency",
  get_table_value(paper_refinement_summary, "Method", "DeepSeekCell-ConfidenceK", "CorrectionEfficiency"),
  7 / 16,
  unit = "fraction"
)
add_numeric_row(
  "DEEPSEEK_FULL_EFFICIENCY_8_8",
  "FullRefined correction efficiency was 8.8%.",
  paper_refinement_summary,
  "Method == 'DeepSeekCell-FullRefined'",
  "CorrectionEfficiency",
  get_table_value(paper_refinement_summary, "Method", "DeepSeekCell-FullRefined", "CorrectionEfficiency"),
  22 / 249,
  unit = "fraction"
)

if (nrow(cross_model_table) > 0) {
  evidence_refined_total <- sum(suppressWarnings(as.numeric(cross_model_table$EvidenceRefined)), na.rm = TRUE)
  wrong_to_correct_total <- sum(suppressWarnings(as.numeric(cross_model_table$WrongToCorrect)), na.rm = TRUE)
  correct_to_wrong_total <- sum(suppressWarnings(as.numeric(cross_model_table$CorrectToWrong)), na.rm = TRUE)
  gpt5_row <- cross_model_table[cross_model_table$Backend == "GPT-5", , drop = FALSE]
} else {
  evidence_refined_total <- wrong_to_correct_total <- correct_to_wrong_total <- NA_real_
  gpt5_row <- data.frame()
}

add_numeric_row(
  "CROSS_MODEL_EVIDENCE_SELECTED_28",
  "Across evaluated backends, Evidence-k selected 28 predictions for refinement.",
  cross_model_main_table,
  "sum(EvidenceRefined)",
  "EvidenceRefined",
  evidence_refined_total,
  28,
  tolerance = 0,
  unit = "count"
)
add_numeric_row(
  "CROSS_MODEL_CORRECTED_25",
  "Across evaluated backends, 25 Evidence-k-selected predictions were corrected.",
  cross_model_main_table,
  "sum(WrongToCorrect)",
  "WrongToCorrect",
  wrong_to_correct_total,
  25,
  tolerance = 0,
  unit = "count"
)
add_numeric_row(
  "CROSS_MODEL_HARMFUL_0",
  "Across evaluated backends, Evidence-k introduced 0 correct-to-wrong changes.",
  cross_model_main_table,
  "sum(CorrectToWrong)",
  "CorrectToWrong",
  correct_to_wrong_total,
  0,
  tolerance = 0,
  unit = "count"
)
add_numeric_row(
  "GPT5_EVIDENCE_SELECTED_4",
  "GPT-5 Evidence-k selected 4 predictions for refinement.",
  cross_model_main_table,
  "Backend == 'GPT-5'",
  "EvidenceRefined",
  if (nrow(gpt5_row) == 1) safe_numeric(gpt5_row$EvidenceRefined) else NA_real_,
  4,
  tolerance = 0,
  unit = "count"
)
add_numeric_row(
  "GPT5_EVIDENCE_CORRECTED_3",
  "GPT-5 Evidence-k corrected 3 predictions.",
  cross_model_main_table,
  "Backend == 'GPT-5'",
  "WrongToCorrect",
  if (nrow(gpt5_row) == 1) safe_numeric(gpt5_row$WrongToCorrect) else NA_real_,
  3,
  tolerance = 0,
  unit = "count"
)
add_numeric_row(
  "GPT5_EVIDENCE_HARMFUL_0",
  "GPT-5 Evidence-k introduced 0 correct-to-wrong changes.",
  cross_model_main_table,
  "Backend == 'GPT-5'",
  "CorrectToWrong",
  if (nrow(gpt5_row) == 1) safe_numeric(gpt5_row$CorrectToWrong) else NA_real_,
  0,
  tolerance = 0,
  unit = "count"
)
add_numeric_row(
  "GPT5_FULL_REFINED_NET_MINUS_2",
  "GPT-5 FullRefined had a net correction count of -2 in the cross-model table.",
  cross_model_main_table,
  "Backend == 'GPT-5'",
  "FullRefinedNet",
  if (nrow(gpt5_row) == 1) safe_numeric(gpt5_row$FullRefinedNet) else NA_real_,
  -2,
  tolerance = 0,
  unit = "count"
)

for (spec in list(
  list("SCOPE_CHECK_SELFREFINED_RECOVERY", "Method == 'DeepSeekCell-SelfRefined'", "RecoveryFraction", 13 / 22),
  list("SCOPE_CHECK_SELFREFINED_BUDGET", "Method == 'DeepSeekCell-SelfRefined'", "RelativeRefinementBudget", 16 / 249),
  list("SCOPE_CHECK_SELFREFINED_EFFICIENCY", "Method == 'DeepSeekCell-SelfRefined'", "CorrectionEfficiency", 13 / 16),
  list("SCOPE_CHECK_RANDOM_EFFICIENCY", "Method == 'DeepSeekCell-RandomK'", "CorrectionEfficiency", 9 / 16),
  list("SCOPE_CHECK_CONFIDENCE_EFFICIENCY", "Method == 'DeepSeekCell-ConfidenceK'", "CorrectionEfficiency", 7 / 16),
  list("SCOPE_CHECK_FULL_EFFICIENCY", "Method == 'DeepSeekCell-FullRefined'", "CorrectionEfficiency", 22 / 249)
)) {
  method <- sub("Method == '([^']+)'.*", "\\1", spec[[2]])
  value <- get_table_value(current_refinement_summary, "Method", method, spec[[3]])
  add_numeric_row(
    spec[[1]],
    "Scope guard: current results/refinement_efficiency_summary.csv should not silently replace the 18-block paper snapshot.",
    current_refinement_summary,
    spec[[2]],
    spec[[3]],
    value,
    spec[[4]],
    tolerance = 1e-6,
    unit = "fraction",
    severity = "WARN"
  )
}

numeric_claim_audit <- if (length(numeric_rows) == 0) data.frame() else do.call(rbind, numeric_rows)
utils::write.csv(
  numeric_claim_audit,
  file.path(results_dir, "submission_numeric_claim_audit.csv"),
  row.names = FALSE
)

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
  file.path(results_dir, "submission_build_audit.csv"),
  file.path(results_dir, "submission_pdf_visual_audit.csv"),
  file.path(results_dir, "submission_numeric_claim_audit.csv"),
  file.path("inst", "extdata", "reliability_spec_v1.0.json"),
  file.path("inst", "extdata", "marker_profiles_v1.0.csv"),
  file.path("inst", "extdata", "marker_aliases_v1.0.csv"),
  file.path(paper_dir, "deepseekcell_bmc_bioinformatics.tex"),
  file.path(paper_dir, "deepseekcell_bmc_refs.bib"),
  file.path(paper_dir, "cross_model_selector_net_efficiency.pdf"),
  file.path(paper_dir, "cross_model_selector_main_table.csv"),
  file.path(paper_dir, "refinement_efficiency_summary.csv"),
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
  paste("Manuscript build status:", build_audit$BuildStatus),
  paste("Manuscript PDF newer than sources:", build_audit$PdfNewerThanSources),
  paste("Manuscript clean build verified:", build_audit$CleanBuildVerified),
  paste("TeX engine available:", build_audit$TeXEngineAvailable),
  paste("PDF visual audit status:", visual_audit$RenderStatus),
  paste("Claim-audit rows:", nrow(claim_audit)),
  paste(
    "Numeric claim-audit PASS/WARN/FAIL:",
    paste(
      as.integer(table(factor(numeric_claim_audit$Status, levels = c("PASS", "WARN", "FAIL")))),
      collapse = "/"
    )
  ),
  paste("OpenSSL SHA256 available:", requireNamespace("openssl", quietly = TRUE))
)
writeLines(summary_lines, file.path(results_dir, "submission_hardening_manifest_summary.txt"))

cat(paste(summary_lines, collapse = "\n"), "\n")
