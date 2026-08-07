# benchmarks/reproduce_submission_bundle.R
#
# Rebuild reviewer-facing benchmark summaries from cached outputs. This script
# intentionally does not call any LLM or rerun expensive annotation workflows.
# It validates the expected cached artifacts, regenerates statistical summaries,
# robustness outputs, universal transfer summaries, and the release manifest,
# then writes a machine-readable reproducibility report.
#
# Usage:
#   Rscript benchmarks/reproduce_submission_bundle.R

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

if (!file.exists(file.path(repo_root, "DESCRIPTION"))) {
  script_path <- tryCatch(sys.frame(1)$ofile, error = function(e) "")
  candidate <- normalizePath(file.path(dirname(script_path %||% "."), ".."), winslash = "/", mustWork = FALSE)
  if (file.exists(file.path(candidate, "DESCRIPTION"))) {
    repo_root <- candidate
    setwd(repo_root)
  }
}

timestamp <- function() {
  format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
}

relative_path <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  root <- normalizePath(repo_root, winslash = "/", mustWork = TRUE)
  sub(paste0("^", gsub("([\\^$.|?*+(){}\\[\\]])", "\\\\\\1", root), "/?"), "", path)
}

file_status <- function(paths) {
  paths <- unique(paths)
  data.frame(
    RelativePath = vapply(paths, relative_path, character(1)),
    Exists = file.exists(paths),
    FileSizeBytes = ifelse(file.exists(paths), as.numeric(file.info(paths)$size), NA_real_),
    MD5 = ifelse(file.exists(paths), unname(tools::md5sum(paths)), NA_character_),
    stringsAsFactors = FALSE
  )
}

run_step <- function(name, command, args, env = character()) {
  started <- timestamp()
  message("Running ", name, "...")
  old_env <- character()
  env_names <- character()
  if (length(env) > 0) {
    split_env <- strsplit(env, "=", fixed = TRUE)
    env_names <- vapply(split_env, `[[`, character(1), 1)
    env_values <- vapply(split_env, function(x) paste(x[-1], collapse = "="), character(1))
    old_env <- Sys.getenv(env_names, unset = NA_character_)
    names(old_env) <- env_names
    do.call(Sys.setenv, stats::setNames(as.list(env_values), env_names))
  }
  on.exit({
    if (length(env_names) > 0) {
      for (env_name in env_names) {
        if (is.na(old_env[[env_name]])) {
          Sys.unsetenv(env_name)
        } else {
          do.call(Sys.setenv, stats::setNames(list(old_env[[env_name]]), env_name))
        }
      }
    }
  }, add = TRUE)
  status <- tryCatch(
    system2(command, args = args, stdout = TRUE, stderr = TRUE),
    error = function(e) structure(conditionMessage(e), status = 1L)
  )
  exit_status <- attr(status, "status")
  if (is.null(exit_status)) {
    exit_status <- 0L
  }
  data.frame(
    Step = name,
    Command = paste(c(command, args), collapse = " "),
    StartedAt = started,
    FinishedAt = timestamp(),
    ExitStatus = as.integer(exit_status),
    Output = paste(status, collapse = "\n"),
    stringsAsFactors = FALSE
  )
}

write_json_report <- function(report, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::write_json(report, path, pretty = TRUE, auto_unbox = TRUE, null = "null")
  } else {
    saveRDS(report, sub("\\.json$", ".rds", path))
  }
}

run_reproducibility_bundle <- function(output_prefix = file.path("results", "submission_reproducibility"),
                                       run_tests = tolower(Sys.getenv("DEEPSEEKCELL_REPRODUCE_RUN_TESTS", "false")) %in% c("1", "true", "yes", "y")) {
  rscript <- file.path(R.home("bin"), "Rscript")
  required_inputs <- c(
    "results/external_validation_results_full.csv",
    "results/external_validation_refinement_behavior.csv",
    "results/external_validation_confidence_quality.csv",
    "results/benchmark_debug",
    "results/reliability_model_v1.1_error.rds",
    "benchmarks/external_validation_manifest.csv",
    "inst/extdata/reliability_spec_v1.0.json",
    "inst/extdata/reliability_model_spec_v1.1.json"
  )
  input_status <- file_status(required_inputs)
  missing <- input_status$RelativePath[!input_status$Exists]
  if (length(missing) > 0) {
    stop(
      "Cannot rebuild submission bundle because required artifact(s) are missing: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  steps <- list()
  if (isTRUE(run_tests)) {
    steps[[length(steps) + 1]] <- run_step(
      "devtools_tests",
      rscript,
      c("-e", "devtools::test()")
    )
  }
  steps[[length(steps) + 1]] <- run_step(
    "external_meta_analysis",
    rscript,
    c(
      "benchmarks/analyse_external_validation_meta_analysis.R",
      "results/external_validation_refinement_behavior.csv",
      "results/external_validation_confidence_quality.csv",
      "results/external_validation_meta_analysis"
    )
  )
  steps[[length(steps) + 1]] <- run_step(
    "external_robustness",
    rscript,
    c(
      "benchmarks/analyse_external_validation_robustness.R",
      "results/external_validation_results_full.csv",
      "results/external_validation_refinement_behavior.csv",
      "results/external_validation_confidence_quality.csv",
      "results/external_validation_robustness"
    )
  )
  steps[[length(steps) + 1]] <- run_step(
    "external_selector_inference",
    rscript,
    c(
      "benchmarks/analyse_selector_inference.R",
      "results/external_validation_refinement_behavior.csv",
      "results/external_validation_selector_inference"
    ),
    env = c(
      "DEEPSEEKCELL_BOOTSTRAP_ITER=5000",
      "DEEPSEEKCELL_PERMUTATION_ITER=5000"
    )
  )
  steps[[length(steps) + 1]] <- run_step(
    "biological_validation",
    rscript,
    c(
      "benchmarks/analyse_biological_validation.R",
      "results/benchmark_debug",
      "benchmarks/external_validation_manifest.csv",
      "results/biological_validation"
    ),
    env = c("DEEPSEEKCELL_ENABLE_GO_ENRICHMENT=false")
  )
  steps[[length(steps) + 1]] <- run_step(
    "universal_reliability_transfer",
    rscript,
    c(
      "benchmarks/analyse_universal_reliability_transfer.R",
      "results/benchmark_debug",
      "results/reliability_model_v1.1_error.rds",
      "results/universal_reliability_transfer",
      "deepseek-chat"
    ),
    env = c(
      "DEEPSEEKCELL_UNIVERSAL_BUDGET_FRACTION=0.2",
      "DEEPSEEKCELL_UNIVERSAL_BUDGET_MIN=1"
    )
  )
  steps[[length(steps) + 1]] <- run_step(
    "benchmark_release_manifest",
    rscript,
    c(
      "benchmarks/build_benchmark_release_manifest.R",
      "results/benchmark_release_manifest.csv",
      "false"
    )
  )

  step_table <- do.call(rbind, steps)
  dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(step_table, paste0(output_prefix, "_steps.csv"), row.names = FALSE)

  expected_outputs <- c(
    "results/external_validation_meta_analysis_selector_pooled_effects.csv",
    "results/external_validation_meta_analysis_selector_contrast_pooled_effects.csv",
    "results/external_validation_robustness_dataset_scorecard.csv",
    "results/external_validation_selector_inference_paired_selector_inference.csv",
    "results/biological_validation_cluster_biological_validation.csv",
    "results/biological_validation_summary_by_method.csv",
    "results/biological_validation_corrected_clusters.csv",
    "results/biological_validation_failure_taxonomy.csv",
    "results/universal_reliability_transfer_risk_transfer_metrics.csv",
    "results/universal_reliability_transfer_selector_summary_by_model.csv",
    "results/benchmark_release_manifest.csv",
    "results/benchmark_release_manifest_summary.csv"
  )
  output_status <- file_status(expected_outputs)
  utils::write.csv(output_status, paste0(output_prefix, "_artifacts.csv"), row.names = FALSE)

  git_commit <- tryCatch(
    system2("git", c("rev-parse", "--short=12", "HEAD"), stdout = TRUE),
    error = function(e) NA_character_
  )
  git_status <- tryCatch(
    system2("git", c("status", "--short"), stdout = TRUE),
    error = function(e) NA_character_
  )

  report <- list(
    created_at = timestamp(),
    repo_root = repo_root,
    git_commit = unname(git_commit),
    git_status_short = git_status,
    r_version = R.version.string,
    platform = R.version$platform,
    run_tests = isTRUE(run_tests),
    all_steps_succeeded = all(step_table$ExitStatus == 0L),
    all_expected_outputs_present = all(output_status$Exists),
    required_inputs = input_status,
    generated_outputs = output_status,
    steps_csv = paste0(output_prefix, "_steps.csv"),
    artifacts_csv = paste0(output_prefix, "_artifacts.csv")
  )
  write_json_report(report, paste0(output_prefix, "_report.json"))

  if (!isTRUE(report$all_steps_succeeded)) {
    failed <- step_table$Step[step_table$ExitStatus != 0L]
    stop("Reproducibility bundle failed at step(s): ", paste(failed, collapse = ", "), call. = FALSE)
  }
  if (!isTRUE(report$all_expected_outputs_present)) {
    missing_outputs <- output_status$RelativePath[!output_status$Exists]
    stop("Expected output(s) were not generated: ", paste(missing_outputs, collapse = ", "), call. = FALSE)
  }

  message("Wrote reproducibility report to ", paste0(output_prefix, "_report.json"))
  invisible(report)
}

args <- commandArgs(trailingOnly = TRUE)
output_prefix <- if (length(args) >= 1) args[[1]] else file.path("results", "submission_reproducibility")
run_reproducibility_bundle(output_prefix)
