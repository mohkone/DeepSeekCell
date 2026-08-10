# benchmarks/write_reliability_spec.R
#
# Regenerate the frozen DeepSeekCell reliability specification and marker audit
# tables from the constants in R/refinement.R. This script is intentionally
# deterministic and should be run before, not after, held-out validation.

script_path <- tryCatch(sys.frame(1)$ofile, error = function(e) NA_character_)
if (is.null(script_path) || is.na(script_path) || !nzchar(script_path)) {
  args <- commandArgs(trailingOnly = FALSE)
  file_args <- args[grepl("^--file=", args)]
  script_path <- if (length(file_args) > 0) sub("^--file=", "", file_args[1]) else ""
}
repo_root <- if (nzchar(script_path)) {
  normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}
if (!file.exists(file.path(repo_root, "DESCRIPTION"))) {
  repo_root <- normalizePath(getwd(), mustWork = TRUE)
}
setwd(repo_root)

source("R/api.R")
source("R/utils.R")
source("R/ontology.R")
source("R/refinement.R")

output_dir <- file.path("inst", "extdata")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

package_version <- read.dcf("DESCRIPTION")[1, "Version"]
git_commit <- tryCatch(
  system2("git", c("rev-parse", "--short=12", "HEAD"), stdout = TRUE),
  error = function(e) NA_character_
)

ontology_path <- file.path("data", "cl.obo")
ontology_hash <- if (file.exists(ontology_path)) {
  unname(tools::md5sum(ontology_path))
} else {
  NA_character_
}

spec <- get_reliability_spec(include_marker_profiles = TRUE)
spec$software <- list(
  package = "DeepSeekCell",
  version = unname(package_version),
  commit_recorded_at_freeze = unname(git_commit),
  release_tag = "v0.1.0"
)
spec$ontology <- list(
  name = "Cell Ontology",
  local_path = "data/cl.obo",
  md5 = ontology_hash,
  terms_reported_by_loader = 3437
)
spec$model_defaults <- list(
  deepseek = list(
    model_id = MODELS$deepseek$model_id,
    temperature = MODELS$deepseek$temperature,
    max_tokens = MODELS$deepseek$max_tokens,
    api_url_env = MODELS$deepseek$api_url_env,
    model_id_env = MODELS$deepseek$model_id_env
  ),
  ollama = list(
    model_id = MODELS$ollama$model_id,
    temperature = MODELS$ollama$temperature,
    max_tokens = MODELS$ollama$max_tokens,
    api_url_env = MODELS$ollama$api_url_env,
    model_id_env = MODELS$ollama$model_id_env
  ),
  openai = list(
    model_id = MODELS$openai$model_id,
    temperature = MODELS$openai$temperature,
    max_tokens = MODELS$openai$max_tokens,
    max_tokens_env = MODELS$openai$max_tokens_env,
    reasoning_effort = MODELS$openai$reasoning_effort,
    reasoning_effort_env = MODELS$openai$reasoning_effort_env,
    text_verbosity = MODELS$openai$text_verbosity,
    text_verbosity_env = MODELS$openai$text_verbosity_env,
    input_cost_per_1m_tokens = 1000 * MODELS$openai$input_cost_per_1k,
    output_cost_per_1m_tokens = 1000 * MODELS$openai$output_cost_per_1k,
    pricing_date = MODELS$openai$pricing_date,
    pricing_source = MODELS$openai$pricing_source,
    api_url_env = MODELS$openai$api_url_env,
    model_id_env = MODELS$openai$model_id_env
  )
)
spec$development_benchmark <- list(
  datasets = c(
    "PBMC",
    "BaronPancreas",
    "MuraroPancreas",
    "TasicBrain",
    "ZeiselBrain",
    "ZilionisLung"
  ),
  role = paste(
    "Development benchmark used before the locked external-validation",
    "protocol. Reliability specification v1.0 must not be changed after",
    "external validation begins."
  )
)
spec$locked_external_validation <- list(
  primary_endpoint = "CorrectionEfficiency = (WrongToCorrect - CorrectToWrong) / NRefined",
  primary_hypothesis = paste(
    "Evidence-k has higher correction efficiency than matched Confidence-k",
    "and Random-k on held-out dataset-replicate blocks."
  ),
  paired_design = paste(
    "Generate one first-pass response per model-dataset-replicate block, cache",
    "and hash it, and apply every selector to that exact response."
  ),
  matched_controls = c(
    "NoRefinement",
    "Random-k",
    "Confidence-k",
    "NoOntology-k",
    "Evidence-k",
    "FullRefined"
  ),
  minimum_dataset_plan = c(
    "independent immune dataset",
    "independent pancreas dataset",
    "independent neural dataset",
    "challenging disease or tumour dataset",
    "preferably one previously unseen tissue"
  )
)
spec$sensitivity_plan <- list(
  weight_sets = c(
    "default",
    "equal",
    "no_ontology",
    "no_tissue",
    "marker_dominant",
    "llm_confidence_dominant",
    "simplex_random_500_to_1000"
  ),
  threshold_grid = seq(0.25, 0.70, by = 0.05),
  rule = paste(
    "Sensitivity analyses are secondary robustness checks; do not select a",
    "new best configuration using held-out validation results."
  )
)

jsonlite::write_json(
  spec,
  path = file.path(output_dir, "reliability_spec_v1.0.json"),
  pretty = TRUE,
  auto_unbox = TRUE,
  null = "null"
)

profile_rows <- do.call(
  rbind,
  lapply(names(MARKER_EVIDENCE_PROFILES), function(scope) {
    profiles <- MARKER_EVIDENCE_PROFILES[[scope]]
    do.call(
      rbind,
      lapply(names(profiles), function(cell_type) {
        data.frame(
          spec_id = RELIABILITY_SPEC_ID,
          tissue_scope = scope,
          canonical_cell_type = cell_type,
          positive_markers = paste(profiles[[cell_type]], collapse = ";"),
          aliases = paste(MARKER_PROFILE_ALIASES[[cell_type]] %||% character(), collapse = ";"),
          source_class = "manual_curated_frozen_v1.0",
          developed_on = "development benchmark only",
          date_frozen = RELIABILITY_SPEC_FROZEN_DATE,
          frozen_before_external_validation = TRUE,
          stringsAsFactors = FALSE
        )
      })
    )
  })
)

utils::write.csv(
  profile_rows,
  file.path(output_dir, "marker_profiles_v1.0.csv"),
  row.names = FALSE
)

alias_rows <- do.call(
  rbind,
  lapply(names(MARKER_PROFILE_ALIASES), function(cell_type) {
    data.frame(
      spec_id = RELIABILITY_SPEC_ID,
      canonical_cell_type = cell_type,
      alias = MARKER_PROFILE_ALIASES[[cell_type]],
      date_frozen = RELIABILITY_SPEC_FROZEN_DATE,
      stringsAsFactors = FALSE
    )
  })
)

utils::write.csv(
  alias_rows,
  file.path(output_dir, "marker_aliases_v1.0.csv"),
  row.names = FALSE
)

message("Wrote frozen reliability specification to ", output_dir)
