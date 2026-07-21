#!/usr/bin/env Rscript

if (!file.exists("benchmarks/run_benchmark.R")) {
  stop("Run this script from the DeepSeekCell repository root.", call. = FALSE)
}

Sys.setenv(
  DEEPSEEKCELL_RUN_BENCHMARK_ON_SOURCE = "false",
  DEEPSEEKCELL_USE_LLM_CACHE = Sys.getenv("DEEPSEEKCELL_USE_LLM_CACHE", unset = "true")
)

source("benchmarks/run_benchmark.R")

args <- commandArgs(trailingOnly = TRUE)
arg_or_env <- function(index, env_var, default) {
  if (length(args) >= index && nzchar(args[[index]])) {
    return(args[[index]])
  }
  Sys.getenv(env_var, unset = default)
}

dataset_name <- arg_or_env(1, "OLLAMA_PILOT_DATASET", "ZilionisLung")
replicate_id <- as.integer(arg_or_env(2, "OLLAMA_PILOT_REPLICATE", "1"))
model_id <- arg_or_env(3, "OLLAMA_MODEL_ID", MODELS$ollama$model_id)
seed_value <- replicate_id * 100
MODELS$ollama$model_id <- model_id
model_slug <- gsub("[^A-Za-z0-9]+", "_", tolower(sub(":latest$", "", model_id)))
cache_slug <- Sys.getenv("OLLAMA_CACHE_SLUG", unset = paste0("ollama_", model_slug))
method_prefix <- paste0("Ollama_", model_slug)
include_full <- include_ollama_full_refinement()

ont_data <- load_benchmark_ontology(ONTOLOGY_FILE)
dataset <- load_benchmark_dataset(dataset_name, seed = seed_value)

result <- run_llm_ablation_wrapper(
  dataset_name = dataset_name,
  data = dataset,
  ont_data = ont_data,
  replicate = replicate_id,
  model_key = "ollama",
  api_key = NULL,
  method_prefix = method_prefix,
  cache_slug = cache_slug,
  include_full_refinement = include_full
)

dir.create("results", showWarnings = FALSE)
prefix <- paste0("ollama_pilot_", tolower(dataset_name), "_", model_slug)
write.csv(result$results, file.path("results", paste0(prefix, "_results.csv")), row.names = FALSE)
write.csv(
  result$refinement_behavior,
  file.path("results", paste0(prefix, "_refinement_behavior.csv")),
  row.names = FALSE
)

print(result$results[, c(
  "Dataset", "Method", "MacroF1", "Accuracy", "CladeAcc",
  "RefinementBudgetK", "SecondPassCalls", "LLMModelID"
)])
print(result$refinement_behavior[, c(
  "Dataset", "Method", "RefinementSelector", "NRefined",
  "WrongToCorrect", "CorrectToWrong", "CorrectionEfficiency"
)])
