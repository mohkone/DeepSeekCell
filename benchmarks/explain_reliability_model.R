# benchmarks/explain_reliability_model.R
#
# Recompute explainability outputs for a trained Risk-k reliability model.
#
# Usage:
#   Rscript benchmarks/explain_reliability_model.R results/reliability_model_v1.1_error.rds results/reliability_model_v1.1_error_training_features.csv

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

source("R/api.R")
source("R/utils.R")
source("R/ontology.R")
source("R/refinement.R")
source("R/reliability_model.R")
Sys.setenv(DEEPSEEKCELL_RUN_RELIABILITY_TRAINING_ON_SOURCE = "false")
source("benchmarks/train_reliability_model.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop(
    "Usage: Rscript benchmarks/explain_reliability_model.R model.rds training_features.csv [output_prefix]",
    call. = FALSE
  )
}

model_path <- args[1]
features_path <- args[2]
output_prefix <- if (length(args) >= 3) args[3] else sub("\\.rds$", "", model_path, ignore.case = TRUE)

if (!file.exists(model_path)) {
  stop("Model RDS not found: ", model_path, call. = FALSE)
}
if (!file.exists(features_path)) {
  stop("Feature CSV not found: ", features_path, call. = FALSE)
}

model <- readRDS(model_path)
features <- utils::read.csv(features_path, stringsAsFactors = FALSE)

write_explainability_outputs(
  model,
  features,
  paste0(output_prefix, ".rds")
)

message("Wrote explainability outputs with prefix ", output_prefix)
