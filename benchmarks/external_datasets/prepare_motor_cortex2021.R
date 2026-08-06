# benchmarks/external_datasets/prepare_motor_cortex2021.R
#
# Adapter entrypoint for a BICCN primary motor cortex atlas dataset. Provide a
# local Seurat/SCE/list RDS, h5ad, 10x directory, or matrix path in registry.csv
# before running.

source(file.path("benchmarks", "external_datasets", "prepare_utils.R"))

prepare_motor_cortex2021 <- function(registry_path = file.path("benchmarks", "external_datasets", "registry.csv"),
                                     output_dir = file.path("data", "external_prepared"),
                                     strict_confirmatory = TRUE) {
  prepare_external_dataset_from_registry(
    "MotorCortex2021",
    registry_path = registry_path,
    output_dir = output_dir,
    strict_confirmatory = strict_confirmatory
  )
}

if (any(grepl("prepare_motor_cortex2021\\.R$", commandArgs(trailingOnly = FALSE)))) {
  args <- commandArgs(trailingOnly = TRUE)
  registry_path <- if (length(args) >= 1) args[1] else file.path("benchmarks", "external_datasets", "registry.csv")
  output_dir <- if (length(args) >= 2) args[2] else file.path("data", "external_prepared")
  prepare_motor_cortex2021(registry_path, output_dir)
}
