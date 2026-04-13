#!/usr/bin/env Rscript
# DeepSeekCell Benchmark – Main Runner

suppressPackageStartupMessages({
  library(Seurat)
  library(mclust)
  library(ggplot2)
  library(dplyr)
  library(reshape2)
})

# Source all modules
source("benchmarks/config.R")
source("benchmarks/datasets.R")
source("benchmarks/methods.R")
source("benchmarks/metrics.R")
source("benchmarks/plot_results.R")

run_benchmark_wrapper <- function(dataset_name, data, method, ont_data, db_file, api_key = NULL) {
  markers <- data$markers
  truth   <- data$truth
  tissue  <- data$tissue
  species <- data$species
  seu     <- data$seurat_obj
  
  pred <- switch(method,
                 DeepSeek = run_llm_annotation(markers, tissue, species, "deepseek", api_key),
                 GPT4o    = run_llm_annotation(markers, tissue, species, "gpt4", api_key),
                 scType   = {
                   scType_tissue <- switch(dataset_name,
                                           PBMC = "Immune system", Pancreas = "Pancreas", Brain = "Brain")
                   run_sctype_custom(seu, scType_tissue, db_file)
                 },
                 SingleR  = run_singler(seu, dataset_name)
  )
  
  if (method %in% c("DeepSeek", "GPT4o")) {
    pred <- setNames(pred, names(markers))
  } else {
    cluster_ids <- as.character(seu$seurat_clusters)
    names(pred) <- colnames(seu)
    cluster_pred <- sapply(unique(cluster_ids), function(cl) {
      cells <- WhichCells(seu, idents = cl)
      if (length(cells) == 0) return("Unknown")
      pred_cl <- pred[cells]
      names(sort(table(pred_cl), decreasing = TRUE))[1]
    })
    names(cluster_pred) <- unique(cluster_ids)
    pred <- cluster_pred[names(markers)]
  }
  truth_aligned <- truth[names(pred)]
  evaluate_metrics(pred, truth_aligned, dataset_name, ont_data)
}

main <- function() {
  ont_data <- load_cell_ontology(ONTOLOGY_FILE)
  sctype_db_path <- find_sctype_db()
  message("Using scType database: ", sctype_db_path)
  
  datasets <- list(
    PBMC = get_pbmc_data(),
    Pancreas = get_pancreas_data(),
    Brain = get_brain_data()
  )
  methods <- c("DeepSeek", "GPT4o", "scType", "SingleR")
  results <- data.frame()
  
  for (ds in names(datasets)) {
    for (m in methods) {
      cat("\n========================================\n")
      cat("Running", m, "on", ds, "...\n")
      api_key <- if (m == "DeepSeek") DEEPSEEK_KEY else if (m == "GPT4o") OFOX_KEY else NULL
      if ((m %in% c("DeepSeek", "GPT4o")) && is.null(api_key)) {
        cat("Skipping", m, "– API key not set.\n")
        next
      }
      metrics <- run_benchmark_wrapper(ds, datasets[[ds]], m, ont_data,
                                       sctype_db_path, api_key)
      results <- rbind(results, data.frame(Dataset = ds, Method = m, t(metrics)))
    }
  }
  
  if (length(unmatched_log) > 0) {
    writeLines(unique(unmatched_log), "unmatched_celltypes.log")
    message("Unmatched cell types written to unmatched_celltypes.log")
  }
  
  if (nrow(results) > 0) {
    write.csv(results, "benchmark_results.csv", row.names = FALSE)
    print(results)
    plot_benchmark()
  } else {
    message("No results produced.")
  }
}

if (interactive()) main()