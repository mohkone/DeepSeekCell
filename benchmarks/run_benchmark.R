#!/usr/bin/env Rscript
# benchmarks/run_benchmark.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

source("benchmarks/config.R")
source("benchmarks/datasets.R")
source("benchmarks/methods.R")
source("benchmarks/metrics.R")
source("benchmarks/statistics.R")

custom_mapping <- list(
  immune = c(
    "naive cd4 positive t cell" = "Naive T cell",
    "naive cd4 t cell" = "Naive T cell",
    "naive cd8 positive t cells" = "Cytotoxic T cell",
    "naive cd8 t cells" = "Cytotoxic T cell",
    "cd8 positive nkt like cells" = "Cytotoxic T cell",
    "cd8 nkt like cells" = "Cytotoxic T cell",
    "nkt like cells" = "Cytotoxic T cell",

    "naive b cells" = "B cell",
    "memory b cells" = "B cell",
    "plasma b cells" = "B cell",

    "classical monocytes" = "Classical monocyte",
    "non classical monocytes" = "Non-classical monocyte",
    "non classical monocyte" = "Non-classical monocyte",

    "natural killer cells" = "Natural killer cell",
    "natural killer  cells" = "Natural killer cell",

    "activated cd4 t cell" = "Naive T cell",
    "cd4 t cell" = "Naive T cell",
    "cd8 t cell" = "Cytotoxic T cell",
    "cytotoxic t cell" = "Cytotoxic T cell",
    "monocyte" = "Classical monocyte",
    "nk cell" = "Natural killer cell",
    "natural killer cell" = "Natural killer cell",
    "b" = "B cell",
    "b cell" = "B cell",
    "dendritic cell" = "Dendritic cell",
    "platelet" = "Platelet",
    "mast cells" = "Mast cell",
    "megakaryocyte" = "Platelet"
  ),
  pancreas = c(
    "alpha" = "Alpha cell",
    "alpha cell" = "Alpha cell",
    "beta" = "Beta cell",
    "beta cell" = "Beta cell",
    "delta" = "Delta cell",
    "delta cell" = "Delta cell",
    "gamma" = "Gamma cell",
    "gamma cell" = "Gamma cell",
    "epsilon" = "Epsilon cell",
    "epsilon cell" = "Epsilon cell",
    "pp cell" = "Gamma cell",
    "acinar" = "Acinar cell",
    "acinar cell" = "Acinar cell",
    "ductal" = "Ductal cell",
    "ductal cell" = "Ductal cell",
    "endothelial" = "Endothelial cell",
    "endothelial cell" = "Endothelial cell",
    "stellate" = "Stellate cell",
    "stellate cell" = "Stellate cell",
    "macrophage" = "Macrophage",
    "immune" = "Immune cell",
    "immune cell" = "Immune cell",
    "gamma pp cells" = "Gamma cell",
    "gamma pp cell" = "Gamma cell",
    "pp cells" = "Gamma cell",
    "pp cell" = "Gamma cell",
    "pancreatic polypeptide cell" = "Gamma cell",
    "mesenchymal cell" = "Stellate cell",
    "mesenchymal cells" = "Stellate cell",
    "pancreatic stellate cells" = "Stellate cell",
    "pancreatic stellate cell" = "Stellate cell"
  ),
  brain = c(
    "excitatory neuron" = "Excitatory neuron",
    "glutamatergic neuron" = "Excitatory neuron",
    "cholinergic neuron" = "Excitatory neuron",
    "inhibitory neuron" = "Inhibitory neuron",
    "gabaergic neuron" = "Inhibitory neuron",
    "gabaergic interneuron" = "Inhibitory neuron",
    "pvalb interneuron" = "Inhibitory neuron",
    "vip interneuron" = "Inhibitory neuron",
    "sst interneuron" = "Inhibitory neuron",
    "astrocyte" = "Astrocyte",
    "microglia" = "Microglia",
    "macrophage" = "Macrophage",
    "perivascular macrophage" = "Macrophage",
    "oligodendrocyte" = "Oligodendrocyte",
    "oligodendrocyte precursor cell" = "Oligodendrocyte precursor cell",
    "opc" = "Oligodendrocyte precursor cell",
    "endothelial" = "Endothelial cell",
    "endothelial cell" = "Endothelial cell",
    "vascular endothelial cell" = "Endothelial cell",
    "pericyte" = "Pericyte",
    "smooth muscle cell" = "Smooth muscle cell",
    "vascular smooth muscle cell" = "Smooth muscle cell",
    "ependymal cell" = "Ependymal cell"
  ),
  lung = c(
    "epithelial cell" = "Epithelial cell",
    "alveolar epithelial cell" = "Epithelial cell",
    "club cell" = "Epithelial cell",
    "ciliated cell" = "Epithelial cell",
    "t cell" = "T cell",
    "b cell" = "B cell",
    "plasma cell" = "Plasma cell",
    "natural killer cell" = "Natural killer cell",
    "nk cell" = "Natural killer cell",
    "macrophage" = "Macrophage",
    "platelet" = "Platelet",
    "platelets" = "Platelet",
    "monocyte" = "Monocyte",
    "neutrophil" = "Neutrophil",
    "dendritic cell" = "Dendritic cell",
    "mast cell" = "Mast cell",
    "fibroblast" = "Fibroblast",
    "endothelial cell" = "Endothelial cell",
    "smooth muscle cell" = "Smooth muscle cell",
    "alveolar macrophage" = "Macrophage",
    "alveolar macrophages" = "Macrophage",
    "cd4 t cell" = "T cell",
    "cd8 t cell" = "T cell",
    "regulatory t cell" = "T cell",
    "tumor cell" = "Epithelial cell",
    "cancer cell" = "Epithelial cell",
    "malignant cell" = "Epithelial cell",
    "cancer stem cells" = "Epithelial cell",
    "cancer stem cell" = "Epithelial cell",
    "epithelial cells" = "Epithelial cell",
    "type i cells" = "Epithelial cell",
    "type ii cells" = "Epithelial cell",
    "type i cell" = "Epithelial cell",
    "type ii cell" = "Epithelial cell",
    "pulmonary alveolar type i cells" = "Epithelial cell",
    "pulmonary alveolar type ii cells" = "Epithelial cell",
    "basal cells airway progenitor cells" = "Epithelial cell",
    "basal cell airway progenitor cell" = "Epithelial cell",
    "airway goblet cells" = "Epithelial cell",
    "airway goblet cell" = "Epithelial cell",
    "endothelial cells" = "Endothelial cell",
    "fibroblasts" = "Fibroblast",
    "mast cells" = "Mast cell",
    "immune system cells" = "Macrophage",
    "immune system cell" = "Macrophage",
    "alveolar macrophages" = "Macrophage",
    "alveolar macrophage" = "Macrophage",
    "secretory cell" = "Epithelial cell",
    "secretory cells" = "Epithelial cell",
    "ionocytes" = "Epithelial cell",
    "ionocyte" = "Epithelial cell",
    "airway epithelial cells" = "Epithelial cell",
    "airway epithelial cell" = "Epithelial cell",
    "endothelial cell" = "Endothelial cell",
    "endothelial cells" = "Endothelial cell",
    "immune system cells" = "Macrophage",
    "immune system cell" = "Macrophage",
    "airway epithelial cells" = "Epithelial cell",
    "airway epithelial cell" = "Epithelial cell",
    "secretory cell" = "Epithelial cell",
    "secretory cells" = "Epithelial cell",
    "ionocytes" = "Epithelial cell",
    "ionocyte" = "Epithelial cell",
    "alveolar macrophages" = "Macrophage",
    "alveolar macrophage" = "Macrophage",
    "endothelial cell" = "Endothelial cell",
    "endothelial cells" = "Endothelial cell",
    "erythrocyte" = "Erythrocyte",
    "erythrocytes" = "Erythrocyte",
    "red blood cell" = "Erythrocyte",
    "red blood cells" = "Erythrocyte"
  )
)

get_domain <- function(tissue) {
  tissue_low <- tolower(tissue)
  if (tissue_low %in% c("pbmc", "blood", "immune system")) return("immune")
  if (grepl("pancreas", tissue_low)) return("pancreas")
  if (grepl("brain", tissue_low)) return("brain")
  if (grepl("lung", tissue_low)) return("lung")
  "general"
}

harmonise_labels <- function(labels, tissue, ont_data = NULL) {
  labs <- normalize_label(labels)
  out <- rep("Unknown", length(labs))

  domain <- get_domain(tissue)

  maps <- c(
    custom_mapping[[domain]] %||% character(),
    custom_mapping$immune,
    custom_mapping$pancreas,
    custom_mapping$brain,
    custom_mapping$lung
  )

  for (i in seq_along(labs)) {
    lab <- labs[i]

    if (.is_unknown_label(lab)) {
      out[i] <- "Unknown"
      next
    }

    if (lab %in% names(maps)) {
      out[i] <- maps[[lab]]
      next
    }

    if (domain == "immune") {
      if (grepl("non classical|nonclassical|fcgr3a", lab)) {
        out[i] <- "Non-classical monocyte"
      } else if (grepl("classical monocyte|monocytes|monocyte|cd14", lab)) {
        out[i] <- "Classical monocyte"
      } else if (grepl("cd8|cytotoxic|nkt", lab)) {
        out[i] <- "Cytotoxic T cell"
      } else if (grepl("cd4|naive|memory|t cell", lab)) {
        out[i] <- "Naive T cell"
      } else if (.is_nk_label(lab)) {
        out[i] <- "Natural killer cell"
      } else if (grepl("b cells|b cell|plasma", lab)) {
        out[i] <- "B cell"
      } else if (grepl("platelet|platelets|megakaryocyte", lab)) {
        out[i] <- "Platelet"
      } else if (grepl("dendritic|dc", lab)) {
        out[i] <- "Dendritic cell"
      } else if (grepl("mast", lab)) {
        out[i] <- "Mast cell"
      }
    } else if (domain == "pancreas") {
      if (grepl("alpha", lab)) out[i] <- "Alpha cell"
      else if (grepl("beta", lab)) out[i] <- "Beta cell"
      else if (grepl("delta", lab)) out[i] <- "Delta cell"
      else if (grepl("gamma|pp", lab)) out[i] <- "Gamma cell"
      else if (grepl("acinar", lab)) out[i] <- "Acinar cell"
      else if (grepl("duct", lab)) out[i] <- "Ductal cell"
      else if (grepl("endothelial", lab)) out[i] <- "Endothelial cell"
      else if (grepl("stellate", lab)) out[i] <- "Stellate cell"
      else if (grepl("macrophage|immune", lab)) out[i] <- "Immune cell"
    } else if (domain == "brain") {
      if (grepl("inhibitory|gaba|interneuron|pvalb|vip|sst", lab)) out[i] <- "Inhibitory neuron"
      else if (grepl("excitatory|glutamatergic|cholinergic|neuron", lab)) out[i] <- "Excitatory neuron"
      else if (grepl("astro", lab)) out[i] <- "Astrocyte"
      else if (grepl("macrophage|perivascular|(^| )pvm", lab)) out[i] <- "Macrophage"
      else if (grepl("micro", lab)) out[i] <- "Microglia"
      else if (grepl("oligo", lab) && grepl("precursor|opc", lab)) out[i] <- "Oligodendrocyte precursor cell"
      else if (grepl("oligo", lab)) out[i] <- "Oligodendrocyte"
      else if (grepl("endothelial|(^| )vend|vascular endothelial", lab)) out[i] <- "Endothelial cell"
      else if (grepl("pericyte", lab)) out[i] <- "Pericyte"
      else if (grepl("smooth muscle|(^| )vsmc", lab)) out[i] <- "Smooth muscle cell"
      else if (grepl("ependymal", lab)) out[i] <- "Ependymal cell"
    } else if (domain == "lung") {
      if (grepl("patient[0-9]+.*specific|epithelial|alveolar|club|ciliated|basal|goblet|ionocyte|secretory|tumor|malignant|cancer|type i|type ii", lab)) out[i] <- "Epithelial cell"
      else if (grepl("t cell", lab)) out[i] <- "T cell"
      else if (grepl("b cell", lab)) out[i] <- "B cell"
      else if (grepl("plasmacytoid|(^| )pdc|dendritic|(^| )dc|dcs", lab)) out[i] <- "Dendritic cell"
      else if (grepl("plasma", lab)) out[i] <- "Plasma cell"
      else if (.is_nk_label(lab)) out[i] <- "Natural killer cell"
      else if (grepl("platelet", lab)) out[i] <- "Platelet"
      else if (grepl("macrophage|(^| )mph($| )", lab)) out[i] <- "Macrophage"
      else if (grepl("monocyte", lab)) out[i] <- "Monocyte"
      else if (grepl("neutrophil", lab)) out[i] <- "Neutrophil"
      else if (grepl("mast", lab)) out[i] <- "Mast cell"
      else if (grepl("fibroblast", lab)) out[i] <- "Fibroblast"
      else if (grepl("endothelial|(^| )ec($| )|capillary|venous|arterial|lymphatic", lab)) out[i] <- "Endothelial cell"
      else if (grepl("smooth muscle", lab)) out[i] <- "Smooth muscle cell"
      else if (grepl("^[bt]?rbc$|erythrocyte|erythroid|red blood cell", lab)) out[i] <- "Erythrocyte"
      else if (grepl("neuroendocrine", lab)) out[i] <- "Epithelial cell"
    }

    if (out[i] == "Unknown") {
      out[i] <- labels[i]
    }
  }

  names(out) <- names(labels)
  trimws(out)
}

cell_to_cluster_prediction <- function(cell_pred, seu, cluster_names) {

  cluster_ids <- as.character(seu$seurat_clusters)

  if (is.null(names(cluster_ids))) {
    names(cluster_ids) <- colnames(seu)
  }

  if (is.null(names(cell_pred))) {
    names(cell_pred) <- colnames(seu)
  }

  cluster_pred <- sapply(cluster_names, function(cl) {

    cell_names <- names(cluster_ids)[cluster_ids == cl]

    if (length(cell_names) == 0) {
      return("Unknown")
    }

    pred_cl <- cell_pred[cell_names]
    pred_cl <- pred_cl[!is.na(pred_cl) & pred_cl != ""]

    if (length(pred_cl) == 0) {
      return("Unknown")
    }

    names(sort(table(pred_cl), decreasing = TRUE))[1]
  })

  names(cluster_pred) <- cluster_names
  cluster_pred
}

run_benchmark_wrapper <- function(dataset_name,
                                  data,
                                  method,
                                  ont_data,
                                  db_file,
                                  replicate = NA_integer_,
                                  api_key = NULL) {

  markers <- data$markers
  truth <- data$truth
  tissue <- data$tissue
  species <- data$species
  seu <- data$seurat_obj

  method_res <- switch(
    method,
    DeepSeek = run_llm_annotation(markers, tissue, species, dataset_name, "deepseek", api_key),
    SingleR = run_singler(seu, dataset_name),
    scType = run_sctype_custom(seu, tissue, db_file, species = species, verbose = TRUE),
    scmap = run_scmap(seu, dataset_name),
    CellTypist = run_celltypist(seu, tissue, species = species, dataset_name = dataset_name),
    stop("Unknown method: ", method, call. = FALSE)
  )

  pred <- method_res$predictions

  if (method %in% c("SingleR", "scType", "scmap")) {
    pred <- cell_to_cluster_prediction(pred, seu, names(markers))
  } else {
    pred <- setNames(pred, names(markers))
  }

  truth_aligned <- truth[names(pred)]

  pred_h <- harmonise_labels(pred, tissue, ont_data)
  truth_h <- harmonise_labels(truth_aligned, tissue, ont_data)

  cluster_purity <- if (!is.null(data$purity)) {
    as.numeric(data$purity[names(pred)])
  } else {
    rep(NA_real_, length(pred))
  }

  debug_df <- data.frame(
    Replicate = replicate,
    Dataset = dataset_name,
    Method = method,
    Cluster = names(pred),
    ClusterPurity = cluster_purity,
    RawPrediction = as.character(pred),
    RawTruth = as.character(truth_aligned),
    HarmonisedPrediction = as.character(pred_h),
    HarmonisedTruth = as.character(truth_h),
    stringsAsFactors = FALSE
  )

  debug_dir <- file.path("results", "benchmark_debug")
  dir.create(debug_dir, recursive = TRUE, showWarnings = FALSE)

  write.csv(
    debug_df,
    file = file.path(
      debug_dir,
      paste0("rep", replicate, "_", dataset_name, "_", method, "_debug.csv")
    ),
    row.names = FALSE
  )

  metrics <- evaluate_metrics(pred_h, truth_h, ont_data, tissue = tissue)

  data.frame(
    Dataset = dataset_name,
    Tissue = tissue,
    Species = species,
    Method = method,
    NClusters = length(markers),
    RuntimeSec = method_res$runtime_sec,
    CostUSD = method_res$cost_usd,
    Tokens = method_res$tokens,
    MeanClusterPurity = mean(cluster_purity, na.rm = TRUE),
    MinClusterPurity = min(cluster_purity, na.rm = TRUE),
    t(metrics),
    stringsAsFactors = FALSE
  )
}

run_replicated_benchmark <- function(n_replicates = 1,
                                     methods_list,
                                     ont_data,
                                     sctype_db_path,
                                     deepseek_key) {

  all_results <- list()
  dataset_summaries <- list()
  cluster_summaries <- list()

  for (rep in seq_len(n_replicates)) {
    cat("\n========== REPLICATE", rep, "of", n_replicates, "==========\n")

    seed_val <- rep * 100
    set.seed(seed_val)

    datasets <- load_benchmark_datasets(seed = seed_val)

    dataset_summaries[[length(dataset_summaries) + 1]] <- bind_rows(
      lapply(datasets, function(x) x$dataset_summary)
    ) %>%
      mutate(Replicate = rep, Seed = seed_val, .before = 1)

    cluster_summaries[[length(cluster_summaries) + 1]] <- bind_rows(
      lapply(datasets, function(x) x$cluster_summary)
    ) %>%
      mutate(Replicate = rep, Seed = seed_val, .before = 1)

    for (ds in names(datasets)) {
      for (m in methods_list) {
        cat("\n--- Replicate", rep, "-", m, "on", ds, "---\n")

        api_key <- if (m == "DeepSeek") {
          deepseek_key
        } else {
          NULL
        }

        if (m == "DeepSeek" && (is.null(api_key) || api_key == "")) {
          cat("Skipping ", m, ": API key not set.\n", sep = "")
          next
        }

        if (m == "CellTypist" && tolower(datasets[[ds]]$species) != "human") {
          cat("Skipping CellTypist on ", ds, ": human-only baseline configuration.\n", sep = "")
          next
        }

        row <- tryCatch(
          run_benchmark_wrapper(
            dataset_name = ds,
            data = datasets[[ds]],
            method = m,
            ont_data = ont_data,
            db_file = sctype_db_path,
            replicate = rep,
            api_key = api_key
          ),
          error = function(e) {
            warning("Failed ", m, " on ", ds, ": ", e$message)
            NULL
          }
        )

        if (!is.null(row)) {
          row$Replicate <- rep
          all_results[[length(all_results) + 1]] <- row
        }
      }
    }
  }

  if (length(all_results) == 0) return(NULL)

  full <- bind_rows(all_results)

  full <- full %>%
    mutate(across(where(is.numeric), ~ ifelse(is.nan(.x), NA_real_, .x)))

  mean_or_na <- function(x) {
    if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  }

  sd_or_na <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) < 2) NA_real_ else stats::sd(x)
  }

  se_or_na <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) < 2) NA_real_ else stats::sd(x) / sqrt(length(x))
  }

  ci95_or_na <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) < 2) NA_real_ else stats::qt(0.975, df = length(x) - 1) * stats::sd(x) / sqrt(length(x))
  }

  summary <- full %>%
    group_by(Dataset, Tissue, Species, Method) %>%
    summarise(
      NReplicates = dplyr::n(),
      SuccessfulRuns = sum(!is.na(MacroF1)),
      across(
        c(
          ARI, MacroF1, Accuracy, BalancedAcc, CladeAcc, UnknownRate,
          RuntimeSec, CostUSD, Tokens, MeanClusterPurity, MinClusterPurity,
          EvaluatedClusters
        ),
        list(
          mean = mean_or_na,
          sd = sd_or_na,
          se = se_or_na,
          ci95 = ci95_or_na
        ),
        .names = "{.col}_{.fn}"
      ),
      NClusters = mean_or_na(NClusters),
      .groups = "drop"
    )

  list(
    full = full,
    summary = summary,
    dataset_summary = bind_rows(dataset_summaries),
    cluster_summary = bind_rows(cluster_summaries)
  )
}

plot_summary_with_sd <- function(summary_df,
                                 metric = "MacroF1",
                                 output_pdf = "benchmark_plot_with_sd.pdf") {

  mean_col <- paste0(metric, "_mean")
  err_col <- paste0(metric, "_ci95")

  if (!err_col %in% names(summary_df)) {
    err_col <- paste0(metric, "_sd")
  }

  bounded_metric <- metric %in% c(
    "ARI", "MacroF1", "Accuracy", "BalancedAcc", "CladeAcc", "UnknownRate"
  )

  err <- if (err_col %in% names(summary_df)) summary_df[[err_col]] else NA_real_
  summary_df$ErrorMin <- summary_df[[mean_col]] - err
  summary_df$ErrorMax <- summary_df[[mean_col]] + err
  summary_df$ErrorMin <- pmax(summary_df$ErrorMin, 0)

  if (bounded_metric) {
    summary_df$ErrorMax <- pmin(summary_df$ErrorMax, 1)
  }

  p <- ggplot(
    summary_df,
    aes(x = Method, y = .data[[mean_col]], fill = Method)
  ) +
    geom_col(position = "dodge") +
    geom_errorbar(
      aes(
        ymin = ErrorMin,
        ymax = ErrorMax
      ),
      width = 0.2
    ) +
    facet_wrap(~ Dataset, scales = "free_y") +
    theme_minimal(base_size = 12) +
    labs(
      title = paste0("DeepSeekCell Benchmark: ", metric, " mean with 95% CI"),
      x = "Method",
      y = metric
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(output_pdf, p, width = 14, height = 8)
  message("Plot saved to ", output_pdf)

  invisible(p)
}

plot_global_metric <- function(summary_df, metric, output_pdf) {
  mean_col <- paste0(metric, "_mean")

  p <- ggplot2::ggplot(
    summary_df,
    ggplot2::aes(
      x = Dataset,
      y = .data[[mean_col]],
      fill = Method
    )
  ) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(
      title = paste("DeepSeekCell benchmark:", metric),
      x = "Dataset",
      y = metric
    )

  ggplot2::ggsave(output_pdf, p, width = 11, height = 7)
  invisible(p)
}

write_benchmark_manifest <- function(results,
                                     methods,
                                     n_replicates,
                                     output_file = "results/benchmark_manifest.txt") {
  lines <- c(
    "DeepSeekCell benchmark manifest",
    "================================",
    paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste("Benchmark mode:", BENCHMARK_MODE),
    paste("Replicates requested:", n_replicates),
    paste("Methods:", paste(methods, collapse = ", ")),
    paste("Top markers per cluster:", TOP_MARKERS),
    paste("Seurat resolution:", SEURAT_RESOLUTION),
    paste("PCA dimensions configured:", N_PCS),
    paste("Cache version:", BENCHMARK_CACHE_VERSION),
    "",
    "Interpretation notes:",
    "- DeepSeekCell is evaluated as a closed-label, marker-guided annotation workflow.",
    "- SingleR, scType, scmap, and CellTypist are evaluated as baseline methods when their dependencies are available.",
    "- scmap is run with scmapCell nearest-neighbor label transfer before cluster-level majority voting.",
    "- CellTypist is evaluated on cluster-average pseudo-profiles using tissue-aware pretrained models unless CELLTYPIST_MODEL overrides the model choice.",
    "- Pairwise Wilcoxon and Friedman test outputs are exploratory because replicate blocks are small and some datasets are deterministic.",
    "- LLM stability is summarized for fixed cluster structures using per-cluster debug prediction files.",
    "- Results summarize cluster-level labels after explicit label harmonization.",
    "- Clade accuracy uses the package Cell Ontology mapper and lineage relation checks.",
    "- Treat datasets with low cluster purity as lower-confidence benchmark evidence.",
    "",
    "Dataset summaries:",
    paste(capture.output(print(results$dataset_summary)), collapse = "\n")
  )

  writeLines(lines, output_file)
}

available_benchmark_methods <- function(requested_methods) {
  selected <- character()

  for (method in requested_methods) {
    if (identical(method, "scmap") && !is_scmap_available()) {
      message("Skipping scmap: Bioconductor package scmap is not installed.")
      next
    }

    if (identical(method, "CellTypist") && !is_celltypist_available()) {
      message(
        "Skipping CellTypist: reticulate or Python modules ",
        "celltypist/anndata/numpy are not available."
      )
      next
    }

    selected <- c(selected, method)
  }

  selected
}

main <- function(n_replicates = 1) {
  deepseek_key <- Sys.getenv("DEEPSEEK_API_KEY")

  ont_data <- load_benchmark_ontology(ONTOLOGY_FILE)
  sctype_db_path <- find_sctype_db()

  methods <- available_benchmark_methods(
    c("DeepSeek", "scType", "SingleR", "scmap", "CellTypist")
  )


  results <- run_replicated_benchmark(
    n_replicates = n_replicates,
    methods_list = methods,
    ont_data = ont_data,
    sctype_db_path = sctype_db_path,
    deepseek_key = deepseek_key
  )

  if (is.null(results)) {
    message("No benchmark results produced.")
    return(invisible(NULL))
  }

  dir.create("results", showWarnings = FALSE)

  write.csv(results$full, "results/benchmark_results_full.csv", row.names = FALSE)
  write.csv(results$summary, "results/benchmark_results_summary.csv", row.names = FALSE)
  write.csv(results$dataset_summary, "results/dataset_summary.csv", row.names = FALSE)
  write.csv(results$cluster_summary, "results/cluster_summary.csv", row.names = FALSE)
  writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
  write_benchmark_manifest(results, methods, n_replicates)

  final_table <- results$summary %>%
    dplyr::select(
      Dataset,
      Tissue,
      Species,
      Method,
      NReplicates,
      SuccessfulRuns,
      ARI_mean,
      ARI_ci95,
      MacroF1_mean,
      MacroF1_ci95,
      Accuracy_mean,
      Accuracy_ci95,
      CladeAcc_mean,
      CladeAcc_ci95,
      UnknownRate_mean,
      EvaluatedClusters_mean,
      NClusters,
      RuntimeSec_mean,
      MeanClusterPurity_mean
    ) %>%
    dplyr::arrange(Dataset, desc(MacroF1_mean))

  write.csv(
    final_table,
    "results/final_benchmark_table.csv",
    row.names = FALSE
  )

  statistical_tests <- run_benchmark_statistical_tests(results$full)
  write.csv(
    statistical_tests$pairwise_wilcoxon,
    "results/benchmark_pairwise_wilcoxon.csv",
    row.names = FALSE
  )
  write.csv(
    statistical_tests$friedman,
    "results/benchmark_friedman_tests.csv",
    row.names = FALSE
  )
  write.csv(
    statistical_tests$llm_stability,
    "results/benchmark_llm_stability.csv",
    row.names = FALSE
  )

  print(results$summary)

  plot_summary_with_sd(results$summary, metric = "MacroF1", output_pdf = "results/benchmark_macroF1.pdf")
  plot_summary_with_sd(results$summary, metric = "Accuracy", output_pdf = "results/benchmark_accuracy.pdf")
  plot_summary_with_sd(results$summary, metric = "CladeAcc", output_pdf = "results/benchmark_clade_accuracy.pdf")
  plot_summary_with_sd(results$summary, metric = "RuntimeSec", output_pdf = "results/benchmark_runtime.pdf")
  plot_global_metric(results$summary, "ARI", "results/global_ARI.pdf")
  plot_global_metric(results$summary, "MacroF1", "results/global_MacroF1.pdf")
  plot_global_metric(results$summary, "Accuracy", "results/global_Accuracy.pdf")
  plot_global_metric(results$summary, "RuntimeSec", "results/global_Runtime.pdf")

  invisible(results)
}

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  n_rep <- if (length(args) > 0) {
    as.integer(args[1])
  } else {
    DEFAULT_BENCHMARK_REPLICATES
  }
  main(n_replicates = n_rep)
}
