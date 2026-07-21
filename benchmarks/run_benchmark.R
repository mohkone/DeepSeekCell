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
source("R/utils.R")
source("benchmarks/ablation.R")

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

empty_native_cell_metrics <- function() {
  c(
    NativeCellARI = NA_real_,
    NativeCellMacroF1 = NA_real_,
    NativeCellAccuracy = NA_real_,
    NativeCellBalancedAcc = NA_real_,
    NativeCellCladeAcc = NA_real_,
    NativeCellUnknownRate = NA_real_,
    NativeCellEvaluatedCells = NA_real_
  )
}

uses_raw_llm_confidence <- function(arm_name) {
  grepl("-Plain$", arm_name) || grepl("Cell-Evidence$", arm_name)
}

refinement_method_family <- function(method) {
  sub("-(RandomK|ConfidenceK|NoOntologyK|SelfRefined|FullRefined)$", "", method)
}

ablation_method_names <- function(method_prefix = "DeepSeek") {
  cell_prefix <- if (identical(method_prefix, "DeepSeek")) {
    "DeepSeekCell"
  } else {
    paste0(method_prefix, "Cell")
  }

  list(
    plain = paste0(method_prefix, "-Plain"),
    evidence = paste0(cell_prefix, "-Evidence"),
    calibrated = paste0(cell_prefix, "-Calibrated"),
    random = paste0(cell_prefix, "-RandomK"),
    confidence = paste0(cell_prefix, "-ConfidenceK"),
    no_ontology = paste0(cell_prefix, "-NoOntologyK"),
    self = paste0(cell_prefix, "-SelfRefined"),
    full = paste0(cell_prefix, "-FullRefined")
  )
}

rename_ablation_arms <- function(arms, method_names) {
  mapping <- c(
    "DeepSeek-Plain" = method_names$plain,
    "DeepSeekCell-Evidence" = method_names$evidence,
    "DeepSeekCell-Calibrated" = method_names$calibrated,
    "DeepSeekCell-SelfRefined" = method_names$self
  )

  out <- list()
  for (old_name in names(arms)) {
    new_name <- mapping[[old_name]] %||% old_name
    arms[[old_name]]$AblationArm <- new_name
    out[[new_name]] <- arms[[old_name]]
  }

  out
}

compute_native_cell_metrics <- function(cell_pred, seu, tissue, ont_data) {
  if (is.null(cell_pred) || length(cell_pred) == 0 || is.null(seu$true_label)) {
    return(empty_native_cell_metrics())
  }

  if (is.null(names(cell_pred))) {
    names(cell_pred) <- colnames(seu)
  }

  cell_names <- intersect(names(cell_pred), colnames(seu))
  if (length(cell_names) < 2) {
    return(empty_native_cell_metrics())
  }

  cell_truth <- as.character(seu$true_label)
  names(cell_truth) <- colnames(seu)

  pred_h <- harmonise_labels(cell_pred[cell_names], tissue, ont_data)
  truth_h <- harmonise_labels(cell_truth[cell_names], tissue, ont_data)
  metrics <- evaluate_metrics(pred_h, truth_h, ont_data, tissue = tissue)

  c(
    NativeCellARI = as.numeric(metrics[["ARI"]]),
    NativeCellMacroF1 = as.numeric(metrics[["MacroF1"]]),
    NativeCellAccuracy = as.numeric(metrics[["Accuracy"]]),
    NativeCellBalancedAcc = as.numeric(metrics[["BalancedAcc"]]),
    NativeCellCladeAcc = as.numeric(metrics[["CladeAcc"]]),
    NativeCellUnknownRate = as.numeric(metrics[["UnknownRate"]]),
    NativeCellEvaluatedCells = as.numeric(metrics[["EvaluatedClusters"]])
  )
}

evaluate_ablation_arm <- function(annotations,
                                  arm_name,
                                  dataset_name,
                                  data,
                                  ont_data,
                                  replicate,
                                  method_res,
                                  first_pass_hash,
                                  llm_backend = NA_character_,
                                  llm_model_id = NA_character_,
                                  refinement_selector = "none",
                                  refinement_budget_k = 0,
                                  selected_clusters = character(),
                                  refinement_runtime_sec = 0,
                                  refinement_cost_usd = 0,
                                  refinement_tokens = 0,
                                  second_pass_calls = 0) {
  markers <- data$markers
  truth <- data$truth
  tissue <- data$tissue
  species <- data$species

  pred <- setNames(as.character(annotations$CellType), annotations$Cluster)
  pred <- pred[names(markers)]
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
    Method = arm_name,
    LLMBackend = llm_backend,
    LLMModelID = llm_model_id,
    RefinementSelector = refinement_selector,
    RefinementBudgetK = refinement_budget_k,
    Cluster = names(pred),
    ClusterPurity = cluster_purity,
    RawPrediction = as.character(pred),
    RawTruth = as.character(truth_aligned),
    HarmonisedPrediction = as.character(pred_h),
    HarmonisedTruth = as.character(truth_h),
    LLMConfidence = annotations$LLMConfidence %||% annotations$Confidence,
    Confidence = annotations$Confidence,
    EvidenceBestCellType = annotations$EvidenceBestCellType %||% NA_character_,
    EvidenceConflict = annotations$EvidenceConflict %||% NA,
    RequiresRefinement = annotations$RequiresRefinement %||% NA,
    FirstPassCellType = annotations$FirstPassCellType %||% annotations$CellType,
    SelectedForRefinement = names(pred) %in% selected_clusters,
    WasFlagged = annotations$WasFlagged %||% FALSE,
    WasRefined = annotations$WasRefined %||% FALSE,
    RefinementChangedLabel = annotations$RefinementChangedLabel %||% FALSE,
    FirstPassResponseHash = first_pass_hash,
    stringsAsFactors = FALSE
  )

  debug_dir <- file.path("results", "benchmark_debug")
  dir.create(debug_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(
    debug_df,
    file = file.path(
      debug_dir,
      paste0("rep", replicate, "_", dataset_name, "_", arm_name, "_debug.csv")
    ),
    row.names = FALSE
  )

  metrics <- evaluate_metrics(pred_h, truth_h, ont_data, tissue = tissue)

  confidence <- annotations$Confidence
  if (uses_raw_llm_confidence(arm_name)) {
    confidence <- annotations$LLMConfidence %||% annotations$Confidence
  }

  correctness <- as.character(pred_h) == as.character(truth_h)
  confidence_quality <- tryCatch(
    compute_confidence_quality(correctness, confidence),
    error = function(e) data.frame()
  )

  reliability <- tryCatch(
    compute_reliability_bins(correctness, confidence) %>%
      dplyr::mutate(
        Dataset = dataset_name,
        Tissue = tissue,
        Species = species,
        Method = arm_name,
        LLMBackend = llm_backend,
        LLMModelID = llm_model_id,
        Replicate = replicate,
        RefinementSelector = refinement_selector,
        RefinementBudgetK = refinement_budget_k,
        .before = 1
      ),
    error = function(e) data.frame()
  )

  runtime <- method_res$runtime_sec
  cost <- method_res$cost_usd
  tokens <- method_res$tokens
  if (isTRUE(second_pass_calls > 0)) {
    runtime <- runtime + refinement_runtime_sec
    cost <- cost + refinement_cost_usd
    tokens <- tokens + refinement_tokens
  }

  result_row <- data.frame(
    Dataset = dataset_name,
    Tissue = tissue,
    Species = species,
    Method = arm_name,
    LLMBackend = llm_backend,
    LLMModelID = llm_model_id,
    RefinementSelector = refinement_selector,
    RefinementBudgetK = refinement_budget_k,
    NClusters = length(markers),
    RuntimeSec = runtime,
    CostUSD = cost,
    Tokens = tokens,
    FirstPassRuntimeSec = method_res$runtime_sec,
    RefinementRuntimeSec = refinement_runtime_sec,
    FirstPassTokens = method_res$tokens,
    RefinementTokens = refinement_tokens,
    SecondPassCalls = second_pass_calls,
    FirstPassResponseHash = first_pass_hash,
    MeanClusterPurity = mean(cluster_purity, na.rm = TRUE),
    MinClusterPurity = min(cluster_purity, na.rm = TRUE),
    t(metrics),
    stringsAsFactors = FALSE
  )

  if (nrow(confidence_quality) > 0) {
    confidence_quality <- confidence_quality %>%
      dplyr::mutate(
        Dataset = dataset_name,
        Tissue = tissue,
        Species = species,
        Method = arm_name,
        LLMBackend = llm_backend,
        LLMModelID = llm_model_id,
        Replicate = replicate,
        RefinementSelector = refinement_selector,
        RefinementBudgetK = refinement_budget_k,
        ConfidenceColumn = ifelse(
          uses_raw_llm_confidence(arm_name),
          "LLMConfidence",
          "Confidence"
        ),
        .before = 1
      )
  }

  list(
    result = result_row,
    debug = debug_df,
    confidence_quality = confidence_quality,
    reliability = reliability,
    pred_h = pred_h,
    truth_h = truth_h
  )
}

safe_usage_number <- function(x) {
  out <- suppressWarnings(as.numeric(x %||% 0))
  if (length(out) == 0 || is.na(out)) {
    return(0)
  }
  out
}

disable_ontology_annotation_context <- function(annotations) {
  if (!is.data.frame(annotations)) {
    return(annotations)
  }
  annotations$CL_ID <- NA_character_
  annotations$OntologyLabel <- NA_character_
  annotations$MatchMethod <- "ontology_disabled"
  annotations$OntologyMatchScore <- NA_real_
  if ("OntologyEvidenceScore" %in% names(annotations)) {
    annotations$OntologyEvidenceScore <- 0.5
  }
  annotations
}

use_llm_response_cache <- function() {
  value <- tolower(Sys.getenv("DEEPSEEKCELL_USE_LLM_CACHE", unset = "true"))
  value %in% c("1", "true", "yes", "y")
}

include_ollama_full_refinement <- function() {
  value <- tolower(Sys.getenv("DEEPSEEKCELL_OLLAMA_FULL_REFINEMENT", unset = "false"))
  value %in% c("1", "true", "yes", "y")
}

read_llm_response_cache <- function(cache_file) {
  if (!use_llm_response_cache() || !file.exists(cache_file)) {
    return(NULL)
  }

  tryCatch(read_first_pass_cache(cache_file), error = function(e) NULL)
}

select_refinement_clusters <- function(first_pass_evidence,
                                       selector,
                                       k,
                                       replicate,
                                       dataset_name) {
  clusters <- as.character(first_pass_evidence$Cluster)
  k <- suppressWarnings(as.integer(k %||% 0))
  if (length(k) == 0 || is.na(k)) {
    k <- 0
  }
  k <- min(max(k, 0), length(clusters))

  if (identical(selector, "full")) {
    return(clusters)
  }

  if (k == 0) {
    return(character())
  }

  if (identical(selector, "evidence-k")) {
    requires_refinement <- as_flag(first_pass_evidence$RequiresRefinement)
    requires_refinement[is.na(requires_refinement)] <- FALSE
    return(clusters[requires_refinement])
  }

  if (identical(selector, "evidence-no-ontology-k")) {
    requires_refinement <- as_flag(first_pass_evidence$RequiresRefinement)
    requires_refinement[is.na(requires_refinement)] <- FALSE
    marker_gap <- suppressWarnings(
      as.numeric(first_pass_evidence$BestMarkerEvidenceScore %||% 0) -
        as.numeric(first_pass_evidence$MarkerEvidenceScore %||% 0)
    )
    marker_gap[is.na(marker_gap)] <- 0
    best_marker <- suppressWarnings(as.numeric(first_pass_evidence$BestMarkerEvidenceScore %||% 0))
    best_marker[is.na(best_marker)] <- 0
    confidence <- as_confidence(first_pass_evidence$LLMConfidence %||% first_pass_evidence$Confidence)
    confidence[is.na(confidence)] <- 0.5

    flagged_ord <- order(
      marker_gap[requires_refinement],
      best_marker[requires_refinement],
      -confidence[requires_refinement],
      decreasing = TRUE
    )
    flagged <- clusters[requires_refinement][flagged_ord]

    remainder_idx <- !requires_refinement
    remainder_ord <- order(
      marker_gap[remainder_idx],
      best_marker[remainder_idx],
      -confidence[remainder_idx],
      decreasing = TRUE
    )
    remainder <- clusters[remainder_idx][remainder_ord]

    return(head(unique(c(flagged, remainder)), k))
  }

  if (identical(selector, "random-k")) {
    seed <- as.integer(replicate * 100000 + sum(utf8ToInt(dataset_name)))
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
      get(".Random.seed", envir = .GlobalEnv)
    } else {
      NULL
    }
    on.exit({
      if (!is.null(old_seed)) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
    return(sort(sample(clusters, k)))
  }

  if (identical(selector, "confidence-k")) {
    confidence <- first_pass_evidence$LLMConfidence %||% first_pass_evidence$Confidence
    confidence <- as_confidence(confidence)
    ord <- order(confidence, decreasing = FALSE, na.last = TRUE)
    return(clusters[head(ord, k)])
  }

  stop("Unknown refinement selector: ", selector, call. = FALSE)
}

run_refinement_control <- function(selector,
                                   selected_annotations,
                                   dataset_name,
                                   replicate,
                                   markers,
                                   tissue,
                                   species,
                                   api_key,
                                   model_key = "deepseek",
                                   cache_slug = model_key) {
  selected_annotations <- selected_annotations[
    !duplicated(selected_annotations$Cluster),
    ,
    drop = FALSE
  ]
  selected_clusters <- as.character(selected_annotations$Cluster)
  empty_refined <- selected_annotations[0, , drop = FALSE]

  out <- list(
    selector = selector,
    selected_clusters = selected_clusters,
    budget_k = length(selected_clusters),
    refined_annotations = empty_refined,
    runtime_sec = 0,
    cost_usd = 0,
    tokens = 0,
    prompt_tokens = 0,
    completion_tokens = 0,
    second_pass_calls = 0,
    response_hash = NA_character_,
    prompt_hash = NA_character_
  )

  if (nrow(selected_annotations) == 0) {
    return(out)
  }

  refinement_dir <- file.path("results", "benchmark_debug", "refinement")
  selector_slug <- gsub("[^A-Za-z0-9]+", "_", tolower(selector))
  refinement_cache <- file.path(
    refinement_dir,
    paste0(
      "rep", replicate, "_", dataset_name, "_", cache_slug,
      "_refinement_", selector_slug, ".rds"
    )
  )
  cached <- read_llm_response_cache(refinement_cache)

  refinement_prompt <- create_refinement_prompt(
    markers = markers,
    annotations = selected_annotations,
    tissue = tissue,
    species = species
  )

  if (!is.null(cached)) {
    content <- cached$response_text %||% ""
    out$runtime_sec <- safe_usage_number(cached$metadata$runtime_sec)
    out$cost_usd <- safe_usage_number(cached$metadata$cost_usd)
    out$tokens <- safe_usage_number(cached$metadata$tokens)
    out$prompt_tokens <- safe_usage_number(cached$metadata$prompt_tokens)
    out$completion_tokens <- safe_usage_number(cached$metadata$completion_tokens)
    out$second_pass_calls <- 1
    out$response_hash <- cached$response_hash %||% hash_text_md5(content)
    out$prompt_hash <- cached$metadata$prompt_hash %||% hash_text_md5(refinement_prompt)
    refinement_success <- nzchar(trimws(content))
  } else {
    refinement_res <- call_llm_api(refinement_prompt, model_key, api_key)
    content <- refinement_res$content %||% ""

    out$runtime_sec <- safe_usage_number(refinement_res$runtime_sec)
    out$cost_usd <- safe_usage_number(refinement_res$cost_usd)
    out$tokens <- safe_usage_number(refinement_res$usage$total_tokens)
    out$prompt_tokens <- safe_usage_number(refinement_res$usage$prompt_tokens)
    out$completion_tokens <- safe_usage_number(refinement_res$usage$completion_tokens)
    out$second_pass_calls <- 1
    out$response_hash <- hash_text_md5(content)
    out$prompt_hash <- hash_text_md5(refinement_prompt)
    refinement_success <- isTRUE(refinement_res$success)

    write_first_pass_cache(
      content,
      refinement_cache,
      metadata = list(
        replicate = replicate,
        dataset = dataset_name,
        selector = selector,
        model_key = model_key,
        model_id = MODELS[[model_key]]$model_id %||% NA_character_,
        n_selected = nrow(selected_annotations),
        prompt_hash = out$prompt_hash,
        runtime_sec = out$runtime_sec,
        tokens = out$tokens,
        prompt_tokens = out$prompt_tokens,
        completion_tokens = out$completion_tokens,
        cost_usd = out$cost_usd
      )
    )
  }

  if (isTRUE(refinement_success)) {
    refined_annotations <- parse_ablation_annotation_response(content)
    refined_match <- match_ablation_clusters(refined_annotations$Cluster, selected_clusters)
    keep <- !is.na(refined_match)
    out$refined_annotations <- refined_annotations[keep, , drop = FALSE]
    if (nrow(out$refined_annotations) > 0) {
      out$refined_annotations$Cluster <- selected_clusters[refined_match[keep]]
    }
  }

  out
}

run_llm_ablation_wrapper <- function(dataset_name,
                                     data,
                                     ont_data,
                                     replicate,
                                     model_key = "deepseek",
                                     api_key = NULL,
                                     method_prefix = MODELS[[model_key]]$name %||% model_key,
                                     cache_slug = model_key,
                                     include_full_refinement = TRUE) {
  markers <- data$markers
  tissue <- data$tissue
  species <- data$species
  model_id <- MODELS[[model_key]]$model_id %||% NA_character_
  method_names <- ablation_method_names(method_prefix)

  first_pass_dir <- file.path("results", "benchmark_debug", "first_pass")
  first_pass_cache <- file.path(
    first_pass_dir,
    paste0("rep", replicate, "_", dataset_name, "_", cache_slug, "_first_pass.rds")
  )
  cached_first_pass <- read_llm_response_cache(first_pass_cache)

  if (!is.null(cached_first_pass)) {
    method_res <- list(
      predictions = NULL,
      response_text = cached_first_pass$response_text %||% "",
      prompt_text = "",
      runtime_sec = safe_usage_number(cached_first_pass$metadata$runtime_sec),
      cost_usd = safe_usage_number(cached_first_pass$metadata$cost_usd),
      tokens = safe_usage_number(cached_first_pass$metadata$tokens),
      prompt_tokens = safe_usage_number(cached_first_pass$metadata$prompt_tokens),
      completion_tokens = safe_usage_number(cached_first_pass$metadata$completion_tokens)
    )
    first_pass_hash <- cached_first_pass$response_hash %||%
      hash_text_md5(method_res$response_text %||% "")
  } else {
    method_res <- run_llm_annotation(
      markers,
      tissue,
      species,
      dataset_name,
      model_key,
      api_key
    )

    first_pass_hash <- hash_text_md5(method_res$response_text %||% "")
    write_first_pass_cache(
      method_res$response_text %||% "",
      first_pass_cache,
      metadata = list(
        replicate = replicate,
        dataset = dataset_name,
        tissue = tissue,
        species = species,
        model_key = model_key,
        model_id = model_id,
        prompt_hash = hash_text_md5(method_res$prompt_text %||% ""),
        runtime_sec = method_res$runtime_sec,
        tokens = method_res$tokens,
        prompt_tokens = method_res$prompt_tokens %||% NA_real_,
        completion_tokens = method_res$completion_tokens %||% NA_real_,
        cost_usd = method_res$cost_usd
      )
    )
  }

  first_pass_annotations <- prepare_first_pass_annotations(
    method_res$response_text %||% "",
    markers = markers,
    tissue = tissue,
    ontology_data = ont_data
  )

  first_pass_evidence <- calibrate_annotation_confidence(
    first_pass_annotations,
    markers,
    tissue = tissue
  )
  first_pass_no_ontology_evidence <- calibrate_annotation_confidence(
    first_pass_annotations,
    markers,
    tissue = tissue,
    use_ontology_evidence = FALSE
  )
  evidence_clusters <- select_refinement_clusters(
    first_pass_evidence,
    selector = "evidence-k",
    k = nrow(first_pass_evidence),
    replicate = replicate,
    dataset_name = dataset_name
  )
  evidence_k <- length(evidence_clusters)
  no_ontology_clusters <- select_refinement_clusters(
    first_pass_no_ontology_evidence,
    selector = "evidence-no-ontology-k",
    k = evidence_k,
    replicate = replicate,
    dataset_name = dataset_name
  )

  selector_clusters <- list()
  selector_clusters[[method_names$random]] <- select_refinement_clusters(
    first_pass_evidence,
    selector = "random-k",
    k = evidence_k,
    replicate = replicate,
    dataset_name = dataset_name
  )
  selector_clusters[[method_names$confidence]] <- select_refinement_clusters(
    first_pass_evidence,
    selector = "confidence-k",
    k = evidence_k,
    replicate = replicate,
    dataset_name = dataset_name
  )
  selector_clusters[[method_names$no_ontology]] <- no_ontology_clusters
  selector_clusters[[method_names$self]] <- evidence_clusters
  if (isTRUE(include_full_refinement)) {
    selector_clusters[[method_names$full]] <- select_refinement_clusters(
      first_pass_evidence,
      selector = "full",
      k = nrow(first_pass_evidence),
      replicate = replicate,
      dataset_name = dataset_name
    )
  }

  selector_names <- stats::setNames(
    c(
      "random-k",
      "confidence-k",
      "evidence-no-ontology-k",
      "evidence-k",
      if (isTRUE(include_full_refinement)) "full" else character()
    ),
    c(
      method_names$random,
      method_names$confidence,
      method_names$no_ontology,
      method_names$self,
      if (isTRUE(include_full_refinement)) method_names$full else character()
    )
  )

  refinement_runs <- lapply(names(selector_clusters), function(arm_name) {
    evidence_source <- if (identical(arm_name, method_names$no_ontology)) {
      disable_ontology_annotation_context(first_pass_no_ontology_evidence)
    } else {
      first_pass_evidence
    }
    selected <- evidence_source[
      evidence_source$Cluster %in% selector_clusters[[arm_name]],
      ,
      drop = FALSE
    ]
    run_refinement_control(
      selector = selector_names[[arm_name]],
      selected_annotations = selected,
      dataset_name = dataset_name,
      replicate = replicate,
      markers = markers,
      tissue = tissue,
      species = species,
      api_key = api_key,
      model_key = model_key,
      cache_slug = cache_slug
    )
  })
  names(refinement_runs) <- names(selector_clusters)

  arms <- create_deepseekcell_ablation_arms(
    first_pass_annotations = first_pass_annotations,
    markers = markers,
    tissue = tissue,
    ontology_data = ont_data,
    refined_annotations = refinement_runs[[method_names$self]]$refined_annotations,
    include_self_refined = TRUE
  )
  arms <- rename_ablation_arms(arms, method_names)
  arms[[method_names$random]] <- create_refined_ablation_arm(
    first_pass_annotations = first_pass_annotations,
    markers = markers,
    tissue = tissue,
    ontology_data = ont_data,
    refined_annotations = refinement_runs[[method_names$random]]$refined_annotations,
    selected_clusters = refinement_runs[[method_names$random]]$selected_clusters,
    arm_name = method_names$random
  )
  arms[[method_names$confidence]] <- create_refined_ablation_arm(
    first_pass_annotations = first_pass_annotations,
    markers = markers,
    tissue = tissue,
    ontology_data = ont_data,
    refined_annotations = refinement_runs[[method_names$confidence]]$refined_annotations,
    selected_clusters = refinement_runs[[method_names$confidence]]$selected_clusters,
    arm_name = method_names$confidence
  )
  arms[[method_names$no_ontology]] <- create_refined_ablation_arm(
    first_pass_annotations = first_pass_annotations,
    markers = markers,
    tissue = tissue,
    ontology_data = ont_data,
    refined_annotations = refinement_runs[[method_names$no_ontology]]$refined_annotations,
    selected_clusters = refinement_runs[[method_names$no_ontology]]$selected_clusters,
    arm_name = method_names$no_ontology
  )
  if (isTRUE(include_full_refinement)) {
    arms[[method_names$full]] <- create_refined_ablation_arm(
      first_pass_annotations = first_pass_annotations,
      markers = markers,
      tissue = tissue,
      ontology_data = ont_data,
      refined_annotations = refinement_runs[[method_names$full]]$refined_annotations,
      selected_clusters = refinement_runs[[method_names$full]]$selected_clusters,
      arm_name = method_names$full
    )
  }

  evaluated <- lapply(names(arms), function(arm_name) {
    refinement_run <- refinement_runs[[arm_name]]
    if (is.null(refinement_run)) {
      refinement_run <- list(
        selector = "none",
        selected_clusters = character(),
        budget_k = 0,
        runtime_sec = 0,
        cost_usd = 0,
        tokens = 0,
        second_pass_calls = 0
      )
    }
    evaluate_ablation_arm(
      annotations = arms[[arm_name]],
      arm_name = arm_name,
      dataset_name = dataset_name,
      data = data,
      ont_data = ont_data,
      replicate = replicate,
      method_res = method_res,
      first_pass_hash = first_pass_hash,
      llm_backend = model_key,
      llm_model_id = model_id,
      refinement_selector = refinement_run$selector,
      refinement_budget_k = refinement_run$budget_k,
      selected_clusters = refinement_run$selected_clusters,
      refinement_runtime_sec = refinement_run$runtime_sec,
      refinement_cost_usd = refinement_run$cost_usd,
      refinement_tokens = refinement_run$tokens,
      second_pass_calls = refinement_run$second_pass_calls
    )
  })
  names(evaluated) <- names(arms)

  plain_behavior <- arms[[method_names$plain]]
  truth <- data$truth[plain_behavior$Cluster]
  truth_h <- harmonise_labels(truth, tissue, ont_data)
  names(truth_h) <- plain_behavior$Cluster
  plain_behavior$CellType <- harmonise_labels(plain_behavior$CellType, tissue, ont_data)

  refinement_behavior <- dplyr::bind_rows(lapply(names(refinement_runs), function(arm_name) {
    refined_behavior <- arms[[arm_name]]
    refined_behavior$CellType <- harmonise_labels(refined_behavior$CellType, tissue, ont_data)
    refinement_run <- refinement_runs[[arm_name]]

    compute_refinement_behavior(
      plain_behavior,
      refined_behavior,
      truth_h
    ) %>%
      dplyr::mutate(
        Dataset = dataset_name,
        Tissue = tissue,
        Species = species,
        Method = arm_name,
        LLMBackend = model_key,
        LLMModelID = model_id,
        RefinementSelector = refinement_run$selector,
        RefinementBudgetK = refinement_run$budget_k,
        Replicate = replicate,
        FirstPassResponseHash = first_pass_hash,
        RefinementResponseHash = refinement_run$response_hash,
        RefinementPromptTokens = refinement_run$prompt_tokens,
        RefinementCompletionTokens = refinement_run$completion_tokens,
        RefinementTokens = refinement_run$tokens,
        RefinementRuntimeSec = refinement_run$runtime_sec,
        RefinementCostUSD = refinement_run$cost_usd,
        SecondPassCalls = refinement_run$second_pass_calls,
        CostPerCorrectedError = safe_divide(
          refinement_run$cost_usd,
          .data$WrongToCorrect
        ),
        .before = 1
      )
  }))

  list(
    results = dplyr::bind_rows(lapply(evaluated, `[[`, "result")),
    debug = dplyr::bind_rows(lapply(evaluated, `[[`, "debug")),
    confidence_quality = dplyr::bind_rows(lapply(evaluated, `[[`, "confidence_quality")),
    reliability = dplyr::bind_rows(lapply(evaluated, `[[`, "reliability")),
    refinement_behavior = refinement_behavior
  )
}

run_deepseek_ablation_wrapper <- function(dataset_name,
                                          data,
                                          ont_data,
                                          replicate,
                                          api_key) {
  run_llm_ablation_wrapper(
    dataset_name = dataset_name,
    data = data,
    ont_data = ont_data,
    replicate = replicate,
    model_key = "deepseek",
    api_key = api_key,
    method_prefix = "DeepSeek",
    cache_slug = "deepseek",
    include_full_refinement = TRUE
  )
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
  native_cell_metrics <- empty_native_cell_metrics()

  if (identical(method, "CellTypist")) {
    native_cell_metrics <- compute_native_cell_metrics(
      pred,
      seu,
      tissue,
      ont_data
    )
  }

  if (method %in% c("SingleR", "scType", "scmap") ||
      identical(method_res$prediction_level %||% NULL, "cell")) {
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
    RefinementSelector = "baseline",
    RefinementBudgetK = 0,
    NClusters = length(markers),
    RuntimeSec = method_res$runtime_sec,
    CostUSD = method_res$cost_usd,
    Tokens = method_res$tokens,
    PredictionLevel = method_res$prediction_level %||% "cluster",
    CellTypistModel = method_res$model %||% NA_character_,
    FirstPassRuntimeSec = NA_real_,
    RefinementRuntimeSec = 0,
    FirstPassTokens = NA_real_,
    RefinementTokens = 0,
    SecondPassCalls = 0,
    FirstPassResponseHash = NA_character_,
    MeanClusterPurity = mean(cluster_purity, na.rm = TRUE),
    MinClusterPurity = min(cluster_purity, na.rm = TRUE),
    t(native_cell_metrics),
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
  ablation_confidence_quality <- list()
  ablation_reliability <- list()
  ablation_refinement_behavior <- list()

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

        is_llm_ablation <- m %in% c("DeepSeek", "Ollama")
        api_key <- if (identical(m, "DeepSeek")) {
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

        if (is_llm_ablation) {
          model_key <- tolower(m)
          include_full_refinement <- !identical(model_key, "ollama") ||
            include_ollama_full_refinement()
          ablation <- tryCatch(
            run_llm_ablation_wrapper(
              dataset_name = ds,
              data = datasets[[ds]],
              ont_data = ont_data,
              replicate = rep,
              model_key = model_key,
              api_key = api_key,
              method_prefix = m,
              cache_slug = model_key,
              include_full_refinement = include_full_refinement
            ),
            error = function(e) {
              warning("Failed ", m, " ablation on ", ds, ": ", e$message)
              NULL
            }
          )

          if (!is.null(ablation)) {
            ablation$results$Replicate <- rep
            all_results[[length(all_results) + 1]] <- ablation$results
            ablation_confidence_quality[[length(ablation_confidence_quality) + 1]] <-
              ablation$confidence_quality
            ablation_reliability[[length(ablation_reliability) + 1]] <-
              ablation$reliability
            ablation_refinement_behavior[[length(ablation_refinement_behavior) + 1]] <-
              ablation$refinement_behavior
          }
        } else {
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
          RuntimeSec, CostUSD, Tokens, FirstPassRuntimeSec,
          RefinementRuntimeSec, FirstPassTokens, RefinementTokens,
          SecondPassCalls, RefinementBudgetK, MeanClusterPurity, MinClusterPurity,
          EvaluatedClusters, NativeCellARI, NativeCellMacroF1,
          NativeCellAccuracy, NativeCellBalancedAcc, NativeCellCladeAcc,
          NativeCellUnknownRate, NativeCellEvaluatedCells
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

  refinement_behavior <- bind_rows(ablation_refinement_behavior)
  if (nrow(refinement_behavior) > 0) {
    refinement_behavior <- refinement_behavior %>%
      dplyr::mutate(
        RefinementFamily = refinement_method_family(.data$Method)
      )

    full_reference <- refinement_behavior %>%
      dplyr::filter(grepl("-FullRefined$", .data$Method)) %>%
      dplyr::select(
        .data$Replicate,
        .data$Dataset,
        .data$RefinementFamily,
        FullWrongToCorrect = .data$WrongToCorrect,
        FullNRefined = .data$NRefined
      )

    refinement_behavior <- refinement_behavior %>%
      dplyr::left_join(
        full_reference,
        by = c("Replicate", "Dataset", "RefinementFamily")
      ) %>%
      dplyr::mutate(
        CorrectionEfficiency = safe_divide(
          .data$WrongToCorrect,
          .data$NRefined
        ),
        HarmRate = safe_divide(
          .data$CorrectToWrong,
          .data$NRefined
        ),
        RecoveryFraction = safe_divide(
          .data$WrongToCorrect,
          .data$FullWrongToCorrect
        ),
        RelativeRefinementBudget = safe_divide(
          .data$NRefined,
          .data$FullNRefined
        )
      )
  }

  list(
    full = full,
    summary = summary,
    dataset_summary = bind_rows(dataset_summaries),
    cluster_summary = bind_rows(cluster_summaries),
    ablation_confidence_quality = bind_rows(ablation_confidence_quality),
    ablation_reliability = bind_rows(ablation_reliability),
    ablation_refinement_behavior = refinement_behavior
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

metric_ci <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(c(NA_real_, NA_real_))
  }
  stats::quantile(x, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE)
}

bootstrap_refinement_efficiency <- function(refinement_behavior,
                                            iterations = 2000,
                                            seed = 20260719) {
  if (is.null(refinement_behavior) ||
      nrow(refinement_behavior) == 0 ||
      iterations <= 0) {
    return(data.frame())
  }

  if (!"RefinementFamily" %in% names(refinement_behavior)) {
    refinement_behavior <- refinement_behavior %>%
      dplyr::mutate(
        RefinementFamily = refinement_method_family(.data$Method)
      )
  }

  block_keys <- unique(refinement_behavior[c("Replicate", "Dataset")])
  block_keys$BlockId <- seq_len(nrow(block_keys))

  behavior <- dplyr::left_join(
    refinement_behavior,
    block_keys,
    by = c("Replicate", "Dataset")
  )

  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
    get(".Random.seed", envir = .GlobalEnv)
  } else {
    NULL
  }
  on.exit({
    if (!is.null(old_seed)) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)

  boot <- dplyr::bind_rows(lapply(seq_len(iterations), function(i) {
    sampled_blocks <- sample(block_keys$BlockId, nrow(block_keys), replace = TRUE)
    sampled <- dplyr::bind_rows(lapply(sampled_blocks, function(block_id) {
      behavior[behavior$BlockId == block_id, , drop = FALSE]
    }))

    full_reference <- sampled %>%
      dplyr::filter(grepl("-FullRefined$", .data$Method)) %>%
      dplyr::group_by(.data$RefinementFamily) %>%
      dplyr::summarise(
        FullWrongToCorrect = sum(.data$WrongToCorrect, na.rm = TRUE),
        FullNRefined = sum(.data$NRefined, na.rm = TRUE),
        .groups = "drop"
      )

    sampled %>%
      dplyr::group_by(.data$RefinementFamily, .data$Method, .data$RefinementSelector) %>%
      dplyr::summarise(
        TotalRefined = sum(.data$NRefined, na.rm = TRUE),
        TotalWrongToCorrect = sum(.data$WrongToCorrect, na.rm = TRUE),
        TotalCorrectToWrong = sum(.data$CorrectToWrong, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::left_join(
        full_reference,
        by = "RefinementFamily"
      ) %>%
      dplyr::mutate(
        BootstrapIteration = i,
        CorrectionEfficiency = safe_divide(
          .data$TotalWrongToCorrect,
          .data$TotalRefined
        ),
        NetCorrectionRate = safe_divide(
          .data$TotalWrongToCorrect - .data$TotalCorrectToWrong,
          .data$TotalRefined
        ),
        RecoveryFraction = safe_divide(
          .data$TotalWrongToCorrect,
          .data$FullWrongToCorrect
        ),
        RelativeRefinementBudget = safe_divide(
          .data$TotalRefined,
          .data$FullNRefined
        )
      )
  }))

  boot %>%
    dplyr::group_by(.data$RefinementFamily, .data$Method, .data$RefinementSelector) %>%
    dplyr::summarise(
      CorrectionEfficiency_ci95_low = metric_ci(.data$CorrectionEfficiency)[1],
      CorrectionEfficiency_ci95_high = metric_ci(.data$CorrectionEfficiency)[2],
      NetCorrectionRate_ci95_low = metric_ci(.data$NetCorrectionRate)[1],
      NetCorrectionRate_ci95_high = metric_ci(.data$NetCorrectionRate)[2],
      RecoveryFraction_ci95_low = metric_ci(.data$RecoveryFraction)[1],
      RecoveryFraction_ci95_high = metric_ci(.data$RecoveryFraction)[2],
      RelativeRefinementBudget_ci95_low = metric_ci(.data$RelativeRefinementBudget)[1],
      RelativeRefinementBudget_ci95_high = metric_ci(.data$RelativeRefinementBudget)[2],
      .groups = "drop"
    )
}

summarise_refinement_efficiency <- function(refinement_behavior,
                                            bootstrap_iterations = 2000) {
  if (is.null(refinement_behavior) || nrow(refinement_behavior) == 0) {
    return(data.frame())
  }

  if (!"RefinementFamily" %in% names(refinement_behavior)) {
    refinement_behavior <- refinement_behavior %>%
      dplyr::mutate(
        RefinementFamily = refinement_method_family(.data$Method)
      )
  }

  full_reference <- refinement_behavior %>%
    dplyr::filter(grepl("-FullRefined$", .data$Method)) %>%
    dplyr::group_by(.data$RefinementFamily) %>%
    dplyr::summarise(
      FullWrongToCorrect = sum(.data$WrongToCorrect, na.rm = TRUE),
      FullNRefined = sum(.data$NRefined, na.rm = TRUE),
      .groups = "drop"
    )

  summary <- refinement_behavior %>%
    dplyr::group_by(.data$RefinementFamily, .data$Method, .data$RefinementSelector) %>%
    dplyr::summarise(
      Blocks = dplyr::n(),
      TotalRefined = sum(.data$NRefined, na.rm = TRUE),
      TotalWrongToCorrect = sum(.data$WrongToCorrect, na.rm = TRUE),
      TotalCorrectToWrong = sum(.data$CorrectToWrong, na.rm = TRUE),
      MeanSelectionPrecision = mean(.data$SelectionPrecision, na.rm = TRUE),
      MeanSelectionRecall = mean(.data$SelectionRecall, na.rm = TRUE),
      MeanSelectionMCC = mean(.data$SelectionMCC, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::left_join(
      full_reference,
      by = "RefinementFamily"
    ) %>%
    dplyr::mutate(
      CorrectionEfficiency = safe_divide(
        .data$TotalWrongToCorrect,
        .data$TotalRefined
      ),
      HarmRate = safe_divide(
        .data$TotalCorrectToWrong,
        .data$TotalRefined
      ),
      NetCorrectionRate = safe_divide(
        .data$TotalWrongToCorrect - .data$TotalCorrectToWrong,
        .data$TotalRefined
      ),
      RecoveryFraction = safe_divide(
        .data$TotalWrongToCorrect,
        .data$FullWrongToCorrect
      ),
      RelativeRefinementBudget = safe_divide(
        .data$TotalRefined,
        .data$FullNRefined
      )
    )

  ci <- bootstrap_refinement_efficiency(
    refinement_behavior,
    iterations = bootstrap_iterations
  )

  if (nrow(ci) > 0) {
    summary <- dplyr::left_join(
      summary,
      ci,
      by = c("RefinementFamily", "Method", "RefinementSelector")
    )
  }

  summary %>%
    dplyr::arrange(dplyr::desc(.data$CorrectionEfficiency))
}

plot_refinement_efficiency <- function(refinement_summary,
                                       output_pdf = "results/refinement_efficiency.pdf") {
  if (is.null(refinement_summary) || nrow(refinement_summary) == 0) {
    return(invisible(NULL))
  }

  percent_label <- function(x) {
    out <- floor((100 * x) * 10 + 0.5) / 10
    ifelse(is.na(x), NA_character_, sprintf("%.1f%%", out))
  }

  if (!"RefinementFamily" %in% names(refinement_summary)) {
    refinement_summary <- refinement_summary %>%
      dplyr::mutate(
        RefinementFamily = refinement_method_family(.data$Method)
      )
  }

  plot_df <- refinement_summary %>%
    dplyr::mutate(
      Selector = dplyr::recode(
        .data$RefinementSelector,
        "full" = "Full",
        "evidence-k" = "Evidence",
        "evidence-no-ontology-k" = "No ontology-k",
        "random-k" = "Random-k",
        "confidence-k" = "Confidence-k",
        .default = .data$RefinementSelector
      ),
      Selector = factor(
        .data$Selector,
        levels = c("Full", "Evidence", "No ontology-k", "Random-k", "Confidence-k")
      )
    )

  long_df <- dplyr::bind_rows(
    plot_df %>%
      dplyr::transmute(
        RefinementFamily,
        Selector,
        Panel = "Wrong-to-correct revisions",
        Value = .data$TotalWrongToCorrect,
        Label = as.character(.data$TotalWrongToCorrect)
      ),
    plot_df %>%
      dplyr::transmute(
        RefinementFamily,
        Selector,
        Panel = "Corrections per refined cluster",
        Value = .data$CorrectionEfficiency,
        Label = percent_label(.data$CorrectionEfficiency)
      ),
    plot_df %>%
      dplyr::transmute(
        RefinementFamily,
        Selector,
        Panel = "Recovery fraction of full refinement",
        Value = .data$RecoveryFraction,
        Label = percent_label(.data$RecoveryFraction)
      ),
    plot_df %>%
      dplyr::transmute(
        RefinementFamily,
        Selector,
        Panel = "Refinement budget",
        Value = .data$TotalRefined,
        Label = as.character(.data$TotalRefined)
      )
  )

  long_df$Panel <- factor(
    long_df$Panel,
    levels = c(
      "Wrong-to-correct revisions",
      "Corrections per refined cluster",
      "Recovery fraction of full refinement",
      "Refinement budget"
    )
  )
  long_df$Facet <- if (dplyr::n_distinct(long_df$RefinementFamily) > 1) {
    paste(long_df$RefinementFamily, long_df$Panel, sep = ": ")
  } else {
    as.character(long_df$Panel)
  }

  p <- ggplot2::ggplot(
    long_df,
    ggplot2::aes(x = Selector, y = Value, fill = Selector)
  ) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::geom_text(
      ggplot2::aes(label = Label),
      vjust = -0.25,
      size = 3.2
    ) +
    ggplot2::facet_wrap(~ Facet, scales = "free_y", ncol = 2) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
      legend.position = "none",
      panel.grid.major.x = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title = "Evidence-guided selective refinement efficiency",
      subtitle = "Aggregated across benchmark dataset-replicate blocks",
      x = "Refinement selector",
      y = NULL
    ) +
    ggplot2::expand_limits(y = 0)

  plot_height <- if (dplyr::n_distinct(long_df$RefinementFamily) > 1) 10 else 7
  ggplot2::ggsave(output_pdf, p, width = 10, height = plot_height)
  message("Plot saved to ", output_pdf)
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
    "- DeepSeek ablations are paired: one first-pass raw response is cached and reused for all DeepSeekCell ablation arms.",
    "- DeepSeekCell-Evidence adds deterministic evidence fields without changing labels or confidence.",
    "- DeepSeekCell-Calibrated changes confidence but not labels.",
    "- DeepSeekCell-RandomK, DeepSeekCell-ConfidenceK, DeepSeekCell-NoOntologyK, DeepSeekCell-SelfRefined, and DeepSeekCell-FullRefined are allowed to change labels through second-pass refinement.",
    "- DeepSeekCell-NoOntologyK disables ontology-derived evidence during selection and refinement prompting while preserving the same per-dataset refinement budget as Evidence-k.",
    "- RandomK, ConfidenceK, and NoOntologyK use the same per-dataset refinement budget as the evidence-conflict selector; FullRefined refines every cluster as an upper-bound control.",
    "- Refinement behavior metrics report both selector quality and label-change benefit/harm, including wrong-to-correct and correct-to-wrong transitions.",
    "- Refinement efficiency is reported as wrong-to-correct revisions per refined cluster; recovery fraction is reported relative to FullRefined wrong-to-correct revisions.",
    "- Selector-pairwise Wilcoxon tests compare Evidence-k with Random-k and Confidence-k across paired dataset-replicate blocks.",
    "- Although the benchmark uses DeepSeek, the evidence-guided selector operates on structured annotation outputs and deterministic biological evidence, so it is model-agnostic.",
    "- Confidence quality uses binary correctness negative log-likelihood, Brier score, ECE, AUROC, and risk-coverage analysis for selected-label correctness.",
    "- SingleR, scType, scmap, and CellTypist are evaluated as baseline methods when their dependencies are available.",
    "- scmap is run with scmapCell nearest-neighbor label transfer before cluster-level majority voting.",
    "- CellTypist is run in its native cell-level setting using tissue-aware pretrained models unless CELLTYPIST_MODEL overrides the model choice; native cell-level metrics are retained, and cluster-level majority voting is used only for the shared cluster-label comparison.",
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

    if (identical(method, "Ollama") && !is_ollama_available()) {
      message("Skipping Ollama: local Ollama server or requested model is not available.")
      next
    }

    selected <- c(selected, method)
  }

  selected
}

is_ollama_available <- function() {
  model <- MODELS$ollama
  tags_url <- sub("/api/generate$", "/api/tags", model$api_url)

  tryCatch({
    resp <- httr2::request(tags_url) |>
      httr2::req_timeout(3) |>
      httr2::req_perform()
    data <- httr2::resp_body_json(resp, simplifyVector = FALSE)
    model_names <- vapply(data$models %||% list(), function(x) x$name %||% "", character(1))

    if (length(model_names) == 0) {
      return(TRUE)
    }

    model$model_id %in% model_names
  }, error = function(e) {
    FALSE
  })
}

should_auto_run_benchmark <- function() {
  value <- tolower(Sys.getenv("DEEPSEEKCELL_RUN_BENCHMARK_ON_SOURCE", unset = "true"))
  value %in% c("1", "true", "yes", "y")
}

main <- function(n_replicates = 1) {
  deepseek_key <- Sys.getenv("DEEPSEEK_API_KEY")

  ont_data <- load_benchmark_ontology(ONTOLOGY_FILE)
  sctype_db_path <- find_sctype_db()

  extra_llm_methods <- trimws(unlist(strsplit(
    Sys.getenv("DEEPSEEKCELL_EXTRA_LLM_METHODS", unset = ""),
    ","
  )))
  extra_llm_methods <- extra_llm_methods[nzchar(extra_llm_methods)]

  methods <- available_benchmark_methods(
    c("DeepSeek", extra_llm_methods, "scType", "SingleR", "scmap", "CellTypist")
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

  celltypist_native <- results$summary %>%
    dplyr::filter(Method == "CellTypist") %>%
    dplyr::select(
      Dataset,
      Tissue,
      Species,
      Method,
      NReplicates,
      SuccessfulRuns,
      NativeCellARI_mean,
      NativeCellARI_ci95,
      NativeCellMacroF1_mean,
      NativeCellMacroF1_ci95,
      NativeCellAccuracy_mean,
      NativeCellAccuracy_ci95,
      NativeCellBalancedAcc_mean,
      NativeCellBalancedAcc_ci95,
      NativeCellCladeAcc_mean,
      NativeCellCladeAcc_ci95,
      NativeCellUnknownRate_mean,
      NativeCellEvaluatedCells_mean,
      RuntimeSec_mean
    )

  write.csv(
    celltypist_native,
    "results/celltypist_native_cell_metrics.csv",
    row.names = FALSE
  )

  write.csv(results$dataset_summary, "results/dataset_summary.csv", row.names = FALSE)
  write.csv(results$cluster_summary, "results/cluster_summary.csv", row.names = FALSE)
  write.csv(
    results$ablation_confidence_quality,
    "results/ablation_confidence_quality.csv",
    row.names = FALSE
  )
  write.csv(
    results$ablation_reliability,
    "results/ablation_reliability_bins.csv",
    row.names = FALSE
  )
  write.csv(
    results$ablation_refinement_behavior,
    "results/ablation_refinement_behavior.csv",
    row.names = FALSE
  )

  refinement_efficiency_summary <- summarise_refinement_efficiency(
    results$ablation_refinement_behavior
  )
  write.csv(
    refinement_efficiency_summary,
    "results/refinement_efficiency_summary.csv",
    row.names = FALSE
  )

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
      NativeCellMacroF1_mean,
      NativeCellAccuracy_mean,
      NativeCellBalancedAcc_mean,
      NativeCellEvaluatedCells_mean,
      FirstPassRuntimeSec_mean,
      RefinementRuntimeSec_mean,
      SecondPassCalls_mean,
      MeanClusterPurity_mean
    ) %>%
    dplyr::arrange(Dataset, desc(MacroF1_mean))

  write.csv(
    final_table,
    "results/final_benchmark_table.csv",
    row.names = FALSE
  )

  statistical_tests <- run_benchmark_statistical_tests(
    results$full,
    refinement_behavior = results$ablation_refinement_behavior
  )
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
    statistical_tests$selector_pairwise,
    "results/refinement_selector_pairwise_wilcoxon.csv",
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
  plot_refinement_efficiency(
    refinement_efficiency_summary,
    output_pdf = "results/refinement_efficiency.pdf"
  )
  plot_global_metric(results$summary, "ARI", "results/global_ARI.pdf")
  plot_global_metric(results$summary, "MacroF1", "results/global_MacroF1.pdf")
  plot_global_metric(results$summary, "Accuracy", "results/global_Accuracy.pdf")
  plot_global_metric(results$summary, "RuntimeSec", "results/global_Runtime.pdf")

  invisible(results)
}

if (!interactive() && should_auto_run_benchmark()) {
  args <- commandArgs(trailingOnly = TRUE)
  n_rep <- if (length(args) > 0) {
    as.integer(args[1])
  } else {
    DEFAULT_BENCHMARK_REPLICATES
  }
  main(n_replicates = n_rep)
}
