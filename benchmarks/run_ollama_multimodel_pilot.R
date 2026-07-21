#!/usr/bin/env Rscript

if (!file.exists("benchmarks/run_benchmark.R")) {
  stop("Run this script from the DeepSeekCell repository root.", call. = FALSE)
}

Sys.setenv(
  DEEPSEEKCELL_RUN_BENCHMARK_ON_SOURCE = "false",
  DEEPSEEKCELL_USE_LLM_CACHE = Sys.getenv("DEEPSEEKCELL_USE_LLM_CACHE", unset = "true")
)

source("benchmarks/run_benchmark.R")

split_csv <- function(x) {
  out <- trimws(unlist(strsplit(x %||% "", ",")))
  out[nzchar(out)]
}

ollama_model_slug <- function(model_id) {
  gsub("[^A-Za-z0-9]+", "_", tolower(sub(":latest$", "", model_id)))
}

list_installed_ollama_models <- function() {
  tags_url <- sub("/api/generate$", "/api/tags", MODELS$ollama$api_url)
  tryCatch({
    resp <- httr2::request(tags_url) |>
      httr2::req_timeout(3) |>
      httr2::req_perform()
    data <- httr2::resp_body_json(resp, simplifyVector = FALSE)
    vapply(data$models %||% list(), function(x) x$name %||% "", character(1))
  }, error = function(e) character())
}

pairwise_cluster_overlap <- function(selected_clusters) {
  evidence <- selected_clusters[
    selected_clusters$RefinementSelector == "evidence-k",
    ,
    drop = FALSE
  ]
  models <- sort(unique(evidence$LLMModelID))
  if (length(models) < 2) {
    return(data.frame())
  }

  pairs <- utils::combn(models, 2, simplify = FALSE)
  rows <- lapply(pairs, function(pair) {
    a <- unique(evidence$Cluster[evidence$LLMModelID == pair[[1]]])
    b <- unique(evidence$Cluster[evidence$LLMModelID == pair[[2]]])
    intersection_n <- length(intersect(a, b))
    union_n <- length(union(a, b))

    data.frame(
      Dataset = unique(evidence$Dataset)[1] %||% NA_character_,
      ModelA = pair[[1]],
      ModelB = pair[[2]],
      NSelectedA = length(a),
      NSelectedB = length(b),
      Intersection = intersection_n,
      Union = union_n,
      Jaccard = safe_divide(intersection_n, union_n),
      OverlapOfSmallerSet = safe_divide(intersection_n, min(length(a), length(b))),
      stringsAsFactors = FALSE
    )
  })

  dplyr::bind_rows(rows)
}

summarise_cluster_overlap <- function(overlap) {
  if (nrow(overlap) == 0) {
    return(data.frame())
  }

  data.frame(
    Dataset = unique(overlap$Dataset)[1] %||% NA_character_,
    NModelPairs = nrow(overlap),
    MeanJaccard = mean(overlap$Jaccard, na.rm = TRUE),
    SDJaccard = stats::sd(overlap$Jaccard, na.rm = TRUE),
    MinJaccard = min(overlap$Jaccard, na.rm = TRUE),
    MaxJaccard = max(overlap$Jaccard, na.rm = TRUE),
    MeanOverlapOfSmallerSet = mean(overlap$OverlapOfSmallerSet, na.rm = TRUE),
    SDOverlapOfSmallerSet = stats::sd(overlap$OverlapOfSmallerSet, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

read_deepseek_selected_clusters <- function(dataset_name, replicate_id) {
  debug_file <- file.path(
    "results",
    "benchmark_debug",
    sprintf("rep%s_%s_DeepSeekCell-SelfRefined_debug.csv", replicate_id, dataset_name)
  )
  if (!file.exists(debug_file)) {
    return(data.frame())
  }

  debug <- read.csv(debug_file, check.names = FALSE)
  if (!"SelectedForRefinement" %in% names(debug)) {
    return(data.frame())
  }

  if (!"LLMModelID" %in% names(debug)) {
    debug$LLMModelID <- "deepseek-chat"
  }
  if (!"RefinementSelector" %in% names(debug)) {
    debug$RefinementSelector <- "evidence-k"
  }
  if (!"FirstPassResponseHash" %in% names(debug)) {
    debug$FirstPassResponseHash <- NA_character_
  }
  if (!"EvidenceBestCellType" %in% names(debug)) {
    debug$EvidenceBestCellType <- NA_character_
  }

  debug %>%
    dplyr::filter(.data$SelectedForRefinement) %>%
    dplyr::mutate(Cluster = as.character(.data$Cluster)) %>%
    dplyr::select(
      dplyr::all_of(c(
        "Dataset", "LLMModelID", "Method", "RefinementSelector", "Cluster",
        "RawPrediction", "RawTruth", "EvidenceBestCellType", "FirstPassResponseHash"
      ))
    )
}

plot_selector_efficiency <- function(selector_summary, output_pdf) {
  plot_data <- selector_summary %>%
    dplyr::filter(.data$RefinementSelector %in% c(
      "evidence-k", "evidence-no-ontology-k", "confidence-k", "random-k"
    )) %>%
    dplyr::mutate(
      LLMModelID = factor(
        .data$LLMModelID,
        levels = c("deepseek-chat", "llama3.2:latest", "mistral:latest", "gemma2:2b")
      ),
      RefinementSelector = factor(
        .data$RefinementSelector,
        levels = c("evidence-k", "evidence-no-ontology-k", "confidence-k", "random-k"),
        labels = c("Evidence-k", "No ontology-k", "Confidence-k", "Random-k")
      )
    )

  if (nrow(plot_data) == 0) {
    return(invisible(NULL))
  }

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = .data$LLMModelID,
      y = .data$CorrectionEfficiency,
      fill = .data$RefinementSelector
    )
  ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(width = 0.75),
      width = 0.68
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = scales::percent(.data$CorrectionEfficiency, accuracy = 0.1)),
      position = ggplot2::position_dodge(width = 0.75),
      vjust = -0.35,
      size = 3
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, 1.08)
    ) +
    ggplot2::labs(
      x = "LLM",
      y = "Correction efficiency",
      fill = "Selector"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 25, hjust = 1),
      legend.position = "top"
    )

  ggplot2::ggsave(output_pdf, p, width = 7.2, height = 4.5)
}

plot_evidence_selection_matrix <- function(selected_clusters, output_pdf) {
  plot_data <- selected_clusters %>%
    dplyr::filter(.data$RefinementSelector == "evidence-k") %>%
    dplyr::mutate(
      Cluster = factor(as.character(.data$Cluster), levels = sort(unique(as.character(.data$Cluster)))),
      LLMModelID = factor(
        .data$LLMModelID,
        levels = c("deepseek-chat", "llama3.2:latest", "mistral:latest", "gemma2:2b")
      )
    )

  if (nrow(plot_data) == 0) {
    return(invisible(NULL))
  }

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data$Cluster, y = .data$LLMModelID)
  ) +
    ggplot2::geom_tile(fill = "#2C7FB8", color = "white", linewidth = 0.35) +
    ggplot2::labs(
      x = "Selected cluster",
      y = "LLM"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5)
    )

  ggplot2::ggsave(output_pdf, p, width = 7.2, height = 3.2)
}

plot_cost_performance <- function(cost_performance, output_pdf) {
  plot_data <- cost_performance %>%
    dplyr::mutate(
      Arm = dplyr::case_when(
        grepl("-Plain$", .data$Method) ~ "Plain",
        grepl("Cell-SelfRefined$", .data$Method) ~ "SelfRefined",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(.data$Arm)) %>%
    dplyr::select(
      dplyr::all_of(c("Dataset", "LLMModelID", "Arm", "MacroF1", "RuntimeSec", "CostUSD"))
    ) %>%
    tidyr::pivot_wider(
      names_from = "Arm",
      values_from = c("MacroF1", "RuntimeSec", "CostUSD")
    ) %>%
    dplyr::mutate(
      DeltaMacroF1 = .data$MacroF1_SelfRefined - .data$MacroF1_Plain,
      AddedRuntimeSec = .data$RuntimeSec_SelfRefined - .data$RuntimeSec_Plain,
      AddedCostUSD = .data$CostUSD_SelfRefined - .data$CostUSD_Plain
    ) %>%
    dplyr::select(
      dplyr::all_of(c("LLMModelID", "DeltaMacroF1", "AddedRuntimeSec", "AddedCostUSD"))
    ) %>%
    tidyr::pivot_longer(
      cols = c("DeltaMacroF1", "AddedRuntimeSec", "AddedCostUSD"),
      names_to = "Metric",
      values_to = "Value"
    ) %>%
    dplyr::mutate(
      Metric = factor(
        .data$Metric,
        levels = c("DeltaMacroF1", "AddedRuntimeSec", "AddedCostUSD"),
        labels = c("Delta Macro-F1", "Added runtime (s)", "Added API cost (USD)")
      ),
      LLMModelID = factor(
        .data$LLMModelID,
        levels = c("deepseek-chat", "llama3.2:latest", "mistral:latest", "gemma2:2b")
      )
    )

  if (nrow(plot_data) == 0) {
    return(invisible(NULL))
  }

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data$LLMModelID, y = .data$Value)
  ) +
    ggplot2::geom_col(fill = "#2C7FB8", width = 0.65) +
    ggplot2::facet_wrap(~Metric, scales = "free_y", ncol = 1) +
    ggplot2::labs(x = "LLM", y = NULL) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))

  ggplot2::ggsave(output_pdf, p, width = 7.2, height = 7.0)
}

args <- commandArgs(trailingOnly = TRUE)
arg_or_env <- function(index, env_var, default) {
  if (length(args) >= index && nzchar(args[[index]])) {
    return(args[[index]])
  }
  Sys.getenv(env_var, unset = default)
}

dataset_name <- arg_or_env(1, "OLLAMA_PILOT_DATASET", "ZilionisLung")
replicate_id <- as.integer(arg_or_env(2, "OLLAMA_PILOT_REPLICATE", "1"))
requested_models <- split_csv(arg_or_env(
  3,
  "OLLAMA_PILOT_MODELS",
  "llama3.2:latest,mistral:latest,gemma2:2b,qwen2.5:latest"
))
seed_value <- replicate_id * 100
include_full <- include_ollama_full_refinement()

installed_models <- list_installed_ollama_models()
availability <- data.frame(
  Model = requested_models,
  Installed = requested_models %in% installed_models,
  stringsAsFactors = FALSE
)
models_to_run <- availability$Model[availability$Installed]

if (length(models_to_run) == 0) {
  print(availability)
  stop("None of the requested Ollama models are installed.", call. = FALSE)
}

ont_data <- load_benchmark_ontology(ONTOLOGY_FILE)
dataset <- load_benchmark_dataset(dataset_name, seed = seed_value)

all_results <- list()
all_behavior <- list()
all_debug <- list()

for (model_id in models_to_run) {
  model_slug <- ollama_model_slug(model_id)
  cache_slug <- paste0("ollama_", model_slug)
  method_prefix <- paste0("Ollama_", model_slug)

  message("\n=== Ollama pilot: ", model_id, " on ", dataset_name, " ===")
  MODELS$ollama$model_id <- model_id

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

  all_results[[model_id]] <- result$results
  all_behavior[[model_id]] <- result$refinement_behavior
  all_debug[[model_id]] <- result$debug
}

results <- dplyr::bind_rows(all_results)
behavior <- dplyr::bind_rows(all_behavior)
debug <- dplyr::bind_rows(all_debug)

selector_summary <- behavior %>%
  dplyr::group_by(.data$Dataset, .data$LLMModelID, .data$RefinementSelector) %>%
  dplyr::summarise(
    NRefined = sum(.data$NRefined, na.rm = TRUE),
    WrongToCorrect = sum(.data$WrongToCorrect, na.rm = TRUE),
    CorrectToWrong = sum(.data$CorrectToWrong, na.rm = TRUE),
    CorrectionEfficiency = safe_divide(WrongToCorrect, NRefined),
    RefinementRuntimeSec = sum(.data$RefinementRuntimeSec, na.rm = TRUE),
    RefinementTokens = sum(.data$RefinementTokens, na.rm = TRUE),
    RefinementCostUSD = sum(.data$RefinementCostUSD, na.rm = TRUE),
    .groups = "drop"
  )

  cost_performance <- results %>%
  dplyr::filter(
    grepl("-Plain$|Cell-SelfRefined$|Cell-RandomK$|Cell-ConfidenceK$|Cell-NoOntologyK$|Cell-FullRefined$", .data$Method)
  ) %>%
  dplyr::select(
    dplyr::all_of(c(
      "Dataset", "LLMModelID", "Method", "MacroF1", "Accuracy", "CladeAcc",
      "RuntimeSec", "CostUSD", "Tokens", "RefinementBudgetK", "SecondPassCalls"
    ))
  )

selected_clusters <- debug %>%
  dplyr::filter(.data$SelectedForRefinement) %>%
  dplyr::mutate(Cluster = as.character(.data$Cluster)) %>%
  dplyr::select(
    dplyr::all_of(c(
      "Dataset", "LLMModelID", "Method", "RefinementSelector", "Cluster",
      "RawPrediction", "RawTruth", "EvidenceBestCellType", "FirstPassResponseHash"
    ))
  )

deepseek_selected_clusters <- read_deepseek_selected_clusters(dataset_name, replicate_id)
combined_selected_clusters <- dplyr::bind_rows(
  deepseek_selected_clusters,
  selected_clusters
)
overlap <- pairwise_cluster_overlap(selected_clusters)
combined_overlap <- pairwise_cluster_overlap(combined_selected_clusters)
overlap_summary <- summarise_cluster_overlap(overlap)
combined_overlap_summary <- summarise_cluster_overlap(combined_overlap)

dir.create("results", showWarnings = FALSE)
prefix <- paste0("ollama_multimodel_", tolower(dataset_name))
write.csv(availability, file.path("results", paste0(prefix, "_model_availability.csv")), row.names = FALSE)
write.csv(results, file.path("results", paste0(prefix, "_results.csv")), row.names = FALSE)
write.csv(behavior, file.path("results", paste0(prefix, "_refinement_behavior.csv")), row.names = FALSE)
write.csv(selector_summary, file.path("results", paste0(prefix, "_selector_summary.csv")), row.names = FALSE)
write.csv(cost_performance, file.path("results", paste0(prefix, "_cost_performance.csv")), row.names = FALSE)
write.csv(selected_clusters, file.path("results", paste0(prefix, "_selected_clusters.csv")), row.names = FALSE)
write.csv(overlap, file.path("results", paste0(prefix, "_evidence_overlap.csv")), row.names = FALSE)
write.csv(overlap_summary, file.path("results", paste0(prefix, "_evidence_overlap_summary.csv")), row.names = FALSE)
write.csv(
  combined_selected_clusters,
  file.path("results", paste0(prefix, "_combined_selected_clusters.csv")),
  row.names = FALSE
)
write.csv(
  combined_overlap,
  file.path("results", paste0(prefix, "_combined_evidence_overlap.csv")),
  row.names = FALSE
)
write.csv(
  combined_overlap_summary,
  file.path("results", paste0(prefix, "_combined_evidence_overlap_summary.csv")),
  row.names = FALSE
)

combined_cost_performance <- cost_performance
combined_selector_summary <- selector_summary

if (file.exists("results/benchmark_results_full.csv")) {
  deepseek_results <- read.csv("results/benchmark_results_full.csv", check.names = FALSE)
  if (!"LLMModelID" %in% names(deepseek_results)) {
    deepseek_results$LLMModelID <- ifelse(
      grepl("^DeepSeek", deepseek_results$Method),
      "deepseek-chat",
      NA_character_
    )
  }
  deepseek_cost_performance <- deepseek_results %>%
    dplyr::filter(
      .data$Dataset == dataset_name,
      .data$Replicate == replicate_id,
      .data$Method %in% c(
        "DeepSeek-Plain",
        "DeepSeekCell-SelfRefined",
        "DeepSeekCell-RandomK",
        "DeepSeekCell-ConfidenceK",
        "DeepSeekCell-NoOntologyK",
        "DeepSeekCell-FullRefined"
      )
    ) %>%
    dplyr::select(
      dplyr::all_of(c(
        "Dataset", "LLMModelID", "Method", "MacroF1", "Accuracy", "CladeAcc",
        "RuntimeSec", "CostUSD", "Tokens", "RefinementBudgetK", "SecondPassCalls"
      ))
    )
  combined_cost_performance <- dplyr::bind_rows(
    deepseek_cost_performance,
    cost_performance
  )
}

if (file.exists("results/ablation_refinement_behavior.csv")) {
  deepseek_behavior <- read.csv("results/ablation_refinement_behavior.csv", check.names = FALSE)
  if (!"LLMModelID" %in% names(deepseek_behavior)) {
    deepseek_behavior$LLMModelID <- ifelse(
      grepl("^DeepSeek", deepseek_behavior$Method),
      "deepseek-chat",
      NA_character_
    )
  }
  deepseek_selector_summary <- deepseek_behavior %>%
    dplyr::filter(
      .data$Dataset == dataset_name,
      .data$Replicate == replicate_id,
      .data$RefinementSelector %in% c("random-k", "confidence-k", "evidence-k", "full")
        | .data$RefinementSelector == "evidence-no-ontology-k"
    ) %>%
    dplyr::group_by(.data$Dataset, .data$LLMModelID, .data$RefinementSelector) %>%
    dplyr::summarise(
      NRefined = sum(.data$NRefined, na.rm = TRUE),
      WrongToCorrect = sum(.data$WrongToCorrect, na.rm = TRUE),
      CorrectToWrong = sum(.data$CorrectToWrong, na.rm = TRUE),
      CorrectionEfficiency = safe_divide(WrongToCorrect, NRefined),
      RefinementRuntimeSec = sum(.data$RefinementRuntimeSec, na.rm = TRUE),
      RefinementTokens = sum(.data$RefinementTokens, na.rm = TRUE),
      RefinementCostUSD = sum(.data$RefinementCostUSD, na.rm = TRUE),
      .groups = "drop"
    )
  combined_selector_summary <- dplyr::bind_rows(
    deepseek_selector_summary,
    selector_summary
  )
}

plain_self_delta <- combined_cost_performance %>%
  dplyr::mutate(
    Arm = dplyr::case_when(
      grepl("-Plain$", .data$Method) ~ "Plain",
      grepl("Cell-SelfRefined$", .data$Method) ~ "SelfRefined",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(.data$Arm)) %>%
  dplyr::select(
    dplyr::all_of(c(
      "Dataset", "LLMModelID", "Arm", "MacroF1", "Accuracy", "CladeAcc",
      "RuntimeSec", "CostUSD", "Tokens", "SecondPassCalls"
    ))
  ) %>%
  tidyr::pivot_wider(
    names_from = Arm,
    values_from = c(
      MacroF1,
      Accuracy,
      CladeAcc,
      RuntimeSec,
      CostUSD,
      Tokens,
      SecondPassCalls
    )
  ) %>%
  dplyr::mutate(
    DeltaMacroF1 = .data$MacroF1_SelfRefined - .data$MacroF1_Plain,
    DeltaAccuracy = .data$Accuracy_SelfRefined - .data$Accuracy_Plain,
    DeltaCladeAcc = .data$CladeAcc_SelfRefined - .data$CladeAcc_Plain
  )

write.csv(
  combined_cost_performance,
  file.path("results", paste0(prefix, "_combined_cost_performance.csv")),
  row.names = FALSE
)
write.csv(
  combined_selector_summary,
  file.path("results", paste0(prefix, "_combined_selector_summary.csv")),
  row.names = FALSE
)
write.csv(
  plain_self_delta,
  file.path("results", paste0(prefix, "_plain_selfrefined_delta.csv")),
  row.names = FALSE
)

plot_selector_efficiency(
  combined_selector_summary,
  file.path("results", paste0(prefix, "_selector_efficiency.pdf"))
)
plot_evidence_selection_matrix(
  combined_selected_clusters,
  file.path("results", paste0(prefix, "_evidence_selection_matrix.pdf"))
)
plot_cost_performance(
  combined_cost_performance,
  file.path("results", paste0(prefix, "_cost_performance.pdf"))
)

cat("\nModel availability:\n")
print(availability)
cat("\nSelector summary:\n")
print(selector_summary)
cat("\nCombined selector summary:\n")
print(combined_selector_summary)
cat("\nCost/performance:\n")
print(cost_performance)
cat("\nPlain to SelfRefined deltas:\n")
print(plain_self_delta)
cat("\nEvidence-k selected-cluster overlap:\n")
print(overlap)
cat("\nCombined Evidence-k selected-cluster overlap:\n")
print(combined_overlap)
cat("\nCombined Evidence-k overlap summary:\n")
print(combined_overlap_summary)
