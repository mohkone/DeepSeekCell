#!/usr/bin/env Rscript

if (!file.exists("results")) {
  stop("Run this script from the DeepSeekCell repository root.", call. = FALSE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(tidyr)
})

safe_divide <- function(x, y) {
  ifelse(is.na(y) | y == 0, NA_real_, x / y)
}

read_csvs <- function(pattern) {
  files <- Sys.glob(file.path("results", pattern))
  files <- files[!grepl("ollama_multimodel_all_", basename(files))]
  if (length(files) == 0) {
    return(data.frame())
  }
  tables <- lapply(files, function(path) {
    if (file.info(path)$size == 0) {
      return(data.frame())
    }
    tryCatch({
      x <- suppressMessages(readr::read_csv(path, show_col_types = FALSE))
      x[, !startsWith(names(x), "..."), drop = FALSE]
    }, error = function(e) data.frame())
  })
  tables <- tables[vapply(tables, nrow, integer(1)) > 0]
  if (length(tables) == 0) {
    return(data.frame())
  }
  dplyr::bind_rows(tables)
}

selector_summary <- read_csvs("ollama_multimodel_*_combined_selector_summary.csv")
overlap_pairwise <- read_csvs("ollama_multimodel_*_combined_evidence_overlap.csv")
plain_self_delta <- read_csvs("ollama_multimodel_*_plain_selfrefined_delta.csv")

if (nrow(selector_summary) == 0) {
  stop("No combined selector summaries found.", call. = FALSE)
}

selector_blocks <- selector_summary %>%
  dplyr::filter(.data$RefinementSelector %in% c(
    "evidence-k", "evidence-no-ontology-k", "confidence-k", "random-k"
  )) %>%
  dplyr::mutate(BlockCorrectionEfficiency = safe_divide(.data$WrongToCorrect, .data$NRefined))

selector_totals <- selector_blocks %>%
  dplyr::group_by(.data$RefinementSelector) %>%
  dplyr::summarise(
    NModelDatasetBlocks = dplyr::n(),
    NInformativeBlocks = sum(.data$NRefined > 0, na.rm = TRUE),
    NRefined = sum(.data$NRefined, na.rm = TRUE),
    WrongToCorrect = sum(.data$WrongToCorrect, na.rm = TRUE),
    CorrectToWrong = sum(.data$CorrectToWrong, na.rm = TRUE),
    CorrectionEfficiency = safe_divide(.data$WrongToCorrect, .data$NRefined),
    MeanBlockEfficiency = mean(.data$BlockCorrectionEfficiency[.data$NRefined > 0], na.rm = TRUE),
    SDBlockEfficiency = stats::sd(.data$BlockCorrectionEfficiency[.data$NRefined > 0], na.rm = TRUE),
    .groups = "drop"
  )

overlap_summary <- if (nrow(overlap_pairwise) > 0) {
  overlap_pairwise %>%
    dplyr::summarise(
      NModelPairs = dplyr::n(),
      MeanJaccard = mean(.data$Jaccard, na.rm = TRUE),
      SDJaccard = stats::sd(.data$Jaccard, na.rm = TRUE),
      MinJaccard = min(.data$Jaccard, na.rm = TRUE),
      MaxJaccard = max(.data$Jaccard, na.rm = TRUE),
      MeanOverlapOfSmallerSet = mean(.data$OverlapOfSmallerSet, na.rm = TRUE),
      SDOverlapOfSmallerSet = stats::sd(.data$OverlapOfSmallerSet, na.rm = TRUE)
    )
} else {
  data.frame()
}

dir.create("results", showWarnings = FALSE)
readr::write_csv(selector_summary, "results/ollama_multimodel_all_selector_summary.csv")
readr::write_csv(selector_totals, "results/ollama_multimodel_all_selector_totals.csv")
readr::write_csv(overlap_pairwise, "results/ollama_multimodel_all_evidence_overlap_pairwise.csv")
readr::write_csv(overlap_summary, "results/ollama_multimodel_all_evidence_overlap_summary.csv")
readr::write_csv(plain_self_delta, "results/ollama_multimodel_all_plain_selfrefined_delta.csv")

selector_plot_data <- selector_blocks %>%
  dplyr::filter(.data$NRefined > 0) %>%
  dplyr::mutate(
    RefinementSelector = factor(
      .data$RefinementSelector,
      levels = c("evidence-k", "evidence-no-ontology-k", "confidence-k", "random-k"),
      labels = c("Evidence-k", "No ontology-k", "Confidence-k", "Random-k")
    )
  )

if (nrow(selector_plot_data) > 0) {
  p <- ggplot2::ggplot(
    selector_plot_data,
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
    ggplot2::facet_wrap(~Dataset, nrow = 1, scales = "free_x") +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, 1.08)
    ) +
    ggplot2::labs(
      x = "LLM",
      y = "Correction efficiency",
      fill = "Selector"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
      legend.position = "top"
    )
  ggplot2::ggsave("results/ollama_multimodel_all_selector_efficiency.pdf", p, width = 9.5, height = 4.3)
}

if (nrow(plain_self_delta) > 0) {
  delta_plot_data <- plain_self_delta %>%
    dplyr::select(dplyr::all_of(c("Dataset", "LLMModelID", "DeltaMacroF1", "DeltaAccuracy"))) %>%
    tidyr::pivot_longer(
      cols = c("DeltaMacroF1", "DeltaAccuracy"),
      names_to = "Metric",
      values_to = "Delta"
    ) %>%
    dplyr::mutate(
      Metric = factor(
        .data$Metric,
        levels = c("DeltaMacroF1", "DeltaAccuracy"),
        labels = c("Delta Macro-F1", "Delta accuracy")
      )
    )

  p <- ggplot2::ggplot(
    delta_plot_data,
    ggplot2::aes(x = .data$LLMModelID, y = .data$Delta)
  ) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
    ggplot2::geom_col(fill = "#2C7FB8", width = 0.65) +
    ggplot2::facet_grid(.data$Metric ~ .data$Dataset, scales = "free_x") +
    ggplot2::labs(x = "LLM", y = "SelfRefined - Plain") +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
  ggplot2::ggsave("results/ollama_multimodel_all_plain_selfrefined_delta.pdf", p, width = 9.5, height = 5.8)
}

cat("\nSelector totals:\n")
print(selector_totals)
cat("\nEvidence-k overlap summary:\n")
print(overlap_summary)
