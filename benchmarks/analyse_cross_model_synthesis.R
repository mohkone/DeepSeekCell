#!/usr/bin/env Rscript

# Cross-model selector synthesis for frozen DeepSeekCell benchmark outputs.
# This script performs no API calls. It only reads completed benchmark CSVs and
# summarizes whether evidence-guided selection behaves consistently across LLM
# backends and evaluation scopes.

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

results_dir <- Sys.getenv("DEEPSEEKCELL_RESULTS_DIR", unset = "results")
output_prefix <- file.path(results_dir, "cross_model_selector")

read_optional_csv <- function(path) {
  if (!file.exists(path)) {
    return(data.frame())
  }
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

safe_num <- function(x) {
  out <- suppressWarnings(as.numeric(x %||% 0))
  out[is.na(out)] <- 0
  out
}

safe_divide <- function(numerator, denominator) {
  ifelse(is.na(denominator) | denominator == 0, NA_real_, numerator / denominator)
}

selector_label <- function(x) {
  labels <- c(
    "evidence-k" = "Evidence-k",
    "confidence-k" = "Confidence-k",
    "random-k" = "Random-k",
    "evidence-no-ontology-k" = "NoOntology-k",
    "full" = "FullRefined",
    "risk-k" = "Risk-k"
  )
  out <- labels[as.character(x)]
  out[is.na(out)] <- as.character(x)[is.na(out)]
  unname(out)
}

model_label <- function(x) {
  labels <- c(
    "deepseek-chat" = "DeepSeek",
    "gpt-5" = "GPT-5",
    "llama3.2:latest" = "Llama 3.2",
    "gemma2:2b" = "Gemma2-2B",
    "mistral:latest" = "Mistral"
  )
  out <- labels[as.character(x)]
  out[is.na(out)] <- as.character(x)[is.na(out)]
  unname(out)
}

standardize_behavior <- function(x,
                                 source,
                                 scope,
                                 include_models = NULL,
                                 exclude_models = NULL) {
  if (nrow(x) == 0) {
    return(data.frame())
  }

  required <- c("LLMModelID", "RefinementSelector", "NRefined",
                "WrongToCorrect", "CorrectToWrong")
  missing <- setdiff(required, names(x))
  if (length(missing) > 0) {
    stop(source, " is missing required columns: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }

  if (!"Dataset" %in% names(x)) x$Dataset <- NA_character_
  if (!"Replicate" %in% names(x)) x$Replicate <- 1
  if (!"RefinementTokens" %in% names(x)) x$RefinementTokens <- 0
  if (!"RefinementRuntimeSec" %in% names(x)) x$RefinementRuntimeSec <- 0
  if (!"RefinementCostUSD" %in% names(x)) x$RefinementCostUSD <- 0
  if (!"SecondPassCalls" %in% names(x)) x$SecondPassCalls <- x$NRefined

  if (!is.null(include_models)) {
    x <- x[x$LLMModelID %in% include_models, , drop = FALSE]
  }
  if (!is.null(exclude_models)) {
    x <- x[!x$LLMModelID %in% exclude_models, , drop = FALSE]
  }
  if (nrow(x) == 0) {
    return(data.frame())
  }

  data.frame(
    Source = source,
    Scope = scope,
    Dataset = as.character(x$Dataset),
    Replicate = as.character(x$Replicate),
    LLMModelID = as.character(x$LLMModelID),
    Model = model_label(x$LLMModelID),
    Selector = selector_label(x$RefinementSelector),
    NRefined = safe_num(x$NRefined),
    WrongToCorrect = safe_num(x$WrongToCorrect),
    CorrectToWrong = safe_num(x$CorrectToWrong),
    RefinementTokens = safe_num(x$RefinementTokens),
    RefinementRuntimeSec = safe_num(x$RefinementRuntimeSec),
    RefinementCostUSD = safe_num(x$RefinementCostUSD),
    SecondPassCalls = safe_num(x$SecondPassCalls),
    stringsAsFactors = FALSE
  )
}

aggregate_selector_behavior <- function(x) {
  if (nrow(x) == 0) {
    return(data.frame())
  }

  split_key <- paste(x$Source, x$Scope, x$LLMModelID, x$Model, x$Selector, sep = "\r")
  rows <- lapply(split(x, split_key), function(d) {
    n_refined <- sum(d$NRefined, na.rm = TRUE)
    wrong_to_correct <- sum(d$WrongToCorrect, na.rm = TRUE)
    correct_to_wrong <- sum(d$CorrectToWrong, na.rm = TRUE)
    data.frame(
      Source = d$Source[[1]],
      Scope = d$Scope[[1]],
      Model = d$Model[[1]],
      LLMModelID = d$LLMModelID[[1]],
      Selector = d$Selector[[1]],
      Blocks = nrow(d),
      Datasets = length(unique(stats::na.omit(d$Dataset))),
      Replicates = length(unique(stats::na.omit(d$Replicate))),
      NRefined = n_refined,
      WrongToCorrect = wrong_to_correct,
      CorrectToWrong = correct_to_wrong,
      NetCorrections = wrong_to_correct - correct_to_wrong,
      CorrectionEfficiency = safe_divide(wrong_to_correct, n_refined),
      NetCorrectionEfficiency = safe_divide(wrong_to_correct - correct_to_wrong, n_refined),
      HarmRate = safe_divide(correct_to_wrong, n_refined),
      RefinementTokens = sum(d$RefinementTokens, na.rm = TRUE),
      RefinementRuntimeSec = sum(d$RefinementRuntimeSec, na.rm = TRUE),
      RefinementCostUSD = sum(d$RefinementCostUSD, na.rm = TRUE),
      SecondPassCalls = sum(d$SecondPassCalls, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  row.names(out) <- NULL
  out[order(out$Model, out$Selector, out$Scope), , drop = FALSE]
}

wide_net_table <- function(summary) {
  evidence <- summary[summary$Selector == "Evidence-k", , drop = FALSE]
  if (nrow(evidence) == 0) {
    return(data.frame())
  }

  selectors <- c("Random-k", "Confidence-k", "NoOntology-k", "FullRefined")
  out <- evidence[, c(
    "Source", "Scope", "Model", "LLMModelID", "Blocks", "Datasets",
    "Replicates", "NRefined", "WrongToCorrect", "CorrectToWrong",
    "NetCorrections", "CorrectionEfficiency", "NetCorrectionEfficiency",
    "HarmRate", "RefinementTokens", "RefinementCostUSD"
  )]
  names(out) <- c(
    "Source", "Scope", "Model", "LLMModelID", "Blocks", "Datasets",
    "Replicates", "EvidenceRefined", "EvidenceWrongToCorrect",
    "EvidenceCorrectToWrong", "EvidenceNet", "EvidenceCorrectionEfficiency",
    "EvidenceNetCorrectionEfficiency", "EvidenceHarmRate",
    "EvidenceRefinementTokens", "EvidenceRefinementCostUSD"
  )

  for (selector in selectors) {
    d <- summary[summary$Selector == selector, c(
      "Source", "Scope", "LLMModelID", "NRefined", "WrongToCorrect",
      "CorrectToWrong", "NetCorrections", "NetCorrectionEfficiency"
    ), drop = FALSE]
    prefix <- gsub("[^A-Za-z0-9]", "", selector)
    names(d)[4:8] <- paste0(prefix, c(
      "Refined", "WrongToCorrect", "CorrectToWrong", "Net",
      "NetCorrectionEfficiency"
    ))
    out <- merge(out, d, by = c("Source", "Scope", "LLMModelID"), all.x = TRUE, sort = FALSE)
  }

  out[order(out$Model, out$Scope), , drop = FALSE]
}

contrast_table <- function(summary) {
  evidence <- summary[summary$Selector == "Evidence-k", , drop = FALSE]
  comparators <- summary[summary$Selector %in% c(
    "Random-k", "Confidence-k", "NoOntology-k", "FullRefined"
  ), , drop = FALSE]

  if (nrow(evidence) == 0 || nrow(comparators) == 0) {
    return(data.frame())
  }

  merged <- merge(
    evidence,
    comparators,
    by = c("Source", "Scope", "LLMModelID", "Model"),
    suffixes = c("_Evidence", "_Comparator"),
    all = FALSE,
    sort = FALSE
  )

  data.frame(
    Source = merged$Source,
    Scope = merged$Scope,
    Model = merged$Model,
    LLMModelID = merged$LLMModelID,
    Comparison = paste("Evidence-k vs", merged$Selector_Comparator),
    DeltaNetCorrections = merged$NetCorrections_Evidence - merged$NetCorrections_Comparator,
    DeltaNetCorrectionEfficiency = merged$NetCorrectionEfficiency_Evidence -
      merged$NetCorrectionEfficiency_Comparator,
    DeltaHarmRate = merged$HarmRate_Evidence - merged$HarmRate_Comparator,
    EvidenceRefined = merged$NRefined_Evidence,
    ComparatorRefined = merged$NRefined_Comparator,
    stringsAsFactors = FALSE
  )
}

core_behavior <- read_optional_csv(file.path(results_dir, "ablation_refinement_behavior.csv"))
ollama_summary <- read_optional_csv(file.path(results_dir, "ollama_multimodel_all_selector_summary.csv"))

records <- list(
  standardize_behavior(
    core_behavior,
    source = "development_benchmark",
    scope = "six development datasets; cached paired ablation",
    include_models = c("deepseek-chat", "gpt-5")
  ),
  standardize_behavior(
    ollama_summary,
    source = "ollama_multimodel_pilot",
    scope = "local multimodel pilot datasets",
    exclude_models = c("deepseek-chat")
  )
)

records <- do.call(rbind, records)
summary <- aggregate_selector_behavior(records)
wide <- wide_net_table(summary)
contrasts <- contrast_table(summary)

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
utils::write.csv(summary, paste0(output_prefix, "_summary.csv"), row.names = FALSE)
utils::write.csv(wide, paste0(output_prefix, "_net_table.csv"), row.names = FALSE)
utils::write.csv(contrasts, paste0(output_prefix, "_contrasts.csv"), row.names = FALSE)

if (requireNamespace("ggplot2", quietly = TRUE) && nrow(summary) > 0) {
  plot_data <- summary[summary$Selector %in% c(
    "Evidence-k", "Random-k", "Confidence-k", "NoOntology-k", "FullRefined"
  ), , drop = FALSE]
  plot_data$Selector <- factor(
    plot_data$Selector,
    levels = c("Evidence-k", "Random-k", "Confidence-k", "NoOntology-k", "FullRefined")
  )
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = Selector, y = NetCorrectionEfficiency, fill = Selector)
  ) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey45") +
    ggplot2::geom_col(width = 0.72, colour = "grey20", linewidth = 0.2) +
    ggplot2::facet_wrap(~ Model, nrow = 1, scales = "free_y") +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(
      x = NULL,
      y = "Net correction efficiency",
      title = "Cross-model selective-refinement synthesis",
      subtitle = "Positive values indicate wrong-to-correct revisions exceed correct-to-wrong revisions"
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
      panel.grid.major.x = ggplot2::element_blank()
    )
  ggplot2::ggsave(
    paste0(output_prefix, "_net_efficiency.pdf"),
    p,
    width = 9.5,
    height = 4.8
  )
}

cat("\nCross-model selector net table:\n")
print(wide, row.names = FALSE)
cat("\nCross-model selector contrasts:\n")
print(contrasts, row.names = FALSE)
