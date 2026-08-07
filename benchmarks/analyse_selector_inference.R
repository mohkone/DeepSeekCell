# benchmarks/analyse_selector_inference.R
#
# Paired statistical inference for fixed-budget refinement selectors. This
# script consumes a refinement-behaviour table, keeps model-dataset-replicate
# blocks paired, and reports bootstrap confidence intervals, permutation tests,
# Wilcoxon tests, and simple effect-size summaries for selector contrasts.
#
# Usage:
#   Rscript benchmarks/analyse_selector_inference.R \
#     results/external_validation_refinement_behavior.csv \
#     results/external_validation_selector_inference

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

read_required_csv <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " not found: ", path, call. = FALSE)
  }
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

safe_num <- function(x) {
  out <- suppressWarnings(as.numeric(x))
  out
}

safe_divide <- function(numerator, denominator) {
  out <- numerator / denominator
  out[!is.finite(out)] <- NA_real_
  out
}

first_existing <- function(x, candidates, default = NA_character_) {
  for (column in candidates) {
    if (column %in% names(x)) {
      value <- x[[column]]
      value[is.na(value) | !nzchar(as.character(value))] <- default
      return(as.character(value))
    }
  }
  rep(default, nrow(x))
}

normalise_selector <- function(method) {
  selector <- as.character(method)
  selector <- sub("^DeepSeekCell-", "", selector)
  selector <- sub("^DeepSeek-", "Plain", selector)
  selector[selector == "SelfRefined"] <- "Evidence-k"
  selector[selector == "RiskK"] <- "Risk-k"
  selector[selector == "RandomK"] <- "Random-k"
  selector[selector == "ConfidenceK"] <- "Confidence-k"
  selector[selector == "NoOntologyK"] <- "NoOntology-k"
  selector[selector == "FullRefined"] <- "FullRefined"
  selector[selector == "Evidence"] <- "Evidence"
  selector[selector == "Calibrated"] <- "Calibrated"
  selector
}

prepare_refinement_table <- function(refinement) {
  required <- c("Dataset", "Replicate", "Method")
  missing <- setdiff(required, names(refinement))
  if (length(missing) > 0) {
    stop(
      "Refinement table is missing required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  refinement$Selector <- normalise_selector(refinement$Method)
  refinement$LLMBackend <- first_existing(refinement, c("LLMBackend", "Backend"), "unknown")
  refinement$LLMModelID <- first_existing(refinement, c("LLMModelID", "ModelID", "Model"), "unknown")
  refinement$Replicate <- safe_num(refinement$Replicate)
  refinement$BlockID <- paste(
    refinement$Dataset,
    refinement$Replicate,
    refinement$LLMBackend,
    refinement$LLMModelID,
    sep = "|"
  )

  numeric_defaults <- list(
    NRefined = 0,
    WrongToCorrect = 0,
    CorrectToWrong = 0,
    RecoveryFraction = NA_real_,
    SelectionPrecision = NA_real_,
    SelectionRecall = NA_real_,
    SelectorMCC = NA_real_,
    TotalRuntimeSec = NA_real_,
    TotalTokens = NA_real_,
    TotalCostUSD = NA_real_
  )

  for (column in names(numeric_defaults)) {
    if (!column %in% names(refinement)) {
      refinement[[column]] <- numeric_defaults[[column]]
    }
    refinement[[column]] <- safe_num(refinement[[column]])
  }

  refinement$NetCorrectionCount <- refinement$WrongToCorrect - refinement$CorrectToWrong
  gross_efficiency <- safe_divide(refinement$WrongToCorrect, refinement$NRefined)
  net_efficiency <- safe_divide(refinement$NetCorrectionCount, refinement$NRefined)

  if (!"CorrectionEfficiency" %in% names(refinement)) {
    refinement$CorrectionEfficiency <- gross_efficiency
  } else {
    refinement$CorrectionEfficiency <- safe_num(refinement$CorrectionEfficiency)
    missing_efficiency <- !is.finite(refinement$CorrectionEfficiency)
    refinement$CorrectionEfficiency[missing_efficiency] <- gross_efficiency[missing_efficiency]
  }

  # Net efficiency is the conservative primary endpoint for fixed-budget
  # refinement because it penalizes correct-to-wrong revisions.
  refinement$NetCorrectionRate <- net_efficiency
  refinement$NetCorrectionEfficiency <- net_efficiency
  refinement
}

paired_contrast_values <- function(refinement, focal, reference, metric) {
  x <- refinement[refinement$Selector %in% c(focal, reference), , drop = FALSE]
  if (nrow(x) == 0 || !metric %in% names(x)) {
    return(data.frame())
  }
  id_columns <- c("BlockID", "Dataset", "Replicate", "LLMBackend", "LLMModelID")
  x <- x[, c(id_columns, "Selector", metric), drop = FALSE]
  names(x)[names(x) == metric] <- "Value"
  x$Value <- safe_num(x$Value)

  focal_df <- x[x$Selector == focal, c(id_columns, "Value"), drop = FALSE]
  reference_df <- x[x$Selector == reference, c(id_columns, "Value"), drop = FALSE]
  names(focal_df)[names(focal_df) == "Value"] <- "FocalValue"
  names(reference_df)[names(reference_df) == "Value"] <- "ReferenceValue"

  wide <- merge(
    focal_df,
    reference_df,
    by = id_columns,
    all = FALSE,
    sort = FALSE
  )
  wide <- wide[is.finite(wide$FocalValue) & is.finite(wide$ReferenceValue), , drop = FALSE]
  if (nrow(wide) == 0) {
    return(data.frame())
  }
  wide$Difference <- wide$FocalValue - wide$ReferenceValue
  wide
}

bootstrap_mean_ci <- function(deltas, n_bootstrap, seed) {
  deltas <- deltas[is.finite(deltas)]
  if (length(deltas) == 0) {
    return(list(mean = NA_real_, lower = NA_real_, upper = NA_real_, draws = numeric()))
  }
  set.seed(seed)
  draws <- replicate(
    n_bootstrap,
    mean(sample(deltas, size = length(deltas), replace = TRUE), na.rm = TRUE)
  )
  list(
    mean = mean(deltas, na.rm = TRUE),
    lower = unname(stats::quantile(draws, 0.025, na.rm = TRUE, names = FALSE)),
    upper = unname(stats::quantile(draws, 0.975, na.rm = TRUE, names = FALSE)),
    draws = draws
  )
}

permutation_p_value <- function(deltas, n_permutations, seed) {
  deltas <- deltas[is.finite(deltas)]
  deltas <- deltas[deltas != 0]
  if (length(deltas) == 0) {
    return(NA_real_)
  }
  observed <- mean(deltas)
  set.seed(seed)
  null <- replicate(
    n_permutations,
    mean(deltas * sample(c(-1, 1), size = length(deltas), replace = TRUE))
  )
  (sum(abs(null) >= abs(observed)) + 1) / (length(null) + 1)
}

wilcoxon_p_value <- function(deltas) {
  deltas <- deltas[is.finite(deltas)]
  if (length(deltas) < 2 || !any(deltas != 0)) {
    return(c(statistic = NA_real_, p_value = NA_real_))
  }
  test <- tryCatch(
    stats::wilcox.test(deltas, mu = 0, exact = FALSE),
    error = function(e) NULL
  )
  if (is.null(test)) {
    return(c(statistic = NA_real_, p_value = NA_real_))
  }
  c(statistic = unname(test$statistic), p_value = unname(test$p.value))
}

tost_equivalence <- function(deltas, margin) {
  deltas <- deltas[is.finite(deltas)]
  if (length(deltas) < 2 || !is.finite(margin) || margin <= 0) {
    return(c(p_lower = NA_real_, p_upper = NA_real_, equivalent = NA))
  }
  mean_delta <- mean(deltas)
  se <- stats::sd(deltas) / sqrt(length(deltas))
  if (!is.finite(se) || se == 0) {
    return(c(p_lower = NA_real_, p_upper = NA_real_, equivalent = NA))
  }
  df <- length(deltas) - 1
  t_lower <- (mean_delta + margin) / se
  t_upper <- (mean_delta - margin) / se
  p_lower <- stats::pt(t_lower, df = df, lower.tail = FALSE)
  p_upper <- stats::pt(t_upper, df = df, lower.tail = TRUE)
  c(
    p_lower = p_lower,
    p_upper = p_upper,
    equivalent = as.logical(max(p_lower, p_upper, na.rm = TRUE) < 0.05)
  )
}

default_contrasts <- function(selectors) {
  contrasts <- data.frame(
    FocalSelector = character(),
    ReferenceSelector = character(),
    stringsAsFactors = FALSE
  )
  add <- function(focal, refs) {
    refs <- intersect(refs, selectors)
    if (!focal %in% selectors || length(refs) == 0) {
      return(data.frame())
    }
    data.frame(
      FocalSelector = focal,
      ReferenceSelector = refs,
      stringsAsFactors = FALSE
    )
  }
  contrasts <- rbind(
    contrasts,
    add("Risk-k", c("Evidence-k", "Random-k", "Confidence-k", "NoOntology-k")),
    add("Evidence-k", c("Random-k", "Confidence-k", "NoOntology-k"))
  )
  unique(contrasts)
}

analyse_selector_inference <- function(refinement_csv,
                                       output_prefix,
                                       n_bootstrap = as.integer(Sys.getenv("DEEPSEEKCELL_BOOTSTRAP_ITER", "5000")),
                                       n_permutations = as.integer(Sys.getenv("DEEPSEEKCELL_PERMUTATION_ITER", "5000")),
                                       seed = as.integer(Sys.getenv("DEEPSEEKCELL_INFERENCE_SEED", "100")),
                                       equivalence_margin = as.numeric(Sys.getenv("DEEPSEEKCELL_EQUIVALENCE_MARGIN", "0.05"))) {
  refinement <- prepare_refinement_table(read_required_csv(refinement_csv, "Refinement behavior table"))
  metrics <- c(
    "NetCorrectionEfficiency",
    "NetCorrectionRate",
    "CorrectionEfficiency",
    "NetCorrectionCount",
    "WrongToCorrect",
    "CorrectToWrong",
    "RecoveryFraction",
    "SelectionPrecision",
    "SelectionRecall",
    "SelectorMCC",
    "TotalRuntimeSec",
    "TotalTokens",
    "TotalCostUSD"
  )
  metrics <- intersect(metrics, names(refinement))
  contrasts <- default_contrasts(sort(unique(refinement$Selector)))

  rows <- list()
  bootstrap_rows <- list()
  counter <- 0L

  for (i in seq_len(nrow(contrasts))) {
    focal <- contrasts$FocalSelector[[i]]
    reference <- contrasts$ReferenceSelector[[i]]
    for (metric in metrics) {
      paired <- paired_contrast_values(refinement, focal, reference, metric)
      deltas <- paired$Difference
      n_blocks <- length(deltas)
      n_nonzero <- sum(deltas != 0, na.rm = TRUE)
      ci <- bootstrap_mean_ci(deltas, n_bootstrap, seed + i + match(metric, metrics))
      permutation_p <- permutation_p_value(deltas, n_permutations, seed + 1000L + i + match(metric, metrics))
      wilcox <- wilcoxon_p_value(deltas)
      tost <- tost_equivalence(deltas, equivalence_margin)
      sd_delta <- stats::sd(deltas, na.rm = TRUE)
      standardized <- if (is.finite(sd_delta) && sd_delta > 0) mean(deltas, na.rm = TRUE) / sd_delta else NA_real_
      focal_mean <- if (n_blocks > 0) mean(paired$FocalValue, na.rm = TRUE) else NA_real_
      reference_mean <- if (n_blocks > 0) mean(paired$ReferenceValue, na.rm = TRUE) else NA_real_
      median_delta <- if (n_blocks > 0) stats::median(deltas, na.rm = TRUE) else NA_real_
      wilcox_statistic <- if ("statistic" %in% names(wilcox)) wilcox[["statistic"]] else NA_real_
      wilcox_p <- if ("p_value" %in% names(wilcox)) wilcox[["p_value"]] else NA_real_
      tost_lower <- if ("p_lower" %in% names(tost)) as.numeric(tost[["p_lower"]]) else NA_real_
      tost_upper <- if ("p_upper" %in% names(tost)) as.numeric(tost[["p_upper"]]) else NA_real_
      tost_equivalent <- if ("equivalent" %in% names(tost)) as.logical(tost[["equivalent"]]) else NA

      rows[[length(rows) + 1]] <- data.frame(
        FocalSelector = focal,
        ReferenceSelector = reference,
        Metric = metric,
        Blocks = n_blocks,
        NonzeroDifferences = n_nonzero,
        FocalMean = focal_mean,
        ReferenceMean = reference_mean,
        MeanDifference = ci$mean,
        BootstrapCI95Lower = ci$lower,
        BootstrapCI95Upper = ci$upper,
        MedianDifference = median_delta,
        StandardizedMeanDifference = standardized,
        ImprovementBlocks = sum(deltas > 0, na.rm = TRUE),
        HarmBlocks = sum(deltas < 0, na.rm = TRUE),
        TieBlocks = sum(deltas == 0, na.rm = TRUE),
        PermutationPValue = permutation_p,
        WilcoxonStatistic = wilcox_statistic,
        WilcoxonPValue = wilcox_p,
        EquivalenceMargin = equivalence_margin,
        TOSTLowerPValue = tost_lower,
        TOSTUpperPValue = tost_upper,
        EquivalentWithinMargin = tost_equivalent,
        stringsAsFactors = FALSE
      )

      if (length(ci$draws) > 0) {
        counter <- counter + 1L
        bootstrap_rows[[length(bootstrap_rows) + 1]] <- data.frame(
          ContrastID = counter,
          FocalSelector = focal,
          ReferenceSelector = reference,
          Metric = metric,
          BootstrapMeanDifference = ci$draws,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  inference <- do.call(rbind, rows)
  if (nrow(inference) > 0) {
    inference$PermutationPValueBH <- ave(
      inference$PermutationPValue,
      inference$Metric,
      FUN = function(p) stats::p.adjust(p, method = "BH")
    )
    inference$WilcoxonPValueBH <- ave(
      inference$WilcoxonPValue,
      inference$Metric,
      FUN = function(p) stats::p.adjust(p, method = "BH")
    )
  }

  dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(
    inference,
    paste0(output_prefix, "_paired_selector_inference.csv"),
    row.names = FALSE
  )
  if (length(bootstrap_rows) > 0) {
    utils::write.csv(
      do.call(rbind, bootstrap_rows),
      paste0(output_prefix, "_bootstrap_distributions.csv"),
      row.names = FALSE
    )
  }

  if (requireNamespace("ggplot2", quietly = TRUE) && nrow(inference) > 0) {
    plot_df <- inference[inference$Metric %in% c("CorrectionEfficiency", "NetCorrectionRate", "NetCorrectionCount"), , drop = FALSE]
    if (nrow(plot_df) > 0) {
      plot_df$Contrast <- paste(plot_df$FocalSelector, "vs", plot_df$ReferenceSelector)
      p <- ggplot2::ggplot(
        plot_df,
        ggplot2::aes(
          x = .data$MeanDifference,
          y = stats::reorder(.data$Contrast, .data$MeanDifference)
        )
      ) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
        ggplot2::geom_errorbar(
          ggplot2::aes(
            xmin = .data$BootstrapCI95Lower,
            xmax = .data$BootstrapCI95Upper
          ),
          orientation = "y",
          width = 0.2
        ) +
        ggplot2::geom_point(size = 2) +
        ggplot2::facet_wrap(~Metric, scales = "free_x") +
        ggplot2::labs(x = "Paired mean difference", y = NULL) +
        ggplot2::theme_minimal(base_size = 11)
      ggplot2::ggsave(
        paste0(output_prefix, "_paired_selector_inference.pdf"),
        p,
        width = 9,
        height = 5
      )
    }
  }

  message("Wrote paired selector inference to ", paste0(output_prefix, "_paired_selector_inference.csv"))
  invisible(inference)
}

args <- commandArgs(trailingOnly = TRUE)
refinement_csv <- if (length(args) >= 1) args[[1]] else file.path("results", "external_validation_refinement_behavior.csv")
output_prefix <- if (length(args) >= 2) args[[2]] else file.path("results", "external_validation_selector_inference")

analyse_selector_inference(refinement_csv, output_prefix)
