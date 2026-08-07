# benchmarks/analyse_external_validation_meta_analysis.R
#
# Study-level meta-analysis for locked external validation. This script consumes
# completed external-validation CSVs and estimates pooled selector efficiency,
# calibration improvement, cost reduction, and heterogeneity. It does not call
# an LLM, tune thresholds, train models, or alter the frozen reliability method.
#
# Usage:
#   Rscript benchmarks/analyse_external_validation_meta_analysis.R \
#     results/external_validation_refinement_behavior.csv \
#     results/external_validation_confidence_quality.csv \
#     results/external_validation_meta_analysis

as_number <- function(x) {
  suppressWarnings(as.numeric(x))
}

as_text <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  trimws(x)
}

read_csv_required <- function(path, description) {
  if (!file.exists(path)) {
    stop(description, " not found: ", path, call. = FALSE)
  }
  utils::read.csv(path, stringsAsFactors = FALSE)
}

ensure_column <- function(x, column, default = NA) {
  if (!column %in% names(x)) {
    x[[column]] <- default
  }
  x
}

block_id <- function(x) {
  for (column in c("Dataset", "Replicate", "LLMBackend", "LLMModelID")) {
    x <- ensure_column(x, column, "")
  }
  paste(x$Dataset, x$Replicate, x$LLMBackend, x$LLMModelID, sep = "|")
}

selector_label <- function(x) {
  method <- if ("Method" %in% names(x)) as_text(x$Method) else rep("", nrow(x))
  selector <- if ("RefinementSelector" %in% names(x)) as_text(x$RefinementSelector) else rep("", nrow(x))
  out <- ifelse(nzchar(selector), selector, method)
  out <- gsub("^DeepSeekCell-", "", out)
  out <- gsub("^DeepSeek-", "", out)
  out <- gsub("_", "-", out)
  out
}

method_family <- function(x) {
  x <- selector_label(x)
  out <- x
  out[x %in% c("evidence-k", "SelfRefined", "self-refined")] <- "Evidence-k"
  out[x %in% c("risk-k", "RiskK")] <- "Risk-k"
  out[x %in% c("confidence-k", "ConfidenceK")] <- "Confidence-k"
  out[x %in% c("random-k", "RandomK")] <- "Random-k"
  out[x %in% c("evidence-no-ontology-k", "NoOntologyK")] <- "NoOntology-k"
  out[x %in% c("full", "FullRefined")] <- "FullRefined"
  out[x %in% c("none", "Plain", "Evidence", "Calibrated")] <- x[x %in% c("none", "Plain", "Evidence", "Calibrated")]
  out
}

safe_divide <- function(num, den) {
  out <- rep(NA_real_, length(num))
  ok <- !is.na(num) & !is.na(den) & den != 0
  out[ok] <- num[ok] / den[ok]
  out
}

bounded_p_se <- function(events, n) {
  events <- as_number(events)
  n <- as_number(n)
  p <- safe_divide(events, n)
  p_clip <- pmin(pmax(p, 0), 1)
  se <- sqrt(p_clip * (1 - p_clip) / pmax(n, 1))
  se[is.na(se) | se == 0] <- sqrt(0.5 * 0.5 / pmax(n[is.na(se) | se == 0], 1))
  se
}

normal_ci <- function(estimate, se, level = 0.95) {
  z <- stats::qnorm(1 - (1 - level) / 2)
  c(estimate - z * se, estimate + z * se)
}

dl_meta <- function(yi, sei) {
  yi <- as_number(yi)
  sei <- as_number(sei)
  ok <- is.finite(yi) & is.finite(sei) & sei > 0
  yi <- yi[ok]
  sei <- sei[ok]
  k <- length(yi)
  if (k == 0) {
    return(list(
      k = 0, fixed = NA_real_, fixed_se = NA_real_,
      random = NA_real_, random_se = NA_real_,
      q = NA_real_, q_df = NA_integer_, q_p = NA_real_,
      tau2 = NA_real_, i2 = NA_real_
    ))
  }

  w <- 1 / sei^2
  fixed <- sum(w * yi) / sum(w)
  fixed_se <- sqrt(1 / sum(w))
  q <- sum(w * (yi - fixed)^2)
  q_df <- max(k - 1, 0)
  c_val <- sum(w) - sum(w^2) / sum(w)
  tau2 <- if (k > 1 && c_val > 0) max(0, (q - q_df) / c_val) else 0
  wr <- 1 / (sei^2 + tau2)
  random <- sum(wr * yi) / sum(wr)
  random_se <- sqrt(1 / sum(wr))
  i2 <- if (k > 1 && q > 0) max(0, (q - q_df) / q) * 100 else 0

  list(
    k = k,
    fixed = fixed,
    fixed_se = fixed_se,
    random = random,
    random_se = random_se,
    q = q,
    q_df = q_df,
    q_p = if (q_df > 0) stats::pchisq(q, df = q_df, lower.tail = FALSE) else NA_real_,
    tau2 = tau2,
    i2 = i2
  )
}

summarise_meta <- function(effects, group_columns) {
  if (nrow(effects) == 0) {
    return(data.frame())
  }

  split_key <- interaction(effects[group_columns], drop = TRUE, lex.order = TRUE)
  rows <- lapply(split(effects, split_key), function(x) {
    meta <- dl_meta(x$Effect, x$StdError)
    ci_fixed <- normal_ci(meta$fixed, meta$fixed_se)
    ci_random <- normal_ci(meta$random, meta$random_se)
    out <- as.data.frame(x[1, group_columns, drop = FALSE], stringsAsFactors = FALSE)
    out$NStudies <- length(unique(x$Dataset))
    out$NBlocks <- length(unique(x$BlockID))
    out$FixedEffectEstimate <- meta$fixed
    out$FixedEffectSE <- meta$fixed_se
    out$FixedEffectCI95Lower <- ci_fixed[[1]]
    out$FixedEffectCI95Upper <- ci_fixed[[2]]
    out$RandomEffectsEstimate <- meta$random
    out$RandomEffectsSE <- meta$random_se
    out$RandomEffectsCI95Lower <- ci_random[[1]]
    out$RandomEffectsCI95Upper <- ci_random[[2]]
    out$CochranQ <- meta$q
    out$Qdf <- meta$q_df
    out$QpValue <- meta$q_p
    out$TauSquared <- meta$tau2
    out$I2Percent <- meta$i2
    out
  })
  do.call(rbind, rows)
}

prepare_refinement_effects <- function(refinement) {
  refinement$BlockID <- block_id(refinement)
  refinement$Selector <- method_family(refinement)
  for (column in c(
    "NRefined", "WrongToCorrect", "CorrectToWrong", "NClusters",
    "RefinementCostUSD", "RefinementTokens", "RefinementRuntimeSec",
    "SecondPassCalls"
  )) {
    refinement <- ensure_column(refinement, column, NA_real_)
    refinement[[column]] <- as_number(refinement[[column]])
  }

  refinement$NetCorrections <- refinement$WrongToCorrect - refinement$CorrectToWrong
  refinement$WrongToCorrectRate <- safe_divide(refinement$WrongToCorrect, refinement$NRefined)
  refinement$CorrectToWrongRate <- safe_divide(refinement$CorrectToWrong, refinement$NRefined)
  refinement$NetCorrectionEfficiency <- safe_divide(refinement$NetCorrections, refinement$NRefined)
  refinement$RefinementBudgetFraction <- safe_divide(refinement$NRefined, refinement$NClusters)
  refinement$WrongToCorrectRateSE <- bounded_p_se(refinement$WrongToCorrect, refinement$NRefined)
  refinement$CorrectToWrongRateSE <- bounded_p_se(refinement$CorrectToWrong, refinement$NRefined)

  pw <- pmin(pmax(refinement$WrongToCorrectRate, 0), 1)
  pc <- pmin(pmax(refinement$CorrectToWrongRate, 0), 1)
  n <- pmax(refinement$NRefined, 1)
  refinement$NetCorrectionEfficiencySE <- sqrt((pw * (1 - pw) + pc * (1 - pc)) / n)
  refinement$NetCorrectionEfficiencySE[
    !is.finite(refinement$NetCorrectionEfficiencySE) |
      refinement$NetCorrectionEfficiencySE == 0
  ] <- sqrt(0.5 / pmax(n[
    !is.finite(refinement$NetCorrectionEfficiencySE) |
      refinement$NetCorrectionEfficiencySE == 0
  ], 1))

  outcomes <- list(
    WrongToCorrectRate = c("WrongToCorrectRate", "WrongToCorrectRateSE"),
    CorrectToWrongRate = c("CorrectToWrongRate", "CorrectToWrongRateSE"),
    NetCorrectionEfficiency = c("NetCorrectionEfficiency", "NetCorrectionEfficiencySE"),
    RefinementBudgetFraction = c("RefinementBudgetFraction", NA_character_)
  )

  rows <- list()
  for (outcome in names(outcomes)) {
    effect_col <- outcomes[[outcome]][[1]]
    se_col <- outcomes[[outcome]][[2]]
    x <- refinement[is.finite(refinement[[effect_col]]) & refinement$NRefined > 0, , drop = FALSE]
    if (nrow(x) == 0) next
    se <- if (!is.na(se_col) && se_col %in% names(x)) x[[se_col]] else 1 / sqrt(pmax(x$NClusters, 1))
    rows[[length(rows) + 1]] <- data.frame(
      Dataset = x$Dataset,
      Tissue = x$Tissue,
      Species = x$Species,
      BlockID = x$BlockID,
      Selector = x$Selector,
      Outcome = outcome,
      Effect = x[[effect_col]],
      StdError = se,
      NRefined = x$NRefined,
      WrongToCorrect = x$WrongToCorrect,
      CorrectToWrong = x$CorrectToWrong,
      NClusters = x$NClusters,
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0) data.frame() else do.call(rbind, rows)
}

prepare_selector_contrasts <- function(refinement,
                                       reference_selectors = c("Random-k", "Confidence-k", "NoOntology-k"),
                                       focal_selectors = c("Evidence-k", "Risk-k")) {
  refinement$BlockID <- block_id(refinement)
  refinement$Selector <- method_family(refinement)
  for (column in c("NRefined", "WrongToCorrect", "CorrectToWrong")) {
    refinement <- ensure_column(refinement, column, NA_real_)
    refinement[[column]] <- as_number(refinement[[column]])
  }
  refinement$NetCorrectionEfficiency <- safe_divide(
    refinement$WrongToCorrect - refinement$CorrectToWrong,
    refinement$NRefined
  )
  refinement$WrongToCorrectRate <- safe_divide(refinement$WrongToCorrect, refinement$NRefined)
  refinement$CorrectToWrongRate <- safe_divide(refinement$CorrectToWrong, refinement$NRefined)
  refinement$NetSE <- sqrt(
    (
      pmin(pmax(refinement$WrongToCorrectRate, 0), 1) *
        (1 - pmin(pmax(refinement$WrongToCorrectRate, 0), 1)) +
        pmin(pmax(refinement$CorrectToWrongRate, 0), 1) *
        (1 - pmin(pmax(refinement$CorrectToWrongRate, 0), 1))
    ) / pmax(refinement$NRefined, 1)
  )
  refinement$NetSE[!is.finite(refinement$NetSE) | refinement$NetSE == 0] <- 1 / sqrt(pmax(refinement$NRefined[!is.finite(refinement$NetSE) | refinement$NetSE == 0], 1))

  rows <- list()
  for (block in unique(refinement$BlockID)) {
    x <- refinement[refinement$BlockID == block, , drop = FALSE]
    for (focal in focal_selectors) {
      f <- x[x$Selector == focal, , drop = FALSE]
      if (nrow(f) == 0) next
      for (reference in reference_selectors) {
        r <- x[x$Selector == reference, , drop = FALSE]
        if (nrow(r) == 0) next
        rows[[length(rows) + 1]] <- data.frame(
          Dataset = f$Dataset[[1]],
          Tissue = f$Tissue[[1]],
          Species = f$Species[[1]],
          BlockID = block,
          FocalSelector = focal,
          ReferenceSelector = reference,
          Outcome = "NetCorrectionEfficiencyDifference",
          Effect = f$NetCorrectionEfficiency[[1]] - r$NetCorrectionEfficiency[[1]],
          StdError = sqrt(f$NetSE[[1]]^2 + r$NetSE[[1]]^2),
          FocalNRefined = f$NRefined[[1]],
          ReferenceNRefined = r$NRefined[[1]],
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(rows) == 0) data.frame() else do.call(rbind, rows)
}

prepare_cost_effects <- function(refinement) {
  refinement$BlockID <- block_id(refinement)
  refinement$Selector <- method_family(refinement)
  for (column in c("NRefined", "RefinementCostUSD", "RefinementTokens", "RefinementRuntimeSec", "SecondPassCalls")) {
    refinement <- ensure_column(refinement, column, NA_real_)
    refinement[[column]] <- as_number(refinement[[column]])
  }

  rows <- list()
  for (block in unique(refinement$BlockID)) {
    x <- refinement[refinement$BlockID == block, , drop = FALSE]
    full <- x[x$Selector == "FullRefined", , drop = FALSE]
    if (nrow(full) == 0) next
    candidates <- x[x$Selector != "FullRefined", , drop = FALSE]
    for (i in seq_len(nrow(candidates))) {
      cnd <- candidates[i, , drop = FALSE]
      metrics <- data.frame(
        Outcome = c("RefinementCallReduction", "CostReduction", "TokenReduction", "RuntimeReduction"),
        Effect = c(
          1 - safe_divide(cnd$NRefined, full$NRefined),
          1 - safe_divide(cnd$RefinementCostUSD, full$RefinementCostUSD),
          1 - safe_divide(cnd$RefinementTokens, full$RefinementTokens),
          1 - safe_divide(cnd$RefinementRuntimeSec, full$RefinementRuntimeSec)
        ),
        stringsAsFactors = FALSE
      )
      metrics <- metrics[is.finite(metrics$Effect), , drop = FALSE]
      if (nrow(metrics) == 0) next
      rows[[length(rows) + 1]] <- data.frame(
        Dataset = cnd$Dataset,
        Tissue = cnd$Tissue,
        Species = cnd$Species,
        BlockID = block,
        Selector = cnd$Selector,
        metrics,
        StdError = 1 / sqrt(pmax(full$NRefined, 1)),
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) == 0) data.frame() else do.call(rbind, rows)
}

prepare_calibration_effects <- function(confidence) {
  confidence$BlockID <- block_id(confidence)
  confidence$Selector <- method_family(confidence)
  for (column in c("Brier", "ECE", "BinaryCorrectnessNLL", "NLL", "AURC", "ClusterCount")) {
    confidence <- ensure_column(confidence, column, NA_real_)
    confidence[[column]] <- as_number(confidence[[column]])
  }
  if (!"BinaryCorrectnessNLL" %in% names(confidence) || all(is.na(confidence$BinaryCorrectnessNLL))) {
    confidence$BinaryCorrectnessNLL <- confidence$NLL
  }

  comparison_methods <- c("Calibrated" = "DeepSeekCell-Calibrated", "SelfRefined" = "DeepSeekCell-SelfRefined")
  rows <- list()
  for (block in unique(confidence$BlockID)) {
    x <- confidence[confidence$BlockID == block, , drop = FALSE]
    plain <- x[x$Method == "DeepSeek-Plain", , drop = FALSE]
    if (nrow(plain) == 0) next
    for (label in names(comparison_methods)) {
      cmp <- x[x$Method == comparison_methods[[label]], , drop = FALSE]
      if (nrow(cmp) == 0) next
      metrics <- c("Brier", "ECE", "BinaryCorrectnessNLL", "AURC")
      for (metric in metrics) {
        if (!is.finite(plain[[metric]][[1]]) || !is.finite(cmp[[metric]][[1]])) next
        # Positive values indicate improvement over raw first-pass confidence.
        effect <- plain[[metric]][[1]] - cmp[[metric]][[1]]
        rows[[length(rows) + 1]] <- data.frame(
          Dataset = cmp$Dataset[[1]],
          Tissue = cmp$Tissue[[1]],
          Species = cmp$Species[[1]],
          BlockID = block,
          Selector = label,
          Outcome = paste0(metric, "Improvement"),
          Effect = effect,
          StdError = 1 / sqrt(pmax(cmp$ClusterCount[[1]], 1)),
          PlainValue = plain[[metric]][[1]],
          ComparedValue = cmp[[metric]][[1]],
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(rows) == 0) data.frame() else do.call(rbind, rows)
}

plot_meta_forest <- function(prefix, pooled) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || nrow(pooled) == 0) {
    return(invisible(FALSE))
  }
  focus <- pooled[
    pooled$Outcome %in% c("NetCorrectionEfficiency", "NetCorrectionEfficiencyDifference") &
      pooled$Selector %in% c("Evidence-k", "Risk-k", "Random-k", "Confidence-k", "NoOntology-k"),
    ,
    drop = FALSE
  ]
  if (nrow(focus) == 0) {
    return(invisible(FALSE))
  }
  focus$Label <- if ("ReferenceSelector" %in% names(focus)) {
    ifelse(
      is.na(focus$ReferenceSelector),
      focus$Selector,
      paste(focus$Selector, "vs", focus$ReferenceSelector)
    )
  } else {
    focus$Selector
  }
  p <- ggplot2::ggplot(
    focus,
    ggplot2::aes(
      x = .data$RandomEffectsEstimate,
      y = stats::reorder(.data$Label, .data$RandomEffectsEstimate)
    )
  ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = .data$RandomEffectsCI95Lower,
        xend = .data$RandomEffectsCI95Upper,
        yend = stats::reorder(.data$Label, .data$RandomEffectsEstimate)
      )
    ) +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_wrap(~Outcome, scales = "free_x") +
    ggplot2::labs(x = "Random-effects estimate", y = NULL) +
    ggplot2::theme_minimal(base_size = 11)
  ggplot2::ggsave(paste0(prefix, "_selector_forest.pdf"), p, width = 8, height = 5)
  invisible(TRUE)
}

run_external_meta_analysis <- function(refinement_csv,
                                       confidence_csv,
                                       output_prefix) {
  refinement <- read_csv_required(refinement_csv, "External refinement behavior")
  confidence <- read_csv_required(confidence_csv, "External confidence-quality table")
  dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)

  selector_effects <- prepare_refinement_effects(refinement)
  selector_pooled <- summarise_meta(selector_effects, c("Selector", "Outcome"))

  selector_contrasts <- prepare_selector_contrasts(refinement)
  selector_contrast_pooled <- summarise_meta(
    selector_contrasts,
    c("FocalSelector", "ReferenceSelector", "Outcome")
  )

  cost_effects <- prepare_cost_effects(refinement)
  cost_pooled <- summarise_meta(cost_effects, c("Selector", "Outcome"))

  calibration_effects <- prepare_calibration_effects(confidence)
  calibration_pooled <- summarise_meta(calibration_effects, c("Selector", "Outcome"))

  utils::write.csv(selector_effects, paste0(output_prefix, "_selector_study_effects.csv"), row.names = FALSE)
  utils::write.csv(selector_pooled, paste0(output_prefix, "_selector_pooled_effects.csv"), row.names = FALSE)
  utils::write.csv(selector_contrasts, paste0(output_prefix, "_selector_contrasts.csv"), row.names = FALSE)
  utils::write.csv(selector_contrast_pooled, paste0(output_prefix, "_selector_contrast_pooled_effects.csv"), row.names = FALSE)
  utils::write.csv(cost_effects, paste0(output_prefix, "_cost_study_effects.csv"), row.names = FALSE)
  utils::write.csv(cost_pooled, paste0(output_prefix, "_cost_pooled_effects.csv"), row.names = FALSE)
  utils::write.csv(calibration_effects, paste0(output_prefix, "_calibration_study_effects.csv"), row.names = FALSE)
  utils::write.csv(calibration_pooled, paste0(output_prefix, "_calibration_pooled_effects.csv"), row.names = FALSE)

  forest_input <- selector_pooled
  if (nrow(selector_contrast_pooled) > 0) {
    contrast_plot <- selector_contrast_pooled
    contrast_plot$Selector <- contrast_plot$FocalSelector
    forest_input <- rbind(
      forest_input[intersect(names(forest_input), names(contrast_plot))],
      contrast_plot[intersect(names(forest_input), names(contrast_plot))]
    )
  }
  plot_meta_forest(output_prefix, forest_input)

  invisible(list(
    selector_effects = selector_effects,
    selector_pooled = selector_pooled,
    selector_contrasts = selector_contrasts,
    selector_contrast_pooled = selector_contrast_pooled,
    cost_effects = cost_effects,
    cost_pooled = cost_pooled,
    calibration_effects = calibration_effects,
    calibration_pooled = calibration_pooled
  ))
}

args <- commandArgs(trailingOnly = TRUE)
refinement_csv <- if (length(args) >= 1) args[[1]] else file.path("results", "external_validation_refinement_behavior.csv")
confidence_csv <- if (length(args) >= 2) args[[2]] else file.path("results", "external_validation_confidence_quality.csv")
output_prefix <- if (length(args) >= 3) args[[3]] else file.path("results", "external_validation_meta_analysis")

run_external_meta_analysis(refinement_csv, confidence_csv, output_prefix)
