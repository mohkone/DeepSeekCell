# benchmarks/ablation.R
#
# Utilities for paired DeepSeekCell ablation experiments. These functions are
# deliberately API-free: the first-pass LLM response should be generated once,
# saved with a hash, parsed into annotations, and then reused by every arm.

prepare_first_pass_annotations <- function(response_text,
                                           markers,
                                           tissue,
                                           ontology_data = NULL) {
  annotations <- parse_ablation_annotation_response(response_text)

  if (nrow(annotations) == 0) {
    annotations <- data.frame(
      Cluster = names(markers),
      CellType = "Unknown",
      Confidence = 0.5,
      CandidateCellTypes = "",
      IsMixed = FALSE,
      PrimaryCellType = "Unknown",
      SecondaryCellType = "",
      TissueConsistency = "unknown",
      Reasoning = "No parseable annotation returned by model.",
      stringsAsFactors = FALSE
    )
  }

  annotations <- annotations[!duplicated(annotations$Cluster), , drop = FALSE]
  annotation_match <- match_ablation_clusters(names(markers), annotations$Cluster)
  missing_clusters <- names(markers)[is.na(annotation_match)]
  if (length(missing_clusters) > 0) {
    annotations <- rbind(
      annotations,
      data.frame(
        Cluster = missing_clusters,
        CellType = "Unknown",
        Confidence = 0.5,
        CandidateCellTypes = "",
        IsMixed = FALSE,
        PrimaryCellType = "Unknown",
        SecondaryCellType = "",
        TissueConsistency = "unknown",
        Reasoning = "No annotation returned by model.",
        stringsAsFactors = FALSE
      )
    )
  }

  annotation_match <- match_ablation_clusters(names(markers), annotations$Cluster)
  annotations <- annotations[annotation_match, , drop = FALSE]
  annotations$Cluster <- names(markers)
  rownames(annotations) <- NULL
  add_ablation_ontology(annotations, ontology_data, tissue)
}

parse_ablation_annotation_response <- function(response_text) {
  if (is.null(response_text) || !nzchar(trimws(response_text))) {
    return(empty_ablation_annotations())
  }

  annotation_key_count <- length(
    gregexpr("\"annotations\"[[:space:]]*:", response_text, perl = TRUE)[[1]]
  )
  if (annotation_key_count > 1) {
    ann <- parse_ablation_annotation_fragments(response_text)
  }

  json_str <- extract_ablation_json(response_text)
  if (exists("ann", inherits = FALSE)) {
    parsed <- NULL
  } else if (is.null(json_str)) {
    ann <- parse_ablation_annotation_fragments(response_text)
    if (!is.data.frame(ann) || nrow(ann) == 0) {
      return(empty_ablation_annotations())
    }
    parsed <- NULL
  } else {
    parsed <- tryCatch(
      jsonlite::fromJSON(json_str, simplifyDataFrame = TRUE),
      error = function(e) NULL
    )
  }

  if (!exists("ann", inherits = FALSE)) {
    if (is.null(parsed)) {
      ann <- parse_ablation_annotation_fragments(response_text)
      if (!is.data.frame(ann) || nrow(ann) == 0) {
        return(empty_ablation_annotations())
      }
    } else {
      ann <- parsed$annotations %||% parsed
    }
  }

  if (is.list(ann) && !is.data.frame(ann)) {
    ann <- tryCatch(as.data.frame(ann, stringsAsFactors = FALSE), error = function(e) NULL)
  }
  if (!is.data.frame(ann) || nrow(ann) == 0) {
    return(empty_ablation_annotations())
  }

  names(ann) <- tolower(names(ann))
  get_col <- function(candidates, default) {
    hit <- candidates[candidates %in% names(ann)]
    if (length(hit) == 0) return(rep(default, nrow(ann)))
    ann[[hit[1]]]
  }

  candidates <- collapse_ablation_candidates(
    get_col(c("candidate_cell_types", "candidatecelltypes", "candidates"), ""),
    nrow(ann)
  )

  out <- data.frame(
    Cluster = as.character(get_col("cluster", paste0("Cluster", seq_len(nrow(ann))))),
    CellType = as.character(get_col(c("cell_type", "celltype"), "Unknown")),
    Confidence = as_confidence(get_col("confidence", 0.5)),
    CandidateCellTypes = candidates,
    IsMixed = as_flag(get_col(c("is_mixed", "ismixed", "mixed"), FALSE)),
    PrimaryCellType = as.character(get_col(c("primary_cell_type", "primarycelltype"), "")),
    SecondaryCellType = as.character(get_col(c("secondary_cell_type", "secondarycelltype"), "")),
    TissueConsistency = as.character(get_col(c("tissue_consistency", "tissueconsistency"), "unknown")),
    Reasoning = as.character(get_col("reasoning", "")),
    stringsAsFactors = FALSE
  )

  out$Cluster <- trimws(out$Cluster)
  out$CellType <- trimws(out$CellType)
  out$CellType[is.na(out$CellType) | !nzchar(out$CellType)] <- "Unknown"
  out$TissueConsistency <- tolower(trimws(out$TissueConsistency))
  out$TissueConsistency[!out$TissueConsistency %in% c(
    "expected", "unexpected", "possible_contamination", "possible_doublet", "unknown"
  )] <- "unknown"
  out <- out[!duplicated(normalise_ablation_cluster_id(out$Cluster)), , drop = FALSE]
  rownames(out) <- NULL
  out
}

parse_ablation_annotation_fragments <- function(response_text) {
  objects <- extract_ablation_json_objects(response_text)
  if (length(objects) == 0) {
    return(NULL)
  }

  rows <- lapply(objects, function(obj) {
    parsed <- tryCatch(
      jsonlite::fromJSON(obj, simplifyDataFrame = TRUE),
      error = function(e) NULL
    )
    if (is.null(parsed)) {
      return(NULL)
    }

    ann <- parsed$annotations %||% parsed
    if (is.list(ann) && !is.data.frame(ann)) {
      ann <- tryCatch(as.data.frame(ann, stringsAsFactors = FALSE), error = function(e) NULL)
    }
    if (!is.data.frame(ann) || nrow(ann) == 0) {
      return(NULL)
    }

    names(ann) <- tolower(names(ann))
    if (!any(c("cell_type", "celltype") %in% names(ann))) {
      return(NULL)
    }

    ann
  })

  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0) {
    return(NULL)
  }

  all_names <- unique(unlist(lapply(rows, names)))
  rows <- lapply(rows, function(x) {
    missing <- setdiff(all_names, names(x))
    for (column in missing) {
      x[[column]] <- NA
    }
    x[all_names]
  })

  do.call(rbind, rows)
}

empty_ablation_annotations <- function() {
  data.frame(
    Cluster = character(),
    CellType = character(),
    Confidence = numeric(),
    CandidateCellTypes = character(),
    IsMixed = logical(),
    PrimaryCellType = character(),
    SecondaryCellType = character(),
    TissueConsistency = character(),
    Reasoning = character(),
    stringsAsFactors = FALSE
  )
}

collapse_ablation_candidates <- function(candidates, n_rows) {
  if (is.null(candidates)) {
    return(rep("", n_rows))
  }

  if (is.data.frame(candidates)) {
    candidates <- as.list(candidates)
  }

  if (is.list(candidates)) {
    out <- vapply(candidates, function(x) {
      x <- unique(trimws(as.character(unlist(x))))
      x <- x[nzchar(x) & !is.na(x)]
      paste(x, collapse = "; ")
    }, character(1))
    if (length(out) == n_rows) return(out)
    return(rep(paste(out, collapse = "; "), n_rows))
  }

  out <- trimws(as.character(candidates))
  out[is.na(out)] <- ""
  if (length(out) == n_rows) return(out)
  rep(out[1] %||% "", n_rows)
}

extract_ablation_json <- function(text) {
  cleaned <- trimws(as.character(text))
  cleaned <- gsub("^```(?:json)?\\s*", "", cleaned, perl = TRUE)
  cleaned <- gsub("\\s*```$", "", cleaned, perl = TRUE)

  start <- regexpr("\\{", cleaned, perl = TRUE)[1]
  if (start == -1) return(NULL)

  chars <- strsplit(cleaned, "", fixed = TRUE)[[1]]
  depth <- 0
  in_string <- FALSE
  escaped <- FALSE

  for (i in seq.int(start, length(chars))) {
    ch <- chars[[i]]
    if (in_string) {
      if (escaped) {
        escaped <- FALSE
      } else if (identical(ch, "\\")) {
        escaped <- TRUE
      } else if (identical(ch, "\"")) {
        in_string <- FALSE
      }
      next
    }
    if (identical(ch, "\"")) {
      in_string <- TRUE
      next
    }
    if (identical(ch, "{")) {
      depth <- depth + 1
    } else if (identical(ch, "}")) {
      depth <- depth - 1
      if (depth == 0) return(paste(chars[start:i], collapse = ""))
    }
  }

  NULL
}

extract_ablation_json_objects <- function(text) {
  cleaned <- trimws(as.character(text))
  cleaned <- gsub("^```(?:json)?\\s*", "", cleaned, perl = TRUE)
  cleaned <- gsub("\\s*```$", "", cleaned, perl = TRUE)

  chars <- strsplit(cleaned, "", fixed = TRUE)[[1]]
  starts <- which(chars == "{")
  if (length(starts) == 0) {
    return(character())
  }

  objects <- character()
  for (start in starts) {
    depth <- 0
    in_string <- FALSE
    escaped <- FALSE

    for (i in seq.int(start, length(chars))) {
      ch <- chars[[i]]
      if (in_string) {
        if (escaped) {
          escaped <- FALSE
        } else if (identical(ch, "\\")) {
          escaped <- TRUE
        } else if (identical(ch, "\"")) {
          in_string <- FALSE
        }
        next
      }
      if (identical(ch, "\"")) {
        in_string <- TRUE
        next
      }
      if (identical(ch, "{")) {
        depth <- depth + 1
      } else if (identical(ch, "}")) {
        depth <- depth - 1
        if (depth == 0) {
          objects <- c(objects, paste(chars[start:i], collapse = ""))
          break
        }
      }
    }
  }

  unique(objects)
}

create_deepseekcell_ablation_arms <- function(first_pass_annotations,
                                              markers,
                                              tissue,
                                              ontology_data = NULL,
                                              refined_annotations = NULL,
                                              include_self_refined = !is.null(refined_annotations)) {
  first_pass_annotations <- add_ablation_ontology(
    first_pass_annotations,
    ontology_data,
    tissue
  )
  raw_confidence <- as_confidence(first_pass_annotations$Confidence)

  plain <- first_pass_annotations
  plain$LLMConfidence <- raw_confidence
  plain$Confidence <- raw_confidence
  plain$AblationArm <- "DeepSeek-Plain"
  plain$ConfidenceMethod <- "llm_raw"

  evidence_scores <- score_annotation_evidence(first_pass_annotations, markers, tissue)
  evidence <- first_pass_annotations
  for (column in setdiff(names(evidence_scores), c("Cluster", "CalibratedConfidence"))) {
    evidence[[column]] <- evidence_scores[[column]]
  }
  evidence$LLMConfidence <- raw_confidence
  evidence$Confidence <- raw_confidence
  evidence$AblationArm <- "DeepSeekCell-Evidence"
  evidence$ConfidenceMethod <- "llm_raw_with_evidence"

  calibrated <- calibrate_annotation_confidence(first_pass_annotations, markers, tissue)
  calibrated$AblationArm <- "DeepSeekCell-Calibrated"

  arms <- list(
    `DeepSeek-Plain` = plain,
    `DeepSeekCell-Evidence` = evidence,
    `DeepSeekCell-Calibrated` = calibrated
  )

  if (isTRUE(include_self_refined)) {
    if (is.null(refined_annotations)) {
      refined_annotations <- first_pass_annotations[0, , drop = FALSE]
    }
    requires_refinement <- as_flag(calibrated$RequiresRefinement)
    requires_refinement[is.na(requires_refinement)] <- FALSE
    selected_clusters <- calibrated$Cluster[requires_refinement]
    self_refined <- create_refined_ablation_arm(
      first_pass_annotations = first_pass_annotations,
      markers = markers,
      tissue = tissue,
      ontology_data = ontology_data,
      refined_annotations = refined_annotations,
      selected_clusters = selected_clusters,
      arm_name = "DeepSeekCell-SelfRefined"
    )
    arms$`DeepSeekCell-SelfRefined` <- self_refined
  }

  arms
}

create_refined_ablation_arm <- function(first_pass_annotations,
                                        markers,
                                        tissue,
                                        ontology_data = NULL,
                                        refined_annotations = NULL,
                                        selected_clusters = character(),
                                        arm_name = "DeepSeekCell-Refined") {
  first_pass_annotations <- add_ablation_ontology(
    first_pass_annotations,
    ontology_data,
    tissue
  )

  if (is.null(refined_annotations)) {
    refined_annotations <- first_pass_annotations[0, , drop = FALSE]
  }

  refined_annotations <- add_ablation_ontology(refined_annotations, ontology_data, tissue)
  refined <- replace_ablation_rows(first_pass_annotations, refined_annotations)
  refined <- add_ablation_ontology(refined, ontology_data, tissue)
  refined <- calibrate_annotation_confidence(refined, markers, tissue)

  change_summary <- summarise_ablation_refinement_changes(
    first_pass_annotations,
    refined_annotations
  )
  refined <- add_ablation_refinement_provenance(
    refined,
    first_pass_annotations,
    flagged_clusters = selected_clusters,
    refined_clusters = refined_annotations$Cluster,
    label_changed_clusters = change_summary$label_changed_clusters,
    confidence_changed_clusters = change_summary$confidence_changed_clusters
  )
  refined$AblationArm <- arm_name
  refined
}

write_first_pass_cache <- function(response_text,
                                   cache_file,
                                   metadata = list()) {
  dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
  hash <- hash_text_md5(response_text)
  payload <- list(
    response_text = response_text,
    response_hash = hash,
    metadata = metadata,
    timestamp = Sys.time()
  )
  saveRDS(payload, cache_file)
  payload
}

read_first_pass_cache <- function(cache_file) {
  readRDS(cache_file)
}

compute_refinement_behavior <- function(plain_annotations,
                                        refined_annotations,
                                        truth) {
  common <- intersect(intersect(plain_annotations$Cluster, refined_annotations$Cluster), names(truth))
  if (length(common) == 0) {
    stop("No overlapping clusters among plain, refined, and truth labels.", call. = FALSE)
  }

  plain <- plain_annotations[match(common, plain_annotations$Cluster), , drop = FALSE]
  refined <- refined_annotations[match(common, refined_annotations$Cluster), , drop = FALSE]
  truth <- truth[common]

  plain_correct <- normalize_cell_type(plain$CellType) == normalize_cell_type(truth)
  refined_correct <- normalize_cell_type(refined$CellType) == normalize_cell_type(truth)
  flagged <- if ("WasFlagged" %in% names(refined)) as_flag(refined$WasFlagged) else rep(FALSE, length(common))
  was_refined <- if ("WasRefined" %in% names(refined)) as_flag(refined$WasRefined) else rep(FALSE, length(common))
  label_changed <- normalize_cell_type(plain$CellType) != normalize_cell_type(refined$CellType)

  wrong_to_correct <- !plain_correct & refined_correct
  correct_to_wrong <- plain_correct & !refined_correct
  true_positive <- sum(flagged & !plain_correct, na.rm = TRUE)
  false_positive <- sum(flagged & plain_correct, na.rm = TRUE)
  true_negative <- sum(!flagged & plain_correct, na.rm = TRUE)
  false_negative <- sum(!flagged & !plain_correct, na.rm = TRUE)
  mcc_denominator <- sqrt(
    (true_positive + false_positive) *
      (true_positive + false_negative) *
      (true_negative + false_positive) *
      (true_negative + false_negative)
  )

  data.frame(
    NClusters = length(common),
    NFlagged = sum(flagged, na.rm = TRUE),
    NRefined = sum(was_refined, na.rm = TRUE),
    FlaggedRate = mean(flagged, na.rm = TRUE),
    NInitiallyIncorrect = sum(!plain_correct, na.rm = TRUE),
    SelectionPrecision = safe_divide(true_positive, true_positive + false_positive),
    SelectionRecall = safe_divide(true_positive, true_positive + false_negative),
    SelectionSpecificity = safe_divide(true_negative, true_negative + false_positive),
    SelectionNPV = safe_divide(true_negative, true_negative + false_negative),
    SelectionMCC = safe_divide(
      true_positive * true_negative - false_positive * false_negative,
      mcc_denominator
    ),
    ConflictPrecision = safe_divide(true_positive, true_positive + false_positive),
    ConflictRecall = safe_divide(true_positive, true_positive + false_negative),
    NLabelChanged = sum(label_changed, na.rm = TRUE),
    WrongToCorrect = sum(wrong_to_correct, na.rm = TRUE),
    CorrectToWrong = sum(correct_to_wrong, na.rm = TRUE),
    NetCorrectionRate = safe_divide(
      sum(wrong_to_correct, na.rm = TRUE) - sum(correct_to_wrong, na.rm = TRUE),
      sum(was_refined, na.rm = TRUE)
    ),
    CorrectionEfficiency = safe_divide(
      sum(wrong_to_correct, na.rm = TRUE),
      sum(was_refined, na.rm = TRUE)
    ),
    HarmRate = safe_divide(
      sum(correct_to_wrong, na.rm = TRUE),
      sum(was_refined, na.rm = TRUE)
    ),
    NetCorrectionPerChangedLabel = safe_divide(
      sum(wrong_to_correct, na.rm = TRUE) - sum(correct_to_wrong, na.rm = TRUE),
      sum(label_changed, na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )
}

compute_confidence_quality <- function(correct,
                                       confidence,
                                       n_bins = 10,
                                       thresholds = c(0.7, 0.8, 0.9)) {
  correct <- as.numeric(correct)
  confidence <- pmin(pmax(as.numeric(confidence), 1e-6), 1 - 1e-6)
  keep <- !is.na(correct) & !is.na(confidence)
  correct <- correct[keep]
  confidence <- confidence[keep]

  if (length(correct) == 0) {
    stop("No complete correctness/confidence pairs.", call. = FALSE)
  }

  bins <- cut(
    confidence,
    breaks = seq(0, 1, length.out = n_bins + 1),
    include.lowest = TRUE,
    labels = FALSE
  )

  ece <- 0
  for (bin in seq_len(n_bins)) {
    idx <- bins == bin
    if (!any(idx)) next
    ece <- ece + mean(idx) * abs(mean(confidence[idx]) - mean(correct[idx]))
  }

  threshold_acc <- vapply(thresholds, function(threshold) {
    idx <- confidence >= threshold
    if (!any(idx)) return(NA_real_)
    mean(correct[idx])
  }, numeric(1))
  names(threshold_acc) <- paste0("AccuracyAt", thresholds)

  data.frame(
    Brier = mean((confidence - correct)^2),
    BinaryCorrectnessNLL = -mean(
      correct * log(confidence) + (1 - correct) * log(1 - confidence)
    ),
    NLL = -mean(correct * log(confidence) + (1 - correct) * log(1 - confidence)),
    ECE = ece,
    AUROC = binary_auc(confidence, correct),
    AURC = area_under_risk_coverage(confidence, correct),
    MeanConfidence = mean(confidence),
    Accuracy = mean(correct),
    CoverageAt0.7 = mean(confidence >= 0.7),
    CoverageAt0.8 = mean(confidence >= 0.8),
    CoverageAt0.9 = mean(confidence >= 0.9),
    as.list(threshold_acc),
    check.names = FALSE
  )
}

compute_reliability_bins <- function(correct,
                                     confidence,
                                     n_bins = 10) {
  correct <- as.numeric(correct)
  confidence <- pmin(pmax(as.numeric(confidence), 0), 1)
  keep <- !is.na(correct) & !is.na(confidence)
  correct <- correct[keep]
  confidence <- confidence[keep]

  bins <- cut(
    confidence,
    breaks = seq(0, 1, length.out = n_bins + 1),
    include.lowest = TRUE,
    labels = FALSE
  )

  rows <- lapply(seq_len(n_bins), function(bin) {
    idx <- bins == bin
    data.frame(
      Bin = bin,
      BinLower = (bin - 1) / n_bins,
      BinUpper = bin / n_bins,
      N = sum(idx),
      Coverage = mean(idx),
      MeanConfidence = if (any(idx)) mean(confidence[idx]) else NA_real_,
      Accuracy = if (any(idx)) mean(correct[idx]) else NA_real_,
      AbsCalibrationError = if (any(idx)) {
        abs(mean(confidence[idx]) - mean(correct[idx]))
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  })

  dplyr::bind_rows(rows)
}

add_ablation_ontology <- function(annotations, ontology_data, tissue) {
  if (is.null(ontology_data)) {
    return(annotations)
  }
  if (nrow(annotations) == 0) {
    annotations$CL_ID <- character()
    annotations$OntologyLabel <- character()
    annotations$MatchMethod <- character()
    annotations$OntologyMatchScore <- numeric()
    return(annotations)
  }

  mappings <- do.call(
    rbind,
    lapply(
      annotations$CellType,
      map_to_cell_ontology,
      ontology = ontology_data,
      tissue = tissue
    )
  )
  annotations$CL_ID <- mappings$CL_ID
  annotations$OntologyLabel <- mappings$OntologyLabel
  annotations$MatchMethod <- mappings$MatchMethod
  annotations$OntologyMatchScore <- mappings$OntologyMatchScore
  annotations
}

replace_ablation_rows <- function(annotations, refined_annotations) {
  cols <- setdiff(intersect(names(annotations), names(refined_annotations)), "Cluster")
  for (i in seq_len(nrow(refined_annotations))) {
    idx <- match_ablation_clusters(refined_annotations$Cluster[[i]], annotations$Cluster)
    if (is.na(idx)) next
    for (col in cols) {
      annotations[[col]][idx] <- refined_annotations[[col]][i]
    }
  }
  annotations
}

summarise_ablation_refinement_changes <- function(first_pass_annotations, refined_annotations) {
  idx <- match_ablation_clusters(refined_annotations$Cluster, first_pass_annotations$Cluster)
  valid <- !is.na(idx)
  refined <- refined_annotations[valid, , drop = FALSE]
  first <- first_pass_annotations[idx[valid], , drop = FALSE]

  label_changed <- normalize_cell_type(refined$CellType) != normalize_cell_type(first$CellType)
  confidence_changed <- abs(as_confidence(refined$Confidence) - as_confidence(first$Confidence)) > 1e-6

  list(
    label_changed_clusters = refined$Cluster[label_changed],
    confidence_changed_clusters = refined$Cluster[confidence_changed]
  )
}

add_ablation_refinement_provenance <- function(annotations,
                                              first_pass_annotations,
                                              flagged_clusters,
                                              refined_clusters,
                                              label_changed_clusters,
                                              confidence_changed_clusters) {
  idx <- match(annotations$Cluster, first_pass_annotations$Cluster)
  annotations$FirstPassCellType <- first_pass_annotations$CellType[idx]
  annotations$FirstPassConfidence <- as_confidence(first_pass_annotations$Confidence[idx])
  annotations$WasFlagged <- annotations$Cluster %in% flagged_clusters
  annotations$WasRefined <- annotations$Cluster %in% refined_clusters
  annotations$RefinementChangedLabel <- annotations$Cluster %in% label_changed_clusters
  annotations$RefinementChangedConfidence <- annotations$Cluster %in% confidence_changed_clusters
  annotations
}

normalise_ablation_cluster_id <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("^cluster[[:space:]_:-]*", "", x, perl = TRUE)
  x <- gsub("[[:space:]]+", "", x, perl = TRUE)
  x
}

match_ablation_clusters <- function(query, reference) {
  exact <- match(query, reference)
  unresolved <- which(is.na(exact))

  if (length(unresolved) == 0) {
    return(exact)
  }

  query_norm <- normalise_ablation_cluster_id(query[unresolved])
  reference_norm <- normalise_ablation_cluster_id(reference)
  norm_match <- match(query_norm, reference_norm)
  exact[unresolved] <- norm_match
  exact
}

hash_text_md5 <- function(text) {
  tmp <- tempfile(fileext = ".txt")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(text, tmp, useBytes = TRUE)
  unname(tools::md5sum(tmp))
}

safe_divide <- function(numerator, denominator) {
  out <- numerator / denominator
  out[is.na(denominator) | denominator == 0] <- NA_real_
  out
}

binary_auc <- function(score, truth) {
  truth <- as.integer(truth)
  if (length(unique(truth)) < 2) {
    return(NA_real_)
  }

  positive <- score[truth == 1]
  negative <- score[truth == 0]
  ranks <- rank(c(positive, negative), ties.method = "average")
  n_pos <- length(positive)
  n_neg <- length(negative)
  (sum(ranks[seq_len(n_pos)]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

area_under_risk_coverage <- function(confidence, correct) {
  ord <- order(confidence, decreasing = TRUE, na.last = NA)
  if (length(ord) == 0) {
    return(NA_real_)
  }

  correct <- as.numeric(correct)[ord]
  coverage <- c(0, seq_along(correct) / length(correct))
  risk <- c(0, 1 - cumsum(correct) / seq_along(correct))

  sum(diff(coverage) * (head(risk, -1) + tail(risk, -1)) / 2)
}
