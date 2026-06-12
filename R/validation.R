#' Validation and Quality Control
#' 
#' Functions for validating annotation results and quality control.

#' Validate annotations for quality control
#' 
#' @param annotations_df Data frame of annotations
#' @param thresholds List of validation thresholds
#' @param markers Optional normalized marker list used to generate annotations.
#' @param metadata Optional run metadata.
#' @return Validation result object
#' @export
validate_annotations <- function(annotations_df, 
                                 thresholds = list(
                                   min_confidence = 0.5,
                                   max_unknown_rate = 0.3,
                                   min_ontology_coverage = 0.5
                                 ),
                                 markers = NULL,
                                 metadata = NULL) {
  default_thresholds <- list(
    min_confidence = 0.5,
    max_unknown_rate = 0.3,
    min_ontology_coverage = 0.5
  )
  thresholds <- utils::modifyList(default_thresholds, thresholds)
  
  if (!is.data.frame(annotations_df) || nrow(annotations_df) == 0) {
    return(list(
      valid = FALSE,
      issues = "No annotations to validate",
      summary = data.frame()
    ))
  }
  
  issues <- c()
  warnings <- c()

  required_cols <- c("Cluster", "CellType", "Confidence")
  missing_cols <- setdiff(required_cols, colnames(annotations_df))

  if (length(missing_cols) > 0) {
    return(list(
      valid = FALSE,
      issues = paste("Annotation table is missing required columns:", paste(missing_cols, collapse = ", ")),
      warnings = character(),
      summary = data.frame(),
      timestamp = Sys.time()
    ))
  }

  annotations_df$Confidence <- as_confidence(annotations_df$Confidence)
  
  # Check confidence scores
  low_conf <- sum(annotations_df$Confidence < thresholds$min_confidence, na.rm = TRUE)
  if (low_conf > 0) {
    warnings <- c(warnings, sprintf("%d clusters have low confidence (< %.2f)", 
                                    low_conf, thresholds$min_confidence))
  }
  
  # Check for unknown annotations
  unknown_pattern <- "unknown|unidentified|not determined"
  unknown_count <- sum(grepl(unknown_pattern, annotations_df$CellType, ignore.case = TRUE))
  unknown_rate <- unknown_count / nrow(annotations_df)
  
  if (unknown_rate > thresholds$max_unknown_rate) {
    issues <- c(issues, sprintf("Unknown rate (%.1f%%) exceeds threshold (%.1f%%)",
                                unknown_rate * 100, thresholds$max_unknown_rate * 100))
  }
  
  # Check ontology coverage
  ontology_coverage <- NA_real_
  if ("CL_ID" %in% colnames(annotations_df)) {
    ontology_missing <- sum(is.na(annotations_df$CL_ID))
    ontology_coverage <- 1 - (ontology_missing / nrow(annotations_df))
    
    if (ontology_coverage < thresholds$min_ontology_coverage) {
      warnings <- c(warnings, sprintf("Low ontology coverage: %.1f%%", 
                                      ontology_coverage * 100))
    }
  }
  
  # Check for duplicate cluster names
  if (any(duplicated(annotations_df$Cluster))) {
    issues <- c(issues, "Duplicate cluster names found")
  }

  if (!is.null(markers)) {
    missing_marker_clusters <- setdiff(annotations_df$Cluster, names(markers))
    missing_annotation_clusters <- setdiff(names(markers), annotations_df$Cluster)

    if (length(missing_marker_clusters) > 0) {
      warnings <- c(
        warnings,
        sprintf(
          "%d annotated clusters were not present in the marker input",
          length(missing_marker_clusters)
        )
      )
    }

    if (length(missing_annotation_clusters) > 0) {
      issues <- c(
        issues,
        sprintf(
          "%d marker clusters are missing annotations",
          length(missing_annotation_clusters)
        )
      )
    }
  }

  mean_confidence <- mean(annotations_df$Confidence, na.rm = TRUE)
  if (is.nan(mean_confidence)) mean_confidence <- NA_real_

  mixed_rate <- if ("IsMixed" %in% colnames(annotations_df)) {
    mean(as_flag(annotations_df$IsMixed), na.rm = TRUE)
  } else {
    NA_real_
  }

  summary <- data.frame(
    n_clusters = nrow(annotations_df),
    mean_confidence = mean_confidence,
    min_confidence = min(annotations_df$Confidence, na.rm = TRUE),
    max_confidence = max(annotations_df$Confidence, na.rm = TRUE),
    unknown_rate = unknown_rate,
    mixed_rate = mixed_rate,
    ontology_coverage = ontology_coverage,
    stringsAsFactors = FALSE
  )

  summary$min_confidence[is.infinite(summary$min_confidence)] <- NA_real_
  summary$max_confidence[is.infinite(summary$max_confidence)] <- NA_real_
  
  result <- list(
    valid = length(issues) == 0,
    issues = issues,
    warnings = warnings,
    summary = summary,
    thresholds = thresholds,
    metadata = metadata,
    timestamp = Sys.time()
  )

  result$quality_score <- calculate_quality_score(result)
  result
}

#' Quality control score for annotations
#' 
#' @param validation_result Output from validate_annotations()
#' @return Quality score between 0 and 100
#' @export
calculate_quality_score <- function(validation_result) {
  
  score <- 100
  
  # Deduct for issues
  score <- score - (length(validation_result$issues) * 20)
  
  # Deduct for warnings (half penalty)
  score <- score - (length(validation_result$warnings) * 5)
  
  # Adjust based on metrics
  summary <- validation_result$summary
  
  if (!is.na(summary$mean_confidence)) {
    score <- score + (summary$mean_confidence - 0.5) * 20
  }
  
  if (!is.na(summary$ontology_coverage)) {
    score <- score + (summary$ontology_coverage - 0.5) * 10
  }
  
  # Clamp to 0-100
  return(min(max(round(score), 0), 100))
}
