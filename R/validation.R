#' Validation and Quality Control
#' 
#' Functions for validating annotation results and quality control.

#' Validate annotations for quality control
#' 
#' @param annotations_df Data frame of annotations
#' @param thresholds List of validation thresholds
#' @return Validation result object
#' @export
validate_annotations <- function(annotations_df, 
                                 thresholds = list(
                                   min_confidence = 0.5,
                                   max_unknown_rate = 0.3,
                                   min_ontology_coverage = 0.5
                                 )) {
  
  if (nrow(annotations_df) == 0) {
    return(list(
      valid = FALSE,
      issues = "No annotations to validate",
      summary = data.frame()
    ))
  }
  
  issues <- c()
  warnings <- c()
  
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
  
  list(
    valid = length(issues) == 0,
    issues = issues,
    warnings = warnings,
    summary = data.frame(
      n_clusters = nrow(annotations_df),
      mean_confidence = mean(annotations_df$Confidence, na.rm = TRUE),
      unknown_rate = unknown_count / nrow(annotations_df),
      ontology_coverage = ifelse("CL_ID" %in% colnames(annotations_df),
                                 sum(!is.na(annotations_df$CL_ID)) / nrow(annotations_df),
                                 NA),
      stringsAsFactors = FALSE
    ),
    timestamp = Sys.time()
  )
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