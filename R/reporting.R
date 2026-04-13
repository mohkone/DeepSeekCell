#' Reporting and Summary Functions
#' 
#' Generates comprehensive reports and summaries from annotation results.

#' Generate annotation report
#' 
#' @param result Annotation result object
#' @return List with summary statistics
#' @export
generate_annotation_report <- function(result) {
  
  if (!result$success) {
    return(list(error = result$error))
  }
  
  df <- result$annotations
  
  list(
    summary = list(
      n_clusters = nrow(df),
      n_unique_cell_types = length(unique(df$CellType)),
      mean_confidence = mean(df$Confidence, na.rm = TRUE),
      sd_confidence = sd(df$Confidence, na.rm = TRUE),
      min_confidence = min(df$Confidence, na.rm = TRUE),
      max_confidence = max(df$Confidence, na.rm = TRUE),
      ontology_coverage = sum(!is.na(df$CL_ID)) / nrow(df),
      total_cost_usd = result$metadata$estimated_cost_usd,
      total_time_sec = result$metadata$total_runtime_sec
    ),
    confidence_distribution = .get_confidence_distribution(df),
    cell_type_frequencies = table(df$CellType),
    ontology_match_summary = .get_ontology_summary(df)
  )
}

.get_confidence_distribution <- function(df) {
  breaks <- seq(0, 1, by = 0.1)
  cut(df$Confidence, breaks = breaks, include.lowest = TRUE)
}

.get_ontology_summary <- function(df) {
  if (!"MatchMethod" %in% colnames(df)) {
    return(NULL)
  }
  table(df$MatchMethod)
}

#' Print formatted report to console
#' 
#' @param result Annotation result object
#' @export
print_annotation_report <- function(result) {
  
  report <- generate_annotation_report(result)
  
  if ("error" %in% names(report)) {
    cat("❌ Annotation failed:", report$error, "\n")
    return(invisible(NULL))
  }
  
  cat("\n", rep("=", 60), "\n", sep = "")
  cat("📊 DeepSeekCell Annotation Report\n")
  cat(rep("=", 60), "\n\n", sep = "")
  
  cat("📈 Summary Statistics:\n")
  cat(sprintf("  • Clusters annotated: %d\n", report$summary$n_clusters))
  cat(sprintf("  • Unique cell types: %d\n", report$summary$n_unique_cell_types))
  cat(sprintf("  • Mean confidence: %.3f ± %.3f\n", 
              report$summary$mean_confidence, report$summary$sd_confidence))
  cat(sprintf("  • Ontology coverage: %.1f%%\n", 
              report$summary$ontology_coverage * 100))
  
  cat("\n💰 Cost & Performance:\n")
  cat(sprintf("  • Estimated cost: $%.4f\n", report$summary$total_cost_usd))
  cat(sprintf("  • Total time: %.2f seconds\n", report$summary$total_time_sec))
  
  cat("\n📋 Cell Type Distribution:\n")
  print(report$cell_type_frequencies)
  
  cat("\n", rep("=", 60), "\n", sep = "")
}

#' Save report to file
#' 
#' @param result Annotation result object
#' @param filename Output filename
#' @param format Report format ("txt", "json", "yaml")
#' @export
save_report <- function(result, filename, format = "txt") {
  
  report <- generate_annotation_report(result)
  
  switch(format,
         txt = .save_txt_report(report, filename),
         json = .save_json_report(report, filename),
         yaml = .save_yaml_report(report, filename),
         stop("Unknown format: ", format)
  )
}

.save_txt_report <- function(report, filename) {
  sink(filename)
  print_annotation_report(report)
  sink()
}

.save_json_report <- function(report, filename) {
  jsonlite::write_json(report, paste0(filename, ".json"), pretty = TRUE)
}

.save_yaml_report <- function(report, filename) {
  yaml::write_yaml(report, paste0(filename, ".yaml"))
}