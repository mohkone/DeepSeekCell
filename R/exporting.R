#' Export Functions for Annotation Results
#' 
#' Supports multiple export formats: CSV, Excel, JSON, RDS, and HTML reports.
#' 
#' @import openxlsx jsonlite logger

#' Export annotations to various formats
#' 
#' @param result Annotation result object from annotate_cell_types()
#' @param filename Base filename (without extension)
#' @param format Export format ("csv", "xlsx", "json", "rds")
#' @return Path to exported file
#' @export
export_annotations <- function(result, filename, format = "csv") {
  
  if (!result$success) {
    stop("Cannot export: annotation failed")
  }
  
  df <- result$annotations
  
  outfile <- switch(format,
                    csv = .export_csv(df, filename),
                    xlsx = .export_xlsx(result, filename),
                    json = .export_json(result, filename),
                    rds = .export_rds(result, filename),
                    stop("Unknown format: ", format)
  )
  
  message("Exported to: {outfile}")
  return(outfile)
}

.export_csv <- function(df, filename) {
  outfile <- paste0(filename, ".csv")
  write.csv(df, outfile, row.names = FALSE)
  return(outfile)
}

.export_xlsx <- function(result, filename) {
  outfile <- paste0(filename, ".xlsx")
  
  wb <- openxlsx::createWorkbook()
  
  # Annotations sheet
  openxlsx::addWorksheet(wb, "Annotations")
  openxlsx::writeData(wb, "Annotations", result$annotations)
  
  # Metadata sheet
  openxlsx::addWorksheet(wb, "Metadata")
  metadata_df <- data.frame(
    Parameter = names(result$metadata),
    Value = as.character(unlist(result$metadata))
  )
  openxlsx::writeData(wb, "Metadata", metadata_df)
  
  # Validation sheet (if available)
  if (!is.null(result$validation)) {
    openxlsx::addWorksheet(wb, "Validation")
    openxlsx::writeData(wb, "Validation", result$validation$summary)
  }
  
  openxlsx::saveWorkbook(wb, outfile, overwrite = TRUE)
  return(outfile)
}

.export_json <- function(result, filename) {
  outfile <- paste0(filename, ".json")
  
  export_list <- list(
    metadata = result$metadata,
    annotations = result$annotations,
    validation = result$validation,
    timestamp = Sys.time()
  )
  
  jsonlite::write_json(export_list, outfile, pretty = TRUE, auto_unbox = TRUE)
  return(outfile)
}

.export_rds <- function(result, filename) {
  outfile <- paste0(filename, ".rds")
  saveRDS(result, outfile)
  return(outfile)
}

#' Generate comprehensive HTML report
#' 
#' @param result Annotation result object
#' @param output_file Output file path (optional)
#' @return Path to generated HTML file
#' @export
generate_html_report <- function(result, output_file = NULL) {
  
  if (is.null(output_file)) {
    output_file <- paste0("annotation_report_", 
                          format(Sys.time(), "%Y%m%d_%H%M%S"), ".html")
  }
  
  # Create summary statistics
  summary_stats <- .create_summary_stats(result)
  
  # Generate HTML
  html_content <- .render_html_template(result, summary_stats)
  
  writeLines(html_content, output_file)
  
  # Save confidence plot
  p <- plot_confidence(result)
  ggplot2::ggsave("confidence_plot.png", p, width = 8, height = 6, dpi = 150)
  
  message("HTML report saved to: {output_file}")
  return(output_file)
}

.create_summary_stats <- function(result) {
  data.frame(
    Metric = c("Tissue", "Species", "Model", "Number of Clusters", 
               "Total Runtime (sec)", "API Latency (sec)", "Tokens Used",
               "Estimated Cost (USD)", "Mean Confidence"),
    Value = c(
      result$metadata$tissue,
      result$metadata$species,
      result$metadata$model,
      result$metadata$n_clusters,
      sprintf("%.2f", result$metadata$total_runtime_sec),
      sprintf("%.2f", result$metadata$api_latency_sec),
      result$metadata$tokens_used,
      sprintf("$%.4f", result$metadata$estimated_cost_usd),
      sprintf("%.3f", mean(result$annotations$Confidence, na.rm = TRUE))
    )
  )
}

.render_html_template <- function(result, summary_stats) {
  sprintf(
    '<!DOCTYPE html>
    <html>
    <head>
      <title>DeepSeekCell Annotation Report</title>
      <style>
        body { font-family: Arial, sans-serif; margin: 40px; line-height: 1.6; }
        h1 { color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }
        h2 { color: #34495e; margin-top: 30px; }
        .summary { background-color: #ecf0f1; padding: 20px; border-radius: 8px; margin: 20px 0; }
        .metric { margin: 8px 0; }
        .metric-label { font-weight: bold; display: inline-block; width: 220px; }
        .success { color: #27ae60; }
        .warning { color: #e67e22; }
        .error { color: #e74c3c; }
        table { border-collapse: collapse; width: 100%%; margin: 20px 0; }
        th, td { border: 1px solid #ddd; padding: 12px; text-align: left; }
        th { background-color: #3498db; color: white; }
        tr:nth-child(even) { background-color: #f2f2f2; }
        .footer { margin-top: 50px; font-size: 12px; color: #7f8c8d; text-align: center; }
      </style>
    </head>
    <body>
      <h1>🧬 DeepSeekCell Annotation Report</h1>
      <p>Generated: %s | Version: 2.0</p>
      
      <div class="summary">
        <h2>📊 Summary Statistics</h2>
        %s
      </div>
      
      <div>
        <h2>📈 Confidence Plot</h2>
        <img src="confidence_plot.png" alt="Confidence Plot" style="max-width: 100%%; border: 1px solid #ddd; border-radius: 4px;">
      </div>
      
      <div>
        <h2>📋 Annotation Results</h2>
        %s
      </div>
      
      <div class="summary">
        <h2>✅ Validation</h2>
        <p>Status: %s</p>
        <p>Issues: %s</p>
      </div>
      
      <div class="footer">
        <p>Generated by DeepSeekCell | <a href="https://github.com/...">Source Code</a></p>
      </div>
    </body>
    </html>',
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    paste(sprintf('<div class="metric"><span class="metric-label">%s:</span> %s</div>', 
                  summary_stats$Metric, summary_stats$Value), collapse = "\n"),
    knitr::kable(result$annotations, format = "html", table.attr = 'class="data-table"'),
    ifelse(result$validation$valid, '<span class="success">✓ Valid</span>', '<span class="warning">⚠ Issues Found</span>'),
    paste(result$validation$issues, collapse = "<br>")
  )
}