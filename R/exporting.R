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
  format <- match.arg(format, c("csv", "xlsx", "json", "rds"))
  
  if (!isTRUE(result$success)) {
    stop("Cannot export: annotation failed.", call. = FALSE)
  }
  
  df <- result$annotations
  
  outfile <- switch(format,
                    csv = .export_csv(df, filename),
                    xlsx = .export_xlsx(result, filename),
                    json = .export_json(result, filename),
                    rds = .export_rds(result, filename),
                    stop("Unknown format: ", format, call. = FALSE)
  )
  
  message("Exported to: ", outfile)
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
  if (!isTRUE(result$success)) {
    stop("Cannot generate report: annotation failed.", call. = FALSE)
  }
  
  if (is.null(output_file)) {
    output_file <- paste0("annotation_report_", 
                          format(Sys.time(), "%Y%m%d_%H%M%S"), ".html")
  }
  
  # Create summary statistics
  summary_stats <- .create_summary_stats(result)
  
  # Generate HTML
  html_content <- .render_html_template(result, summary_stats)
  
  writeLines(html_content, output_file, useBytes = TRUE)
  
  message("HTML report saved to: ", output_file)
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
  
  validation_status <- if (is.null(result$validation)) {
    '<span class="warning">Not performed</span>'
  } else if (isTRUE(result$validation$valid)) {
    '<span class="success">Valid</span>'
  } else {
    '<span class="warning">Issues found</span>'
  }
  
  validation_issues <- if (is.null(result$validation)) {
    "Validation was not performed."
  } else if (length(result$validation$issues) == 0) {
    "None"
  } else {
    paste(html_escape(result$validation$issues), collapse = "<br>")
  }

  validation_warnings <- if (is.null(result$validation) || length(result$validation$warnings) == 0) {
    "None"
  } else {
    paste(html_escape(result$validation$warnings), collapse = "<br>")
  }

  summary_html <- paste(
    sprintf(
      '<div class="metric"><span class="metric-label">%s:</span> %s</div>',
      html_escape(summary_stats$Metric),
      html_escape(summary_stats$Value)
    ),
    collapse = "\n"
  )

  annotations_html <- .render_html_table(result$annotations)
  
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
        table { border-collapse: collapse; width: 100%%; margin: 20px 0; }
        th, td { border: 1px solid #ddd; padding: 12px; text-align: left; }
        th { background-color: #3498db; color: white; }
        tr:nth-child(even) { background-color: #f2f2f2; }
        .footer { margin-top: 50px; font-size: 12px; color: #7f8c8d; text-align: center; }
      </style>
    </head>
    <body>
      <h1>DeepSeekCell Annotation Report</h1>
      <p>Generated: %s | Schema version: %s</p>

      <div class="summary">
        <h2>Summary Statistics</h2>
        %s
      </div>

      <div>
        <h2>Annotation Results</h2>
        %s
      </div>

      <div class="summary">
        <h2>Validation</h2>
        <p>Status: %s</p>
        <p>Issues: %s</p>
        <p>Warnings: %s</p>
      </div>

      <div class="footer">
        <p>Generated by DeepSeekCell</p>
      </div>
    </body>
    </html>',
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    html_escape(result$metadata$schema_version %||% deepseekcell_version()),
    summary_html,
    annotations_html,
    validation_status,
    validation_issues,
    validation_warnings
  )
}

.render_html_table <- function(df) {
  if (!is.data.frame(df) || nrow(df) == 0) {
    return("<p>No annotation rows available.</p>")
  }

  header <- paste(sprintf("<th>%s</th>", html_escape(names(df))), collapse = "")
  rows <- apply(df, 1, function(row) {
    paste0(
      "<tr>",
      paste(sprintf("<td>%s</td>", html_escape(row)), collapse = ""),
      "</tr>"
    )
  })

  paste0(
    '<table class="data-table"><thead><tr>',
    header,
    "</tr></thead><tbody>",
    paste(rows, collapse = "\n"),
    "</tbody></table>"
  )
}
