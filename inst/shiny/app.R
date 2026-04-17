#!/usr/bin/env Rscript
# DeepSeekCell Shiny Application – Publication Version
suppressPackageStartupMessages({
  library(shiny)
  library(shinythemes)
  library(shinycssloaders)
  library(DT)
  library(ggplot2)
  library(plotly)
  library(dplyr)
})

# Source R functions (exclude benchmark scripts)
r_dir <- file.path("..", "..", "R")
if (!dir.exists(r_dir)) r_dir <- "R"
r_files <- list.files(r_dir, pattern = "\\.R$", full.names = TRUE)
r_files <- r_files[!grepl("benchmark|^main\\.R$", basename(r_files))]
for (f in r_files) source(f, local = TRUE)

# Define %||% operator if not already present
`%||%` <- function(x, y) if (is.null(x)) y else x

ui <- fluidPage(
  theme = shinytheme("flatly"),
  tags$head(tags$style(HTML("
    .live-badge { background-color: #27ae60; color: white; padding: 5px 12px;
                  border-radius: 20px; font-size: 14px; font-weight: bold; }
    .confidence-high { color: #27ae60; font-weight: bold; }
    .confidence-mid { color: #f39c12; }
    .confidence-low { color: #e74c3c; }
    .metric-box {
      background: white; border: 1px solid #dee2e6; border-radius: 8px;
      padding: 15px; text-align: center; margin: 10px 0; box-shadow: 0 2px 4px rgba(0,0,0,0.1);
    }
    .metric-value { font-size: 24px; font-weight: bold; color: #2c3e50; }
    .metric-label { font-size: 12px; color: #7f8c8d; margin-top: 5px; }
  "))),
  div(class = "text-center", style = "background-color: #2c3e50; color: white; padding: 20px;",
      h1("🧬 DeepSeekCell", span(class = "live-badge", "LIVE MODE")),
      h4("AI-Powered Cell Type Annotation with Ontology-Aware Reasoning"),
      p("Version 2.0 | For research use only")),
  sidebarLayout(
    sidebarPanel(width = 4,
                 textInput("api_key", "API Key", placeholder = "sk-..."),
                 selectInput("model", "Model", choices = c("DeepSeek" = "deepseek", "GPT-4o" = "gpt4")),
                 textInput("tissue", "Tissue", value = "PBMC"),
                 selectInput("species", "Species", choices = c("Human", "Mouse", "Rat")),
                 checkboxInput("use_ontology", "Map to Cell Ontology", value = TRUE),
                 textAreaInput("c1", "Cluster 1", rows = 2, placeholder = "CD3D, CD3E, CD8A"),
                 textAreaInput("c2", "Cluster 2", rows = 2, placeholder = "CD14, LYZ, FCGR3A"),
                 textAreaInput("c3", "Cluster 3", rows = 2, placeholder = "CD79A, MS4A1, CD19"),
                 textAreaInput("c4", "Cluster 4", rows = 2, placeholder = "NCAM1, KLRB1, NKG7"),
                 textAreaInput("c5", "Cluster 5", rows = 2, placeholder = "PPBP, PF4, ITGA2B"),
                 actionButton("run", "Annotate", class = "btn-primary btn-lg", width = "100%"),
                 hr(),
                 fluidRow(
                   column(6, actionButton("example_pbmc", "PBMC", class = "btn-info btn-sm", width = "100%")),
                   column(6, actionButton("example_pancreas", "Pancreas", class = "btn-info btn-sm", width = "100%"))
                 ),
                 br(),
                 actionButton("clear", "Clear All", class = "btn-warning btn-sm", width = "100%")
    ),
    mainPanel(width = 8,
              tabsetPanel(
                tabPanel("Results",
                         br(),
                         DTOutput("results_table") %>% withSpinner(),
                         br(),
                         fluidRow(
                           column(4, downloadButton("download_csv", "CSV", class = "btn-success", width = "100%")),
                           column(4, downloadButton("download_xlsx", "Excel", class = "btn-success", width = "100%")),
                           column(4, downloadButton("download_report", "HTML Report", class = "btn-info", width = "100%"))
                         )),
                tabPanel("Confidence", br(), plotlyOutput("confidence_plot") %>% withSpinner()),
                tabPanel("Cost & Performance", br(),
                         fluidRow(
                           column(4, uiOutput("cost_metric")),
                           column(4, uiOutput("time_metric")),
                           column(4, uiOutput("tokens_metric"))
                         ),
                         hr(),
                         h4("Detailed Metrics"),
                         tableOutput("performance_table"),
                         br(),
                         h4("Validation Report"),
                         verbatimTextOutput("validation_report")),
                tabPanel("Metadata", br(), verbatimTextOutput("metadata")),
                tabPanel("Help", br(),
                         div(class = "well",
                             h4("How to Use DeepSeekCell"),
                             tags$ol(
                               tags$li("Enter your API key for DeepSeek or GPT-4o."),
                               tags$li("Specify tissue and species."),
                               tags$li("Input marker genes for up to 5 clusters (comma‑separated)."),
                               tags$li("Click 'Annotate' and wait 15–30 seconds."),
                               tags$li("Download results as CSV, Excel, or HTML report.")
                             ),
                             h4("Marker Gene Tips"),
                             tags$ul(
                               tags$li("Use 5–15 markers per cluster for best results."),
                               tags$li("Avoid mitochondrial (MT-), ribosomal (RP), and pseudogenes.")
                             ),
                             h4("Citation"),
                             p("If you use DeepSeekCell, please cite:"),
                             p(tags$em("DeepSeekCell: Benchmarking Large Language Model‑Powered Cell Type Annotation with Ontology‑Aware Evaluation and an Interactive Shiny Application. Computers in Biology and Medicine, 2026."))
                         ))
              )
    )
  )
)

server <- function(input, output, session) {
  values <- reactiveValues(result = NULL)
  
  observeEvent(input$example_pbmc, {
    updateTextAreaInput(session, "c1", value = "CD3D, CD3E, CD8A, CD4, CD247")
    updateTextAreaInput(session, "c2", value = "CD14, LYZ, FCGR3A, MS4A7, CST3")
    updateTextAreaInput(session, "c3", value = "CD79A, MS4A1, CD19, BANK1, CD22")
    updateTextAreaInput(session, "c4", value = "NCAM1, KLRB1, KLRD1, NKG7, GNLY")
    updateTextAreaInput(session, "c5", value = "PPBP, PF4, ITGA2B, GP1BB, TUBB1")
    updateTextInput(session, "tissue", value = "PBMC")
  })
  
  observeEvent(input$example_pancreas, {
    updateTextAreaInput(session, "c1", value = "GCG, GC, ARX, SLC38A4")
    updateTextAreaInput(session, "c2", value = "INS, IAPP, MAFA, PDX1")
    updateTextAreaInput(session, "c3", value = "SST, HHEX, RBP4, PCSK2")
    updateTextAreaInput(session, "c4", value = "PRSS1, CPA1, CTRB1, AMY2A")
    updateTextAreaInput(session, "c5", value = "KRT19, SOX9, HNF1B, TFF2")
    updateTextInput(session, "tissue", value = "Pancreas")
  })
  
  observeEvent(input$clear, {
    for (i in 1:5) updateTextAreaInput(session, paste0("c", i), value = "")
  })
  
  collect_markers <- reactive({
    markers <- list()
    for (i in 1:5) {
      txt <- input[[paste0("c", i)]]
      if (!is.null(txt) && nchar(trimws(txt)) > 0) {
        genes <- trimws(unlist(strsplit(txt, ",")))
        genes <- genes[genes != ""]
        if (length(genes) > 0) markers[[paste0("Cluster", i)]] <- genes
      }
    }
    if (length(markers) == 0) {
      showNotification("Enter at least one cluster.", type = "error")
      return(NULL)
    }
    markers
  })
  
  observeEvent(input$run, {
    markers <- collect_markers()
    if (is.null(markers)) return()
    
    if (is.null(input$api_key) || input$api_key == "") {
      showNotification("API key is required.", type = "error")
      return()
    }
    
    showNotification("Calling API... This may take 15-30 seconds.",
                     type = "message", duration = NULL, id = "status")
    
    res <- tryCatch(
      annotate_cell_types(
        markers = markers,
        tissue = input$tissue,
        species = input$species,
        model_name = input$model,          # FIXED: was input$44model
        api_key = input$api_key,
        use_ontology = input$use_ontology,
        validate = TRUE
      ),
      error = function(e) { list(success = FALSE, error = e$message) }
    )
    removeNotification(id = "status")
    values$result <- res
  })
  
  # Metric boxes
  render_metric_box <- function(value, subtitle, color = "blue") {
    colors <- list(green = "#27ae60", blue = "#3498db", red = "#e74c3c", yellow = "#f39c12")
    bg <- colors[[color]] %||% colors$blue
    div(class = "metric-box", style = paste0("border-top: 4px solid ", bg, ";"),
        div(class = "metric-value", value),
        div(class = "metric-label", subtitle))
  }
  
  output$cost_metric <- renderUI({
    req(values$result$success)
    cost <- values$result$metadata$estimated_cost_usd
    color <- ifelse(cost < 0.01, "green", ifelse(cost < 0.05, "yellow", "red"))
    render_metric_box(sprintf("$%.4f", cost), "Estimated Cost", color)
  })
  
  output$time_metric <- renderUI({
    req(values$result$success)
    time <- values$result$metadata$total_runtime_sec
    color <- ifelse(time < 20, "green", ifelse(time < 40, "yellow", "red"))
    render_metric_box(sprintf("%.1f sec", time), "Total Runtime", color)
  })
  
  output$tokens_metric <- renderUI({
    req(values$result$success)
    tokens <- values$result$metadata$tokens_used
    render_metric_box(prettyNum(tokens, big.mark = ","), "Tokens Used", "blue")
  })
  
  output$results_table <- renderDT({
    req(values$result$success)
    df <- values$result$annotations
    df$Confidence_Display <- sprintf('<span class="%s">%.3f</span>',
                                     ifelse(df$Confidence >= 0.7, "confidence-high",
                                            ifelse(df$Confidence >= 0.4, "confidence-mid", "confidence-low")),
                                     df$Confidence)
    display_cols <- c("Cluster", "CellType", "Confidence_Display")
    if ("CL_ID" %in% names(df)) display_cols <- c(display_cols, "CL_ID", "OntologyLabel")
    datatable(df[, display_cols, drop = FALSE], escape = FALSE,
              options = list(pageLength = 10), rownames = FALSE)
  })
  
  output$confidence_plot <- renderPlotly({
    req(values$result$success)
    df <- values$result$annotations
    df$Cluster <- factor(df$Cluster, levels = rev(unique(df$Cluster)))
    p <- ggplot(df, aes(x = Cluster, y = Confidence, fill = Confidence,
                        text = paste("Cluster:", Cluster, "<br>Cell Type:", CellType,
                                     "<br>Confidence:", round(Confidence, 3)))) +
      geom_col(width = 0.7) +
      scale_fill_gradient2(low = "#FF6B6B", mid = "#FFE66D", high = "#4ECDC4", midpoint = 0.5) +
      coord_flip(ylim = c(0, 1)) +
      labs(title = "Annotation Confidence", x = "Cluster", y = "Confidence") +
      theme_minimal()
    ggplotly(p, tooltip = "text")
  })
  
  output$performance_table <- renderTable({
    req(values$result$success)
    m <- values$result$metadata
    data.frame(
      Metric = c("API Latency (sec)", "Tokens per Second", "Cost per Token", "Cost per Cluster"),
      Value = c(
        sprintf("%.2f", m$api_latency_sec),
        sprintf("%.0f", m$tokens_used / max(m$api_latency_sec, 1)),
        sprintf("$%.8f", m$estimated_cost_usd / max(m$tokens_used, 1)),
        sprintf("$%.4f", m$estimated_cost_usd / m$n_clusters)
      )
    )
  })
  
  output$validation_report <- renderPrint({
    req(values$result$success)
    v <- values$result$validation
    if (!is.null(v)) {
      cat("=== Validation Report ===\n\n")
      if (v$valid) cat("✅ All checks passed.\n\n") else cat("❌ Issues found:\n", paste("•", v$issues), "\n")
      if (length(v$warnings)) cat("⚠️ Warnings:\n", paste("•", v$warnings), "\n")
      cat("Summary Statistics:\n")
      print(v$summary)
    } else {
      cat("No validation data available")
    }
  })
  
  output$metadata <- renderPrint({
    req(values$result$success)
    cat("=== Annotation Metadata ===\n\n")
    for (n in names(values$result$metadata)) {
      val <- values$result$metadata[[n]]
      if (is.numeric(val)) val <- round(val, 4)
      cat(sprintf("%s: %s\n", n, val))
    }
  })
  
  # Download handlers
  output$download_csv <- downloadHandler(
    filename = function() paste0("annotation_", Sys.Date(), ".csv"),
    content = function(file) write.csv(values$result$annotations, file, row.names = FALSE)
  )
  
  output$download_xlsx <- downloadHandler(
    filename = function() paste0("annotation_", Sys.Date(), ".xlsx"),
    content = function(file) {
      req(values$result$success)
      if (requireNamespace("openxlsx", quietly = TRUE)) {
        wb <- openxlsx::createWorkbook()
        openxlsx::addWorksheet(wb, "Annotations")
        openxlsx::writeData(wb, "Annotations", values$result$annotations)
        openxlsx::addWorksheet(wb, "Metadata")
        metadata_df <- data.frame(Parameter = names(values$result$metadata),
                                  Value = as.character(unlist(values$result$metadata)))
        openxlsx::writeData(wb, "Metadata", metadata_df)
        if (!is.null(values$result$validation)) {
          openxlsx::addWorksheet(wb, "Validation")
          openxlsx::writeData(wb, "Validation", values$result$validation$summary)
        }
        openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
      } else {
        write.csv(values$result$annotations, file, row.names = FALSE)
      }
    }
  )
  
  output$download_report <- downloadHandler(
    filename = function() paste0("annotation_report_", Sys.Date(), ".html"),
    content = function(file) {
      req(values$result$success)
      m <- values$result$metadata
      summary_stats <- data.frame(
        Metric = c("Tissue", "Species", "Model", "Number of Clusters",
                   "Total Runtime (sec)", "API Latency (sec)", "Tokens Used",
                   "Estimated Cost (USD)", "Mean Confidence"),
        Value = c(
          m$tissue, m$species, m$model, m$n_clusters,
          sprintf("%.2f", m$total_runtime_sec),
          sprintf("%.2f", m$api_latency_sec),
          m$tokens_used,
          sprintf("$%.4f", m$estimated_cost_usd),
          sprintf("%.3f", mean(values$result$annotations$Confidence, na.rm = TRUE))
        )
      )
      html_content <- sprintf(
        '<!DOCTYPE html><html><head><title>DeepSeekCell Report</title>
        <style>body{font-family:Arial;margin:40px}h1{color:#2c3e50}table{border-collapse:collapse;width:100%%}th,td{border:1px solid #ddd;padding:8px}th{background:#3498db;color:#fff}</style>
        </head><body><h1>🧬 DeepSeekCell Annotation Report</h1><p>%s</p>
        <h2>Summary</h2><div>%s</div><h2>Annotations</h2>%s</body></html>',
        Sys.time(),
        paste(sprintf("<b>%s:</b> %s<br>", summary_stats$Metric, summary_stats$Value), collapse = "\n"),
        knitr::kable(values$result$annotations, format = "html")
      )
      writeLines(html_content, file)
    }
  )
}

shinyApp(ui, server)