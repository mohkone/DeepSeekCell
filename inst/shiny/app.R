#!/usr/bin/env Rscript
# DeepSeekCell Shiny Application - Publication Version
suppressPackageStartupMessages({
  library(shiny)
  library(shinythemes)
  library(shinycssloaders)
  library(DT)
  library(ggplot2)
  library(plotly)
})

source_deepseekcell_functions <- function() {
  candidate_dirs <- c(
    file.path("..", "..", "R"),
    "R",
    normalizePath(file.path(getwd(), "..", "..", "R"), mustWork = FALSE)
  )

  r_dir <- candidate_dirs[dir.exists(candidate_dirs)][1]
  if (is.na(r_dir)) {
    stop("Could not locate the package R/ directory.", call. = FALSE)
  }

  r_files <- list.files(r_dir, pattern = "\\.R$", full.names = TRUE)
  r_files <- r_files[!grepl("benchmark|^main\\.R$", basename(r_files))]

  priority <- file.path(r_dir, c("utils.R", "api.R"))
  ordered_files <- c(priority[file.exists(priority)], setdiff(r_files, priority))
  invisible(lapply(ordered_files, source, local = FALSE))
}

source_deepseekcell_functions()

format_display_label <- function(x) {
  x <- gsub("_", " ", trimws(as.character(x)), fixed = TRUE)
  x <- tolower(x)
  tools::toTitleCase(x)
}

format_tissue_badge <- function(x) {
  x <- tolower(trimws(as.character(x)))

  css_class <- if (x == "expected") {
    "badge-expected"
  } else if (x %in% c("possible_contamination", "unexpected")) {
    "badge-contamination"
  } else if (x %in% c("possible_doublet", "doublet")) {
    "badge-doublet"
  } else {
    "badge-unknown"
  }

  sprintf('<span class="%s">%s</span>', css_class, format_display_label(x))
}

format_mixed_badge <- function(x) {
  ifelse(
    as.logical(x),
    '<span class="badge-doublet">Mixed</span>',
    '<span class="badge-expected">Single</span>'
  )
}

format_column_names <- function(x) {
  labels <- c(
    Cluster = "Cluster",
    CellType = "Cell Type",
    Confidence = "Confidence",
    TissueConsistency = "Tissue Consistency",
    IsMixed = "Mixed",
    CL_ID = "CL ID",
    OntologyLabel = "Ontology Label",
    MatchMethod = "Match Method",
    OntologyMatchScore = "Ontology Match Score",
    Markers = "Marker Genes",
    Reasoning = "Reasoning"
  )

  unname(ifelse(x %in% names(labels), labels[x], x))
}

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
    .badge-expected {
  background-color: #27ae60;
  color: white;
  padding: 4px 8px;
  border-radius: 12px;
  font-weight: bold;
}

.badge-contamination {
  background-color: #e67e22;
  color: white;
  padding: 4px 8px;
  border-radius: 12px;
  font-weight: bold;
}

.badge-doublet {
  background-color: #8e44ad;
  color: white;
  padding: 4px 8px;
  border-radius: 12px;
  font-weight: bold;
}

.badge-unknown {
  background-color: #7f8c8d;
  color: white;
  padding: 4px 8px;
  border-radius: 12px;
  font-weight: bold;
}
    .reasoning-cell {
      max-width: 420px;
      white-space: normal;
      line-height: 1.35;
    }
    .ontology-link {
      font-family: monospace;
      white-space: nowrap;
    }
    .marker-cell,
    .reasoning-cell {
      max-width: 520px;
      white-space: normal;
      line-height: 1.35;
    }
    table.dataTable td {
      vertical-align: middle;
    }
  "))),
  div(class = "text-center", style = "background-color: #2c3e50; color: white; padding: 20px;",
      h1("DeepSeekCell", span(class = "live-badge", "LIVE MODE")),
      h4("Ontology-Guided LLM Framework for Explainable Cell Type Annotation"),
      p(paste("Version", deepseekcell_version(), "| For research use only"))),
  sidebarLayout(
    sidebarPanel(width = 4,
                 passwordInput("api_key", "API Key", placeholder = "sk-..."),
                 selectInput("model", "Model", choices = c("DeepSeek" = "deepseek", "Ollama (local)" = "ollama")),
                 # In UI, after selectInput, add:
                 conditionalPanel(
                   condition = "input.model == 'ollama'",
                   helpText("Ollama must be running locally (http://localhost:11434). No API key needed.")
                 ),
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
                         withSpinner(DTOutput("results_table")),
                         br(),
                         fluidRow(
                           column(4, downloadButton("download_csv", "CSV", class = "btn-success", width = "100%")),
                           column(4, downloadButton("download_xlsx", "Excel", class = "btn-success", width = "100%")),
                           column(4, downloadButton("download_report", "HTML Report", class = "btn-info", width = "100%"))
                         )),
                tabPanel("Explainability", br(), withSpinner(DTOutput("explainability_table"))),
                tabPanel("Confidence", br(), withSpinner(plotlyOutput("confidence_plot"))),
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
                                tags$li("Enter your DeepSeek API key, or choose Ollama for local annotation."),
                               tags$li("Specify tissue and species."),
                               tags$li("Input marker genes for up to 5 clusters (comma-, semicolon-, or newline-separated)."),
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
                             p(tags$em("An Ontology-Guided Large Language Model Framework and Interactive Shiny Platform for Explainable Cell Type Annotation in Single-Cell RNA Sequencing."))
                         ))
              )
    )
  )
)

server <- function(input, output, session) {
  values <- reactiveValues(result = NULL)

  collect_markers <- reactive({
    markers <- list()
    for (i in 1:5) {
      txt <- input[[paste0("c", i)]]
      if (!is.null(txt) && nchar(trimws(txt)) > 0) {
        genes <- split_marker_text(txt)
        if (length(genes) > 0) markers[[paste0("Cluster", i)]] <- genes
      }
    }
    if (length(markers) == 0) {
      showNotification("Enter at least one cluster.", type = "error")
      return(NULL)
    }
    markers
  })
  
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
  
  observeEvent(input$run, {
    markers <- collect_markers()
    if (is.null(markers)) return()
    
    model_config <- get_model_config(input$model)
    api_key_to_use <- resolve_api_key(model_config, input$api_key)

    if (isTRUE(model_config$requires_api_key) && is.null(api_key_to_use)) {
      showNotification("API key is required for this model.", type = "error")
      return()
    }
    
    showNotification("Calling API... This may take 15-30 seconds.",
                     type = "message", duration = NULL, id = "status")
    
    res <- tryCatch(
      annotate_cell_types(
        markers = markers,
        tissue = input$tissue,
        species = input$species,
        model_name = input$model,
        api_key = api_key_to_use,
        use_ontology = input$use_ontology,
        validate = TRUE
      ),
      error = function(e) { list(success = FALSE, error = e$message) }
    )
    removeNotification(id = "status")
    values$result <- res
    if (isTRUE(res$success)) {
      showNotification("Annotation completed successfully.", type = "message")
    } else {
      showNotification(paste("Annotation failed:", res$error), type = "error")
    }
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

    text_cols <- names(df)[vapply(df, is.character, logical(1))]
    df[text_cols] <- lapply(df[text_cols], html_escape)
    
    df$Confidence <- sprintf(
      '<span class="%s">%.3f</span>',
      ifelse(
        df$Confidence >= 0.7,
        "confidence-high",
        ifelse(df$Confidence >= 0.4, "confidence-mid", "confidence-low")
      ),
      as.numeric(df$Confidence)
    )
    
    if ("TissueConsistency" %in% names(df)) {
      df$TissueConsistency <- vapply(
        df$TissueConsistency,
        format_tissue_badge,
        character(1)
      )
    }
    
    if ("IsMixed" %in% names(df)) {
      df$IsMixed <- format_mixed_badge(df$IsMixed)
    }

    if ("CL_ID" %in% names(df)) {
      df$CL_ID <- vapply(
        df$CL_ID,
        function(cl_id) {
          cl_id <- trimws(as.character(cl_id))
          if (!nzchar(cl_id)) {
            return("")
          }

          url <- paste0(
            "https://www.ebi.ac.uk/ols4/search?q=",
            utils::URLencode(cl_id, reserved = TRUE),
            "&ontology=cl"
          )

          sprintf(
            '<a class="ontology-link" href="%s" target="_blank" rel="noopener noreferrer">%s</a>',
            url,
            cl_id
          )
        },
        character(1)
      )
    }

    if ("Reasoning" %in% names(df)) {
      df$Reasoning <- ifelse(
        nzchar(trimws(df$Reasoning)),
        sprintf('<div class="reasoning-cell">%s</div>', df$Reasoning),
        ""
      )
    }

    if ("OntologyMatchScore" %in% names(df)) {
      score <- suppressWarnings(as.numeric(df$OntologyMatchScore))
      df$OntologyMatchScore <- ifelse(is.na(score), "", sprintf("%.3f", score))
    }
    
    display_cols <- c("Cluster", "CellType", "Confidence")
    
    if ("TissueConsistency" %in% names(df)) {
      display_cols <- c(display_cols, "TissueConsistency")
    }
    
    if ("IsMixed" %in% names(df)) {
      display_cols <- c(display_cols, "IsMixed")
    }
    
    if ("CL_ID" %in% names(df)) {
      display_cols <- c(display_cols, "CL_ID", "OntologyLabel")
    }

    dt_options <- list(
      pageLength = 10,
      autoWidth = TRUE,
      scrollX = TRUE
    )
    
    DT::datatable(
      df[, display_cols, drop = FALSE],
      escape = FALSE,
      rownames = FALSE,
      colnames = format_column_names(display_cols),
      options = dt_options
    )
  })

  output$explainability_table <- renderDT({
    req(values$result$success)

    df <- values$result$annotations

    if (!is.null(values$result$markers)) {
      marker_text <- vapply(
        df$Cluster,
        function(cluster_name) {
          genes <- values$result$markers[[cluster_name]]
          if (is.null(genes) || length(genes) == 0) {
            return("")
          }
          paste(genes, collapse = ", ")
        },
        character(1)
      )
      df$Markers <- marker_text
    }

    keep_cols <- c(
      "Cluster", "CellType", "TissueConsistency", "IsMixed", "Markers",
      "Reasoning", "CL_ID", "OntologyLabel", "MatchMethod", "OntologyMatchScore"
    )
    keep_cols <- intersect(keep_cols, names(df))
    df <- df[, keep_cols, drop = FALSE]

    text_cols <- names(df)[vapply(df, is.character, logical(1))]
    df[text_cols] <- lapply(df[text_cols], html_escape)

    if ("TissueConsistency" %in% names(df)) {
      df$TissueConsistency <- vapply(
        df$TissueConsistency,
        format_tissue_badge,
        character(1)
      )
    }

    if ("IsMixed" %in% names(df)) {
      df$IsMixed <- format_mixed_badge(df$IsMixed)
    }

    if ("CL_ID" %in% names(df)) {
      df$CL_ID <- vapply(
        df$CL_ID,
        function(cl_id) {
          cl_id <- trimws(as.character(cl_id))
          if (!nzchar(cl_id)) {
            return("")
          }

          url <- paste0(
            "https://www.ebi.ac.uk/ols4/search?q=",
            utils::URLencode(cl_id, reserved = TRUE),
            "&ontology=cl"
          )

          sprintf(
            '<a class="ontology-link" href="%s" target="_blank" rel="noopener noreferrer">%s</a>',
            url,
            cl_id
          )
        },
        character(1)
      )
    }

    if ("Markers" %in% names(df)) {
      df$Markers <- ifelse(
        nzchar(trimws(df$Markers)),
        sprintf('<div class="marker-cell">%s</div>', df$Markers),
        ""
      )
    }

    if ("Reasoning" %in% names(df)) {
      df$Reasoning <- ifelse(
        nzchar(trimws(df$Reasoning)),
        sprintf('<div class="reasoning-cell">%s</div>', df$Reasoning),
        ""
      )
    }

    if ("OntologyMatchScore" %in% names(df)) {
      score <- suppressWarnings(as.numeric(df$OntologyMatchScore))
      df$OntologyMatchScore <- ifelse(is.na(score), "", sprintf("%.3f", score))
    }

    width_targets <- which(names(df) %in% c("Markers", "Reasoning")) - 1
    dt_options <- list(
      pageLength = 5,
      autoWidth = TRUE,
      scrollX = TRUE
    )

    if (length(width_targets) > 0) {
      dt_options$columnDefs <- list(
        list(width = "520px", targets = width_targets)
      )
    }

    DT::datatable(
      df,
      escape = FALSE,
      rownames = FALSE,
      colnames = format_column_names(names(df)),
      options = dt_options
    )
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
      if (v$valid) cat("All checks passed.\n\n") else cat("Issues found:\n", paste("-", v$issues), "\n")
      if (length(v$warnings)) cat("Warnings:\n", paste("-", v$warnings), "\n")
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
    content = function(file) {
      req(values$result$success)
      write.csv(values$result$annotations, file, row.names = FALSE)
    }
  )
  
  output$download_xlsx <- downloadHandler(
    filename = function() paste0("annotation_", Sys.Date(), ".xlsx"),
    content = function(file) {
      req(values$result$success)
      if (!requireNamespace("openxlsx", quietly = TRUE)) {
        stop("Package 'openxlsx' is required for Excel export.", call. = FALSE)
      }

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
    }
  )
  
  output$download_report <- downloadHandler(
    filename = function() paste0("annotation_report_", Sys.Date(), ".html"),
    content = function(file) {
      req(values$result$success)
      generate_html_report(values$result, output_file = file)
    }
  )
}

shinyApp(ui, server)
