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

ui <- fluidPage(
  theme = shinytheme("flatly"),
  tags$head(tags$style(HTML("
    .live-badge { background-color: #27ae60; color: white; padding: 5px 12px;
                  border-radius: 20px; font-size: 14px; font-weight: bold; }
    .confidence-high { color: #27ae60; font-weight: bold; }
    .confidence-mid { color: #f39c12; }
    .confidence-low { color: #e74c3c; }
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
                 actionButton("run", "Annotate", class = "btn-primary btn-lg", width = "100%"),
                 hr(),
                 actionButton("example_pbmc", "Load PBMC Example", class = "btn-info"),
                 actionButton("clear", "Clear", class = "btn-warning")
    ),
    mainPanel(width = 8,
              tabsetPanel(
                tabPanel("Results", DTOutput("results_table") %>% withSpinner(),
                         downloadButton("download_csv", "CSV")),
                tabPanel("Confidence", plotlyOutput("confidence_plot") %>% withSpinner()),
                tabPanel("Performance", verbatimTextOutput("metadata"))
              )
    )
  )
)

server <- function(input, output, session) {
  values <- reactiveValues(result = NULL)
  
  observeEvent(input$example_pbmc, {
    updateTextAreaInput(session, "c1", value = "CD3D, CD3E, CD8A, CD4")
    updateTextAreaInput(session, "c2", value = "CD14, LYZ, FCGR3A, MS4A7")
    updateTextAreaInput(session, "c3", value = "CD79A, MS4A1, CD19")
    updateTextInput(session, "tissue", value = "PBMC")
  })
  
  observeEvent(input$clear, {
    for (i in 1:3) updateTextAreaInput(session, paste0("c", i), value = "")
  })
  
  observeEvent(input$run, {
    markers <- list()
    for (i in 1:3) {
      txt <- input[[paste0("c", i)]]
      if (nzchar(txt)) markers[[paste0("Cluster", i)]] <- trimws(strsplit(txt, ",")[[1]])
    }
    if (length(markers) == 0) {
      showNotification("Enter at least one cluster.", type = "error")
      return()
    }
    
    showNotification("Calling API...", duration = NULL, id = "status")
    res <- tryCatch(
      annotate_cell_types(markers, input$tissue, input$species,
                          input$model, input$api_key, input$use_ontology),
      error = function(e) { list(success = FALSE, error = e$message) }
    )
    removeNotification(id = "status")
    values$result <- res
  })
  
  output$results_table <- renderDT({
    req(values$result$success)
    datatable(values$result$annotations, options = list(pageLength = 5))
  })
  
  output$confidence_plot <- renderPlotly({
    req(values$result$success)
    df <- values$result$annotations
    df$Cluster <- factor(df$Cluster, levels = rev(unique(df$Cluster)))
    p <- ggplot(df, aes(x = Cluster, y = Confidence, fill = Confidence)) +
      geom_col() + coord_flip(ylim = c(0,1)) + theme_minimal()
    ggplotly(p)
  })
  
  output$metadata <- renderPrint({
    req(values$result$success)
    str(values$result$metadata)
  })
  
  output$download_csv <- downloadHandler(
    filename = "annotations.csv",
    content = function(file) write.csv(values$result$annotations, file, row.names = FALSE)
  )
}

shinyApp(ui, server)