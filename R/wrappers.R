#' User-Friendly Wrappers
#' 
#' Simplified interfaces for common use cases.

#' Quick annotation with minimal arguments
#' 
#' @param markers Character vector or list of marker genes
#' @param tissue Tissue name
#' @param api_key API key (optional, uses env var if not provided)
#' @param model Model to use ("deepseek" or "gpt4")
#' @return Annotation result object
#' @export
quick_annotate <- function(markers, tissue, api_key = NULL, model = "deepseek") {
  
  # Convert single vector to list
  if (is.character(markers) && !is.list(markers)) {
    markers <- list(Cluster1 = markers)
  }
  
  # Get API key from environment if not provided
  if (is.null(api_key)) {
    api_key <- Sys.getenv(paste0(toupper(model), "_API_KEY"))
    if (api_key == "") {
      stop("API key not found. Provide api_key or set ", 
           toupper(model), "_API_KEY environment variable")
    }
  }
  
  annotate_cell_types(
    markers = markers,
    tissue = tissue,
    model_name = model,
    api_key = api_key,
    use_ontology = TRUE,
    validate = TRUE
  )
}

#' Annotate from Seurat object
#' 
#' @param seurat_obj Seurat object with clusters
#' @param tissue Tissue name
#' @param api_key API key
#' @param model Model to use
#' @param top_markers Number of top markers per cluster
#' @return Annotation result object
#' @export
annotate_seurat <- function(seurat_obj, tissue, api_key = NULL, 
                            model = "deepseek", top_markers = 15) {
  
  # Check Seurat is installed
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package is required for this function")
  }
  
  # Find all markers
  markers <- Seurat::FindAllMarkers(seurat_obj, only.pos = TRUE)
  
  # Convert to list format
  markers_list <- seurat_markers_to_list(markers, top_n = top_markers)
  
  # Annotate
  quick_annotate(markers_list, tissue, api_key, model)
}

#' Batch annotate multiple samples
#' 
#' @param samples List of samples, each with markers and tissue
#' @param api_key API key
#' @param model Model to use
#' @return List of annotation results
#' @export
batch_annotate <- function(samples, api_key = NULL, model = "deepseek") {
  
  results <- list()
  
  for (sample_name in names(samples)) {
    cat("Processing sample:", sample_name, "\n")
    
    sample_data <- samples[[sample_name]]
    
    result <- quick_annotate(
      markers = sample_data$markers,
      tissue = sample_data$tissue %||% "Unknown",
      api_key = api_key,
      model = model
    )
    
    results[[sample_name]] <- result
  }
  
  return(results)
}

#' Run annotation with progress tracking
#' 
#' @param markers List of markers
#' @param tissue Tissue name
#' @param api_key API key
#' @param model Model to use
#' @return Annotation result
#' @export
annotate_with_progress <- function(markers, tissue, api_key = NULL, model = "deepseek") {
  
  if (requireNamespace("progress", quietly = TRUE)) {
    pb <- progress::progress_bar$new(
      format = "  Annotating [:bar] :percent in :elapsed",
      total = length(markers),
      clear = FALSE
    )
    
    pb$tick(0)
    result <- quick_annotate(markers, tissue, api_key, model)
    pb$tick(length(markers))
  } else {
    result <- quick_annotate(markers, tissue, api_key, model)
  }
  
  return(result)
}