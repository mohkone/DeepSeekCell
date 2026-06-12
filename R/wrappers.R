#' User-Friendly Wrappers
#' 
#' Simplified interfaces for common use cases.

#' Quick annotation with minimal arguments
#' 
#' @param markers Character vector or list of marker genes
#' @param tissue Tissue name
#' @param api_key API key (optional, uses env var if not provided)
#' @param model Model to use ("deepseek" or "ollama")
#' @param species Species name.
#' @return Annotation result object
#' @export
quick_annotate <- function(markers, tissue, api_key = NULL, model = "deepseek", species = "Human") {
  
  # Convert single vector to list
  if (is.character(markers) && !is.list(markers)) {
    markers <- list(Cluster1 = markers)
  }
  
  model_config <- get_model_config(model)
  api_key <- resolve_api_key(model_config, api_key)

  if (isTRUE(model_config$requires_api_key) && (is.null(api_key) || !nzchar(api_key))) {
    env_hint <- paste(model_config$api_key_env %||% character(), collapse = ", ")
    stop(
      "API key not found. Provide api_key or set one of: ",
      env_hint,
      call. = FALSE
    )
  }
  
  annotate_cell_types(
    markers = markers,
    tissue = tissue,
    species = species,
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
  quick_annotate(markers_list, tissue, api_key, model, species = "Human")
}

#' Batch annotate multiple samples
#' 
#' @param samples List of samples, each with markers and tissue
#' @param api_key API key
#' @param model Model to use
#' @return List of annotation results
#' @export
batch_annotate <- function(samples, api_key = NULL, model = "deepseek") {
  if (!is.list(samples) || length(samples) == 0) {
    stop("samples must be a non-empty list.", call. = FALSE)
  }
  
  results <- list()
  
  sample_names <- names(samples) %||% paste0("sample", seq_along(samples))

  for (sample_name in sample_names) {
    cat("Processing sample:", sample_name, "\n")
    
    sample_data <- samples[[sample_name]]
    if (is.null(sample_data)) {
      sample_data <- samples[[which(sample_names == sample_name)[1]]]
    }
    
    result <- quick_annotate(
      markers = sample_data$markers,
      tissue = sample_data$tissue %||% "Unknown",
      api_key = api_key,
      model = model,
      species = sample_data$species %||% "Human"
    )
    
    results[[sample_name]] <- result
  }
  
  results
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
