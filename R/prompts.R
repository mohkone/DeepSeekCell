#' Prompt Engineering for LLM Annotation
#' 
#' Constructs optimized prompts for cell type annotation tasks.

#' Create annotation prompt from marker genes
#' 
#' @param markers Named list of marker genes per cluster
#' @param tissue Tissue name
#' @param species Species (Human/Mouse/Rat)
#' @param include_reasoning Whether to request reasoning in output
#' @return Formatted prompt string
#' @export
create_annotation_prompt <- function(markers, tissue, species, 
                                     include_reasoning = TRUE) {
  
  cluster_text <- .format_cluster_markers(markers)
  
  prompt_parts <- c(
    .create_header(tissue, species),
    "",
    "Marker genes per cluster:",
    cluster_text,
    "",
    .create_instructions(include_reasoning)
  )
  
  paste(prompt_parts, collapse = "\n")
}

.format_cluster_markers <- function(markers) {
  paste(
    sapply(names(markers), function(cluster_name) {
      markers_str <- paste(markers[[cluster_name]][1:min(15, length(markers[[cluster_name]]))], 
                           collapse = ", ")
      sprintf("  %s: %s", cluster_name, markers_str)
    }),
    collapse = "\n"
  )
}

.create_header <- function(tissue, species) {
  sprintf("Tissue: %s\nSpecies: %s", tissue, species)
}

.create_instructions <- function(include_reasoning) {
  if (include_reasoning) {
    paste(
      "Please annotate each cluster with a cell type and confidence score.",
      "For each annotation, provide brief reasoning based on the marker genes.",
      "",
      "Return ONLY valid JSON in this format:",
      '{"annotations": [',
      '  {"cluster": "Cluster1", "cell_type": "T cell", "confidence": 0.95, "reasoning": "CD3D and CD3E are T cell markers"},',
      '  ...',
      ']}',
      sep = "\n"
    )
  } else {
    paste(
      "Please annotate each cluster with a cell type and confidence score.",
      "",
      "Return ONLY valid JSON in this format:",
      '{"annotations": [',
      '  {"cluster": "Cluster1", "cell_type": "T cell", "confidence": 0.95},',
      '  ...',
      ']}',
      sep = "\n"
    )
  }
}

#' Create batch annotation prompt for multiple samples
#' 
#' @param samples List of sample annotations with markers
#' @return Batch prompt string
#' @export
create_batch_prompt <- function(samples) {
  
  batch_text <- paste(
    sapply(names(samples), function(sample_name) {
      markers <- samples[[sample_name]]
      cluster_text <- .format_cluster_markers(markers)
      sprintf("Sample: %s\n%s", sample_name, cluster_text)
    }),
    collapse = "\n\n---\n\n"
  )
  
  paste(
    "You are annotating multiple samples. Process each independently.",
    "",
    batch_text,
    "",
    "For each sample, provide annotations for all clusters.",
    "Return JSON with 'sample' as top-level key.",
    sep = "\n"
  )
}