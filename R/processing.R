#' Input Processing and Normalization
#' 
#' Handles various input formats and prepares marker genes for annotation.

#' Process cell data from various input formats
#' 
#' @param input List, data.frame, or character vector of markers
#' @param top_genes Number of top genes to use per cluster (default 15)
#' @return Processed markers list
#' @export
process_cell_data <- function(input, top_genes = 15) {
  
  markers <- .parse_input(input, top_genes)
  
  if (length(markers) == 0) {
    stop("No valid marker genes found in input")
  }
  
  # Filter low-quality markers
  markers <- lapply(markers, .filter_markers)
  
  # Remove empty clusters
  markers <- markers[lengths(markers) > 0]
  
  message("Processed {length(markers)} clusters with average ",
           "{round(mean(sapply(markers, length)))} markers per cluster")
  
  list(
    markers = markers,
    n_clusters = length(markers),
    total_markers = sum(sapply(markers, length))
  )
}

.parse_input <- function(input, top_genes) {
  
  if (is.list(input)) {
    # Named list of marker vectors
    markers <- lapply(input, function(x) {
      if (is.character(x)) {
        genes <- unlist(strsplit(x, ","))
        genes <- trimws(genes)
        return(head(genes[genes != ""], top_genes))
      }
      return(head(x, top_genes))
    })
    
  } else if (is.data.frame(input)) {
    # Data frame with 'cluster' and 'gene' columns
    if (!all(c("cluster", "gene") %in% colnames(input))) {
      stop("Data frame must contain 'cluster' and 'gene' columns")
    }
    
    markers <- split(input$gene, input$cluster)
    markers <- lapply(markers, function(x) head(unique(x), top_genes))
    
  } else if (is.character(input)) {
    # Single vector of markers
    markers <- list(Cluster1 = head(trimws(unlist(strsplit(input, ","))), top_genes))
    
  } else {
    stop("Unsupported input type. Use list, data.frame, or character vector.")
  }
  
  return(markers)
}

.filter_markers <- function(genes) {
  # Remove common problematic genes
  problematic_patterns <- c("^MT-", "^MTRNR", "^RP[LS]", "^MALAT", 
                            "^XIST", "^RPS", "^RPL", "^HB[AB]", "^HBA")
  
  for (pattern in problematic_patterns) {
    genes <- genes[!grepl(pattern, genes, ignore.case = TRUE)]
  }
  
  # Remove duplicates while preserving order
  genes <- unique(genes)
  
  # Remove empty strings
  genes <- genes[genes != ""]
  
  return(genes)
}

#' Convert Seurat object markers to DeepSeekCell format
#' 
#' @param seurat_markers Output from Seurat::FindAllMarkers()
#' @param top_n Number of top markers per cluster
#' @return List formatted for annotate_cell_types()
#' @export
seurat_markers_to_list <- function(seurat_markers, top_n = 15) {
  
  if (!all(c("cluster", "gene") %in% colnames(seurat_markers))) {
    stop("Seurat markers must have 'cluster' and 'gene' columns")
  }
  
  # Get top markers by avg_log2FC
  seurat_markers <- seurat_markers[order(seurat_markers$avg_log2FC, decreasing = TRUE), ]
  
  markers_list <- split(seurat_markers$gene, seurat_markers$cluster)
  markers_list <- lapply(markers_list, head, top_n)
  
  return(markers_list)
}