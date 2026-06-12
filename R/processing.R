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
  if (!is.numeric(top_genes) || length(top_genes) != 1 || is.na(top_genes) || top_genes < 1) {
    stop("top_genes must be a positive numeric scalar.", call. = FALSE)
  }
  
  markers <- .parse_input(input, top_genes)
  
  if (length(markers) == 0) {
    stop("No valid marker genes found in input")
  }
  
  # Filter low-quality markers
  markers <- lapply(markers, .filter_markers)
  
  # Remove empty clusters
  markers <- markers[lengths(markers) > 0]

  if (length(markers) == 0) {
    stop("No valid marker genes remained after filtering.", call. = FALSE)
  }
  
  message(
    "Processed ", length(markers), " clusters with average ",
    round(mean(lengths(markers))), " markers per cluster"
  )
  
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
      genes <- if (is.character(x) && length(x) == 1) {
        split_marker_text(x)
      } else {
        unique(trimws(as.character(x)))
      }

      head(genes[!is.na(genes) & nzchar(genes)], top_genes)
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
    markers <- list(Cluster1 = head(split_marker_text(input), top_genes))
    
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
  
  effect_col <- intersect(c("avg_log2FC", "avg_logFC"), colnames(seurat_markers))[1]
  if (is.na(effect_col)) {
    stop("Seurat markers must contain avg_log2FC or avg_logFC.", call. = FALSE)
  }

  if ("p_val_adj" %in% colnames(seurat_markers)) {
    seurat_markers <- seurat_markers[is.na(seurat_markers$p_val_adj) | seurat_markers$p_val_adj < 0.05, ]
  }

  seurat_markers <- seurat_markers[order(seurat_markers$cluster, -seurat_markers[[effect_col]]), ]
  
  markers_list <- split(seurat_markers$gene, seurat_markers$cluster)
  markers_list <- lapply(markers_list, head, top_n)
  
  return(markers_list)
}
