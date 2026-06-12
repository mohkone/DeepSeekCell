# R/prompts.R
# Enhanced prompt engineering for LLM-based cell type annotation

#' Create annotation prompt from marker genes
#'
#' This prompt uses marker-based annotation,
#' extended with confidence scoring, reasoning, contamination awareness,
#' mixture detection, and publication-ready JSON output.
#'
#' @param markers Named list of marker genes per cluster.
#' @param tissue Tissue name.
#' @param species Species name, e.g. "Human", "Mouse", "Rat".
#' @param include_reasoning Whether to request concise marker-based reasoning.
#' @return Character prompt.
#' @export
create_annotation_prompt <- function(markers,
                                     tissue,
                                     species,
                                     include_reasoning = TRUE) {
  
  if (!is.list(markers) || length(markers) == 0) {
    stop("markers must be a non-empty named list.", call. = FALSE)
  }
  
  if (is.null(names(markers)) || any(names(markers) == "")) {
    names(markers) <- paste0("Cluster", seq_along(markers))
  }
  
  cluster_text <- .format_cluster_markers(markers)
  tissue_key <- trimws(as.character(tissue))
  
  tissue_prior <- .get_tissue_prior(tissue_key)
  
  reasoning_rule <- if (isTRUE(include_reasoning)) {
    '"reasoning": "brief marker-based explanation, including whether the cell type is expected or unexpected in the stated tissue"'
  } else {
    '"reasoning": ""'
  }
  
  sprintf(
    'Species: %s
Tissue: %s

You are an expert single-cell transcriptomics annotator.

Your task:
For each cluster, infer the most likely biological cell type from the marker genes.

Important rules:
1. Identify the most likely biological cell type first, even if it is unexpected for the stated tissue.
2. Do NOT label a cluster as Unknown only because it is unexpected in the tissue.
3. If markers suggest contamination, ambient RNA, doublets, or non-native cells, still report the most likely biological cell type and mention this in the reasoning.
4. Use Unknown only when no biologically coherent cell type can be inferred from the marker genes.
5. Use standardized, publication-ready cell type names.
6. Prefer the most specific biologically justified label.
7. If markers indicate more than one cell type, set is_mixed = true and describe the mixture.
8. Confidence must be numeric between 0 and 1.
9. Return only valid JSON. Do not include markdown fences or extra text.

Tissue prior:
%s

Marker genes:
%s

Required JSON schema:
{
  "annotations": [
    {
      "cluster": "Cluster1",
      "cell_type": "T cell",
      "confidence": 0.95,
      "is_mixed": false,
      "primary_cell_type": "T cell",
      "secondary_cell_type": "",
      "tissue_consistency": "expected",
      %s
    }
  ]
}

Allowed values for tissue_consistency:
- "expected"
- "unexpected"
- "possible_contamination"
- "possible_doublet"
- "unknown"',
    species,
    tissue,
    tissue_prior,
    cluster_text,
    reasoning_rule
  )
}


#' Format cluster marker genes for prompt
#' @keywords internal
.format_cluster_markers <- function(markers, max_genes = 30) {
  
  paste(
    vapply(
      names(markers),
      function(cluster_name) {
        genes <- markers[[cluster_name]]
        genes <- unique(trimws(as.character(genes)))
        genes <- genes[nzchar(genes)]
        genes <- head(genes, max_genes)
        
        sprintf("%s: %s", cluster_name, paste(genes, collapse = ", "))
      },
      character(1)
    ),
    collapse = "\n"
  )
}


#' Tissue-specific biological prior
#' @keywords internal
.get_tissue_prior <- function(tissue) {
  
  tissue_key <- tolower(trimws(tissue))

  switch(
    tissue_key,
    pbmc = paste(
      "Common PBMC cell types include T cells, CD4+ T cells, CD8+ T cells,",
      "B cells, plasma cells, NK cells, monocytes, dendritic cells, platelets,",
      "basophils, and hematopoietic progenitors.",
      "However, if markers clearly indicate non-PBMC cells, annotate the true biological cell type",
      "and mark tissue_consistency as possible_contamination or unexpected."
    ),
    
    pancreas = paste(
      "Common pancreatic cell types include acinar cells, ductal cells, alpha cells, beta cells,",
      "delta cells, pancreatic stellate cells, endothelial cells, macrophages, and other immune cells.",
      "If markers indicate non-pancreatic cells, annotate the true biological cell type",
      "and mark tissue_consistency as possible_contamination or unexpected."
    ),
    
    brain = paste(
      "Common brain cell types include neurons, excitatory neurons, inhibitory neurons, astrocytes,",
      "oligodendrocytes, microglia, endothelial cells, pericytes, ependymal cells,",
      "and oligodendrocyte precursor cells.",
      "If markers indicate non-brain cells, annotate the true biological cell type",
      "and mark tissue_consistency as possible_contamination or unexpected."
    ),
    
    prostate = paste(
      "Common prostate cell types include luminal epithelial cells, basal epithelial cells,",
      "fibroblasts, smooth muscle cells, endothelial cells, lymphatic endothelial cells,",
      "pericytes, immune cells, mast cells, macrophages, and T cells."
    ),
    
    paste0(
      "Use the tissue as biological context, but do not force all annotations to belong to ",
      tissue,
      ". Always identify the most likely biological cell type from the markers first."
    )
  )
}


#' Create batch annotation prompt for multiple samples
#'
#' @param samples List of samples, each containing marker lists.
#' @return Batch prompt string.
#' @export
create_batch_prompt <- function(samples) {
  
  if (!is.list(samples) || length(samples) == 0) {
    stop("samples must be a non-empty list.", call. = FALSE)
  }
  
  batch_text <- paste(
    vapply(
      names(samples),
      function(sample_name) {
        markers <- samples[[sample_name]]
        cluster_text <- .format_cluster_markers(markers)
        sprintf("Sample: %s\n%s", sample_name, cluster_text)
      },
      character(1)
    ),
    collapse = "\n\n---\n\n"
  )
  
  paste(
    "You are annotating multiple single-cell RNA-seq samples.",
    "Process each sample independently.",
    "Identify the most likely biological cell type for each cluster.",
    "Report possible contamination, doublets, or unexpected tissue origin when relevant.",
    "",
    batch_text,
    "",
    "Return only valid JSON with one top-level field named 'samples'.",
    sep = "\n"
  )
}
