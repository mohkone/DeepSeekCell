# R/ontology.R

#' Cell Ontology Mapping Functions
#'
#' Provides robust mapping between predicted cell types and Cell Ontology
#' identifiers, with synonym normalization, conservative fuzzy matching,
#' and lineage-based evaluation.
#'
#' @import ontologyIndex stringdist logger


# -------------------------------------------------------------------------
# Cell type synonym normalization dictionary
# -------------------------------------------------------------------------

CELL_TYPE_SYNONYMS <- c(
  "nk cell" = "natural killer cell",
  "natural killer cell" = "natural killer cell",
  "treg" = "regulatory t cell",
  "regulatory t-cell" = "regulatory t cell",
  "cd4 t cell" = "cd4 positive alpha beta t cell",
  "cd8 t cell" = "cd8 positive alpha beta t cell",
  "cytotoxic t cell" = "cd8 positive alpha beta t cell",
  "b lymphocyte" = "b cell",
  "t lymphocyte" = "t cell",
  "monocytes" = "monocyte",
  "macrophages" = "macrophage",
  "dendritic cells" = "dendritic cell",
  "platelets" = "platelet",
  "mast cell basophil" = "mast cell",
  "basophil mast cell" = "mast cell",
  "basophil/mast cell" = "mast cell",
  "mast cell/basophil" = "mast cell",
  "microglia" = "microglial cell",
  "endothelial" = "endothelial cell",
  "fibroblasts" = "fibroblast",
  "smooth muscle" = "smooth muscle cell",
  "pericytes" = "pericyte"
)

TISSUE_CONTEXT_ALIASES <- list(
  pancreas = list(
    "alpha cell" = c(
      "pancreatic A cell",
      "pancreatic alpha cell",
      "alpha cell of islet of Langerhans"
    ),
    "pancreatic alpha cell" = c(
      "pancreatic A cell",
      "pancreatic alpha cell",
      "alpha cell of islet of Langerhans"
    ),
    "islet alpha cell" = c(
      "pancreatic A cell",
      "pancreatic alpha cell",
      "alpha cell of islet of Langerhans"
    ),
    "beta cell" = c(
      "type B pancreatic cell",
      "pancreatic beta cell",
      "beta cell of islet of Langerhans"
    ),
    "pancreatic beta cell" = c(
      "type B pancreatic cell",
      "pancreatic beta cell",
      "beta cell of islet of Langerhans"
    ),
    "islet beta cell" = c(
      "type B pancreatic cell",
      "pancreatic beta cell",
      "beta cell of islet of Langerhans"
    ),
    "delta cell" = c(
      "pancreatic D cell",
      "pancreatic delta cell",
      "delta cell of pancreatic islet"
    ),
    "pancreatic delta cell" = c(
      "pancreatic D cell",
      "pancreatic delta cell",
      "delta cell of pancreatic islet"
    ),
    "islet delta cell" = c(
      "pancreatic D cell",
      "pancreatic delta cell",
      "delta cell of pancreatic islet"
    ),
    "acinar cell" = c(
      "pancreatic acinar cell",
      "acinar cell of pancreas"
    ),
    "pancreatic acinar cell" = c(
      "pancreatic acinar cell",
      "acinar cell of pancreas"
    ),
    "ductal cell" = c(
      "pancreatic ductal cell",
      "duct epithelial cell of pancreas"
    ),
    "pancreatic ductal cell" = c(
      "pancreatic ductal cell",
      "duct epithelial cell of pancreas"
    )
  ),
  brain = list(
    "neuron" = c(
      "central nervous system neuron",
      "neuron"
    ),
    "neuronal cell" = c(
      "central nervous system neuron",
      "neuron"
    ),
    "glial cell" = "glial cell",
    "glia" = "glial cell",
    "astrocyte" = "astrocyte",
    "astrocytes" = "astrocyte",
    "oligodendrocyte" = "oligodendrocyte",
    "oligodendrocytes" = "oligodendrocyte",
    "microglia" = "microglial cell",
    "microglial cell" = "microglial cell",
    "opc" = "oligodendrocyte precursor cell",
    "opcs" = "oligodendrocyte precursor cell",
    "oligodendrocyte precursor" = "oligodendrocyte precursor cell",
    "oligodendrocyte precursor cell" = "oligodendrocyte precursor cell",
    "ng2 cell" = "oligodendrocyte precursor cell",
    "ng2 glia" = "oligodendrocyte precursor cell",
    "ependymal" = "ependymal cell",
    "ependymal cell" = "ependymal cell",
    "interneuron" = "interneuron",
    "interneurons" = "interneuron",
    "cortical interneuron" = "cortical interneuron",
    "gabaergic interneuron" = "GABAergic interneuron",
    "inhibitory interneuron" = "GABAergic interneuron",
    "inhibitory neuron" = "GABAergic neuron",
    "gabaergic neuron" = "GABAergic neuron",
    "gaba ergic neuron" = "GABAergic neuron",
    "excitatory neuron" = "glutamatergic neuron",
    "glutamatergic neuron" = "glutamatergic neuron",
    "dopaminergic neuron" = "dopaminergic neuron",
    "pyramidal cell" = "pyramidal neuron",
    "pyramidal neuron" = "pyramidal neuron",
    "medium spiny neuron" = "medium spiny neuron",
    "msn" = "medium spiny neuron",
    "neural stem cell" = "neural stem cell",
    "radial glia" = "radial glial cell",
    "radial glial cell" = "radial glial cell",
    "neuroblast" = "neuroblast (sensu Vertebrata)"
  ),
  lung = list(
    "lung epithelial cell" = "epithelial cell of lung",
    "epithelial cell" = "epithelial cell of lung",
    "alveolar epithelial cell" = "pulmonary alveolar epithelial cell",
    "pneumocyte" = "pulmonary alveolar epithelial cell",
    "alveolar type 1 cell" = "pulmonary alveolar type 1 cell",
    "alveolar type i cell" = "pulmonary alveolar type 1 cell",
    "type 1 pneumocyte" = "pulmonary alveolar type 1 cell",
    "type i pneumocyte" = "pulmonary alveolar type 1 cell",
    "at1 cell" = "pulmonary alveolar type 1 cell",
    "at1" = "pulmonary alveolar type 1 cell",
    "ati cell" = "pulmonary alveolar type 1 cell",
    "at i cell" = "pulmonary alveolar type 1 cell",
    "alveolar type 2 cell" = "pulmonary alveolar type 2 cell",
    "alveolar type ii cell" = "pulmonary alveolar type 2 cell",
    "type 2 pneumocyte" = "pulmonary alveolar type 2 cell",
    "type ii pneumocyte" = "pulmonary alveolar type 2 cell",
    "at2 cell" = "pulmonary alveolar type 2 cell",
    "at2" = "pulmonary alveolar type 2 cell",
    "atii cell" = "pulmonary alveolar type 2 cell",
    "at ii cell" = "pulmonary alveolar type 2 cell",
    "club cell" = "club cell",
    "clara cell" = "club cell",
    "ciliated cell" = "lung multiciliated epithelial cell",
    "ciliated epithelial cell" = "lung multiciliated epithelial cell",
    "multiciliated cell" = "lung multiciliated epithelial cell",
    "multiciliated epithelial cell" = "lung multiciliated epithelial cell",
    "tracheobronchial ciliated cell" =
      "multiciliated columnar cell of tracheobronchial tree",
    "basal cell" = "respiratory basal cell",
    "respiratory basal cell" = "respiratory basal cell",
    "basal epithelial cell" = "basal epithelial cell of tracheobronchial tree",
    "airway basal cell" = "basal epithelial cell of tracheobronchial tree",
    "goblet cell" = "lung goblet cell",
    "lung goblet cell" = "lung goblet cell",
    "respiratory goblet cell" = "respiratory tract goblet cell",
    "secretory cell" = "lung secretory cell",
    "lung secretory cell" = "lung secretory cell",
    "airway secretory cell" = "respiratory airway secretory cell",
    "rasc" = "respiratory airway secretory cell",
    "respiratory airway secretory cell" = "respiratory airway secretory cell",
    "ionocyte" = "pulmonary ionocyte",
    "pulmonary ionocyte" = "pulmonary ionocyte",
    "neuroendocrine cell" = "pulmonary neuroendocrine cell",
    "pulmonary neuroendocrine cell" = "pulmonary neuroendocrine cell",
    "tuft cell" = "brush cell of tracheobronchial tree",
    "brush cell" = "brush cell of tracheobronchial tree",
    "bronchial epithelial cell" = "bronchial epithelial cell",
    "endothelial cell" = "lung endothelial cell",
    "lung endothelial cell" = "lung endothelial cell",
    "microvascular endothelial cell" = "lung microvascular endothelial cell",
    "capillary endothelial cell" = "pulmonary capillary endothelial cell",
    "pulmonary capillary endothelial cell" =
      "pulmonary capillary endothelial cell",
    "fibroblast" = "fibroblast of lung",
    "fibroblasts" = "fibroblast of lung",
    "lung fibroblast" = "fibroblast of lung",
    "interstitial fibroblast" = "pulmonary interstitial fibroblast",
    "pulmonary interstitial fibroblast" = "pulmonary interstitial fibroblast",
    "pericyte" = "lung pericyte",
    "lung pericyte" = "lung pericyte",
    "macrophage" = "lung macrophage",
    "lung macrophage" = "lung macrophage",
    "alveolar macrophage" = "alveolar macrophage",
    "interstitial macrophage" = "lung interstitial macrophage",
    "lung interstitial macrophage" = "lung interstitial macrophage",
    "smooth muscle cell" = "tracheobronchial smooth muscle cell",
    "airway smooth muscle cell" = "tracheobronchial smooth muscle cell",
    "bronchial smooth muscle cell" = "bronchial smooth muscle cell"
  )
)

TISSUE_CONTEXT_PATTERNS <- c(
  pancreas = "pancreas|pancreatic|islet|langerhans",
  brain = paste(
    "brain|cerebral|cortex|cortical|hippocampus|hippocampal",
    "cerebellum|cerebellar|striatum|striatal|midbrain",
    "forebrain|hindbrain|spinal cord|central nervous system|\\bcns\\b|neural",
    sep = "|"
  ),
  lung = paste(
    "lung|pulmonary|alveolar|airway|bronchus|bronchial",
    "bronchiole|bronchiolar|tracheobronchial|respiratory tract",
    "pneumocyte",
    sep = "|"
  )
)


#' Map cell type to Cell Ontology ID
#'
#' @param cell_type Character string of predicted cell type.
#' @param ontology List returned by load_cell_ontology().
#' @param tissue Optional tissue context used to disambiguate generic labels.
#' @return Data frame with CellType, CL_ID, OntologyLabel, and MatchMethod.
#' @export
map_to_cell_ontology <- function(cell_type, ontology, tissue = NULL) {
  
  if (.is_unknown_cell_type(cell_type)) {
    return(.create_empty_mapping(cell_type))
  }
  
  if (is.null(ontology) || is.null(ontology$dataframe)) {
    warning("Invalid ontology object supplied.")
    return(.create_empty_mapping(cell_type))
  }
  
  ont_df <- ontology$dataframe
  
  required_cols <- c("cl_id", "name", "name_lower")
  if (!all(required_cols %in% colnames(ont_df))) {
    warning("Ontology dataframe is missing required columns.")
    return(.create_empty_mapping(cell_type))
  }
  
  normalized <- normalize_cell_type(cell_type)
  
  if (!nzchar(normalized)) {
    return(.create_empty_mapping(cell_type))
  }

  # Strategy 0: tissue-aware aliases for biologically ambiguous labels.
  # Example: "alpha cell" in pancreas should be pancreatic A cell, not a
  # retinal ganglion alpha cell synonym.
  context_match <- .find_contextual_match(normalized, ont_df, tissue)

  if (!is.null(context_match)) {
    return(
      .create_mapping(
        cell_type,
        ont_df[context_match$idx, ],
        context_match$method,
        context_match$score
      )
    )
  }
  
  # Strategy 1: exact Cell Ontology name match
  exact_idx <- which(ont_df$name_lower == normalized)
  
  if (length(exact_idx) > 0) {
    return(.create_mapping(cell_type, ont_df[exact_idx[1], ], "exact", 1))
  }
  
  # Strategy 2: synonym match
  syn_match <- .find_synonym_match(normalized, ont_df)
  
  if (!is.null(syn_match)) {
    return(.create_mapping(cell_type, ont_df[syn_match, ], "synonym", 0.95))
  }
  
  # Strategy 3: conservative fuzzy match
  fuzzy_match <- .find_fuzzy_match(normalized, ont_df)
  
  if (!is.null(fuzzy_match)) {
    return(
      .create_mapping(
        cell_type,
        ont_df[fuzzy_match$idx, ],
        paste0("fuzzy_jw_", fuzzy_match$dist),
        round(1 - fuzzy_match$dist, 3)
      )
    )
  }
  
  .create_empty_mapping(cell_type)
}


# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

.is_unknown_cell_type <- function(cell_type) {
  is.null(cell_type) ||
    length(cell_type) == 0 ||
    is.na(cell_type) ||
    !nzchar(trimws(cell_type)) ||
    grepl(
      "^(unknown|unidentified|not determined|undetermined|ambiguous)$",
      trimws(cell_type),
      ignore.case = TRUE
    )
}


.create_mapping <- function(cell_type, row, method, score) {
  data.frame(
    CellType = as.character(cell_type),
    CL_ID = as.character(row$cl_id),
    OntologyLabel = as.character(row$name),
    MatchMethod = as.character(method),
    OntologyMatchScore = as.numeric(score),
    stringsAsFactors = FALSE
  )
}


.create_empty_mapping <- function(cell_type) {
  data.frame(
    CellType = as.character(cell_type %||% NA_character_),
    CL_ID = NA_character_,
    OntologyLabel = NA_character_,
    MatchMethod = "none",
    OntologyMatchScore = NA_real_,
    stringsAsFactors = FALSE
  )
}


.find_contextual_match <- function(normalized, ont_df, tissue) {
  aliases <- .context_aliases_for_tissue(tissue)

  if (length(aliases) == 0) {
    return(NULL)
  }

  alias_names <- normalize_cell_type(names(aliases))
  alias_idx <- match(normalized, alias_names)

  if (is.na(alias_idx)) {
    return(NULL)
  }

  preferred_terms <- normalize_cell_type(aliases[[alias_idx]])
  preferred_terms <- preferred_terms[nzchar(preferred_terms)]

  if (length(preferred_terms) == 0) {
    return(NULL)
  }

  exact_idx <- .find_exact_preferred_match(preferred_terms, ont_df)

  if (!is.null(exact_idx)) {
    return(list(idx = exact_idx, method = "context_exact", score = 1))
  }

  syn_idx <- .find_synonym_match(preferred_terms, ont_df)

  if (!is.null(syn_idx)) {
    return(list(idx = syn_idx, method = "context_synonym", score = 0.98))
  }

  NULL
}


.find_exact_preferred_match <- function(preferred_terms, ont_df) {
  for (term in preferred_terms) {
    idx <- which(ont_df$name_lower == term)

    if (length(idx) > 0) {
      return(idx[1])
    }
  }

  NULL
}


.context_aliases_for_tissue <- function(tissue) {
  if (is.null(tissue) || length(tissue) == 0 || is.na(tissue[1])) {
    return(list())
  }

  tissue_norm <- normalize_cell_type(tissue[1])

  if (!nzchar(tissue_norm)) {
    return(list())
  }

  matches <- names(TISSUE_CONTEXT_PATTERNS)[
    vapply(
      TISSUE_CONTEXT_PATTERNS,
      grepl,
      logical(1),
      x = tissue_norm,
      perl = TRUE
    )
  ]

  if (length(matches) == 0) {
    return(list())
  }

  aliases <- do.call(c, unname(TISSUE_CONTEXT_ALIASES[matches]))
  aliases[!duplicated(names(aliases))]
}


.find_synonym_match <- function(normalized, ont_df) {
  
  if (!"synonyms" %in% colnames(ont_df)) {
    return(NULL)
  }
  
  synonyms <- ont_df$synonyms
  
  if (is.null(synonyms) || all(is.na(synonyms))) {
    return(NULL)
  }
  
  synonym_norm <- lapply(
    synonyms,
    function(s) {
      if (is.na(s) || !nzchar(s)) {
        return(NA_character_)
      }
      
      parts <- .extract_synonym_terms(s)
      parts <- normalize_cell_type(parts)
      parts[nzchar(parts)]
    }
  )

  normalized <- unique(normalize_cell_type(normalized))
  normalized <- normalized[nzchar(normalized)]

  if (length(normalized) == 0) {
    return(NULL)
  }
  
  syn_matches <- which(
    vapply(
      synonym_norm,
      function(s) !all(is.na(s)) && any(normalized %in% s),
      logical(1)
    )
  )
  
  if (length(syn_matches) > 0) {
    return(syn_matches[1])
  }
  
  NULL
}

.extract_synonym_terms <- function(s) {
  quoted <- regmatches(s, gregexpr('"[^"]+"', s, perl = TRUE))[[1]]
  if (length(quoted) > 0 && !identical(quoted, character(0))) {
    return(gsub('"', "", quoted, fixed = TRUE))
  }

  unlist(strsplit(s, ";", fixed = TRUE))
}


.find_fuzzy_match <- function(normalized, ont_df) {
  
  if (is.null(normalized) || is.na(normalized) || nchar(normalized) < 6) {
    return(NULL)
  }
  
  candidates <- ont_df$name_lower
  
  valid <- !is.na(candidates) & nzchar(candidates)
  
  if (!any(valid)) {
    return(NULL)
  }
  
  candidates_valid <- candidates[valid]
  valid_indices <- which(valid)
  
  distances <- stringdist::stringdist(
    normalized,
    candidates_valid,
    method = "jw"
  )
  
  best_local_idx <- which.min(distances)
  best_dist <- distances[best_local_idx]
  best_global_idx <- valid_indices[best_local_idx]
  
  best_label <- candidates[best_global_idx]
  
  # Conservative biological safety rule:
  # Do not fuzzy-map very different short labels, e.g. mast cell -> t cell.
  token_overlap <- .token_overlap(normalized, best_label)
  
  if (best_dist <= 0.10 && token_overlap > 0) {
    return(list(
      idx = best_global_idx,
      dist = round(best_dist, 3)
    ))
  }
  
  NULL
}


.token_overlap <- function(a, b) {
  ta <- unique(unlist(strsplit(a, "\\s+")))
  tb <- unique(unlist(strsplit(b, "\\s+")))
  
  ta <- ta[nzchar(ta)]
  tb <- tb[nzchar(tb)]
  
  if (length(ta) == 0 || length(tb) == 0) {
    return(0)
  }
  
  length(intersect(ta, tb)) / length(union(ta, tb))
}


#' Normalize cell type string for ontology matching
#'
#' @param x Character vector of cell type names.
#' @return Normalized character vector.
#' @export
normalize_cell_type <- function(x) {
  
  if (is.null(x)) {
    return(character())
  }
  
  x <- as.character(x)
  x <- tolower(trimws(x))
  
  x <- gsub("\\(.*?\\)", "", x)
  x <- gsub("\\b(likely|probable|possible|putative)\\b", "", x)
  x <- gsub("\\+", " positive ", x)
  x <- gsub("/", " ", x)
  x <- gsub("-", " ", x)
  x <- gsub("[^a-z0-9\\s]", " ", x)
  x <- gsub("\\s+", " ", x)
  x <- trimws(x)
  
  x <- vapply(
    x,
    function(xx) {
      if (xx %in% names(CELL_TYPE_SYNONYMS)) {
        return(unname(CELL_TYPE_SYNONYMS[[xx]]))
      }
      xx
    },
    character(1)
  )
  
  unname(x)
}


# -------------------------------------------------------------------------
# Lineage-based ontology evaluation
# -------------------------------------------------------------------------

#' Check whether two Cell Ontology IDs are lineage-related
#'
#' @param cl_id_1 First CL identifier.
#' @param cl_id_2 Second CL identifier.
#' @param ontology Ontology object from load_cell_ontology().
#' @return Logical.
#' @export
are_related_terms <- function(cl_id_1, cl_id_2, ontology) {
  
  if (is.na(cl_id_1) || is.na(cl_id_2)) {
    return(FALSE)
  }
  
  if (identical(cl_id_1, cl_id_2)) {
    return(TRUE)
  }
  
  if (is.null(ontology$ancestor_cache)) {
    return(FALSE)
  }
  
  ancestors_1 <- ontology$ancestor_cache[[cl_id_1]]
  ancestors_2 <- ontology$ancestor_cache[[cl_id_2]]
  
  if (is.null(ancestors_1) || is.null(ancestors_2)) {
    return(FALSE)
  }
  
  cl_id_2 %in% ancestors_1 || cl_id_1 %in% ancestors_2
}


#' Calculate clade accuracy
#'
#' @param predicted_cl Vector of predicted CL IDs.
#' @param ground_truth_cl Vector of ground-truth CL IDs.
#' @param ontology Ontology object from load_cell_ontology().
#' @return Proportion of lineage-related predictions.
#' @export
calculate_clade_accuracy <- function(predicted_cl, ground_truth_cl, ontology) {
  
  if (length(predicted_cl) != length(ground_truth_cl)) {
    stop("predicted_cl and ground_truth_cl must have the same length.", call. = FALSE)
  }
  
  valid_idx <- !is.na(predicted_cl) & !is.na(ground_truth_cl)
  
  if (sum(valid_idx) == 0) {
    return(NA_real_)
  }
  
  related <- mapply(
    are_related_terms,
    predicted_cl[valid_idx],
    ground_truth_cl[valid_idx],
    MoreArgs = list(ontology = ontology)
  )
  
  mean(related)
}
