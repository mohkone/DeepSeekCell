#' Cell Ontology Mapping Functions
#' 
#' Provides intelligent mapping between predicted cell types and
#' Cell Ontology (CL) identifiers with lineage-based evaluation.
#' 
#' @import ontologyIndex stringdist logger

#' Map cell type to Cell Ontology ID with intelligent matching
#' 
#' @param cell_type Character string of predicted cell type
#' @param ontology List from load_cell_ontology()
#' @return Data frame with CL_ID, OntologyLabel, and MatchMethod
#' @export
map_to_cell_ontology <- function(cell_type, ontology) {
  
  if (is.null(cell_type) || is.na(cell_type) || cell_type == "" || cell_type == "Unknown") {
    return(.create_empty_mapping(cell_type))
  }
  
  normalized <- normalize_cell_type(cell_type)
  ont_df <- ontology$dataframe
  
  # Strategy 1: Exact match
  exact_idx <- which(ont_df$name_lower == normalized)
  if (length(exact_idx) > 0) {
    return(.create_mapping(cell_type, ont_df[exact_idx[1], ], "exact"))
  }
  
  # Strategy 2: Synonym match
  syn_match <- .find_synonym_match(normalized, ont_df)
  if (!is.null(syn_match)) {
    return(.create_mapping(cell_type, ont_df[syn_match, ], "synonym"))
  }
  
  # Strategy 3: Fuzzy matching with adaptive threshold
  fuzzy_match <- .find_fuzzy_match(normalized, ont_df)
  if (!is.null(fuzzy_match)) {
    return(.create_mapping(cell_type, ont_df[fuzzy_match$idx, ], 
                           paste0("fuzzy_", fuzzy_match$dist)))
  }
  
  # No match found
  return(.create_empty_mapping(cell_type))
}

.create_mapping <- function(cell_type, row, method) {
  data.frame(
    CellType = cell_type,
    CL_ID = row$cl_id,
    OntologyLabel = row$name,
    MatchMethod = method,
    stringsAsFactors = FALSE
  )
}

.create_empty_mapping <- function(cell_type) {
  data.frame(
    CellType = cell_type,
    CL_ID = NA_character_,
    OntologyLabel = NA_character_,
    MatchMethod = "none",
    stringsAsFactors = FALSE
  )
}

.find_synonym_match <- function(normalized, ont_df) {
  syn_matches <- which(sapply(ont_df$synonyms, function(s) {
    !is.na(s) && grepl(normalized, s, fixed = TRUE)
  }))
  
  if (length(syn_matches) > 0) return(syn_matches[1])
  return(NULL)
}

.find_fuzzy_match <- function(normalized, ont_df) {
  distances <- stringdist::stringdist(normalized, ont_df$name_lower, method = "lv")
  min_dist <- min(distances)
  adaptive_threshold <- max(3, nchar(normalized) * 0.3)
  
  if (min_dist <= adaptive_threshold) {
    return(list(idx = which.min(distances), dist = min_dist))
  }
  return(NULL)
}

#' Normalize cell type string for matching
#' 
#' @param x Character string to normalize
#' @return Normalized string
#' @export
normalize_cell_type <- function(x) {
  x <- tolower(trimws(x))
  x <- gsub("\\(.*?\\)", "", x)  # Remove parentheses content
  x <- gsub("likely|probable|possible|putative", "", x)
  x <- gsub("\\+", " positive ", x)
  x <- gsub("[^a-z0-9\\s]", " ", x)  # Replace punctuation with space
  x <- gsub("\\s+", " ", x)  # Collapse multiple spaces
  x <- trimws(x)
  return(x)
}

#' Check if two CL IDs are related (lineage-based evaluation)
#' 
#' @param cl_id_1 First CL identifier
#' @param cl_id_2 Second CL identifier
#' @param ontology Ontology object from load_cell_ontology()
#' @return Logical indicating if terms are related
#' @export
are_related_terms <- function(cl_id_1, cl_id_2, ontology) {
  if (is.na(cl_id_1) || is.na(cl_id_2)) return(FALSE)
  if (cl_id_1 == cl_id_2) return(TRUE)
  
  ancestors_1 <- ontology$ancestor_cache[[cl_id_1]]
  ancestors_2 <- ontology$ancestor_cache[[cl_id_2]]
  
  if (is.null(ancestors_1) || is.null(ancestors_2)) return(FALSE)
  
  return(cl_id_2 %in% ancestors_1 || cl_id_1 %in% ancestors_2)
}

#' Calculate clade accuracy (lineage-based evaluation)
#' 
#' @param predicted_cl Vector of predicted CL IDs
#' @param ground_truth_cl Vector of ground truth CL IDs
#' @param ontology Ontology object
#' @return Proportion of related predictions
#' @export
calculate_clade_accuracy <- function(predicted_cl, ground_truth_cl, ontology) {
  valid_idx <- !is.na(predicted_cl) & !is.na(ground_truth_cl)
  if (sum(valid_idx) == 0) return(NA)
  
  related <- mapply(are_related_terms, 
                    predicted_cl[valid_idx], 
                    ground_truth_cl[valid_idx],
                    MoreArgs = list(ontology = ontology))
  
  return(mean(related))
}