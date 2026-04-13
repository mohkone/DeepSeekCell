#' Cell Ontology Loader
#' 
#' Handles loading, caching, and preprocessing of the Cell Ontology OBO file.
#' 
#' @import ontologyIndex

.ontology_cache <- NULL

init_ontology_cache <- function(max_size = 100 * 1024^2) {
  .ontology_cache <<- cachem::cache_mem(max_size = max_size)
  return(.ontology_cache)
}

#' Load Cell Ontology with caching and validation
#' 
#' @param ontology_path Path to cl.obo file (optional)
#' @param force_reload Force reload from disk
#' @return List containing ontology data frame and ancestor mapping
#' @export
load_cell_ontology <- function(ontology_path = NULL, force_reload = FALSE) {
  
  if (is.null(ontology_path)) {
    ontology_path <- get_default_ontology_path()
  }
  
  if (is.null(.ontology_cache)) init_ontology_cache()
  
  cache_key <- "cell_ontology"
  if (!force_reload && .ontology_cache$exists(cache_key)) {
    message("Returning cached ontology")
    return(.ontology_cache$get(cache_key))
  }
  
  if (!file.exists(ontology_path)) {
    warning("Cell Ontology file not found at: ", ontology_path)
    return(create_fallback_ontology())
  }
  
  message("Loading Cell Ontology from: ", ontology_path)
  
  # Use try() to catch *all* conditions, including the vectorised one
  ont <- try(
    ontologyIndex::get_ontology(
      ontology_path,
      extract_tags = c("name", "synonym", "is_a", "relationship", "def")
    ),
    silent = TRUE
  )
  
  if (inherits(ont, "try-error") || is.null(ont) || length(ont$id) == 0) {
    message("Ontology loading failed, using fallback")
    return(create_fallback_ontology())
  }
  
  ontology_df <- .build_ontology_dataframe(ont)
  ancestor_cache <- .build_ancestor_cache(ont)
  
  result <- list(
    dataframe = ontology_df,
    ancestor_cache = ancestor_cache,
    ontology_object = ont,
    load_time = Sys.time(),
    is_fallback = FALSE
  )
  
  .ontology_cache$set(cache_key, result)
  message("Ontology loaded successfully with ", nrow(ontology_df), " terms")
  return(result)
}

#' Get default ontology path
#' 
#' @return Path to cl.obo file
#' @export
get_default_ontology_path <- function() {
  possible_paths <- c(
    "data/cl.obo",
    "../data/cl.obo",
    "../../data/cl.obo",
    system.file("extdata", "cl.obo", package = "DeepSeekCell")
  )
  for (path in possible_paths) {
    if (file.exists(path)) return(normalizePath(path))
  }
  return("data/cl.obo")
}

.build_ontology_dataframe <- function(ont) {
  data.frame(
    cl_id = ont$id,
    name = ont$name,
    name_lower = tolower(ont$name),
    definition = sapply(ont$id, function(id) ont$def[[id]] %||% NA_character_),
    stringsAsFactors = FALSE
  )
}

.build_ancestor_cache <- function(ont) {
  ancestors <- list()
  for (term_id in ont$id) {
    term_ancestors <- .traverse_ancestors(term_id, ont)
    ancestors[[term_id]] <- unique(term_ancestors)
  }
  return(ancestors)
}

.traverse_ancestors <- function(term_id, ont, visited = NULL) {
  if (is.null(visited)) visited <- character()
  if (term_id %in% visited) return(character())
  visited <- c(visited, term_id)
  parents <- ont$is_a[[term_id]]
  if (is.null(parents) || length(parents) == 0) return(character())
  all_ancestors <- parents
  for (parent in parents) {
    all_ancestors <- c(all_ancestors, .traverse_ancestors(parent, ont, visited))
  }
  return(unique(all_ancestors))
}

#' Create fallback ontology when OBO file unavailable
#' 
#' @return Fallback ontology list
#' @export
create_fallback_ontology <- function() {
  message("Creating fallback ontology mapping")
  fallback_df <- data.frame(
    cl_id = c("CL:0000084", "CL:0000236", "CL:0000576", "CL:0000623",
              "CL:0000451", "CL:0000233", "CL:0000540", "CL:0000127",
              "CL:0000128", "CL:0000129", "CL:0000502", "CL:0000169",
              "CL:0000503", "CL:0000168", "CL:0000167"),
    name = c("T cell", "B cell", "monocyte", "natural killer cell",
             "dendritic cell", "platelet", "neuron", "astrocyte",
             "oligodendrocyte", "microglial cell", "alpha cell", 
             "beta cell", "delta cell", "acinar cell", "ductal cell"),
    name_lower = tolower(c("T cell", "B cell", "monocyte", "natural killer cell",
                           "dendritic cell", "platelet", "neuron", "astrocyte",
                           "oligodendrocyte", "microglial cell", "alpha cell", 
                           "beta cell", "delta cell", "acinar cell", "ductal cell")),
    definition = NA_character_,
    stringsAsFactors = FALSE
  )
  fallback_df$synonyms <- NA_character_
  list(
    dataframe = fallback_df,
    ancestor_cache = list(),
    load_time = Sys.time(),
    is_fallback = TRUE
  )
}