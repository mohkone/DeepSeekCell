# Cell Ontology Loader
#
# Handles loading, caching, and preprocessing of the Cell Ontology OBO file.

.ontology_cache_env <- new.env(parent = emptyenv())
.ontology_cache_env$cache <- NULL

init_ontology_cache <- function(max_size = 1024 * 1024^2) {
  .ontology_cache_env$cache <- cachem::cache_mem(max_size = max_size)
  return(.ontology_cache_env$cache)
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
  
  if (is.null(.ontology_cache_env$cache)) init_ontology_cache()
  ontology_cache <- .ontology_cache_env$cache
  
  ontology_path_norm <- normalizePath(ontology_path, winslash = "/", mustWork = FALSE)
  cache_key <- .ontology_cache_key(ontology_path_norm, ontology_path)
  if (!force_reload && ontology_cache$exists(cache_key)) {
    return(ontology_cache$get(cache_key))
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
      extract_tags = "everything"
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
    ontology_path = ontology_path_norm,
    load_time = Sys.time(),
    is_fallback = FALSE
  )
  
  ontology_cache$set(cache_key, result)
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
  
  get_value <- function(x, id) {
    val <- x[[id]]
    if (is.null(val) || length(val) == 0) return(NA_character_)
    paste(val, collapse = "; ")
  }
  
  cl_ids <- ont$id[grepl("^CL:", ont$id)]
  
  names_vec <- unname(ont$name[cl_ids])
  keep <- !is.na(names_vec) & nzchar(names_vec)
  
  cl_ids <- cl_ids[keep]
  names_vec <- names_vec[keep]
  
  data.frame(
    cl_id = cl_ids,
    name = names_vec,
    name_lower = normalize_cell_type(names_vec),
    synonyms = vapply(cl_ids, function(id) get_value(ont$synonym, id), character(1)),
    definition = vapply(cl_ids, function(id) get_value(ont$def, id), character(1)),
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
    cl_id = c(
      "CL:0000084", "CL:0000236", "CL:0000576", "CL:0000623",
      "CL:0000451", "CL:0000233", "CL:0000540", "CL:0000127",
      "CL:0000128", "CL:0000129", "CL:0000097", "CL:0000669",
      "CL:0000192", "CL:0000171", "CL:0000169", "CL:0000173",
      "CL:0002064", "CL:0002079", "CL:0000622", "CL:2000029",
      "CL:0000125", "CL:0002453", "CL:0000065", "CL:0000099",
      "CL:0008031", "CL:0011005", "CL:0000617", "CL:0000679",
      "CL:0000700", "CL:0000598", "CL:1001474", "CL:0000047",
      "CL:0000681", "CL:0000031"
    ),
    name = c(
      "T cell", "B cell", "monocyte", "natural killer cell",
      "dendritic cell", "platelet", "neuron", "astrocyte",
      "oligodendrocyte", "microglial cell", "mast cell", "pericyte",
      "smooth muscle cell", "pancreatic A cell", "type B pancreatic cell",
      "pancreatic D cell", "pancreatic acinar cell",
      "pancreatic ductal cell", "acinar cell",
      "central nervous system neuron", "glial cell",
      "oligodendrocyte precursor cell", "ependymal cell",
      "interneuron", "cortical interneuron", "GABAergic interneuron",
      "GABAergic neuron", "glutamatergic neuron", "dopaminergic neuron",
      "pyramidal neuron", "medium spiny neuron", "neural stem cell",
      "radial glial cell", "neuroblast (sensu Vertebrata)"
    ),
    synonyms = c(
      NA_character_, "B lymphocyte", NA_character_, "NK cell",
      NA_character_, NA_character_, NA_character_, NA_character_,
      NA_character_, "microglia", NA_character_, NA_character_,
      "smooth muscle", paste(
        '"pancreatic alpha cell"',
        '"alpha cell of islet of Langerhans"',
        '"alpha cell"',
        sep = "; "
      ), paste(
        '"pancreatic beta cell"',
        '"beta cell of islet of Langerhans"',
        '"beta cell"',
        sep = "; "
      ), paste(
        '"pancreatic delta cell"',
        '"delta cell of pancreatic islet"',
        '"delta cell"',
        sep = "; "
      ), paste(
        '"acinar cell of pancreas"',
        '"acinar cell"',
        sep = "; "
      ), paste(
        '"duct epithelial cell of pancreas"',
        '"ductal cell"',
        sep = "; "
      ), NA_character_, paste(
        '"brain neuron"',
        '"CNS neuron"',
        '"neuronal cell"',
        sep = "; "
      ), paste(
        '"glia"',
        '"neuroglia"',
        sep = "; "
      ), paste(
        '"OPC"',
        '"oligodendrocyte precursor"',
        '"NG2 cell"',
        '"NG2 glia"',
        sep = "; "
      ), '"ependymal"', NA_character_, NA_character_, paste(
        '"inhibitory interneuron"',
        '"GABA-ergic interneuron"',
        sep = "; "
      ), paste(
        '"inhibitory neuron"',
        '"GABA-ergic neuron"',
        sep = "; "
      ), '"excitatory neuron"', '"dopaminergic cell"', paste(
        '"pyramidal cell"',
        '"projection neuron"',
        sep = "; "
      ), '"MSN"', NA_character_, '"radial glia"', '"neuroblast"'
    ),
    definition = NA_character_,
    stringsAsFactors = FALSE
  )

  lung_df <- data.frame(
    cl_id = c(
      "CL:0000082", "CL:0000322", "CL:0002062", "CL:0002063",
      "CL:0000158", "CL:1000271", "CL:0002145", "CL:0002633",
      "CL:0002329", "CL:1000143", "CL:0002370", "CL:1000272",
      "CL:4052031", "CL:0017000", "CL:1000223", "CL:0002075",
      "CL:0002328", "CL:1001567", "CL:2000016", "CL:4028001",
      "CL:0002553", "CL:0002241", "CL:0009089", "CL:1001603",
      "CL:0000583", "CL:4033043", "CL:0019019", "CL:0002598"
    ),
    name = c(
      "epithelial cell of lung", "pulmonary alveolar epithelial cell",
      "pulmonary alveolar type 1 cell", "pulmonary alveolar type 2 cell",
      "club cell", "lung multiciliated epithelial cell",
      "multiciliated columnar cell of tracheobronchial tree",
      "respiratory basal cell",
      "basal epithelial cell of tracheobronchial tree",
      "lung goblet cell", "respiratory tract goblet cell",
      "lung secretory cell", "respiratory airway secretory cell",
      "pulmonary ionocyte", "pulmonary neuroendocrine cell",
      "brush cell of tracheobronchial tree", "bronchial epithelial cell",
      "lung endothelial cell", "lung microvascular endothelial cell",
      "pulmonary capillary endothelial cell", "fibroblast of lung",
      "pulmonary interstitial fibroblast", "lung pericyte",
      "lung macrophage", "alveolar macrophage",
      "lung interstitial macrophage",
      "tracheobronchial smooth muscle cell",
      "bronchial smooth muscle cell"
    ),
    synonyms = c(
      '"lung epithelial cell"', paste(
        '"alveolar epithelial cell"',
        '"pneumocyte"',
        sep = "; "
      ), paste(
        '"AT1 cell"',
        '"AT1"',
        '"ATI cell"',
        '"type I pneumocyte"',
        '"type 1 pneumocyte"',
        '"alveolar type 1 cell"',
        sep = "; "
      ), paste(
        '"AT2 cell"',
        '"AT2"',
        '"ATII cell"',
        '"type II pneumocyte"',
        '"type 2 pneumocyte"',
        '"alveolar type 2 cell"',
        sep = "; "
      ), '"Clara cell"', paste(
        '"ciliated cell"',
        '"ciliated epithelial cell"',
        '"multiciliated cell"',
        sep = "; "
      ), '"tracheobronchial ciliated cell"', '"basal cell"',
      '"airway basal cell"', '"goblet cell"',
      '"respiratory goblet cell"', '"secretory cell"', paste(
        '"airway secretory cell"',
        '"RASC"',
        sep = "; "
      ), '"ionocyte"', '"neuroendocrine cell"', '"tuft cell"',
      NA_character_, '"endothelial cell"',
      '"microvascular endothelial cell"', '"capillary endothelial cell"',
      paste(
        '"fibroblast"',
        '"lung fibroblast"',
        sep = "; "
      ), '"interstitial fibroblast"', '"pericyte"', '"macrophage"',
      NA_character_, '"interstitial macrophage"',
      paste(
        '"smooth muscle cell"',
        '"airway smooth muscle cell"',
        sep = "; "
      ), NA_character_
    ),
    definition = NA_character_,
    stringsAsFactors = FALSE
  )

  fallback_df <- rbind(fallback_df, lung_df)

  fallback_df$name_lower <- normalize_cell_type(fallback_df$name)
  fallback_df <- fallback_df[, c(
    "cl_id", "name", "name_lower", "synonyms", "definition"
  )]

  list(
    dataframe = fallback_df,
    ancestor_cache = list(),
    ontology_path = NA_character_,
    load_time = Sys.time(),
    is_fallback = TRUE
  )
}

.ontology_cache_key <- function(ontology_path_norm, ontology_path) {
  mtime <- if (file.exists(ontology_path)) {
    format(file.info(ontology_path)$mtime, "%Y%m%d%H%M%S")
  } else {
    "missing"
  }

  raw_key <- paste0("cellontology", ontology_path_norm, mtime)
  safe_key <- gsub("[^a-z0-9]", "", tolower(raw_key), perl = TRUE)

  if (!nzchar(safe_key)) {
    return("cellontology")
  }

  safe_key
}
