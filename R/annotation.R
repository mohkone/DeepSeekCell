# R/annotation.R
#' Core annotation orchestration
#'
#' Coordinates the entire cell type annotation pipeline from marker genes
#' to final annotated results with ontology mapping.
#'
#' @param markers Named list of marker genes per cluster
#' @param tissue Tissue name
#' @param species Species (Human/Mouse/Rat)
#' @param model_name Model to use ("deepseek" or "ollama")
#' @param api_key API key for the selected model
#' @param use_ontology Whether to map to Cell Ontology
#' @param validate Whether to perform validation
#' @param max_genes_per_cluster Maximum marker genes to include per cluster.
#' @param ontology_path Optional path to a Cell Ontology OBO file.
#' @param return_prompt Whether to include the submitted prompt in the result.
#' @return Comprehensive result object
#' @export

annotate_cell_types <- function(markers,
                                tissue,
                                species = "Human",
                                model_name = "deepseek",
                                api_key = NULL,
                                use_ontology = TRUE,
                                validate = TRUE,
                                max_genes_per_cluster = 30,
                                ontology_path = NULL,
                                return_prompt = FALSE) {
  start_time <- Sys.time()
  
  if (is_blank(tissue)) {
    stop("tissue must be a non-empty character value.", call. = FALSE)
  }

  if (is_blank(species)) {
    stop("species must be a non-empty character value.", call. = FALSE)
  }

  markers <- normalize_marker_list(markers, max_genes = max_genes_per_cluster)
  
  if (length(markers) == 0) {
    stop("No valid marker genes provided after filtering.", call. = FALSE)
  }
  
  model_config <- get_model_config(model_name)
  api_key <- resolve_api_key(model_config, api_key)
  
  if (isTRUE(model_config$requires_api_key) && (is.null(api_key) || !nzchar(api_key))) {
    env_hint <- paste(model_config$api_key_env %||% character(), collapse = ", ")
    stop("API key is required for model: ", model_config$name,
         if (nzchar(env_hint)) paste0(". Set one of: ", env_hint) else "",
         call. = FALSE)
  }
  
  prompt <- create_annotation_prompt(
    markers = markers,
    tissue = tissue,
    species = species,
    include_reasoning = TRUE
  )
  
  api_result <- call_llm_api(prompt, model_config, api_key)
  
  if (!isTRUE(api_result$success)) {
    return(list(
      success = FALSE,
      error = api_result$error,
      markers = markers,
      metadata = list(
        tissue = tissue,
        species = species,
        model = model_config$name,
        model_id = model_config$model_id,
        schema_version = deepseekcell_version(),
        timestamp = Sys.time()
      )
    ))
  }
  
  annotations <- parse_annotation_response(api_result$content)
  
  if (nrow(annotations) == 0) {
    annotations <- .unknown_annotation_rows(names(markers), "No parseable annotation returned by model.")
  }
  
  annotations <- .complete_annotation_rows(annotations, names(markers))
  ontology_is_fallback <- NA
  
  if (isTRUE(use_ontology)) {
    ontology_data <- load_cell_ontology(ontology_path = ontology_path)
    ontology_is_fallback <- isTRUE(ontology_data$is_fallback)
    mappings <- do.call(
      rbind,
      lapply(
        annotations$CellType,
        map_to_cell_ontology,
        ontology = ontology_data,
        tissue = tissue
      )
    )
    
    annotations$CL_ID <- mappings$CL_ID
    annotations$OntologyLabel <- mappings$OntologyLabel
    annotations$MatchMethod <- mappings$MatchMethod
    annotations$OntologyMatchScore <- mappings$OntologyMatchScore
  }
  
  validation <- if (isTRUE(validate)) {
    validate_annotations(
      annotations,
      markers = markers,
      metadata = list(tissue = tissue, species = species, model = model_config$name)
    )
  } else {
    NULL
  }
  
  total_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  tokens <- api_result$usage$total_tokens %||% 0
  if (is.na(tokens)) tokens <- 0
  cost <- (tokens / 1000) * (model_config$cost_per_1k_tokens %||% 0)

  result <- list(
    success = TRUE,
    annotations = annotations,
    markers = markers,
    metadata = list(
      tissue = tissue,
      species = species,
      model = model_config$name,
      model_id = model_config$model_id,
      n_clusters = length(markers),
      n_annotations = nrow(annotations),
      api_latency_sec = api_result$latency_sec,
      total_runtime_sec = total_time,
      tokens_used = tokens,
      estimated_cost_usd = cost,
      ontology_enabled = isTRUE(use_ontology),
      ontology_is_fallback = ontology_is_fallback,
      schema_version = deepseekcell_version(),
      timestamp = Sys.time()
    ),
    validation = validation
  )

  if (isTRUE(return_prompt)) {
    result$prompt <- prompt
  }

  result
}

.unknown_annotation_rows <- function(cluster_names, reason) {
  data.frame(
    Cluster = cluster_names,
    CellType = "Unknown",
    Confidence = 0.5,
    IsMixed = FALSE,
    PrimaryCellType = "Unknown",
    SecondaryCellType = "",
    TissueConsistency = "unknown",
    Reasoning = reason,
    stringsAsFactors = FALSE
  )
}

.complete_annotation_rows <- function(annotations, cluster_names) {
  annotations <- annotations[!duplicated(annotations$Cluster), , drop = FALSE]

  missing_clusters <- setdiff(cluster_names, annotations$Cluster)
  if (length(missing_clusters) > 0) {
    annotations <- rbind(
      annotations,
      .unknown_annotation_rows(missing_clusters, "No annotation returned by model.")
    )
  }

  annotations <- annotations[match(cluster_names, annotations$Cluster), , drop = FALSE]
  rownames(annotations) <- NULL
  annotations
}
