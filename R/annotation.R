#' Core annotation orchestration
#'
#' Coordinates the entire cell type annotation pipeline from marker genes
#' to final annotated results with ontology mapping.
#'
#' @param markers Named list of marker genes per cluster
#' @param tissue Tissue name
#' @param species Species (Human/Mouse/Rat)
#' @param model_name Model to use ("deepseek" or "gpt4")
#' @param api_key API key for the selected model
#' @param use_ontology Whether to map to Cell Ontology
#' @param validate Whether to perform validation
#' @return Comprehensive result object
#' @export
annotate_cell_types <- function(markers,
                                tissue,
                                species = "Human",
                                model_name = "deepseek",
                                api_key = NULL,
                                use_ontology = TRUE,
                                validate = TRUE) {
  
  start_time <- Sys.time()
  message("Starting annotation: tissue=", tissue, ", species=", species,
          ", model=", model_name)
  
  # Input validation
  if (length(markers) == 0) stop("No marker genes provided.")
  if (is.null(api_key) || api_key == "") stop("API key is required.")
  
  model_config <- get_model_config(model_name)
  if (is.null(model_config)) stop("Unknown model: ", model_name)
  
  # Build prompt
  prompt <- create_annotation_prompt(markers, tissue, species,
                                     include_reasoning = FALSE)
  
  # Call LLM API
  api_result <- call_llm_api(prompt, model_config, api_key)
  if (!api_result$success) {
    return(list(success = FALSE, error = api_result$error,
                metadata = list(timestamp = Sys.time())))
  }
  
  # Parse response
  annotations <- parse_annotation_response(api_result$content)
  if (nrow(annotations) == 0) {
    annotations <- data.frame(Cluster = names(markers),
                              CellType = "Unknown",
                              Confidence = 0.5,
                              stringsAsFactors = FALSE)
  }
  
  # Ontology mapping
  if (use_ontology) {
    if (!exists("ONTOLOGY_DATA")) {
      ONTOLOGY_DATA <<- load_cell_ontology()
    }
    mappings <- lapply(annotations$CellType,
                       function(ct) map_to_cell_ontology(ct, ONTOLOGY_DATA))
    onto_df <- do.call(rbind, mappings)
    annotations <- cbind(annotations,
                         onto_df[, c("CL_ID", "OntologyLabel", "MatchMethod")])
  }
  
  # Validation
  validation <- NULL
  if (validate) {
    validation <- validate_annotations(annotations)
  }
  
  # Metadata
  total_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  tokens <- api_result$usage$total_tokens %||% 0
  cost <- tokens / 1000 * model_config$cost_per_1k_tokens
  
  list(
    success = TRUE,
    annotations = annotations,
    metadata = list(
      tissue = tissue,
      species = species,
      model = model_config$name,
      n_clusters = length(markers),
      n_annotations = nrow(annotations),
      api_latency_sec = api_result$latency_sec,
      total_runtime_sec = total_time,
      tokens_used = tokens,
      estimated_cost_usd = cost,
      timestamp = Sys.time()
    ),
    validation = validation
  )
}