# R/annotation.R
#' Core annotation orchestration
#'
#' Coordinates the marker-based LLM annotation pipeline from first-pass
#' candidate generation to evidence-guided reliability assessment and optional
#' fixed-budget refinement.
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
#' @param calibrate_confidence Whether to replace raw LLM confidence with
#' evidence-adjusted confidence while preserving LLMConfidence.
#' @param self_refine Whether to run a selective second-pass refinement call.
#' @param refinement_strategy Strategy for selecting refined clusters. Use
#' `"evidence"` for Evidence-k, `"confidence"` for Confidence-k, `"random"` for
#' Random-k, `"full"` for FullRefined, or `"none"` to disable selection.
#' @param refinement_budget Maximum number of clusters to refine. If `NULL`,
#' Evidence-k refines all evidence-conflicted clusters, Confidence-k and
#' Random-k use the same matched budget, and FullRefined refines all clusters.
#' @param use_ontology_evidence Whether ontology mapping contributes to
#' evidence-adjusted confidence and conflict detection.
#' @param refinement_seed Optional random seed for Random-k selection.
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
                                calibrate_confidence = TRUE,
                                self_refine = FALSE,
                                refinement_strategy = c(
                                  "evidence",
                                  "confidence",
                                  "random",
                                  "full",
                                  "none"
                                ),
                                refinement_budget = NULL,
                                use_ontology_evidence = TRUE,
                                refinement_seed = NULL,
                                return_prompt = FALSE) {
  start_time <- Sys.time()
  refinement_strategy <- match.arg(refinement_strategy)
  
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
    annotations <- .add_ontology_mappings(annotations, ontology_data, tissue)
  }
  first_pass_annotations <- annotations

  refinement <- list(
    enabled = isTRUE(self_refine),
    strategy = refinement_strategy,
    budget = refinement_budget,
    attempted = FALSE,
    success = NA,
    n_flagged = 0,
    n_candidates = 0,
    n_reviewed = 0,
    n_returned = 0,
    n_updated = 0,
    n_label_changed = 0,
    n_confidence_changed = 0,
    error = NULL
  )
  flagged_clusters <- character()
  refined_clusters <- character()
  label_changed_clusters <- character()
  confidence_changed_clusters <- character()
  selected_clusters <- character()
  selected_scores <- numeric()

  preliminary_annotations <- if (isTRUE(calibrate_confidence) || isTRUE(self_refine)) {
    calibrate_annotation_confidence(
      annotations,
      markers,
      tissue = tissue,
      use_ontology_evidence = use_ontology_evidence
    )
  } else {
    annotations
  }

  requires_refinement <- preliminary_annotations$RequiresRefinement %||%
    rep(FALSE, nrow(preliminary_annotations))
  requires_refinement <- as_flag(requires_refinement)
  requires_refinement[is.na(requires_refinement)] <- FALSE
  flagged_clusters <- preliminary_annotations$Cluster[requires_refinement]
  refinement$n_flagged <- length(flagged_clusters)

  if (isTRUE(self_refine)) {
    refinement_candidates <- select_refinement_candidates(
      preliminary_annotations,
      strategy = refinement_strategy,
      budget = refinement_budget,
      seed = refinement_seed
    )
    selected_clusters <- refinement_candidates$Cluster %||% character()
    selected_scores <- refinement_candidates$SelectionScore %||% numeric()
    refinement$n_candidates <- length(flagged_clusters)
    refinement$n_reviewed <- nrow(refinement_candidates)
    refinement$n_selected <- nrow(refinement_candidates)
    refinement$selected_clusters <- selected_clusters
    refinement$selection_scores <- selected_scores

    if (nrow(refinement_candidates) > 0) {
      refinement$attempted <- TRUE
      refinement_prompt <- create_refinement_prompt(
        markers = markers,
        annotations = refinement_candidates,
        tissue = tissue,
        species = species
      )
      refinement_result <- call_llm_api(refinement_prompt, model_config, api_key)
      refinement$success <- isTRUE(refinement_result$success)

      if (isTRUE(refinement_result$success)) {
        refined_annotations <- parse_annotation_response(refinement_result$content)
        refined_annotations <- refined_annotations[
          refined_annotations$Cluster %in% annotations$Cluster,
          ,
          drop = FALSE
        ]

        if (nrow(refined_annotations) > 0) {
          change_summary <- .summarise_refinement_changes(
            first_pass_annotations,
            refined_annotations
          )
          annotations <- .replace_annotation_rows(annotations, refined_annotations)
          if (isTRUE(use_ontology)) {
            annotations <- .add_ontology_mappings(annotations, ontology_data, tissue)
          }
          refined_clusters <- refined_annotations$Cluster
          label_changed_clusters <- change_summary$label_changed_clusters
          confidence_changed_clusters <- change_summary$confidence_changed_clusters
          refinement$n_returned <- nrow(refined_annotations)
          refinement$n_updated <- change_summary$n_label_changed
          refinement$n_label_changed <- change_summary$n_label_changed
          refinement$n_confidence_changed <- change_summary$n_confidence_changed
          refinement$label_changed_clusters <- label_changed_clusters
          refinement$confidence_changed_clusters <- confidence_changed_clusters
          refinement$usage <- refinement_result$usage
          refinement$latency_sec <- refinement_result$latency_sec
          if (isTRUE(return_prompt)) {
            refinement$prompt <- refinement_prompt
          }
        }
      } else {
        refinement$error <- refinement_result$error %||% "Refinement call failed."
      }
    }
  }

  if (isTRUE(calibrate_confidence)) {
    annotations <- calibrate_annotation_confidence(
      annotations,
      markers,
      tissue = tissue,
      use_ontology_evidence = use_ontology_evidence
    )
  }
  annotations <- .add_refinement_provenance(
    annotations = annotations,
    first_pass_annotations = first_pass_annotations,
    flagged_clusters = flagged_clusters,
    refined_clusters = refined_clusters,
    label_changed_clusters = label_changed_clusters,
    confidence_changed_clusters = confidence_changed_clusters
  )
  
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
  first_tokens <- .safe_token_count(api_result$usage$total_tokens %||% 0)
  refine_tokens <- .safe_token_count(refinement$usage$total_tokens %||% 0)
  tokens <- first_tokens + refine_tokens
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
      confidence_calibrated = isTRUE(calibrate_confidence),
      ontology_evidence_enabled = isTRUE(use_ontology_evidence),
      self_refinement_enabled = isTRUE(self_refine),
      self_refinement_strategy = refinement_strategy,
      self_refinement_budget = refinement_budget,
      self_refinement_flagged = refinement$n_flagged,
      self_refinement_attempted = isTRUE(refinement$attempted),
      self_refinement_updates = refinement$n_updated,
      self_refinement_reviewed = refinement$n_reviewed,
      self_refinement_label_changes = refinement$n_label_changed,
      self_refinement_confidence_changes = refinement$n_confidence_changed,
      schema_version = deepseekcell_version(),
      timestamp = Sys.time()
    ),
    validation = validation,
    refinement = refinement
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
    CandidateCellTypes = "",
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

.add_ontology_mappings <- function(annotations, ontology_data, tissue) {
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
  annotations
}

.replace_annotation_rows <- function(annotations, refined_annotations) {
  base_cols <- c(
    "Cluster", "CellType", "Confidence", "CandidateCellTypes", "IsMixed",
    "PrimaryCellType", "SecondaryCellType", "TissueConsistency", "Reasoning"
  )
  shared_cols <- intersect(base_cols, names(refined_annotations))

  for (i in seq_len(nrow(refined_annotations))) {
    cluster <- refined_annotations$Cluster[[i]]
    idx <- match(cluster, annotations$Cluster)
    if (is.na(idx)) next

    for (column in shared_cols) {
      annotations[[column]][idx] <- refined_annotations[[column]][i]
    }
  }

  annotations
}

.summarise_refinement_changes <- function(first_pass_annotations, refined_annotations) {
  if (!is.data.frame(refined_annotations) || nrow(refined_annotations) == 0) {
    return(list(
      n_label_changed = 0,
      n_confidence_changed = 0,
      label_changed_clusters = character(),
      confidence_changed_clusters = character()
    ))
  }

  idx <- match(refined_annotations$Cluster, first_pass_annotations$Cluster)
  valid <- !is.na(idx)
  if (!any(valid)) {
    return(list(
      n_label_changed = 0,
      n_confidence_changed = 0,
      label_changed_clusters = character(),
      confidence_changed_clusters = character()
    ))
  }

  refined <- refined_annotations[valid, , drop = FALSE]
  first <- first_pass_annotations[idx[valid], , drop = FALSE]
  label_changed <- normalize_cell_type(refined$CellType) !=
    normalize_cell_type(first$CellType)
  confidence_changed <- abs(
    as_confidence(refined$Confidence) - as_confidence(first$Confidence)
  ) > 1e-6

  list(
    n_label_changed = sum(label_changed, na.rm = TRUE),
    n_confidence_changed = sum(confidence_changed, na.rm = TRUE),
    label_changed_clusters = refined$Cluster[label_changed],
    confidence_changed_clusters = refined$Cluster[confidence_changed]
  )
}

.add_refinement_provenance <- function(annotations,
                                       first_pass_annotations,
                                       flagged_clusters = character(),
                                       refined_clusters = character(),
                                       label_changed_clusters = character(),
                                       confidence_changed_clusters = character()) {
  idx <- match(annotations$Cluster, first_pass_annotations$Cluster)

  annotations$FirstPassCellType <- first_pass_annotations$CellType[idx]
  annotations$FirstPassConfidence <- as_confidence(first_pass_annotations$Confidence[idx])
  annotations$WasFlagged <- annotations$Cluster %in% flagged_clusters
  annotations$WasRefined <- annotations$Cluster %in% refined_clusters
  annotations$RefinementChangedLabel <- annotations$Cluster %in% label_changed_clusters
  annotations$RefinementChangedConfidence <- annotations$Cluster %in%
    confidence_changed_clusters
  annotations$RefinementReason <- ifelse(
    annotations$WasRefined,
    ifelse(
      annotations$RefinementChangedLabel,
      "refinement_changed_label",
      "refinement_reviewed_retained_label"
    ),
    ifelse(annotations$WasFlagged, "flagged_not_refined", "")
  )

  annotations
}

.safe_token_count <- function(x) {
  tokens <- suppressWarnings(as.numeric(x))
  if (length(tokens) == 0 || all(is.na(tokens))) {
    return(0)
  }

  sum(tokens[!is.na(tokens)])
}
