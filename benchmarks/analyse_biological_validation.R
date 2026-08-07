# benchmarks/analyse_biological_validation.R
#
# Biological plausibility checks for refined DeepSeekCell annotations. The
# script evaluates whether label changes are supported by canonical markers,
# Cell Ontology hierarchy, marker-derived differential-expression agreement,
# lineage consistency, and optional GO biological-process enrichment.
#
# Usage:
#   Rscript benchmarks/analyse_biological_validation.R \
#     results/benchmark_debug \
#     benchmarks/external_validation_manifest.csv \
#     results/biological_validation

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

source("R/utils.R")
source("R/ontology.R")
source("R/ontology_loader.R")
source("R/refinement.R")

if (file.exists("benchmarks/metrics.R")) {
  source("benchmarks/metrics.R")
}

as_text <- function(x, default = "") {
  out <- as.character(x %||% default)
  out[is.na(out)] <- default
  out
}

as_number <- function(x, default = NA_real_) {
  out <- suppressWarnings(as.numeric(x))
  out[!is.finite(out)] <- default
  out
}

as_bool <- function(x) {
  value <- tolower(trimws(as.character(x %||% FALSE)))
  value %in% c("true", "1", "yes", "y")
}

safe_divide <- function(numerator, denominator) {
  out <- numerator / denominator
  out[!is.finite(out)] <- NA_real_
  out
}

normalise_key <- function(x) {
  normalize_cell_type(as_text(x))
}

label_exact <- function(a, b) {
  normalise_key(a) == normalise_key(b)
}

.bio_validation_cache <- new.env(parent = emptyenv())
.bio_validation_cache$label_mapping <- new.env(parent = emptyenv())
.bio_validation_cache$ontology_relation <- new.env(parent = emptyenv())
.bio_validation_cache$go_available <- new.env(parent = emptyenv())
.bio_validation_cache$go_terms <- new.env(parent = emptyenv())

cache_key <- function(...) {
  paste(vapply(list(...), function(x) paste(as_text(x), collapse = ";"), character(1)), collapse = "||")
}

read_csv_required <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " not found: ", path, call. = FALSE)
  }
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

load_ontology_for_biology <- function(path = "data/cl.obo") {
  if (file.exists(path)) {
    load_cell_ontology(path)
  } else {
    create_fallback_ontology()
  }
}

load_manifest_datasets <- function(manifest_path) {
  manifest <- read_csv_required(manifest_path, "Validation manifest")
  required <- c("Dataset", "PreparedRdsPath")
  missing <- setdiff(required, names(manifest))
  if (length(missing) > 0) {
    stop("Manifest is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  out <- list()
  for (i in seq_len(nrow(manifest))) {
    dataset_id <- as.character(manifest$Dataset[[i]])
    path <- as.character(manifest$PreparedRdsPath[[i]])
    if (!nzchar(dataset_id) || !file.exists(path)) {
      next
    }
    x <- readRDS(path)
    x$tissue <- x$tissue %||% manifest$Tissue[[i]] %||% NA_character_
    x$species <- x$species %||% manifest$Species[[i]] %||% NA_character_
    x$manifest_row <- as.list(manifest[i, , drop = FALSE])
    out[[dataset_id]] <- x
  }
  out
}

debug_csv_files <- function(debug_dir, datasets = NULL) {
  files <- list.files(debug_dir, pattern = "_debug[.]csv$", full.names = TRUE)
  if (length(files) == 0 || is.null(datasets)) {
    return(files)
  }
  keep <- vapply(files, function(path) {
    any(vapply(datasets, function(dataset) {
      grepl(paste0("_", dataset, "_"), basename(path), fixed = TRUE)
    }, logical(1)))
  }, logical(1))
  files[keep]
}

read_debug_tables <- function(debug_dir, datasets) {
  files <- debug_csv_files(debug_dir, names(datasets))
  if (length(files) == 0) {
    stop("No cluster-level debug CSV files found in ", debug_dir, call. = FALSE)
  }

  rows <- lapply(files, function(path) {
    x <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    x$SourceDebugFile <- normalizePath(path, winslash = "/", mustWork = FALSE)
    x
  })
  debug <- do.call(rbind, rows)
  debug[debug$Dataset %in% names(datasets), , drop = FALSE]
}

cluster_markers <- function(dataset, cluster) {
  markers <- dataset$markers %||% list()
  cluster <- as.character(cluster)
  direct <- markers[[cluster]]
  if (!is.null(direct)) {
    return(as.character(direct))
  }
  alt <- paste0("Cluster", cluster)
  markers[[alt]] %||% character()
}

marker_support <- function(label, markers, tissue) {
  profiles <- .get_marker_evidence_profiles(tissue)
  match <- .score_marker_profiles(markers, profiles, label)
  predicted_idx <- .match_profile_name(label, names(profiles))
  matched <- character()
  if (length(predicted_idx) > 0) {
    scores <- lapply(profiles[predicted_idx], .marker_profile_score, markers = markers)
    best <- which.max(vapply(scores, `[[`, numeric(1), "score"))
    matched <- scores[[best]]$matched_markers
  }

  list(
    score = as.numeric(match$predicted_score),
    best_score = as.numeric(match$best_score),
    best_label = as.character(match$best_cell_type),
    matched_markers = matched,
    best_matched_markers = match$best_markers
  )
}

map_label <- function(label, ontology, tissue) {
  key <- cache_key(normalise_key(label), normalise_key(tissue))
  if (exists(key, envir = .bio_validation_cache$label_mapping, inherits = FALSE)) {
    return(get(key, envir = .bio_validation_cache$label_mapping, inherits = FALSE))
  }
  tryCatch(
    {
      out <- map_to_cell_ontology(label, ontology, tissue = tissue)
      assign(key, out, envir = .bio_validation_cache$label_mapping)
      out
    },
    error = function(e) data.frame(
      CellType = as.character(label),
      CL_ID = NA_character_,
      OntologyLabel = NA_character_,
      MatchMethod = "error",
      OntologyMatchScore = NA_real_,
      stringsAsFactors = FALSE
    )
  )
}

ontology_distance_to_truth <- function(label, truth, ontology, tissue) {
  key <- cache_key(normalise_key(label), normalise_key(truth), normalise_key(tissue))
  if (exists(key, envir = .bio_validation_cache$ontology_relation, inherits = FALSE)) {
    return(get(key, envir = .bio_validation_cache$ontology_relation, inherits = FALSE))
  }
  label_map <- map_label(label, ontology, tissue)
  truth_map <- map_label(truth, ontology, tissue)
  label_id <- as.character(label_map$CL_ID[[1]])
  truth_id <- as.character(truth_map$CL_ID[[1]])
  if (is.na(label_id) || is.na(truth_id) || !nzchar(label_id) || !nzchar(truth_id)) {
    out <- list(
      cl_id = label_id,
      truth_cl_id = truth_id,
      exact = FALSE,
      related = FALSE,
      relation = "unmapped",
      score = 0,
      ontology_label = as.character(label_map$OntologyLabel[[1]] %||% NA_character_),
      match_method = as.character(label_map$MatchMethod[[1]] %||% NA_character_)
    )
    assign(key, out, envir = .bio_validation_cache$ontology_relation)
    return(out)
  }
  exact <- identical(label_id, truth_id)
  related <- are_related_terms(label_id, truth_id, ontology)
  relation <- if (exact) {
    "exact"
  } else if (related) {
    "ancestor_descendant"
  } else {
    "unrelated"
  }
  out <- list(
    cl_id = label_id,
    truth_cl_id = truth_id,
    exact = exact,
    related = related,
    relation = relation,
    score = if (exact) 1 else if (related) 0.5 else 0,
    ontology_label = as.character(label_map$OntologyLabel[[1]] %||% NA_character_),
    match_method = as.character(label_map$MatchMethod[[1]] %||% NA_character_)
  )
  assign(key, out, envir = .bio_validation_cache$ontology_relation)
  out
}

lineage_from_label <- function(label, markers = character()) {
  text <- normalise_key(label)
  genes <- toupper(markers)
  has_any <- function(patterns) any(patterns %in% genes)

  if (grepl("t cell|lymphocyte|cd4|cd8|natural killer|nk cell", text) ||
      has_any(c("CD3D", "CD3E", "NKG7", "GNLY", "CD4", "CD8A"))) return("lymphoid")
  if (grepl("b cell|plasma|lymphoid progenitor", text) ||
      has_any(c("MS4A1", "CD79A", "CD79B", "JCHAIN", "MZB1"))) return("lymphoid")
  if (grepl("monocyte|macrophage|dendritic|granulocyte|myeloid", text) ||
      has_any(c("LYZ", "LST1", "S100A8", "S100A9", "C1QA", "MPO"))) return("myeloid")
  if (grepl("hematopoietic|erythroid|megakaryocyte|progenitor|hspc|stem", text) ||
      has_any(c("CD34", "GATA1", "HBB", "PF4", "PROM1"))) return("hematopoietic_progenitor")
  if (grepl("alpha cell|beta cell|delta cell|gamma cell|epsilon cell|islet|endocrine", text) ||
      has_any(c("INS", "IAPP", "GCG", "SST", "PPY", "GHRL"))) return("endocrine")
  if (grepl("acinar|ductal|epithelial|basal|club|ciliated|goblet|ionocyte|alveolar|tubule|podocyte|intercalated|collecting duct", text) ||
      has_any(c("KRT19", "KRT5", "SFTPC", "SCGB1A1", "AQP2", "SLC34A1", "NPHS1"))) return("epithelial")
  if (grepl("neuron|glia|astrocyte|oligodendrocyte|microglial|radial|neural|ependymal", text) ||
      has_any(c("SNAP25", "SLC17A7", "GAD1", "GFAP", "MBP", "PDGFRA", "SOX2"))) return("neural")
  if (grepl("endothelial", text) || has_any(c("PECAM1", "VWF", "KDR", "CLDN5"))) return("endothelial")
  if (grepl("fibroblast|stellate|smooth muscle|pericyte|mesangial|stromal", text) ||
      has_any(c("COL1A1", "DCN", "ACTA2", "PDGFRB", "RGS5"))) return("stromal")
  if (grepl("unknown|ambiguous|unidentified", text)) return("unknown")
  "other"
}

profile_program_terms <- function(label, tissue) {
  profiles <- .get_marker_evidence_profiles(tissue)
  idx <- .match_profile_name(label, names(profiles))
  if (length(idx) == 0) {
    return(character())
  }
  names(profiles)[idx]
}

program_consistency <- function(label, markers, tissue) {
  terms <- profile_program_terms(label, tissue)
  if (length(terms) == 0) {
    return(list(score = NA_real_, terms = character()))
  }
  profiles <- .get_marker_evidence_profiles(tissue)
  scores <- vapply(terms, function(term) {
    .marker_profile_score(profiles[[term]], markers)$score
  }, numeric(1))
  list(score = max(scores, na.rm = TRUE), terms = terms[scores == max(scores, na.rm = TRUE)])
}

go_available <- function(species) {
  enable_go <- tolower(Sys.getenv("DEEPSEEKCELL_ENABLE_GO_ENRICHMENT", "false")) %in%
    c("1", "true", "yes", "y")
  if (!isTRUE(enable_go)) {
    return(FALSE)
  }
  species <- tolower(as.character(species %||% ""))
  key <- cache_key(species)
  if (exists(key, envir = .bio_validation_cache$go_available, inherits = FALSE)) {
    return(get(key, envir = .bio_validation_cache$go_available, inherits = FALSE))
  }
  annotation_package <- if (grepl("mouse", species)) "org.Mm.eg.db" else "org.Hs.eg.db"
  out <- all(vapply(c("clusterProfiler", "AnnotationDbi", annotation_package), requireNamespace, logical(1), quietly = TRUE))
  assign(key, out, envir = .bio_validation_cache$go_available)
  out
}

optional_go_terms <- function(genes, species, top_n = 10) {
  genes <- unique(toupper(as.character(genes)))
  genes <- genes[nzchar(genes) & !is.na(genes)]
  key <- cache_key(tolower(as.character(species %||% "")), paste(sort(genes), collapse = ";"), top_n)
  if (exists(key, envir = .bio_validation_cache$go_terms, inherits = FALSE)) {
    return(get(key, envir = .bio_validation_cache$go_terms, inherits = FALSE))
  }
  if (length(genes) < 3 || !go_available(species)) {
    assign(key, character(), envir = .bio_validation_cache$go_terms)
    return(character())
  }

  annotation_package <- if (grepl("mouse", tolower(species))) "org.Mm.eg.db" else "org.Hs.eg.db"
  orgdb <- getExportedValue(annotation_package, annotation_package)
  entrez <- tryCatch(
    AnnotationDbi::mapIds(
      orgdb,
      keys = genes,
      keytype = "SYMBOL",
      column = "ENTREZID",
      multiVals = "first"
    ),
    error = function(e) NULL
  )
  entrez <- unique(as.character(stats::na.omit(entrez)))
  if (length(entrez) < 3) {
    assign(key, character(), envir = .bio_validation_cache$go_terms)
    return(character())
  }

  enrich <- tryCatch(
    clusterProfiler::enrichGO(
      gene = entrez,
      OrgDb = orgdb,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      readable = FALSE,
      qvalueCutoff = 0.2
    ),
    error = function(e) NULL
  )
  if (is.null(enrich)) {
    assign(key, character(), envir = .bio_validation_cache$go_terms)
    return(character())
  }
  tab <- as.data.frame(enrich)
  if (!is.data.frame(tab) || nrow(tab) == 0 || !"ID" %in% names(tab)) {
    assign(key, character(), envir = .bio_validation_cache$go_terms)
    return(character())
  }
  out <- head(as.character(tab$ID), top_n)
  assign(key, out, envir = .bio_validation_cache$go_terms)
  out
}

go_consistency <- function(cluster_markers, first_label, final_label, truth_label, species, tissue) {
  cluster_terms <- optional_go_terms(cluster_markers, species)
  if (length(cluster_terms) == 0) {
    return(list(
      available = FALSE,
      first_score = NA_real_,
      final_score = NA_real_,
      truth_score = NA_real_,
      final_delta = NA_real_
    ))
  }

  profiles <- .get_marker_evidence_profiles(tissue)
  profile_genes <- function(label) {
    idx <- .match_profile_name(label, names(profiles))
    unique(unlist(profiles[idx], use.names = FALSE))
  }
  score_label <- function(label) {
    terms <- optional_go_terms(profile_genes(label), species)
    if (length(terms) == 0) return(NA_real_)
    length(intersect(cluster_terms, terms)) / length(union(cluster_terms, terms))
  }
  first_score <- score_label(first_label)
  final_score <- score_label(final_label)
  truth_score <- score_label(truth_label)
  list(
    available = TRUE,
    first_score = first_score,
    final_score = final_score,
    truth_score = truth_score,
    final_delta = final_score - first_score
  )
}

classify_failure <- function(row) {
  if (isTRUE(row$WrongToCorrect)) return("corrected")
  if (isTRUE(row$CorrectToWrong)) return("harmful_revision")
  if (!isTRUE(row$InitiallyCorrect) && !isTRUE(row$FinallyCorrect)) {
    if (!is.na(row$FinalMarkerScore) && row$FinalMarkerScore < 0.2) return("insufficient_canonical_marker_support")
    if (!is.na(row$FinalOntologyScore) && row$FinalOntologyScore == 0) return("ontology_mismatch_or_unmapped_label")
    if (!identical(row$FinalLineage, row$TruthLineage)) return("lineage_confusion")
    if (!is.na(row$ClusterPurity) && row$ClusterPurity < 0.75) return("low_purity_or_mixed_cluster")
    return("uncorrected_related_or_granularity_error")
  }
  if (isTRUE(row$LabelChanged) && isTRUE(row$FinallyCorrect)) return("changed_but_still_correct")
  if (isTRUE(row$LabelChanged)) return("neutral_label_change")
  "unchanged"
}

build_validation_rows <- function(debug, datasets, ontology) {
  rows <- lapply(seq_len(nrow(debug)), function(i) {
    row <- debug[i, , drop = FALSE]
    dataset <- datasets[[as.character(row$Dataset)]]
    if (is.null(dataset)) {
      return(NULL)
    }

    tissue <- as.character(row$Tissue %||% dataset$tissue %||% NA_character_)
    species <- as.character(row$Species %||% dataset$species %||% NA_character_)
    markers <- cluster_markers(dataset, row$Cluster)
    first_label <- as.character(row$FirstPassCellType %||% row$RawPrediction)
    final_label <- as.character(row$RawPrediction)
    truth_label <- as.character(row$RawTruth)
    first_eval <- as.character(row$FirstPassCellType %||% row$RawPrediction)
    final_eval <- as.character(row$HarmonisedPrediction %||% row$RawPrediction)
    truth_eval <- as.character(row$HarmonisedTruth %||% row$RawTruth)

    first_marker <- marker_support(first_label, markers, tissue)
    final_marker <- marker_support(final_label, markers, tissue)
    truth_marker <- marker_support(truth_label, markers, tissue)
    first_ont <- ontology_distance_to_truth(first_label, truth_label, ontology, tissue)
    final_ont <- ontology_distance_to_truth(final_label, truth_label, ontology, tissue)
    first_program <- program_consistency(first_label, markers, tissue)
    final_program <- program_consistency(final_label, markers, tissue)
    truth_program <- program_consistency(truth_label, markers, tissue)
    go <- go_consistency(markers, first_label, final_label, truth_label, species, tissue)

    first_lineage <- lineage_from_label(first_label, markers)
    final_lineage <- lineage_from_label(final_label, markers)
    truth_lineage <- lineage_from_label(truth_label, markers)
    initially_correct <- label_exact(first_eval, truth_eval)
    finally_correct <- label_exact(final_eval, truth_eval)
    label_changed <- !label_exact(first_label, final_label)

    out <- data.frame(
      Replicate = as_number(row$Replicate %||% 1),
      Dataset = as.character(row$Dataset),
      Tissue = tissue,
      Species = species,
      Method = as.character(row$Method),
      RefinementSelector = as.character(row$RefinementSelector %||% ""),
      Cluster = as.character(row$Cluster),
      ClusterPurity = as_number(row$ClusterPurity %||% NA_real_),
      NMarkers = length(markers),
      MarkerGenes = paste(head(markers, 50), collapse = ";"),
      FirstPassLabel = first_label,
      FinalLabel = final_label,
      TruthLabel = truth_label,
      FirstPassEvalLabel = first_eval,
      FinalEvalLabel = final_eval,
      TruthEvalLabel = truth_eval,
      InitiallyCorrect = initially_correct,
      FinallyCorrect = finally_correct,
      LabelChanged = label_changed,
      WrongToCorrect = !initially_correct && finally_correct,
      CorrectToWrong = initially_correct && !finally_correct,
      WasRefined = as_bool(row$WasRefined %||% row$SelectedForRefinement %||% FALSE),
      SelectedForRefinement = as_bool(row$SelectedForRefinement %||% FALSE),
      FirstMarkerScore = first_marker$score,
      FinalMarkerScore = final_marker$score,
      TruthMarkerScore = truth_marker$score,
      MarkerScoreDelta = final_marker$score - first_marker$score,
      FinalVsTruthMarkerGap = abs(final_marker$score - truth_marker$score),
      FirstBestMarkerLabel = first_marker$best_label,
      FinalBestMarkerLabel = final_marker$best_label,
      TruthBestMarkerLabel = truth_marker$best_label,
      FinalMatchedMarkers = paste(final_marker$matched_markers, collapse = ";"),
      TruthMatchedMarkers = paste(truth_marker$matched_markers, collapse = ";"),
      BestMatchedMarkers = paste(final_marker$best_matched_markers, collapse = ";"),
      FirstCLID = first_ont$cl_id,
      FinalCLID = final_ont$cl_id,
      TruthCLID = final_ont$truth_cl_id,
      FirstOntologyRelation = first_ont$relation,
      FinalOntologyRelation = final_ont$relation,
      FirstOntologyScore = first_ont$score,
      FinalOntologyScore = final_ont$score,
      OntologyScoreDelta = final_ont$score - first_ont$score,
      FinalOntologyLabel = final_ont$ontology_label,
      FinalOntologyMatchMethod = final_ont$match_method,
      FirstLineage = first_lineage,
      FinalLineage = final_lineage,
      TruthLineage = truth_lineage,
      FirstLineageMatchesTruth = identical(first_lineage, truth_lineage),
      FinalLineageMatchesTruth = identical(final_lineage, truth_lineage),
      LineageImproved = !identical(first_lineage, truth_lineage) && identical(final_lineage, truth_lineage),
      FirstProgramScore = first_program$score,
      FinalProgramScore = final_program$score,
      TruthProgramScore = truth_program$score,
      ProgramScoreDelta = final_program$score - first_program$score,
      FinalProgramTerms = paste(final_program$terms, collapse = ";"),
      PathwayEnrichmentAvailable = go$available,
      FirstGOConsistency = go$first_score,
      FinalGOConsistency = go$final_score,
      TruthGOConsistency = go$truth_score,
      GOConsistencyDelta = go$final_delta,
      SourceDebugFile = as.character(row$SourceDebugFile),
      stringsAsFactors = FALSE
    )
    out$FailureCategory <- classify_failure(out[1, , drop = FALSE])
    out
  })
  do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
}

summarise_biological_validation <- function(rows) {
  if (nrow(rows) == 0) return(data.frame())
  aggregate_rows <- split(rows, list(rows$Method, rows$RefinementSelector), drop = TRUE)
  out <- lapply(aggregate_rows, function(x) {
    data.frame(
      Method = x$Method[[1]],
      RefinementSelector = x$RefinementSelector[[1]],
      NClusters = nrow(x),
      NRefined = sum(x$WasRefined, na.rm = TRUE),
      NLabelChanged = sum(x$LabelChanged, na.rm = TRUE),
      WrongToCorrect = sum(x$WrongToCorrect, na.rm = TRUE),
      CorrectToWrong = sum(x$CorrectToWrong, na.rm = TRUE),
      MeanMarkerScoreDeltaChanged = mean(x$MarkerScoreDelta[x$LabelChanged], na.rm = TRUE),
      MedianMarkerScoreDeltaChanged = stats::median(x$MarkerScoreDelta[x$LabelChanged], na.rm = TRUE),
      MarkerImprovedAmongChanged = mean(x$MarkerScoreDelta[x$LabelChanged] > 0, na.rm = TRUE),
      MeanOntologyScoreDeltaChanged = mean(x$OntologyScoreDelta[x$LabelChanged], na.rm = TRUE),
      OntologyImprovedAmongChanged = mean(x$OntologyScoreDelta[x$LabelChanged] > 0, na.rm = TRUE),
      LineageImprovedAmongChanged = mean(x$LineageImproved[x$LabelChanged], na.rm = TRUE),
      FinalLineageConsistency = mean(x$FinalLineageMatchesTruth, na.rm = TRUE),
      ProgramImprovedAmongChanged = mean(x$ProgramScoreDelta[x$LabelChanged] > 0, na.rm = TRUE),
      GOAvailableClusters = sum(x$PathwayEnrichmentAvailable, na.rm = TRUE),
      MeanGOConsistencyDeltaChanged = mean(x$GOConsistencyDelta[x$LabelChanged], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

case_studies <- function(rows, n = 50) {
  x <- rows[rows$WrongToCorrect | rows$CorrectToWrong | rows$LabelChanged, , drop = FALSE]
  if (nrow(x) == 0) return(x)
  x$CasePriority <- with(
    x,
    2 * as.numeric(WrongToCorrect) -
      2 * as.numeric(CorrectToWrong) +
      as_number(MarkerScoreDelta, 0) +
      as_number(OntologyScoreDelta, 0) +
      as.numeric(LineageImproved)
  )
  x <- x[order(x$WrongToCorrect, x$CasePriority, decreasing = TRUE), , drop = FALSE]
  head(x, n)
}

failure_taxonomy <- function(rows) {
  if (nrow(rows) == 0) return(data.frame())
  tab <- as.data.frame(table(rows$Method, rows$FailureCategory), stringsAsFactors = FALSE)
  names(tab) <- c("Method", "FailureCategory", "NClusters")
  totals <- stats::aggregate(NClusters ~ Method, tab, sum)
  names(totals)[2] <- "MethodTotal"
  tab <- merge(tab, totals, by = "Method", all.x = TRUE)
  tab$Fraction <- safe_divide(tab$NClusters, tab$MethodTotal)
  tab[order(tab$Method, -tab$NClusters), , drop = FALSE]
}

plot_biological_outputs <- function(prefix, rows, summary) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || nrow(rows) == 0) {
    return(invisible(FALSE))
  }
  changed <- rows[rows$LabelChanged, , drop = FALSE]
  if (nrow(changed) > 0) {
    p <- ggplot2::ggplot(
      changed,
      ggplot2::aes(x = .data$Method, y = .data$MarkerScoreDelta, fill = .data$WrongToCorrect)
    ) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
      ggplot2::geom_boxplot(outlier.alpha = 0.4) +
      ggplot2::coord_flip() +
      ggplot2::labs(x = NULL, y = "Final - first-pass canonical marker support", fill = "Corrected") +
      ggplot2::theme_minimal(base_size = 11)
    ggplot2::ggsave(paste0(prefix, "_marker_recovery_delta.pdf"), p, width = 8, height = 6)

    q <- ggplot2::ggplot(
      changed,
      ggplot2::aes(x = .data$OntologyScoreDelta, y = .data$MarkerScoreDelta, color = .data$FailureCategory)
    ) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
      ggplot2::geom_point(alpha = 0.8) +
      ggplot2::labs(x = "Ontology score delta", y = "Marker support delta", color = "Outcome") +
      ggplot2::theme_minimal(base_size = 11)
    ggplot2::ggsave(paste0(prefix, "_ontology_marker_delta.pdf"), q, width = 7, height = 5)
  }

  if (nrow(summary) > 0) {
    r <- ggplot2::ggplot(
      summary,
      ggplot2::aes(x = .data$Method, y = .data$FinalLineageConsistency)
    ) +
      ggplot2::geom_col(fill = "#4B5563") +
      ggplot2::coord_flip() +
      ggplot2::labs(x = NULL, y = "Final lineage consistency with truth") +
      ggplot2::theme_minimal(base_size = 11)
    ggplot2::ggsave(paste0(prefix, "_lineage_consistency.pdf"), r, width = 8, height = 6)
  }
  invisible(TRUE)
}

run_biological_validation <- function(debug_dir,
                                      manifest_path,
                                      output_prefix,
                                      ontology_path = "data/cl.obo") {
  datasets <- load_manifest_datasets(manifest_path)
  ontology <- load_ontology_for_biology(ontology_path)
  debug <- read_debug_tables(debug_dir, datasets)
  rows <- build_validation_rows(debug, datasets, ontology)
  summary <- summarise_biological_validation(rows)
  corrected <- rows[rows$WrongToCorrect, , drop = FALSE]
  cases <- case_studies(rows)
  taxonomy <- failure_taxonomy(rows)

  dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(rows, paste0(output_prefix, "_cluster_biological_validation.csv"), row.names = FALSE)
  utils::write.csv(corrected, paste0(output_prefix, "_corrected_clusters.csv"), row.names = FALSE)
  utils::write.csv(summary, paste0(output_prefix, "_summary_by_method.csv"), row.names = FALSE)
  utils::write.csv(cases, paste0(output_prefix, "_case_studies.csv"), row.names = FALSE)
  utils::write.csv(taxonomy, paste0(output_prefix, "_failure_taxonomy.csv"), row.names = FALSE)
  plot_biological_outputs(output_prefix, rows, summary)

  message("Wrote biological validation tables with ", nrow(rows), " cluster-method rows.")
  invisible(list(
    cluster_validation = rows,
    corrected = corrected,
    summary = summary,
    case_studies = cases,
    failure_taxonomy = taxonomy
  ))
}

args <- commandArgs(trailingOnly = TRUE)
debug_dir <- if (length(args) >= 1) args[[1]] else file.path("results", "benchmark_debug")
manifest_path <- if (length(args) >= 2) args[[2]] else file.path("benchmarks", "external_validation_manifest.csv")
output_prefix <- if (length(args) >= 3) args[[3]] else file.path("results", "biological_validation")

run_biological_validation(debug_dir, manifest_path, output_prefix)
