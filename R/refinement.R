# R/refinement.R
# Deterministic evidence scoring, confidence calibration, and self-refinement
# prompts for ontology-guided annotation.

MARKER_EVIDENCE_PROFILES <- list(
  generic = list(
    "T cell" = c("CD3D", "CD3E", "CD2", "TRAC", "TRBC1", "TRBC2"),
    "CD4 T cell" = c("CD3D", "CD3E", "CD4", "IL7R", "CCR7", "LTB"),
    "CD8 T cell" = c("CD3D", "CD3E", "CD8A", "CD8B", "GZMK", "NKG7"),
    "B cell" = c("MS4A1", "CD79A", "CD79B", "CD74", "BANK1", "CD19"),
    "plasma cell" = c("MZB1", "JCHAIN", "XBP1", "SDC1", "IGHG1", "IGKC"),
    "natural killer cell" = c("NKG7", "GNLY", "KLRD1", "KLRF1", "PRF1", "GZMB"),
    "monocyte" = c("LYZ", "LST1", "S100A8", "S100A9", "FCN1", "CTSS"),
    "macrophage" = c("C1QA", "C1QB", "C1QC", "CD68", "MARCO", "MRC1"),
    "dendritic cell" = c("FCER1A", "CLEC10A", "CD1C", "ITGAX", "LILRA4", "IRF8"),
    "endothelial cell" = c("PECAM1", "VWF", "KDR", "CLDN5", "FLT1", "EMCN"),
    "fibroblast" = c("COL1A1", "COL1A2", "DCN", "LUM", "PDGFRA", "THY1"),
    "smooth muscle cell" = c("ACTA2", "TAGLN", "MYH11", "CNN1", "MYL9", "TPM2"),
    "pericyte" = c("RGS5", "PDGFRB", "MCAM", "CSPG4", "ACTA2", "NOTCH3"),
    "mast cell" = c("TPSAB1", "TPSB2", "CPA3", "KIT", "MS4A2", "HDC")
  ),
  pancreas = list(
    "alpha cell" = c("GCG", "ARX", "SLC38A4", "TTR", "GC", "IRX2"),
    "beta cell" = c("INS", "IAPP", "MAFA", "PDX1", "NKX6-1", "SLC30A8"),
    "delta cell" = c("SST", "HHEX", "RBP4", "PCSK2", "LEPR", "RGS2"),
    "gamma cell" = c("PPY", "PNOC", "PAX6", "PCSK1", "CHGB"),
    "epsilon cell" = c("GHRL", "MLN", "PCSK1", "CHGA"),
    "acinar cell" = c("PRSS1", "CPA1", "CTRB1", "AMY2A", "CELA3A", "REG1A"),
    "ductal cell" = c("KRT19", "SOX9", "HNF1B", "TFF2", "MUC1", "SLC4A4"),
    "pancreatic stellate cell" = c("COL1A1", "COL3A1", "DCN", "LUM", "RGS5", "ACTA2")
  ),
  brain = list(
    "neuron" = c("SNAP25", "SYT1", "RBFOX3", "TUBB3", "MAP2", "STMN2"),
    "glutamatergic neuron" = c("SLC17A7", "SLC17A6", "CAMK2A", "SATB2", "TBR1"),
    "GABAergic neuron" = c("GAD1", "GAD2", "SLC32A1", "DLX1", "DLX2"),
    "astrocyte" = c("GFAP", "AQP4", "ALDH1L1", "SLC1A3", "GJA1", "S100B"),
    "oligodendrocyte" = c("MBP", "MOG", "PLP1", "MOBP", "MAG", "CLDN11"),
    "oligodendrocyte precursor cell" = c("PDGFRA", "CSPG4", "VCAN", "SOX10", "OLIG1", "OLIG2"),
    "microglial cell" = c("CX3CR1", "P2RY12", "TMEM119", "AIF1", "CSF1R", "TYROBP"),
    "ependymal cell" = c("FOXJ1", "PIFO", "TPPP3", "DNAH5", "CFAP43")
  ),
  lung = list(
    "alveolar macrophage" = c("MARCO", "PPARG", "FABP4", "MRC1", "C1QA", "C1QB"),
    "AT1 cell" = c("AGER", "PDPN", "CAV1", "HOPX", "RTKN2", "AQP5"),
    "AT2 cell" = c("SFTPC", "SFTPA1", "SFTPA2", "SFTPB", "ABCA3", "LAMP3"),
    "club cell" = c("SCGB1A1", "SCGB3A2", "CYP2F1", "BPIFB1", "MUC1"),
    "ciliated cell" = c("FOXJ1", "PIFO", "TPPP3", "DNAH5", "RSPH1", "CFAP43"),
    "basal cell" = c("KRT5", "KRT15", "TP63", "KRT17", "NGFR"),
    "goblet cell" = c("MUC5AC", "MUC5B", "SPDEF", "AGR2", "TFF3"),
    "ionocyte" = c("FOXI1", "CFTR", "ASCL3", "ATP6V1B1", "ATP6V0D2"),
    "lung endothelial cell" = c("PECAM1", "VWF", "CLDN5", "EMCN", "CA4", "RGCC"),
    "lung fibroblast" = c("COL1A1", "COL1A2", "DCN", "LUM", "PDGFRA", "COL3A1")
  ),
  prostate = list(
    "prostate epithelial cell" = c("KLK3", "KLK2", "MSMB", "ACPP", "KRT8", "KRT18"),
    "basal epithelial cell" = c("KRT5", "KRT14", "TP63", "KRT15", "DST"),
    "luminal epithelial cell" = c("KLK3", "KLK2", "NKX3-1", "KRT8", "KRT18", "MSMB"),
    "lymphatic endothelial cell" = c("PROX1", "LYVE1", "FLT4", "PDPN", "CCL21", "MMRN1")
  )
)

MARKER_PROFILE_ALIASES <- list(
  "T cell" = c("t lymphocyte", "t-cell", "t cells"),
  "CD4 T cell" = c("cd4 positive alpha beta t cell", "cd4-positive t cell"),
  "CD8 T cell" = c("cd8 positive alpha beta t cell", "cytotoxic t cell"),
  "B cell" = c("b lymphocyte", "b-cell", "b cells"),
  "natural killer cell" = c("nk cell", "nk cells"),
  "alpha cell" = c("pancreatic alpha cell", "pancreatic a cell", "islet alpha cell"),
  "beta cell" = c("pancreatic beta cell", "type b pancreatic cell", "islet beta cell"),
  "delta cell" = c("pancreatic delta cell", "pancreatic d cell", "islet delta cell"),
  "acinar cell" = c("pancreatic acinar cell"),
  "ductal cell" = c("pancreatic ductal cell", "duct epithelial cell"),
  "glutamatergic neuron" = c("excitatory neuron"),
  "GABAergic neuron" = c("inhibitory neuron", "gabaergic interneuron"),
  "oligodendrocyte precursor cell" = c("opc", "ng2 cell", "ng2 glia"),
  "microglial cell" = c("microglia"),
  "AT1 cell" = c("alveolar type 1 cell", "type 1 pneumocyte", "ati cell"),
  "AT2 cell" = c("alveolar type 2 cell", "type 2 pneumocyte", "atii cell"),
  "club cell" = c("clara cell"),
  "lung fibroblast" = c("fibroblast of lung", "pulmonary interstitial fibroblast"),
  "lung endothelial cell" = c("endothelial cell of lung", "pulmonary endothelial cell")
)

#' Score annotation evidence for marker-driven cell type calls
#'
#' @param annotations Annotation data frame.
#' @param markers Named marker list used for annotation.
#' @param tissue Tissue context.
#' @param use_ontology_evidence Whether Cell Ontology mapping should contribute
#' to evidence scoring and conflict detection. When `FALSE`, ontology evidence is
#' neutralized and ontology-specific conflict rules are disabled.
#' @return Data frame with deterministic evidence components.
#' @export
score_annotation_evidence <- function(annotations,
                                      markers,
                                      tissue = NULL,
                                      use_ontology_evidence = TRUE) {
  if (!is.data.frame(annotations) || nrow(annotations) == 0) {
    return(data.frame())
  }

  profiles <- .get_marker_evidence_profiles(tissue)

  rows <- lapply(seq_len(nrow(annotations)), function(i) {
    row <- annotations[i, , drop = FALSE]
    cluster <- as.character(row$Cluster %||% paste0("Cluster", i))
    cluster_markers <- markers[[cluster]] %||% character()
    marker_match <- .score_marker_profiles(cluster_markers, profiles, row$CellType)
    ontology_score <- if (isTRUE(use_ontology_evidence)) {
      .score_ontology_component(row)
    } else {
      0.5
    }
    tissue_score <- .score_tissue_component(row$TissueConsistency %||% "unknown")
    llm_confidence <- as_confidence(row$Confidence %||% 0.5)

    marker_component <- marker_match$predicted_score
    if (marker_match$best_score < 0.2) {
      marker_component <- 0.5
    }

    conflict <- .has_evidence_conflict(
      predicted_label = row$CellType,
      best_label = marker_match$best_cell_type,
      best_score = marker_match$best_score,
      predicted_score = marker_match$predicted_score,
      ontology_score = ontology_score,
      use_ontology_evidence = isTRUE(use_ontology_evidence)
    )

    consensus_score <- if (isTRUE(conflict)) 0.25 else 1
    mixed_factor <- if (isTRUE(row$IsMixed %||% FALSE)) 0.9 else 1

    calibrated <- (
      0.35 * llm_confidence +
        0.25 * ontology_score +
        0.25 * marker_component +
        0.10 * tissue_score +
        0.05 * consensus_score
    ) * mixed_factor

    calibrated <- pmin(pmax(calibrated, 0), 1)

    data.frame(
      Cluster = cluster,
      LLMConfidence = llm_confidence,
      OntologyEvidenceScore = ontology_score,
      MarkerEvidenceScore = marker_component,
      BestMarkerEvidenceScore = marker_match$best_score,
      TissueEvidenceScore = tissue_score,
      ConsensusEvidenceScore = consensus_score,
      EvidenceBestCellType = marker_match$best_cell_type,
      EvidenceMatchedMarkers = paste(marker_match$best_markers, collapse = ", "),
      EvidenceConflict = isTRUE(conflict),
      RequiresRefinement = isTRUE(conflict) && marker_match$best_score >= 0.45,
      CalibratedConfidence = round(calibrated, 3),
      CalibrationDelta = round(calibrated - llm_confidence, 3),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

#' Attach calibrated confidence and evidence columns to annotations
#'
#' The original model confidence is retained as `LLMConfidence`, while
#' `Confidence` is replaced by a deterministic calibrated confidence.
#'
#' @param annotations Annotation data frame.
#' @param markers Named marker list used for annotation.
#' @param tissue Tissue context.
#' @param use_ontology_evidence Whether ontology evidence contributes to the
#' composite confidence and conflict-detection fields.
#' @return Annotation data frame with evidence and calibrated confidence.
#' @export
calibrate_annotation_confidence <- function(annotations,
                                           markers,
                                           tissue = NULL,
                                           use_ontology_evidence = TRUE) {
  evidence <- score_annotation_evidence(
    annotations,
    markers,
    tissue,
    use_ontology_evidence = use_ontology_evidence
  )
  if (nrow(evidence) == 0) {
    return(annotations)
  }

  keep <- setdiff(
    names(evidence),
    c("Cluster", "CalibratedConfidence")
  )

  annotations$LLMConfidence <- evidence$LLMConfidence
  annotations$Confidence <- evidence$CalibratedConfidence

  for (column in keep) {
    annotations[[column]] <- evidence[[column]]
  }

  annotations$ConfidenceMethod <- if (isTRUE(use_ontology_evidence)) {
    "ontology_marker_calibrated"
  } else {
    "marker_tissue_calibrated_no_ontology"
  }
  annotations
}

#' Create an ontology-guided self-refinement prompt
#'
#' @param markers Named marker list.
#' @param annotations Annotation data frame with evidence columns.
#' @param tissue Tissue context.
#' @param species Species name.
#' @return Character prompt.
#' @export
create_refinement_prompt <- function(markers, annotations, tissue, species = "Human") {
  if (!is.data.frame(annotations) || nrow(annotations) == 0) {
    stop("annotations must contain at least one row.", call. = FALSE)
  }

  cluster_text <- paste(
    vapply(seq_len(nrow(annotations)), function(i) {
      row <- annotations[i, , drop = FALSE]
      cluster <- as.character(row$Cluster)
      marker_text <- paste(markers[[cluster]] %||% character(), collapse = ", ")
      sprintf(
        paste(
          "Cluster: %s",
          "Marker genes: %s",
          "First-pass cell type: %s",
          "First-pass confidence: %.3f",
          "Cell Ontology mapping: %s / %s / %s",
          "Marker-evidence best cell type: %s",
          "Matched marker genes: %s",
          "Evidence conflict: %s",
          "First-pass reasoning: %s",
          sep = "\n"
        ),
        cluster,
        marker_text,
        row$CellType %||% "Unknown",
        as_confidence(row$LLMConfidence %||% row$Confidence %||% 0.5),
        row$CL_ID %||% "",
        row$OntologyLabel %||% "",
        row$MatchMethod %||% "",
        row$EvidenceBestCellType %||% "",
        row$EvidenceMatchedMarkers %||% "",
        row$EvidenceConflict %||% FALSE,
        row$Reasoning %||% ""
      )
    }, character(1)),
    collapse = "\n\n---\n\n"
  )

  sprintf(
    paste(
      "Species: %s",
      "Tissue: %s",
      "",
      "You are performing ontology-guided self-refinement of single-cell cluster annotations.",
      "Review only the clusters listed below. The first-pass LLM label, Cell Ontology mapping,",
      "marker-evidence score, and tissue context may disagree.",
      "",
      "Decision rules:",
      "1. Keep the first-pass label if marker evidence and ontology context support it.",
      "2. Revise the label only when marker evidence strongly supports a better biological identity.",
      "3. Prefer Cell Ontology-compatible, publication-ready cell type names.",
      "4. Do not invent Cell Ontology identifiers.",
      "5. Return only valid JSON with the same annotation schema.",
      "",
      "Clusters requiring review:",
      "%s",
      "",
      "Required JSON schema:",
      "{\"annotations\":[{\"cluster\":\"Cluster1\",\"cell_type\":\"beta cell\",\"confidence\":0.95,",
      "\"is_mixed\":false,\"primary_cell_type\":\"beta cell\",\"secondary_cell_type\":\"\",",
      "\"tissue_consistency\":\"expected\",\"reasoning\":\"revised because INS and IAPP support beta cell identity\"}]}",
      sep = "\n"
    ),
    species,
    tissue,
    cluster_text
  )
}

.get_marker_evidence_profiles <- function(tissue = NULL) {
  tissue_key <- if (is.null(tissue)) "" else normalize_cell_type(tissue[1])
  profiles <- MARKER_EVIDENCE_PROFILES$generic

  for (name in names(MARKER_EVIDENCE_PROFILES)) {
    if (identical(name, "generic")) next
    if (grepl(name, tissue_key, fixed = TRUE)) {
      profiles <- c(profiles, MARKER_EVIDENCE_PROFILES[[name]])
    }
  }

  profiles[!duplicated(names(profiles))]
}

.score_marker_profiles <- function(markers, profiles, predicted_label) {
  if (length(profiles) == 0) {
    return(.empty_marker_match())
  }

  marker_scores <- lapply(profiles, .marker_profile_score, markers = markers)
  scores <- vapply(marker_scores, `[[`, numeric(1), "score")
  best_idx <- which.max(scores)

  predicted_idx <- .match_profile_name(predicted_label, names(profiles))
  predicted_score <- if (length(predicted_idx) > 0) max(scores[predicted_idx]) else 0

  best_markers <- marker_scores[[best_idx]]$matched_markers

  list(
    best_cell_type = names(profiles)[best_idx],
    best_score = unname(scores[best_idx]),
    predicted_score = unname(predicted_score),
    best_markers = best_markers
  )
}

.marker_profile_score <- function(profile_genes, markers) {
  marker_genes <- .normalize_gene_symbols(markers)
  profile_genes <- .normalize_gene_symbols(profile_genes)
  matched <- intersect(profile_genes, marker_genes)

  if (length(marker_genes) == 0 || length(profile_genes) == 0) {
    score <- 0
  } else {
    denominator <- min(length(profile_genes), max(3, length(marker_genes)))
    score <- length(matched) / denominator
  }

  list(
    score = pmin(score, 1),
    matched_markers = matched
  )
}

.match_profile_name <- function(label, profile_names) {
  label_norm <- normalize_cell_type(label %||% "")
  if (!nzchar(label_norm)) {
    return(integer())
  }

  profile_norm <- normalize_cell_type(profile_names)
  profile_in_label <- vapply(
    profile_norm,
    function(pattern) grepl(pattern, label_norm, fixed = TRUE),
    logical(1)
  )
  label_in_profile <- vapply(
    profile_norm,
    function(value) grepl(label_norm, value, fixed = TRUE),
    logical(1)
  )
  direct <- which(profile_norm == label_norm | profile_in_label | label_in_profile)

  alias_hits <- unlist(lapply(names(MARKER_PROFILE_ALIASES), function(profile) {
    aliases <- normalize_cell_type(c(profile, MARKER_PROFILE_ALIASES[[profile]]))
    if (label_norm %in% aliases) {
      which(profile_norm == normalize_cell_type(profile))
    } else {
      integer()
    }
  }))

  unique(c(direct, alias_hits))
}

.has_evidence_conflict <- function(predicted_label,
                                   best_label,
                                   best_score,
                                   predicted_score,
                                   ontology_score,
                                   use_ontology_evidence = TRUE) {
  if (is.null(best_label) || is.na(best_label) || !nzchar(best_label)) {
    return(FALSE)
  }

  predicted_matches_best <- length(.match_profile_name(predicted_label, best_label)) > 0
  strong_marker_disagreement <- best_score >= 0.5 &&
    predicted_score <= max(0.1, best_score - 0.25) &&
    !predicted_matches_best

  unsupported_ontology <- isTRUE(use_ontology_evidence) &&
    best_score >= 0.6 &&
    ontology_score < 0.45

  isTRUE(strong_marker_disagreement || unsupported_ontology)
}

.score_ontology_component <- function(row) {
  if (!"CL_ID" %in% names(row) || is.na(row$CL_ID) || !nzchar(row$CL_ID)) {
    return(0.5)
  }

  score <- suppressWarnings(as.numeric(row$OntologyMatchScore %||% NA_real_))
  if (!is.na(score)) {
    return(pmin(pmax(score, 0), 1))
  }

  method <- tolower(as.character(row$MatchMethod %||% ""))
  if (grepl("^context_exact|^exact", method)) return(1)
  if (grepl("^synonym", method)) return(0.95)
  if (grepl("^fuzzy", method)) return(0.75)
  0.7
}

.score_tissue_component <- function(tissue_consistency) {
  value <- tolower(trimws(as.character(tissue_consistency %||% "unknown")))
  switch(
    value,
    expected = 1,
    possible_contamination = 0.65,
    possible_doublet = 0.55,
    unexpected = 0.35,
    unknown = 0.6,
    0.6
  )
}

.normalize_gene_symbols <- function(x) {
  x <- as.character(x %||% character())
  x <- trimws(x)
  x <- x[nzchar(x) & !is.na(x)]
  toupper(x)
}

.empty_marker_match <- function() {
  list(
    best_cell_type = NA_character_,
    best_score = 0,
    predicted_score = 0,
    best_markers = character()
  )
}
