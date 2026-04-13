# benchmarks/methods.R
# Annotation methods: LLM, SingleR, scType

suppressPackageStartupMessages({
  library(httr2)
  library(jsonlite)
  library(SingleR)
  library(celldex)
  library(SingleCellExperiment)
  library(dplyr)
  library(openxlsx)
  library(HGNChelper)
  library(scales)
  library(Seurat)
})

# -----------------------------------------------------------------------------
# LLM helpers
# -----------------------------------------------------------------------------

call_llm_api <- function(prompt, model_key, api_key) {
  model <- MODELS[[model_key]]
  cat("Calling", model_key, "API...\n")
  
  req <- request(model$api_url) |>
    req_headers(Authorization = paste("Bearer", api_key),
                "Content-Type" = "application/json") |>
    req_body_json(list(
      model = model$model_id,
      messages = list(
        list(role = "system",
             content = "You are a bioinformatics expert. Output ONLY valid JSON. Do not include any other text."),
        list(role = "user", content = prompt)
      ),
      temperature = model$temperature,
      max_tokens = model$max_tokens
    ))
  
  resp <- tryCatch(
    req_perform(req),
    error = function(e) {
      cat("API request error:", e$message, "\n")
      return(NULL)
    }
  )
  
  if (is.null(resp)) return(NULL)
  
  # Check HTTP status
  if (resp_status(resp) != 200) {
    cat("HTTP error:", resp_status(resp), "\n")
    cat("Response body:\n", resp_body_string(resp), "\n")
    return(NULL)
  }
  
  data <- resp_body_json(resp)
  return(data$choices[[1]]$message$content)
}

parse_llm_response <- function(response_text) {
  cleaned <- gsub("^```(?:json)?\\s*\\n?", "", response_text)
  cleaned <- gsub("\\s*```$", "", cleaned)
  first <- regexpr("\\{", cleaned)
  last  <- regexpr("\\}[^}]*$", cleaned, perl = TRUE)
  if (first == -1 || last == -1) return(NULL)
  json_str <- substr(cleaned, first, last + attr(last, "match.length") - 1)
  parsed <- tryCatch(fromJSON(json_str), error = function(e) NULL)
  if (is.null(parsed)) return(NULL)
  if (is.list(parsed) && !is.null(names(parsed))) return(unlist(parsed))
  if (!is.null(parsed$annotations)) {
    df <- as.data.frame(parsed$annotations, stringsAsFactors = FALSE)
    names(df) <- tolower(names(df))
    if ("cluster" %in% names(df) && "cell_type" %in% names(df))
      return(setNames(df$cell_type, df$cluster))
  }
  return(NULL)
}

run_llm_annotation <- function(markers_list, tissue, species, model_key, api_key) {
  clusters_text <- paste(sprintf("%s: %s", names(markers_list),
                                 sapply(markers_list, paste, collapse=", ")),
                         collapse = "\n")
  prompt <- sprintf(
    "Tissue: %s\nSpecies: %s\nMarker genes per cluster:\n%s\n\nAnnotate each cluster with cell type. Output JSON as a dictionary mapping cluster ID to cell type, e.g., {\"0\": \"T cell\", \"1\": \"B cell\"}. Use standard cell type names.",
    tissue, species, clusters_text
  )
  cat("Calling", model_key, "API...\n")
  response <- tryCatch(call_llm_api(prompt, model_key, api_key), error = function(e) NULL)
  if (is.null(response)) {
    cat("No response, using Unknown for all clusters.\n")
    return(setNames(rep("Unknown", length(markers_list)), names(markers_list)))
  }
  cat("Response (first 300 chars):\n", substr(response, 1, 300), "\n")
  pred_dict <- parse_llm_response(response)
  if (is.null(pred_dict)) {
    cat("Parsing failed, using Unknown for all clusters.\n")
    return(setNames(rep("Unknown", length(markers_list)), names(markers_list)))
  }
  expected <- names(markers_list)
  pred <- sapply(expected, function(cl) ifelse(cl %in% names(pred_dict), pred_dict[[cl]], "Unknown"))
  cat("Final predictions:\n"); print(pred)
  return(pred)
}

# -----------------------------------------------------------------------------
# SingleR
# -----------------------------------------------------------------------------

run_singler <- function(seu, dataset_name) {
  message("    Running SingleR for ", dataset_name, "...")
  sce <- as.SingleCellExperiment(seu)
  ref <- NULL
  if (dataset_name == "PBMC") {
    ref <- tryCatch(celldex::MonacoImmuneData(), error = function(e) NULL)
  } else if (dataset_name == "Pancreas") {
    ref <- tryCatch(celldex::BaronPancreasData(), error = function(e) NULL)
    if (is.null(ref)) {
      message("      BaronPancreasData failed, trying HumanPrimaryCellAtlasData...")
      ref <- tryCatch(celldex::HumanPrimaryCellAtlasData(), error = function(e) NULL)
    }
  } else if (dataset_name == "Brain") {
    ref <- tryCatch(celldex::MouseRNAseqData(), error = function(e) NULL)
  }
  if (is.null(ref)) {
    message("      Reference loading failed, returning 'Unknown' for all cells.")
    return(rep("Unknown", ncol(seu)))
  }
  if (!"logcounts" %in% assayNames(ref)) ref <- scuttle::logNormCounts(ref)
  keep <- !is.na(ref$label.main)
  ref <- ref[, keep]
  if (ncol(ref) == 0) {
    message("      Reference has no valid cells after NA removal")
    return(rep("Unknown", ncol(seu)))
  }
  pred <- tryCatch({
    SingleR(test = sce, ref = ref, labels = ref$label.main,
            assay.type.test = "logcounts")
  }, error = function(e) {
    message("      SingleR failed: ", e$message)
    return(NULL)
  })
  if (is.null(pred)) return(rep("Unknown", ncol(seu)))
  pred_labels <- as.character(pred$pruned.labels)
  pred_labels[is.na(pred_labels)] <- "Unknown"
  return(pred_labels)
}

# -----------------------------------------------------------------------------
# scType (custom implementation)
# -----------------------------------------------------------------------------

gene_sets_prepare <- function(path_to_db_file, cell_type) {
  cell_markers <- openxlsx::read.xlsx(path_to_db_file)
  cell_markers <- cell_markers[cell_markers$tissueType == cell_type, ]
  cell_markers$geneSymbolmore1 <- gsub(" ", "", cell_markers$geneSymbolmore1)
  cell_markers$geneSymbolmore2 <- gsub(" ", "", cell_markers$geneSymbolmore2)
  
  cell_markers$geneSymbolmore1 <- sapply(1:nrow(cell_markers), function(i) {
    markers_all <- gsub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore1[i], ",")))
    markers_all <- toupper(markers_all[markers_all != "NA" & markers_all != ""])
    markers_all <- sort(markers_all)
    if (length(markers_all) > 0) {
      suppressMessages({
        markers_all <- unique(na.omit(checkGeneSymbols(markers_all)$Suggested.Symbol))
      })
      paste0(markers_all, collapse = ",")
    } else ""
  })
  cell_markers$geneSymbolmore2 <- sapply(1:nrow(cell_markers), function(i) {
    markers_all <- gsub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore2[i], ",")))
    markers_all <- toupper(markers_all[markers_all != "NA" & markers_all != ""])
    markers_all <- sort(markers_all)
    if (length(markers_all) > 0) {
      suppressMessages({
        markers_all <- unique(na.omit(checkGeneSymbols(markers_all)$Suggested.Symbol))
      })
      paste0(markers_all, collapse = ",")
    } else ""
  })
  cell_markers$geneSymbolmore1 <- gsub("///", ",", cell_markers$geneSymbolmore1)
  cell_markers$geneSymbolmore2 <- gsub("///", ",", cell_markers$geneSymbolmore2)
  gs <- lapply(1:nrow(cell_markers), function(j) {
    gsub(" ", "", unlist(strsplit(toString(cell_markers$geneSymbolmore1[j]), ",")))
  })
  names(gs) <- cell_markers$cellName
  gs2 <- lapply(1:nrow(cell_markers), function(j) {
    gsub(" ", "", unlist(strsplit(toString(cell_markers$geneSymbolmore2[j]), ",")))
  })
  names(gs2) <- cell_markers$cellName
  list(gs_positive = gs, gs_negative = gs2)
}

sctype_score <- function(scRNAseqData, scaled = TRUE, gs, gs2 = NULL) {
  if (!is.matrix(scRNAseqData)) warning("scRNAseqData doesn't seem to be a matrix")
  marker_stat <- sort(table(unlist(gs)), decreasing = TRUE)
  marker_sensitivity <- data.frame(
    score_marker_sensitivity = scales::rescale(as.numeric(marker_stat), to = c(0, 1), from = c(length(gs), 1)),
    gene_ = names(marker_stat),
    stringsAsFactors = FALSE
  )
  rownames(scRNAseqData) <- toupper(rownames(scRNAseqData))
  gs <- lapply(gs, function(x) rownames(scRNAseqData)[rownames(scRNAseqData) %in% x])
  gs2 <- lapply(gs2, function(x) rownames(scRNAseqData)[rownames(scRNAseqData) %in% x])
  cell_markers_genes_score <- marker_sensitivity[marker_sensitivity$gene_ %in% unique(unlist(gs)), ]
  if (!scaled) Z <- t(scale(t(scRNAseqData))) else Z <- scRNAseqData
  for (jj in 1:nrow(cell_markers_genes_score)) {
    Z[cell_markers_genes_score[jj, "gene_"], ] <- Z[cell_markers_genes_score[jj, "gene_"], ] *
      cell_markers_genes_score[jj, "score_marker_sensitivity"]
  }
  Z <- Z[unique(c(unlist(gs), unlist(gs2))), ]
  es <- do.call("rbind", lapply(names(gs), function(gss_) {
    sapply(1:ncol(Z), function(j) {
      gs_z <- Z[gs[[gss_]], j]
      gz_2 <- Z[gs2[[gss_]], j] * -1
      sum_t1 <- sum(gs_z) / sqrt(length(gs_z))
      sum_t2 <- if (length(gz_2) > 0) sum(gz_2) / sqrt(length(gz_2)) else 0
      sum_t1 + sum_t2
    })
  }))
  dimnames(es) <- list(names(gs), colnames(Z))
  es <- es[!apply(is.na(es) | es == "", 1, all), ]
  es
}

run_sctype_custom <- function(seu, tissue, db_file) {
  gs_list <- gene_sets_prepare(db_file, tissue)
  expr_matrix <- GetAssayData(seu, layer = "data")
  expr_matrix <- as.matrix(expr_matrix)
  es_max <- sctype_score(scRNAseqData = expr_matrix, scaled = FALSE,
                         gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
  if (is.null(es_max) || nrow(es_max) == 0 || ncol(es_max) == 0) {
    warning("No scores generated for tissue: ", tissue)
    return(rep("unknown", ncol(seu)))
  }
  if (!is.matrix(es_max)) es_max <- as.matrix(es_max)
  if (!is.numeric(es_max)) storage.mode(es_max) <- "numeric"
  cluster_ids <- as.character(seu$seurat_clusters)
  cL_results <- do.call(rbind, lapply(unique(cluster_ids), function(cl) {
    cells_in_cluster <- which(cluster_ids == cl)
    if (length(cells_in_cluster) == 0) return(NULL)
    scores <- rowMeans(es_max[, cells_in_cluster, drop = FALSE])
    scores <- as.numeric(scores)
    names(scores) <- rownames(es_max)
    best_idx <- which.max(scores)
    data.frame(cluster = cl, type = names(scores)[best_idx],
               score = scores[best_idx], stringsAsFactors = FALSE)
  }))
  if (is.null(cL_results) || nrow(cL_results) == 0) {
    return(rep("unknown", ncol(seu)))
  }
  map_dict <- setNames(cL_results$type, cL_results$cluster)
  result <- unname(map_dict[as.character(seu$seurat_clusters)])
  result[is.na(result)] <- "unknown"
  return(result)
}

find_sctype_db <- function() {
  possible_paths <- c("scType/db.xlsx", "scType/ScTypeDB_full.xlsx")
  for (p in possible_paths) {
    if (file.exists(p)) return(normalizePath(p, mustWork = FALSE))
  }
  stop("scType database not found.")
}