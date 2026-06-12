# benchmarks/methods.R

suppressPackageStartupMessages({
  library(httr2)
  library(jsonlite)
  library(SingleR)
  library(celldex)
  library(SingleCellExperiment)
  library(openxlsx)
  library(HGNChelper)
  library(scales)
  library(Seurat)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

# -----------------------------------------------------------------------------
# LLM API
# -----------------------------------------------------------------------------

call_llm_api <- function(prompt, model_key, api_key, max_retries = 3) {
  model <- MODELS[[model_key]]
  start <- Sys.time()
  last_error <- NULL
  
  for (attempt in seq_len(max_retries)) {
    res <- tryCatch({
      body <- list(
        model = model$model_id,
        messages = list(
          list(
            role = "system",
            content = "You are a bioinformatics expert. Output ONLY valid JSON. Do not include markdown or extra text."
          ),
          list(role = "user", content = prompt)
        ),
        temperature = model$temperature %||% 0,
        max_tokens = model$max_tokens %||% 2000
      )
      
      if (!is.null(model$top_p)) {
        body$top_p <- model$top_p
      }
      
      req <- httr2::request(model$api_url) |>
        httr2::req_timeout(90) |>
        httr2::req_headers(
          Authorization = paste("Bearer", api_key),
          "Content-Type" = "application/json"
        ) |>
        httr2::req_body_json(body)
      
      resp <- httr2::req_perform(req)
      data <- httr2::resp_body_json(resp, simplifyVector = FALSE)
      
      content <- data$choices[[1]]$message$content %||% ""
      usage <- data$usage %||% list(
        prompt_tokens = NA_real_,
        completion_tokens = NA_real_,
        total_tokens = NA_real_
      )
      
      elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
      
      cost <- NA_real_
      if (!is.na(usage$prompt_tokens) && !is.na(usage$completion_tokens)) {
        cost <- (usage$prompt_tokens / 1000) * (model$input_cost_per_1k %||% 0) +
          (usage$completion_tokens / 1000) * (model$output_cost_per_1k %||% 0)
      }
      
      list(
        success = TRUE,
        content = content,
        usage = usage,
        runtime_sec = elapsed,
        cost_usd = cost
      )
    }, error = function(e) {
      last_error <<- conditionMessage(e)
      NULL
    })
    
    if (!is.null(res)) return(res)
    if (attempt < max_retries) Sys.sleep(min(2^attempt, 10))
  }
  
  list(
    success = FALSE,
    content = NULL,
    error = last_error %||% "Unknown API error",
    usage = list(prompt_tokens = 0, completion_tokens = 0, total_tokens = 0),
    runtime_sec = as.numeric(difftime(Sys.time(), start, units = "secs")),
    cost_usd = NA_real_
  )
}

parse_llm_response <- function(response_text) {
  if (is.null(response_text) || !nzchar(trimws(response_text))) return(NULL)
  
  cleaned <- trimws(response_text)
  cleaned <- gsub("^```(?:json)?\\s*", "", cleaned, perl = TRUE)
  cleaned <- gsub("\\s*```$", "", cleaned, perl = TRUE)
  
  first <- regexpr("\\{", cleaned, perl = TRUE)
  if (first[1] == -1) return(NULL)
  
  json_str <- substr(cleaned, first[1], nchar(cleaned))
  
  parsed <- tryCatch(
    jsonlite::fromJSON(json_str, simplifyDataFrame = TRUE),
    error = function(e) NULL
  )
  
  if (is.null(parsed)) return(NULL)
  
  if (!is.null(parsed$annotations)) {
    df <- as.data.frame(parsed$annotations, stringsAsFactors = FALSE)
    names(df) <- tolower(names(df))
    
    if ("cluster" %in% names(df) && "cell_type" %in% names(df)) {
      return(setNames(as.character(df$cell_type), as.character(df$cluster)))
    }
  }
  
  if (is.list(parsed) && !is.data.frame(parsed)) {
    out <- unlist(parsed)
    return(as.character(out))
  }
  
  NULL
}

run_llm_annotation <- function(markers_list,
                               tissue,
                               species,
                               dataset_name,
                               model_key,
                               api_key,
                               save_debug = TRUE) {
  clusters_text <- paste(
    sprintf(
      "%s: %s",
      names(markers_list),
      vapply(markers_list, paste, collapse = ", ", character(1))
    ),
    collapse = "\n"
  )
  
  allowed_labels <- get_allowed_labels(tissue)
  
  label_guide <- switch(
    tolower(tissue),
    pbmc = paste(
      "PBMC label guide:",
      "Naive T cell: CD3D, CD3E, IL7R, CCR7, LTB, LEF1, TCF7. Do NOT use this label if CD8A, NKG7, GNLY, GZMB, PRF1 are dominant.",
      "Cytotoxic T cell: CD8A, CD8B, NKG7, GZMB, GZMH, PRF1, CCL5. Use this label for CD8/NKG7/GZMB-positive T cells.",
      "Classical monocyte: CD14, LST1, LYZ, S100A8, S100A9, FCN1, VCAN.",
      "Non-classical monocyte: FCGR3A, MS4A7, LST1, IFITM3, RHOC.",
      "B cell: MS4A1, CD79A, CD79B, CD74, HLA-DRA, BANK1.",
      "Natural killer cell: GNLY, NKG7, KLRD1, KLRB1, FCGR3A, PRF1. Use NK cell when CD3D/CD3E are absent or weak.",
      "Dendritic cell: FCER1A, CST3, CLEC10A, LILRA4, IRF7.",
      "Platelet: PPBP, PF4, GP9, ITGA2B, TUBB1.",
      sep = "\n"
    ),
    pancreas = paste(
      "Pancreas label guide:",
      "Alpha cell: GCG, TTR, IRX2.",
      "Beta cell: INS, IAPP, MAFA, PDX1.",
      "Delta cell: SST, HHEX.",
      "Gamma cell: PPY.",
      "Acinar cell: PRSS1, CPA1, CTRB1, AMY2A.",
      "Ductal cell: KRT19, SOX9, SPP1.",
      "Endothelial cell: PECAM1, VWF, KDR.",
      "Stellate cell: COL1A1, COL3A1, DCN.",
      "Macrophage: CD68, LYZ, C1QA.",
      sep = "\n"
    ),
    brain = paste(
      "Brain label guide:",
      "Excitatory neuron: SLC17A7, CAMK2A, SNAP25.",
      "Inhibitory neuron: GAD1, GAD2, SLC6A1.",
      "Astrocyte: AQP4, GFAP, ALDH1L1.",
      "Microglia: CX3CR1, C1QA, P2RY12.",
      "Oligodendrocyte: MBP, PLP1, MOG.",
      "Oligodendrocyte precursor cell: PDGFRA, CSPG4.",
      "Endothelial cell: PECAM1, VWF, CLDN5.",
      "Pericyte: PDGFRB, RGS5, CSPG4.",
      "Ependymal cell: FOXJ1, TPPP3.",
      sep = "\n"
    ),
    lung = paste(
      "Lung label guide:",
      "T cell: CD3D, CD3E, CD4, CD8A, IL7R.",
      "B cell: MS4A1, CD79A, CD79B, CD74.",
      "Plasma cell: MZB1, JCHAIN, IGKC, XBP1.",
      "Natural killer cell: GNLY, NKG7, KLRD1, PRF1.",
      "Macrophage: MARCO, APOE, C1QA, LYZ, CD68.",
      "Monocyte: LST1, FCN1, S100A8, S100A9.",
      "Neutrophil: S100A8, S100A9, FCGR3B, CSF3R.",
      "Dendritic cell: FCER1A, CLEC10A, CST3.",
      "Mast cell: TPSAB1, TPSB2, KIT.",
      "Epithelial cell: EPCAM, KRT8, KRT18, KRT19, SFTPC.",
      "Endothelial cell: PECAM1, VWF, CLDN5.",
      "Fibroblast: COL1A1, COL1A2, DCN, LUM.",
      sep = "\n"
    ),
    ""
  )
  
  prompt <- sprintf(
    paste0(
      "You are annotating single-cell RNA-seq clusters.\n\n",
      "Species: %s\n",
      "Tissue: %s\n\n",
      "For each cluster, choose exactly ONE cell type from this allowed label list:\n",
      "%s\n\n",
      "%s\n\n",
      "Rules:\n",
      "1. Use only labels from the allowed list.\n",
      "2. Do not invent new labels.\n",
      "3. Use Unknown only if the marker list is random or contains no recognizable canonical marker genes.\n",
      "4. Use the marker guide to distinguish similar cell types.\n",
      "5. Prefer a biologically plausible allowed label over Unknown.\n",
      "6. Return ONLY valid JSON.\n\n",
      "Marker genes per cluster:\n",
      "%s\n\n",
      "Required JSON format:\n",
      "{\"annotations\":[{\"cluster\":\"0\",\"cell_type\":\"Naive T cell\"}]}"
    ),
    species,
    tissue,
    paste(allowed_labels, collapse = ", "),
    label_guide,
    clusters_text
  )
  
  if (save_debug) {
    dir.create("benchmark_debug/prompts", recursive = TRUE, showWarnings = FALSE)
    writeLines(
      prompt,
      file.path("benchmark_debug/prompts", paste0(dataset_name, "_", model_key, "_prompt.txt"))
    )
  }
  
  message("Calling ", model_key, " API for ", dataset_name, "...")
  
  api_res <- call_llm_api(prompt, model_key, api_key)
  
  if (save_debug && !is.null(api_res$content)) {
    dir.create("benchmark_debug/responses", recursive = TRUE, showWarnings = FALSE)
    writeLines(
      api_res$content,
      file.path("benchmark_debug/responses", paste0(dataset_name, "_", model_key, "_response.txt"))
    )
  }
  
  if (!isTRUE(api_res$success)) {
    warning("LLM call failed for ", dataset_name, ": ", api_res$error)
    
    pred <- setNames(rep("Unknown", length(markers_list)), names(markers_list))
    
    return(list(
      predictions = pred,
      runtime_sec = api_res$runtime_sec,
      cost_usd = api_res$cost_usd,
      tokens = 0
    ))
  }
  
  pred_dict <- parse_llm_response(api_res$content)
  
  if (is.null(pred_dict)) {
    warning("Could not parse LLM response for ", dataset_name)
    pred <- setNames(rep("Unknown", length(markers_list)), names(markers_list))
  } else {
    pred <- sapply(names(markers_list), function(cl) {
      p <- pred_dict[[cl]]
      if (is.null(p) || !nzchar(p)) "Unknown" else p
    })
  }
  
  list(
    predictions = pred,
    runtime_sec = api_res$runtime_sec,
    cost_usd = api_res$cost_usd,
    tokens = api_res$usage$total_tokens %||% NA_real_
  )
}

# -----------------------------------------------------------------------------
# SingleR
# -----------------------------------------------------------------------------

run_singler <- function(seu, dataset_name) {
  start <- Sys.time()
  sce <- as.SingleCellExperiment(seu)
  
  ref <- switch(
    dataset_name,
    PBMC = celldex::MonacoImmuneData(),
    BaronPancreas = celldex::HumanPrimaryCellAtlasData(),
    MuraroPancreas = celldex::HumanPrimaryCellAtlasData(),
    LawlorPancreas = celldex::HumanPrimaryCellAtlasData(),
    SegerstolpePancreas = celldex::HumanPrimaryCellAtlasData(),
    TasicBrain = celldex::MouseRNAseqData(),
    ZeiselBrain = celldex::MouseRNAseqData(),
    RomanovBrain = celldex::MouseRNAseqData(),
    ZilionisLung = celldex::BlueprintEncodeData(),
    TabulaMurisLung = celldex::MouseRNAseqData(),
    celldex::HumanPrimaryCellAtlasData()
  )
  
  if (!"label.main" %in% colnames(colData(ref))) {
    colData(ref)$label.main <- colData(ref)$label
  }
  
  if (!"logcounts" %in% assayNames(ref)) {
    ref <- scuttle::logNormCounts(ref)
  }
  
  if (!"logcounts" %in% assayNames(sce)) {
    sce <- scuttle::logNormCounts(sce)
  }
  
  keep <- !is.na(ref$label.main)
  ref <- ref[, keep]
  
  common_genes <- intersect(rownames(sce), rownames(ref))
  
  if (length(common_genes) < 100) {
    warning("Too few common genes for SingleR on ", dataset_name)
    pred <- rep("Unknown", ncol(seu))
  } else {
    sce <- sce[common_genes, ]
    ref <- ref[common_genes, ]
    
    singler_res <- tryCatch(
      SingleR(
        test = sce,
        ref = ref,
        labels = ref$label.main,
        assay.type.test = "logcounts",
        assay.type.ref = "logcounts"
      ),
      error = function(e) NULL
    )
    
    if (is.null(singler_res)) {
      pred <- rep("Unknown", ncol(seu))
    } else {
      pred <- as.character(singler_res$pruned.labels)
      pred[is.na(pred)] <- "Unknown"
    }
  }
  
  names(pred) <- colnames(seu)
  
  list(
    predictions = pred,
    runtime_sec = as.numeric(difftime(Sys.time(), start, units = "secs")),
    cost_usd = 0,
    tokens = NA_real_
  )
}

# -----------------------------------------------------------------------------
# scType
# -----------------------------------------------------------------------------

.clean_gene_symbols <- function(genes, species = "Human") {
  genes <- trimws(genes)
  genes <- genes[genes != "" & genes != "NA"]
  
  if (length(genes) == 0) return(character())
  
  if (tolower(species) == "human") {
    genes <- toupper(genes)
    
    checked <- suppressMessages(
      suppressWarnings(
        HGNChelper::checkGeneSymbols(genes)
      )
    )
    
    cleaned <- checked$Suggested.Symbol
    cleaned <- cleaned[!is.na(cleaned) & cleaned != ""]
    
    return(unique(cleaned))
  }
  
  unique(genes)
}

gene_sets_prepare <- function(path_to_db_file, cell_type, species = "Human") {
  cell_markers <- openxlsx::read.xlsx(path_to_db_file)
  cell_markers <- cell_markers[cell_markers$tissueType == cell_type, ]
  
  if (nrow(cell_markers) == 0) {
    stop("No scType markers found for tissue: ", cell_type, call. = FALSE)
  }
  
  clean_gene_set <- function(x) {
    x <- gsub(" ", "", x)
    x <- gsub("///", ",", x)
    genes <- unlist(strsplit(x, ","))
    .clean_gene_symbols(genes, species = species)
  }
  
  gs <- lapply(cell_markers$geneSymbolmore1, clean_gene_set)
  names(gs) <- cell_markers$cellName
  
  gs2 <- lapply(cell_markers$geneSymbolmore2, clean_gene_set)
  names(gs2) <- cell_markers$cellName
  
  list(gs_positive = gs, gs_negative = gs2)
}

sctype_score <- function(scRNAseqData, scaled = FALSE, gs, gs2 = NULL) {
  
  
  rownames(scRNAseqData) <- toupper(rownames(scRNAseqData))
  
  gs <- lapply(gs, function(x) intersect(toupper(x), rownames(scRNAseqData)))
  
  if (is.null(gs2)) {
    gs2 <- lapply(gs, function(x) character())
  } else {
    gs2 <- lapply(gs2, function(x) intersect(toupper(x), rownames(scRNAseqData)))
  }
  
  valid <- lengths(gs) > 0
  gs <- gs[valid]
  gs2 <- gs2[valid]
  
  if (length(gs) == 0) {
    stop("No marker genes overlap with expression matrix.", call. = FALSE)
  }
  
  all_marker_genes <- unique(c(unlist(gs), unlist(gs2)))
  all_marker_genes <- intersect(all_marker_genes, rownames(scRNAseqData))
  
  Z <- scRNAseqData[all_marker_genes, , drop = FALSE]
  
  if (!scaled) {
    Z <- as.matrix(Z)
    Z <- t(scale(t(Z)))
    Z[is.na(Z)] <- 0
  } else {
    Z <- as.matrix(Z)
  }
  
  es <- do.call("rbind", lapply(names(gs), function(ct) {
    pos <- gs[[ct]]
    neg <- gs2[[ct]]
    
    pos_score <- colSums(Z[pos, , drop = FALSE]) / sqrt(max(length(pos), 1))
    
    if (length(neg) > 0) {
      neg_score <- colSums(Z[neg, , drop = FALSE]) / sqrt(length(neg))
      pos_score - neg_score
    } else {
      pos_score
    }
  }))
  
  rownames(es) <- names(gs)
  colnames(es) <- colnames(scRNAseqData)
  
  es
}

find_sctype_db <- function() {
  possible_paths <- c("scType/db.xlsx", "scType/ScTypeDB_full.xlsx", SCTYPE_DB)
  
  for (p in possible_paths) {
    if (file.exists(p)) return(normalizePath(p, mustWork = FALSE))
  }
  
  stop("scType database not found.", call. = FALSE)
}

get_expression_matrix_safe <- function(seu) {
  assays <- names(seu@assays)
  
  preferred <- assays[!grepl("ERCC|Spike", assays, ignore.case = TRUE)]
  
  if (length(preferred) == 0) {
    preferred <- assays
  }
  
  for (assay_name in preferred) {
    expr <- tryCatch(
      Seurat::GetAssayData(seu, assay = assay_name, layer = "data"),
      error = function(e) NULL
    )
    
    if (is.null(expr)) {
      expr <- tryCatch(
        Seurat::GetAssayData(seu, assay = assay_name, layer = "counts"),
        error = function(e) NULL
      )
    }
    
    if (!is.null(expr) && nrow(expr) > 500) {
      return(expr)
    }
  }
  
  stop("Could not extract valid non-ERCC expression matrix.", call. = FALSE)
}

run_sctype_custom <- function(seu, tissue, db_file, species = "Human", verbose = TRUE) {
  start <- Sys.time()
  
  tissue_for_sctype <- switch(
    tissue,
    PBMC = "Immune system",
    Blood = "Immune system",
    Pancreas = "Pancreas",
    Brain = "Brain",
    Lung = "Lung",
    tissue
  )
  
  if (verbose) {
    message("scType tissue: ", tissue_for_sctype)
    message("scType DB file: ", db_file)
  }
  
  gs_list <- tryCatch(
    gene_sets_prepare(db_file, tissue_for_sctype, species = species),
    error = function(e) {
      warning("scType marker preparation failed: ", e$message)
      NULL
    }
  )
  
  if (is.null(gs_list)) {
    pred <- rep("Unknown", ncol(seu))
    names(pred) <- colnames(seu)
    
    return(list(
      predictions = pred,
      runtime_sec = as.numeric(difftime(Sys.time(), start, units = "secs")),
      cost_usd = 0,
      tokens = NA_real_
    ))
  }
  
  expr_matrix <- get_expression_matrix_safe(seu)
  
  if (!inherits(expr_matrix, "matrix")) {
    expr_matrix <- as(expr_matrix, "dgCMatrix")
  }
  
  if (verbose) {
    message("scType positive gene sets: ", length(gs_list$gs_positive))
    message("scType negative gene sets: ", length(gs_list$gs_negative))
    message("First scType cell types: ", paste(head(names(gs_list$gs_positive)), collapse = ", "))
    message("Expression matrix dim: ", paste(dim(expr_matrix), collapse = " x "))
    message("First expression genes: ", paste(head(rownames(expr_matrix)), collapse = ", "))
    
    common_genes <- intersect(
      toupper(rownames(expr_matrix)),
      unique(toupper(unlist(gs_list$gs_positive)))
    )
    
    message("Common positive marker genes with expression matrix: ", length(common_genes))
  }
  
  es_max <- tryCatch(
    sctype_score(
      scRNAseqData = expr_matrix,
      scaled = FALSE,
      gs = gs_list$gs_positive,
      gs2 = gs_list$gs_negative
    ),
    error = function(e) {
      warning("scType scoring failed: ", e$message)
      NULL
    }
  )
  
  if (is.null(es_max) || nrow(es_max) == 0 || ncol(es_max) == 0) {
    pred <- rep("Unknown", ncol(seu))
    names(pred) <- colnames(seu)
    
    return(list(
      predictions = pred,
      runtime_sec = as.numeric(difftime(Sys.time(), start, units = "secs")),
      cost_usd = 0,
      tokens = NA_real_
    ))
  }
  
  cluster_ids <- as.character(seu$seurat_clusters)
  
  cluster_pred <- sapply(unique(cluster_ids), function(cl) {
    cells <- which(cluster_ids == cl)
    
    if (length(cells) == 0) return("Unknown")
    
    scores <- rowMeans(es_max[, cells, drop = FALSE], na.rm = TRUE)
    
    if (all(is.na(scores))) return("Unknown")
    
    names(scores)[which.max(scores)]
  })
  
  if (verbose) {
    message("scType cluster predictions:")
    print(cluster_pred)
  }
  
  pred <- unname(cluster_pred[cluster_ids])
  pred[is.na(pred)] <- "Unknown"
  names(pred) <- colnames(seu)
  
  list(
    predictions = pred,
    runtime_sec = as.numeric(difftime(Sys.time(), start, units = "secs")),
    cost_usd = 0,
    tokens = NA_real_
  )
}

# -----------------------------------------------------------------------------
# Allowed labels for closed-label LLM benchmark
# -----------------------------------------------------------------------------

get_allowed_labels <- function(tissue) {
  tissue <- tolower(tissue)
  
  if (tissue == "pbmc") {
    return(c(
      "Naive T cell",
      "Cytotoxic T cell",
      "Classical monocyte",
      "Non-classical monocyte",
      "B cell",
      "Natural killer cell",
      "Dendritic cell",
      "Platelet",
      "Unknown"
    ))
  }
  
  if (tissue == "pancreas") {
    return(c(
      "Alpha cell",
      "Beta cell",
      "Delta cell",
      "Gamma cell",
      "Acinar cell",
      "Ductal cell",
      "Endothelial cell",
      "Stellate cell",
      "Macrophage",
      "Immune cell",
      "Unknown"
    ))
  }
  
  if (tissue == "brain") {
    return(c(
      "Excitatory neuron",
      "Inhibitory neuron",
      "Astrocyte",
      "Microglia",
      "Oligodendrocyte",
      "Oligodendrocyte precursor cell",
      "Endothelial cell",
      "Pericyte",
      "Ependymal cell",
      "Unknown"
    ))
  }
  
  if (tissue == "lung") {
    return(c(
      "T cell",
      "B cell",
      "Plasma cell",
      "Natural killer cell",
      "Macrophage",
      "Monocyte",
      "Neutrophil",
      "Dendritic cell",
      "Mast cell",
      "Epithelial cell",
      "Endothelial cell",
      "Fibroblast",
      "Smooth muscle cell",
      "Unknown"
    ))
  }
  
  "Unknown"
}