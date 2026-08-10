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

  if (is.null(model)) {
    stop("Unknown model key: ", model_key, call. = FALSE)
  }

  if (isTRUE(model$is_ollama)) {
    return(call_ollama_api(prompt, model, max_retries = max_retries))
  }
  api_format <- model$api_format %||% "chat"

  for (attempt in seq_len(max_retries)) {
    res <- tryCatch({
      body <- benchmark_llm_request_body(prompt, model, api_format)

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

      content <- benchmark_llm_response_text(data, api_format)
      benchmark_validate_llm_response_text(content, data, api_format)
      usage <- benchmark_standardise_usage(data$usage %||% NULL, api_format)

      elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))

      cost <- benchmark_estimate_cost_usd(usage, model)
      pricing <- benchmark_pricing_metadata(model)

      list(
        success = TRUE,
        content = content,
        usage = usage,
        runtime_sec = elapsed,
        cost_usd = cost,
        pricing_date = pricing$PricingDate,
        pricing_source = pricing$PricingSource,
        input_cost_per_1m_tokens = pricing$InputCostPer1MTokens,
        output_cost_per_1m_tokens = pricing$OutputCostPer1MTokens
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
    cost_usd = NA_real_,
    pricing_date = model$pricing_date %||% NA_character_,
    pricing_source = model$pricing_source %||% NA_character_,
    input_cost_per_1m_tokens = 1000 * suppressWarnings(as.numeric(model$input_cost_per_1k %||% NA_real_)),
    output_cost_per_1m_tokens = 1000 * suppressWarnings(as.numeric(model$output_cost_per_1k %||% NA_real_))
  )
}

benchmark_llm_request_body <- function(prompt, model, api_format = "chat") {
  system_prompt <- "You are a bioinformatics expert. Output ONLY valid JSON. Do not include markdown or extra text."

  if (identical(api_format, "responses")) {
    body <- list(
      model = model$model_id,
      instructions = system_prompt,
      input = prompt,
      max_output_tokens = model$max_tokens %||% 2000
    )

    if (!is.null(model$reasoning_effort) && nzchar(model$reasoning_effort)) {
      body$reasoning <- list(effort = model$reasoning_effort)
    }
    if (!is.null(model$text_verbosity) && nzchar(model$text_verbosity)) {
      body$text <- list(verbosity = model$text_verbosity)
    }
  } else {
    body <- list(
      model = model$model_id,
      messages = list(
        list(role = "system", content = system_prompt),
        list(role = "user", content = prompt)
      ),
      temperature = model$temperature %||% 0,
      max_tokens = model$max_tokens %||% 2000
    )
  }

  body[!vapply(body, is.null, logical(1))]
}

benchmark_llm_response_text <- function(data, api_format = "chat") {
  if (identical(api_format, "responses")) {
    if (!is.null(data$output_text) && nzchar(data$output_text)) {
      return(data$output_text)
    }

    pieces <- character()
    for (item in data$output %||% list()) {
      for (part in item$content %||% list()) {
        text <- part$text %||% part$output_text %||% ""
        if (nzchar(text)) pieces <- c(pieces, text)
      }
    }
    return(paste(pieces, collapse = "\n"))
  }

  data$choices[[1]]$message$content %||% ""
}

benchmark_response_is_blank <- function(x) {
  is.null(x) ||
    length(x) == 0 ||
    all(!nzchar(trimws(as.character(x))))
}

benchmark_validate_llm_response_text <- function(content, data, api_format = "chat") {
  status <- data$status %||% ""
  reason <- data$incomplete_details$reason %||% ""

  if (identical(api_format, "responses") &&
      identical(status, "incomplete") &&
      benchmark_response_is_blank(content)) {
    stop(
      "OpenAI Responses API returned an incomplete response without text",
      if (nzchar(reason)) paste0(" (reason: ", reason, ")") else "",
      call. = FALSE
    )
  }

  if (benchmark_response_is_blank(content)) {
    stop("LLM API returned empty response text.", call. = FALSE)
  }

  invisible(TRUE)
}

benchmark_standardise_usage <- function(usage, api_format = "chat") {
  if (is.null(usage)) {
    return(list(
      prompt_tokens = NA_real_,
      completion_tokens = NA_real_,
      total_tokens = NA_real_
    ))
  }

  if (identical(api_format, "responses")) {
    prompt_tokens <- usage$input_tokens %||% usage$prompt_tokens %||% NA_real_
    completion_tokens <- usage$output_tokens %||% usage$completion_tokens %||% NA_real_
  } else {
    prompt_tokens <- usage$prompt_tokens %||% NA_real_
    completion_tokens <- usage$completion_tokens %||% NA_real_
  }
  total_tokens <- usage$total_tokens %||% NA_real_
  if (is.na(suppressWarnings(as.numeric(total_tokens)))) {
    total_tokens <- suppressWarnings(as.numeric(prompt_tokens)) +
      suppressWarnings(as.numeric(completion_tokens))
  }

  list(
    prompt_tokens = prompt_tokens,
    completion_tokens = completion_tokens,
    total_tokens = total_tokens
  )
}

benchmark_safe_number <- function(x, default = 0) {
  out <- suppressWarnings(as.numeric(x %||% default))
  if (length(out) == 0 || is.na(out)) {
    return(default)
  }
  out
}

benchmark_estimate_cost_usd <- function(usage, model) {
  prompt_tokens <- benchmark_safe_number(usage$prompt_tokens %||% 0)
  completion_tokens <- benchmark_safe_number(usage$completion_tokens %||% 0)
  input_cost <- suppressWarnings(as.numeric(model$input_cost_per_1k %||% NA_real_))
  output_cost <- suppressWarnings(as.numeric(model$output_cost_per_1k %||% NA_real_))

  if (!is.na(input_cost) || !is.na(output_cost)) {
    input_cost[is.na(input_cost)] <- 0
    output_cost[is.na(output_cost)] <- 0
    return((prompt_tokens / 1000) * input_cost + (completion_tokens / 1000) * output_cost)
  }

  NA_real_
}

benchmark_pricing_metadata <- function(model) {
  input_cost <- suppressWarnings(as.numeric(model$input_cost_per_1k %||% NA_real_))
  output_cost <- suppressWarnings(as.numeric(model$output_cost_per_1k %||% NA_real_))

  list(
    PricingDate = model$pricing_date %||% NA_character_,
    PricingSource = model$pricing_source %||% NA_character_,
    InputCostPer1MTokens = ifelse(is.na(input_cost), NA_real_, input_cost * 1000),
    OutputCostPer1MTokens = ifelse(is.na(output_cost), NA_real_, output_cost * 1000)
  )
}

call_ollama_api <- function(prompt, model, max_retries = 3, timeout_sec = 300) {
  start <- Sys.time()
  last_error <- NULL
  full_prompt <- paste(
    "You are a bioinformatics expert. Output ONLY valid JSON. Do not include markdown or extra text.",
    prompt,
    sep = "\n\n"
  )

  for (attempt in seq_len(max_retries)) {
    res <- tryCatch({
      body <- list(
        model = model$model_id,
        prompt = full_prompt,
        stream = FALSE,
        options = list(
          temperature = model$temperature %||% 0,
          num_predict = model$max_tokens %||% 2000
        )
      )

      req <- httr2::request(model$api_url) |>
        httr2::req_timeout(timeout_sec) |>
        httr2::req_body_json(body)

      resp <- httr2::req_perform(req)
      data <- httr2::resp_body_json(resp, simplifyVector = FALSE)

      content <- data$response %||% ""
      prompt_tokens <- data$prompt_eval_count %||% ceiling(nchar(full_prompt) / 4)
      completion_tokens <- data$eval_count %||% ceiling(nchar(content) / 4)
      total_tokens <- prompt_tokens + completion_tokens
      elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))

      list(
        success = TRUE,
        content = content,
        usage = list(
          prompt_tokens = prompt_tokens,
          completion_tokens = completion_tokens,
          total_tokens = total_tokens
        ),
        runtime_sec = elapsed,
        cost_usd = 0
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
    error = last_error %||% "Unknown Ollama error",
    usage = list(prompt_tokens = 0, completion_tokens = 0, total_tokens = 0),
    runtime_sec = as.numeric(difftime(Sys.time(), start, units = "secs")),
    cost_usd = 0
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

normalise_llm_cluster_id <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("^cluster[[:space:]_:-]*", "", x, perl = TRUE)
  x <- gsub("[[:space:]]+", "", x, perl = TRUE)
  x
}

build_llm_annotation_prompt <- function(markers_list,
                                        tissue,
                                        species) {
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
      "Epsilon cell: GHRL, MLN.",
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
      "Radial glia: HOPX, VIM, SOX2, PAX6, NES.",
      "Neural progenitor cell: SOX2, NES, MKI67, TOP2A, PAX6.",
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
    `hematopoietic system` = paste(
      "Hematopoietic label guide:",
      "Hematopoietic stem/progenitor cell: CD34, PROM1, SPINK2, HLF, MEIS1, GATA2.",
      "Common lymphoid progenitor: IL7R, DNTT, RAG1, RAG2, VPREB1, CD79A.",
      "Granulocyte-monocyte progenitor: MPO, ELANE, AZU1, CTSG, PRTN3, LYZ.",
      "Megakaryocyte-erythroid progenitor: GATA1, KLF1, GYPA, ALAS2, HBB, HBA1.",
      "Erythroid progenitor cell: GATA1, KLF1, ALAS2, HBB, HBA1, HBA2, AHSP.",
      sep = "\n"
    ),
    kidney = paste(
      "Kidney label guide:",
      "Proximal tubule cell: SLC34A1, SLC5A2, LRP2, CUBN, ALDOB, ACSM2, AQP1.",
      "Distal convoluted tubule cell: SLC12A3, TRPM6, CALB1, PVALB, FXYD2.",
      "Connecting tubule cell: CALB1, SCNN1G, TRPV5, FXYD4, KCNJ1.",
      "Loop of Henle cell: SLC12A1, UMOD, CLDN10, KCNJ1, FXYD2.",
      "Collecting duct principal cell: AQP2, AQP3, AQP4, SCNN1G, AVPR2.",
      "Intercalated cell: ATP6V1B1, ATP6V0D2, SLC4A1, SLC26A4, FOXI1.",
      "Podocyte: NPHS1, NPHS2, PODXL, WT1, SYNPO.",
      "Mesangial cell: PDGFRB, MEIS1, POSTN, TAGLN, ACTA2.",
      "Endothelial cell: PECAM1, VWF, KDR, EMCN, CLDN5.",
      "Immune cell: LYZ, LAPTM5, PTPRC, CD68, C1QA.",
      sep = "\n"
    ),
    ""
  )

  sprintf(
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
      "{\"annotations\":[{\"cluster\":\"0\",\"cell_type\":\"Naive T cell\",",
      "\"confidence\":0.95,\"candidate_cell_types\":[\"Naive T cell\",\"Cytotoxic T cell\"],",
      "\"is_mixed\":false,\"primary_cell_type\":\"Naive T cell\",",
      "\"secondary_cell_type\":\"\",\"tissue_consistency\":\"expected\",",
      "\"reasoning\":\"CD3D, CD3E, IL7R, and CCR7 support a naive T cell identity\"}]}"
    ),
    species,
    tissue,
    paste(allowed_labels, collapse = ", "),
    label_guide,
    clusters_text
  )
}

run_llm_annotation <- function(markers_list,
                               tissue,
                               species,
                               dataset_name,
                               model_key,
                               api_key,
                               save_debug = TRUE) {
  prompt <- build_llm_annotation_prompt(markers_list, tissue, species)

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
      tokens = 0,
      prompt_tokens = 0,
      completion_tokens = 0,
      pricing_date = api_res$pricing_date %||% NA_character_,
      pricing_source = api_res$pricing_source %||% NA_character_,
      input_cost_per_1m_tokens = api_res$input_cost_per_1m_tokens %||% NA_real_,
      output_cost_per_1m_tokens = api_res$output_cost_per_1m_tokens %||% NA_real_
    ))
  }

  pred_dict <- NULL
  if (exists("parse_ablation_annotation_response", mode = "function")) {
    parsed_annotations <- parse_ablation_annotation_response(api_res$content)
    if (is.data.frame(parsed_annotations) && nrow(parsed_annotations) > 0) {
      pred_dict <- stats::setNames(
        as.character(parsed_annotations$CellType),
        as.character(parsed_annotations$Cluster)
      )
    }
  }
  if (is.null(pred_dict)) {
    pred_dict <- parse_llm_response(api_res$content)
  }

  if (is.null(pred_dict)) {
    warning("Could not parse LLM response for ", dataset_name)
    pred <- setNames(rep("Unknown", length(markers_list)), names(markers_list))
  } else {
    pred_names <- names(pred_dict)
    pred_name_norm <- normalise_llm_cluster_id(pred_names)
    pred <- sapply(names(markers_list), function(cl) {
      idx <- match(normalise_llm_cluster_id(cl), pred_name_norm)
      p <- if (is.na(idx)) NULL else pred_dict[[idx]]
      if (is.null(p) || !nzchar(p)) "Unknown" else p
    })
  }

  list(
    predictions = pred,
    runtime_sec = api_res$runtime_sec,
    cost_usd = api_res$cost_usd,
    tokens = api_res$usage$total_tokens %||% NA_real_,
    prompt_tokens = api_res$usage$prompt_tokens %||% NA_real_,
    completion_tokens = api_res$usage$completion_tokens %||% NA_real_,
    pricing_date = api_res$pricing_date %||% NA_character_,
    pricing_source = api_res$pricing_source %||% NA_character_,
    input_cost_per_1m_tokens = api_res$input_cost_per_1m_tokens %||% NA_real_,
    output_cost_per_1m_tokens = api_res$output_cost_per_1m_tokens %||% NA_real_,
    response_text = api_res$content %||% "",
    prompt_text = prompt
  )
}

# -----------------------------------------------------------------------------
# SingleR
# -----------------------------------------------------------------------------

get_benchmark_reference <- function(dataset_name) {
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

  ref[, !is.na(ref$label.main)]
}

prepare_sce_for_reference <- function(seu, ref) {
  sce <- as.SingleCellExperiment(seu)

  if (!"logcounts" %in% assayNames(ref)) {
    ref <- scuttle::logNormCounts(ref)
  }

  if (!"logcounts" %in% assayNames(sce)) {
    sce <- scuttle::logNormCounts(sce)
  }

  common_genes <- intersect(rownames(sce), rownames(ref))

  list(
    sce = sce,
    ref = ref,
    common_genes = common_genes
  )
}

as_sce_for_scmap <- function(x) {
  if (inherits(x, "SingleCellExperiment")) {
    return(x)
  }

  coerced <- tryCatch(
    methods::as(x, "SingleCellExperiment"),
    error = function(e) NULL
  )

  if (!is.null(coerced)) {
    return(coerced)
  }

  SingleCellExperiment::SingleCellExperiment(
    assays = SummarizedExperiment::assays(x),
    rowData = SummarizedExperiment::rowData(x),
    colData = SummarizedExperiment::colData(x),
    metadata = S4Vectors::metadata(x)
  )
}

run_singler <- function(seu, dataset_name) {
  start <- Sys.time()
  ref_data <- prepare_sce_for_reference(
    seu = seu,
    ref = get_benchmark_reference(dataset_name)
  )
  sce <- as_sce_for_scmap(ref_data$sce)
  ref <- as_sce_for_scmap(ref_data$ref)
  common_genes <- ref_data$common_genes

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
# scmap
# -----------------------------------------------------------------------------

is_scmap_available <- function() {
  requireNamespace("scmap", quietly = TRUE)
}

extract_scmap_labels <- function(scmap_res, n_cells) {
  labels <- NULL

  if (!is.null(scmap_res$combined_labs)) {
    labels <- scmap_res$combined_labs
  } else if (!is.null(scmap_res$scmap_cluster_labs)) {
    labels <- scmap_res$scmap_cluster_labs
    if (is.matrix(labels) || is.data.frame(labels)) {
      if (ncol(labels) == n_cells) {
        labels <- labels[1, , drop = TRUE]
      } else if (nrow(labels) == n_cells) {
        labels <- labels[, 1, drop = TRUE]
      } else {
        labels <- as.vector(as.matrix(labels))
      }
    }
  }

  if (is.null(labels)) {
    return(rep("Unknown", n_cells))
  }

  labels <- as.character(labels)
  if (length(labels) != n_cells) {
    return(rep("Unknown", n_cells))
  }

  labels[is.na(labels) | labels == "" | labels == "unassigned"] <- "Unknown"
  labels
}

extract_scmap_cell_labels <- function(scmap_res, reference_labels, n_cells) {
  if (is.null(scmap_res) || length(scmap_res) == 0) {
    return(rep("Unknown", n_cells))
  }

  res <- scmap_res[[1]]
  cells <- res$cells

  if (is.null(cells)) {
    return(rep("Unknown", n_cells))
  }

  if (is.vector(cells)) {
    cells <- matrix(cells, ncol = length(cells))
  }

  if (ncol(cells) != n_cells && nrow(cells) == n_cells) {
    cells <- t(cells)
  }

  if (ncol(cells) != n_cells) {
    return(rep("Unknown", n_cells))
  }

  vapply(seq_len(ncol(cells)), function(i) {
    idx <- as.integer(cells[, i])
    idx <- idx[!is.na(idx) & idx > 0 & idx <= length(reference_labels)]

    if (length(idx) == 0) {
      return("Unknown")
    }

    labs <- as.character(reference_labels[idx])
    labs <- labs[!is.na(labs) & labs != ""]

    if (length(labs) == 0) {
      return("Unknown")
    }

    counts <- sort(table(labs), decreasing = TRUE)
    winners <- names(counts)[counts == max(counts)]
    labs[match(TRUE, labs %in% winners)]
  }, character(1))
}

select_scmap_features_robust <- function(ref, dataset_name, n_features = 500) {
  selected_ref <- tryCatch(
    scmap::selectFeatures(ref, n_features = n_features, suppress_plot = TRUE),
    error = function(e) {
      warning(
        "scmap selectFeatures failed on ", dataset_name,
        "; using high-variance reference genes. Reason: ", e$message
      )
      NULL
    }
  )

  if (!is.null(selected_ref) &&
      "scmap_features" %in% colnames(rowData(selected_ref)) &&
      sum(rowData(selected_ref)$scmap_features, na.rm = TRUE) >= 10) {
    return(selected_ref)
  }

  expr <- as.matrix(logcounts(ref))
  vars <- apply(expr, 1, stats::var, na.rm = TRUE)
  vars[!is.finite(vars)] <- NA_real_
  vars <- vars[!is.na(vars) & vars > 0]

  if (length(vars) < 10) {
    stop("Fewer than 10 variable genes available for scmap on ", dataset_name, call. = FALSE)
  }

  selected <- names(sort(vars, decreasing = TRUE))[seq_len(min(n_features, length(vars)))]
  rowData(ref)$scmap_features <- rownames(ref) %in% selected
  rowData(ref)$scmap_scores <- NA_real_
  ref
}

run_scmap <- function(seu, dataset_name) {
  start <- Sys.time()

  if (!is_scmap_available()) {
    stop("scmap is not installed. Install it with BiocManager::install('scmap').", call. = FALSE)
  }

  ref_data <- prepare_sce_for_reference(
    seu = seu,
    ref = get_benchmark_reference(dataset_name)
  )
  sce <- as_sce_for_scmap(ref_data$sce)
  ref <- as_sce_for_scmap(ref_data$ref)
  common_genes <- ref_data$common_genes

  if (length(common_genes) < 100) {
    warning("Too few common genes for scmap on ", dataset_name)
    pred <- rep("Unknown", ncol(seu))
  } else {
    sce <- sce[common_genes, ]
    ref <- ref[common_genes, ]

    rowData(ref)$feature_symbol <- rownames(ref)
    rowData(sce)$feature_symbol <- rownames(sce)
    colData(ref)$scmap_label <- ref$label.main

    scmap_res <- tryCatch({
      ref <- select_scmap_features_robust(ref, dataset_name)
      ref <- scmap::indexCell(ref)
      scmap::scmapCell(
        projection = sce,
        index_list = list(reference = metadata(ref)$scmap_cell_index),
        w = 10
      )
    }, error = function(e) {
      warning("scmap failed on ", dataset_name, ": ", e$message)
      NULL
    })

    pred <- if (is.null(scmap_res)) {
      rep("Unknown", ncol(seu))
    } else {
      extract_scmap_cell_labels(scmap_res, ref$scmap_label, ncol(seu))
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
# CellTypist
# -----------------------------------------------------------------------------

is_celltypist_available <- function() {
  requireNamespace("reticulate", quietly = TRUE) &&
    reticulate::py_module_available("celltypist") &&
    reticulate::py_module_available("anndata") &&
    reticulate::py_module_available("numpy")
}

get_expression_matrix_safe <- function(seu) {
  assay <- Seurat::DefaultAssay(seu)

  expr <- tryCatch(
    Seurat::GetAssayData(seu, assay = assay, layer = "data"),
    error = function(e) NULL
  )

  if (is.null(expr) || nrow(expr) == 0 || ncol(expr) == 0) {
    expr <- Seurat::GetAssayData(seu, assay = assay, slot = "data")
  }

  expr
}

cluster_average_expression <- function(seu) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Matrix is required for cluster-average expression.", call. = FALSE)
  }

  expr <- get_expression_matrix_safe(seu)
  clusters <- as.character(Seurat::Idents(seu))
  cluster_names <- sort(unique(clusters))

  avg <- vapply(
    cluster_names,
    function(cl) Matrix::rowMeans(expr[, clusters == cl, drop = FALSE]),
    numeric(nrow(expr))
  )

  avg <- t(avg)
  rownames(avg) <- cluster_names
  colnames(avg) <- rownames(expr)
  avg
}

parse_celltypist_labels <- function(predicted_labels, sample_names) {
  labels <- reticulate::py_to_r(predicted_labels)

  if (is.data.frame(labels)) {
    label_col <- intersect(c("majority_voting", "predicted_labels"), names(labels))[1]
    if (!is.na(label_col)) {
      labels <- labels[[label_col]]
    } else {
      labels <- labels[[1]]
    }
  }

  labels <- as.character(labels)
  if (length(labels) != length(sample_names)) {
    stop("CellTypist returned ", length(labels), " labels for ",
         length(sample_names), " cell profiles.", call. = FALSE)
  }

  labels[is.na(labels) | labels == ""] <- "Unknown"
  stats::setNames(labels, sample_names)
}

get_celltypist_model <- function(tissue, dataset_name = NULL) {
  if (nzchar(CELLTYPIST_MODEL)) {
    return(CELLTYPIST_MODEL)
  }

  tissue_key <- tolower(tissue)
  dataset_key <- tolower(dataset_name %||% "")

  if (tissue_key %in% c("pbmc", "immune system") || grepl("pbmc", dataset_key)) {
    return("Immune_All_Low.pkl")
  }

  if (tissue_key == "pancreas" || grepl("pancreas", dataset_key)) {
    return("Adult_Human_PancreaticIslet.pkl")
  }

  if (tissue_key == "lung" || grepl("lung", dataset_key)) {
    return("Human_Lung_Atlas.pkl")
  }

  "Immune_All_Low.pkl"
}

run_celltypist <- function(seu, tissue, species = "Human", dataset_name = NULL) {
  start <- Sys.time()

  if (!is_celltypist_available()) {
    stop(
      "CellTypist Python dependencies are not available. Install reticulate plus ",
      "Python modules celltypist, anndata, and numpy.",
      call. = FALSE
    )
  }

  if (tolower(species) != "human") {
    stop("CellTypist is enabled only for human datasets in this benchmark.", call. = FALSE)
  }

  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Matrix is required for native CellTypist evaluation.", call. = FALSE)
  }

  expr <- get_expression_matrix_safe(seu)
  if (!inherits(expr, "sparseMatrix")) {
    expr <- Matrix::Matrix(expr, sparse = TRUE)
  }

  cell_names <- colnames(expr)
  gene_names <- rownames(expr)

  celltypist <- reticulate::import("celltypist", convert = FALSE)
  anndata <- reticulate::import("anndata", convert = FALSE)

  x_py <- reticulate::r_to_py(t(expr), convert = FALSE)
  x_py <- x_py$astype("float32")

  adata <- anndata$AnnData(X = x_py)
  adata$obs_names <- cell_names
  adata$var_names <- gene_names

  model_arg <- get_celltypist_model(tissue, dataset_name)
  message("CellTypist model: ", model_arg)

  annotation <- tryCatch({
    celltypist$annotate(adata, model = model_arg, majority_voting = FALSE)
  }, error = function(e) {
    stop("CellTypist failed for ", tissue, ": ", e$message, call. = FALSE)
  })

  pred <- parse_celltypist_labels(annotation$predicted_labels, cell_names)

  list(
    predictions = pred,
    prediction_level = "cell",
    model = model_arg,
    n_cells = length(cell_names),
    n_genes = length(gene_names),
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
      "Epsilon cell",
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
      "Radial glia",
      "Neural progenitor cell",
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

  if (tissue %in% c("hematopoietic system", "hematopoietic", "hspc")) {
    return(c(
      "Hematopoietic stem/progenitor cell",
      "Common lymphoid progenitor",
      "Granulocyte-monocyte progenitor",
      "Megakaryocyte-erythroid progenitor",
      "Erythroid progenitor cell",
      "Unknown"
    ))
  }

  if (tissue == "kidney") {
    return(c(
      "Proximal tubule cell",
      "Distal convoluted tubule cell",
      "Connecting tubule cell",
      "Loop of Henle cell",
      "Collecting duct principal cell",
      "Intercalated cell",
      "Podocyte",
      "Mesangial cell",
      "Endothelial cell",
      "Immune cell",
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
