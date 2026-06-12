# R/api.R
#' Unified API Handler for LLM Cell Type Annotation
#'
#' Supports DeepSeek and local Ollama endpoints.
#'
#' @importFrom httr2 request req_headers req_body_json req_perform resp_body_json req_timeout
#' @importFrom jsonlite fromJSON

MODELS <- list(
  deepseek = list(
    name = "DeepSeek-Chat",
    api_url = "https://api.deepseek.com/v1/chat/completions",
    api_url_env = "DEEPSEEK_API_URL",
    model_id = "deepseek-chat",
    model_id_env = "DEEPSEEK_MODEL_ID",
    max_tokens = 2000,
    temperature = 0.1,
    cost_per_1k_tokens = 0.00014,
    requires_api_key = TRUE,
    api_key_env = "DEEPSEEK_API_KEY",
    is_ollama = FALSE
  ),
  ollama = list(
    name = "Ollama local",
    api_url = "http://localhost:11434/api/generate",
    api_url_env = "OLLAMA_API_URL",
    model_id = "llama3.2:latest",
    model_id_env = "OLLAMA_MODEL_ID",
    max_tokens = 2000,
    temperature = 0.1,
    cost_per_1k_tokens = 0,
    requires_api_key = FALSE,
    api_key_env = character(),
    is_ollama = TRUE
  )
)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

#' Get model configuration
#' @param model_name Model identifier.
#' @return Model configuration list.
#' @export
get_model_config <- function(model_name) {
  if (!model_name %in% names(MODELS)) {
    stop(
      "Unknown model: ", model_name,
      ". Available models: ", paste(names(MODELS), collapse = ", "),
      call. = FALSE
    )
  }
  
  model <- MODELS[[model_name]]
  api_url_override <- first_env(model$api_url_env %||% character())
  model_id_override <- first_env(model$model_id_env %||% character())

  if (nzchar(api_url_override)) {
    if (grepl("/$", api_url_override) && !grepl("/chat/completions/?$", api_url_override)) {
      api_url_override <- paste0(api_url_override, "chat/completions")
    }
    model$api_url <- api_url_override
  }

  if (nzchar(model_id_override)) {
    model$model_id <- model_id_override
  }

  model
}

#' Resolve API key from explicit argument or documented environment variables
#' @param model Model configuration from get_model_config().
#' @param api_key Optional explicit key.
#' @return API key string or NULL.
#' @keywords internal
resolve_api_key <- function(model, api_key = NULL) {
  if (!is.null(api_key) && nzchar(api_key)) {
    return(api_key)
  }

  env_key <- first_env(model$api_key_env %||% character())
  if (nzchar(env_key)) {
    return(env_key)
  }

  NULL
}

#' Call LLM API with retry logic
#' @param prompt User prompt.
#' @param model Model configuration.
#' @param api_key API key. Not required for Ollama.
#' @param max_retries Maximum retry attempts.
#' @param timeout_sec Request timeout in seconds.
#' @return Standardized API response list.
#' @export
call_llm_api <- function(prompt,
                         model,
                         api_key = NULL,
                         max_retries = 3,
                         timeout_sec = 60) {
  stopifnot(is.character(prompt), length(prompt) == 1)
  
  if (isTRUE(model$is_ollama)) {
    return(call_ollama_api(prompt, model, max_retries, timeout_sec))
  }

  api_key <- resolve_api_key(model, api_key)
  
  if (isTRUE(model$requires_api_key) && (is.null(api_key) || identical(api_key, ""))) {
    env_hint <- paste(model$api_key_env %||% character(), collapse = ", ")
    stop("API key is required for model: ", model$name,
         if (nzchar(env_hint)) paste0(". Set one of: ", env_hint) else "",
         call. = FALSE)
  }
  
  start_time <- Sys.time()
  last_error <- NULL
  
  for (attempt in seq_len(max_retries)) {
    response <- tryCatch({
      req <- httr2::request(model$api_url) |>
        httr2::req_timeout(timeout_sec) |>
        httr2::req_headers(
          Authorization = paste("Bearer", api_key),
          "Content-Type" = "application/json"
        ) |>
        httr2::req_body_json(list(
          model = model$model_id,
          messages = list(
            list(role = "system", content = create_system_prompt()),
            list(role = "user", content = prompt)
          ),
          temperature = model$temperature,
          max_tokens = model$max_tokens
        ))
      
      resp <- httr2::req_perform(req)
      data <- httr2::resp_body_json(resp, simplifyVector = FALSE)
      
      content <- data$choices[[1]]$message$content %||% ""
      usage <- data$usage %||% list(
        prompt_tokens = NA_integer_,
        completion_tokens = NA_integer_,
        total_tokens = NA_integer_
      )
      
      elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      
      list(
        success = TRUE,
        content = content,
        usage = usage,
        model = model$name,
        latency_sec = elapsed,
        attempt = attempt
      )
    }, error = function(e) {
      last_error <<- conditionMessage(e)
      NULL
    })
    
    if (!is.null(response)) return(response)
    
    if (attempt < max_retries) {
      Sys.sleep(min(2^attempt, 10))
    }
  }
  
  list(
    success = FALSE,
    content = NULL,
    error = last_error %||% "Unknown API error",
    usage = list(total_tokens = 0),
    model = model$name,
    latency_sec = as.numeric(difftime(Sys.time(), start_time, units = "secs")),
    attempt = max_retries
  )
}

#' Call local Ollama API
#' @keywords internal
call_ollama_api <- function(prompt,
                            model,
                            max_retries = 3,
                            timeout_sec = 120) {
  start_time <- Sys.time()
  full_prompt <- paste(create_system_prompt(), prompt, sep = "\n\n")
  last_error <- NULL
  
  for (attempt in seq_len(max_retries)) {
    response <- tryCatch({
      req <- httr2::request(model$api_url) |>
        httr2::req_timeout(timeout_sec) |>
        httr2::req_body_json(list(
          model = model$model_id,
          prompt = full_prompt,
          stream = FALSE,
          options = list(
            temperature = model$temperature,
            num_predict = model$max_tokens
          )
        ))
      
      resp <- httr2::req_perform(req)
      data <- httr2::resp_body_json(resp, simplifyVector = FALSE)
      
      content <- data$response %||% ""
      total_tokens <- ceiling((nchar(full_prompt) + nchar(content)) / 4)
      
      list(
        success = TRUE,
        content = content,
        usage = list(
          prompt_tokens = ceiling(nchar(full_prompt) / 4),
          completion_tokens = ceiling(nchar(content) / 4),
          total_tokens = total_tokens
        ),
        model = model$name,
        latency_sec = as.numeric(difftime(Sys.time(), start_time, units = "secs")),
        attempt = attempt
      )
    }, error = function(e) {
      last_error <<- conditionMessage(e)
      NULL
    })
    
    if (!is.null(response)) return(response)
    if (attempt < max_retries) Sys.sleep(min(2^attempt, 10))
  }
  
  list(
    success = FALSE,
    content = NULL,
    error = last_error %||% "Unknown Ollama error",
    usage = list(total_tokens = 0),
    model = model$name,
    latency_sec = as.numeric(difftime(Sys.time(), start_time, units = "secs")),
    attempt = max_retries
  )
}

#' Parse an LLM annotation response
#'
#' Extracts a JSON payload from a model response and normalizes it to the
#' annotation schema used by DeepSeekCell. A simple line-by-line fallback is
#' used for legacy responses.
#'
#' @param response_text Raw response text from an LLM endpoint.
#' @return Data frame with annotation columns.
#' @export
parse_annotation_response <- function(response_text) {
  
  if (is.null(response_text) || !nzchar(trimws(response_text))) {
    return(.empty_annotation_result())
  }
  
  json_str <- extract_json_payload(response_text)

  if (is.null(json_str)) {
    warning("No JSON object found in response.")
    return(.parse_line_by_line(response_text))
  }
  
  parsed <- tryCatch(
    jsonlite::fromJSON(json_str, simplifyDataFrame = TRUE),
    error = function(e) {
      warning("JSON parsing failed: ", conditionMessage(e))
      NULL
    }
  )
  
  if (is.null(parsed)) {
    return(.parse_line_by_line(response_text))
  }
  
  ann <- parsed$annotations %||% parsed

  if (is.list(ann) && !is.data.frame(ann)) {
    ann <- tryCatch(
      as.data.frame(ann, stringsAsFactors = FALSE),
      error = function(e) NULL
    )
  }

  if (!is.data.frame(ann) || nrow(ann) == 0) {
    return(.empty_annotation_result())
  }

  .normalise_annotation_dataframe(ann)
}

.empty_annotation_result <- function() {
  data.frame(
    Cluster = character(),
    CellType = character(),
    Confidence = numeric(),
    IsMixed = logical(),
    PrimaryCellType = character(),
    SecondaryCellType = character(),
    TissueConsistency = character(),
    Reasoning = character(),
    stringsAsFactors = FALSE
  )
}

.normalise_annotation_dataframe <- function(ann) {
  names(ann) <- tolower(names(ann))
  
  get_col <- function(df, candidates, default) {
    hit <- candidates[candidates %in% names(df)]
    if (length(hit) == 0) {
      if (length(default) == nrow(df)) {
        return(default)
      }
      return(rep(default, nrow(df)))
    }
    df[[hit[1]]]
  }
  
  confidence <- as_confidence(get_col(ann, "confidence", 0.5))
  is_mixed <- as_flag(get_col(ann, c("is_mixed", "ismixed", "mixed"), FALSE))
  
  result <- data.frame(
    Cluster = as.character(get_col(ann, "cluster", paste0("Cluster", seq_len(nrow(ann))))),
    CellType = as.character(get_col(ann, c("cell_type", "celltype"), "Unknown")),
    Confidence = confidence,
    IsMixed = is_mixed,
    PrimaryCellType = as.character(get_col(ann, c("primary_cell_type", "primarycelltype"), "")),
    SecondaryCellType = as.character(get_col(ann, c("secondary_cell_type", "secondarycelltype"), "")),
    TissueConsistency = as.character(get_col(ann, c("tissue_consistency", "tissueconsistency"), "unknown")),
    Reasoning = as.character(get_col(ann, "reasoning", NA_character_)),
    stringsAsFactors = FALSE
  )
  
  result$Cluster <- trimws(result$Cluster)
  result$CellType <- trimws(result$CellType)
  result$PrimaryCellType <- trimws(result$PrimaryCellType)
  result$SecondaryCellType <- trimws(result$SecondaryCellType)
  result$TissueConsistency <- trimws(tolower(result$TissueConsistency))
  result$TissueConsistency[!result$TissueConsistency %in% c(
    "expected", "unexpected", "possible_contamination", "possible_doublet", "unknown"
  )] <- "unknown"
  result$CellType[!nzchar(result$CellType) | is.na(result$CellType)] <- "Unknown"
  
  result
}

#' Create system prompt for annotation task
#' 
#' @return System prompt string
#' @export
create_system_prompt <- function() {
  paste(
    "You are a senior bioinformatics scientist specializing in single-cell RNA-seq analysis.",
    "Annotate cell clusters from marker genes using rigorous marker-gene reasoning.",
    "",
    "Core principles:",
    "1. Identify the most likely biological cell type from the markers.",
    "2. Do not return Unknown only because the cell type is unexpected in the stated tissue.",
    "3. If markers indicate contamination, ambient RNA, or doublets, report the likely biological identity and flag it.",
    "4. Use Unknown only when marker evidence is biologically incoherent or insufficient.",
    "5. Use Cell Ontology-compatible names when possible, but do not invent ontology IDs.",
    "6. Use standardized cell type names suitable for publication.",
    "7. Return confidence scores between 0 and 1.",
    "8. Return only valid JSON.",
    "",
    "Required JSON fields:",
    "cluster, cell_type, confidence, is_mixed, primary_cell_type, secondary_cell_type, tissue_consistency, reasoning.",
    sep = "\n"
  )
}


#' Internal: parse line-by-line format (fallback)
#' @param response_text Raw response text.
#' @keywords internal
.parse_line_by_line <- function(response_text) {
  lines <- strsplit(response_text, "\n")[[1]]
  results <- list()
  
  for (line in lines) {
    # Pattern: "ClusterX: Cell Type [0.95]" or "ClusterX: Cell Type (0.95)"
    pattern <- "^(Cluster\\s*\\d+|[^:]+):\\s*([^\\[(]+?)\\s*[\\[(]([0-9.]+)[\\])]"
    matches <- regmatches(line, regexec(pattern, line))
    
    if (length(matches[[1]]) >= 4) {
      results[[length(results) + 1]] <- data.frame(
        Cluster = trimws(matches[[1]][2]),
        CellType = trimws(matches[[1]][3]),
        Confidence = as_confidence(matches[[1]][4]),
        IsMixed = FALSE,
        PrimaryCellType = trimws(matches[[1]][3]),
        SecondaryCellType = "",
        TissueConsistency = "unknown",
        Reasoning = NA_character_,
        stringsAsFactors = FALSE
      )
    }
  }
  
  if (length(results) > 0) {
    return(.normalise_annotation_dataframe(do.call(rbind, results)))
  }
  
  warning("Could not parse response: ", substr(response_text, 1, 200))
  .empty_annotation_result()
}
