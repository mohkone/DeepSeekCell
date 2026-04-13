#' Unified API Handler for Multiple LLM Providers
#' 
#' Supports DeepSeek and GPT-4o (via Ofox.ai proxy for China)
#' 
#' @import httr2 jsonlite

# Model configurations
MODELS <- list(
  deepseek = list(
    name = "DeepSeek-Chat",
    api_url = "https://api.deepseek.com/v1/chat/completions",
    model_id = "deepseek-chat",
    max_tokens = 2000,   # Increased from 500
    temperature = 0.1,
    cost_per_1k_tokens = 0.00014
  ),
  gpt4 = list(
    name = "GPT-4o",
    api_url = "https://api.ofox.ai/v1/chat/completions",
    model_id = "openai/gpt-4o",
    max_tokens = 2000,   # Increased from 500
    temperature = 0.1,
    cost_per_1k_tokens = 0.0025
  )
)

#' Get model configuration by name
#' 
#' @param model_name Model identifier ("deepseek" or "gpt4")
#' @return Model configuration list or NULL if not found
#' @export
get_model_config <- function(model_name) {
  return(MODELS[[model_name]])
}

#' Call LLM API with retry logic and error handling
#' 
#' @param prompt User prompt for annotation
#' @param model Model configuration list
#' @param api_key API key for the service
#' @param max_retries Maximum number of retry attempts
#' @return List containing response content and metadata
#' @export
call_llm_api <- function(prompt, model, api_key, max_retries = 3) {
  
  message("Calling ", model$name, " API...")
  start_time <- Sys.time()
  
  for (attempt in seq_len(max_retries)) {
    
    tryCatch({
      
      req <- httr2::request(model$api_url) |>
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
      data <- httr2::resp_body_json(resp)
      
      content <- data$choices[[1]]$message$content
      usage <- data$usage
      
      elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      
      message("API call successful in ", elapsed, " seconds. ",
              "Tokens used: ", usage$total_tokens)
      
      return(list(
        content = content,
        usage = usage,
        model = model$name,
        latency_sec = elapsed,
        attempt = attempt,
        success = TRUE
      ))
      
    }, error = function(e) {
      warning("API attempt ", attempt, " failed: ", e$message)
      
      if (attempt < max_retries) {
        Sys.sleep(2 * attempt)  # Exponential backoff
      } else {
        message("All ", max_retries, " API attempts failed")
        return(list(
          content = NULL,
          error = e$message,
          success = FALSE,
          model = model$name
        ))
      }
    })
  }
}

#' Create system prompt for annotation task
#' 
#' @return System prompt string
#' @export
create_system_prompt <- function() {
  paste(
    "You are a Senior Bioinformatics Scientist specializing in single-cell RNA-seq analysis.",
    "Your task is to annotate cell clusters based on their marker genes.",
    "",
    "Rules:",
    "1. Provide reasoning for each annotation",
    "2. Output confidence scores between 0 and 1",
    "3. Use standardized cell type names",
    "4. If uncertain, suggest 'Unknown' with low confidence",
    "",
    "Output format (JSON):",
    '{"annotations": [',
    '  {"cluster": "Cluster1", "cell_type": "T cell", "confidence": 0.95, "reasoning": "CD3D and CD3E are T cell markers"},',
    '  ...',
    ']}',
    sep = "\n"
  )
}

#' Parse LLM response to extract annotations
#' 
#' @param response_text Raw response from API
#' @return Data frame with cluster, cell_type, confidence, reasoning
#' @export
parse_annotation_response <- function(response_text) {
  
  if (is.null(response_text) || response_text == "") {
    return(data.frame())
  }
  
  # Step 1: Remove markdown code fences (```json ... ```)
  # Remove opening ```json or ``` (with optional newline)
  cleaned <- gsub("^```(?:json)?\\s*\n?", "", response_text, perl = TRUE)
  # Remove closing ``` at the end (with optional whitespace before)
  cleaned <- gsub("\n?\\s*```$", "", cleaned, perl = TRUE)
  
  # Step 2: Find the first '{' and the last '}' in the cleaned string
  first_brace <- regexpr("\\{", cleaned, perl = TRUE)
  last_brace <- regexpr("\\}[^}]*$", cleaned, perl = TRUE)  # position of last '}'
  
  if (first_brace == -1 || last_brace == -1) {
    warning("No JSON object found in response")
    return(.parse_line_by_line(response_text))
  }
  
  # Extract JSON substring
  json_str <- substr(cleaned, first_brace, last_brace + attr(last_brace, "match.length") - 1)
  
  # Step 3: Parse JSON
  tryCatch({
    parsed <- jsonlite::fromJSON(json_str, simplifyVector = TRUE)
    
    # Extract annotations array
    annotations_list <- if (!is.null(parsed$annotations)) parsed$annotations else parsed
    
    # Convert to data frame
    if (is.data.frame(annotations_list)) {
      df <- annotations_list
    } else if (is.list(annotations_list) && length(annotations_list) > 0) {
      # Flatten list of lists into data frame
      df <- do.call(rbind, lapply(annotations_list, function(x) {
        as.data.frame(t(unlist(x)), stringsAsFactors = FALSE)
      }))
    } else {
      df <- NULL
    }
    
    if (!is.null(df) && nrow(df) > 0) {
      # Standardize column names
      names(df) <- tolower(names(df))
      
      # Map expected columns
      result <- data.frame(
        Cluster = if ("cluster" %in% names(df)) df$cluster else paste0("Cluster", 1:nrow(df)),
        CellType = if ("cell_type" %in% names(df)) df$cell_type else 
          if ("celltype" %in% names(df)) df$celltype else "Unknown",
        Confidence = if ("confidence" %in% names(df)) as.numeric(df$confidence) else 0.5,
        Reasoning = if ("reasoning" %in% names(df)) df$reasoning else NA,
        stringsAsFactors = FALSE
      )
      return(result)
    }
  }, error = function(e) {
    warning("JSON parsing failed: ", e$message)
  })
  
  # Fallback: line-by-line parsing
  .parse_line_by_line(response_text)
}

#' Internal: parse line-by-line format (fallback)
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
        Confidence = as.numeric(matches[[1]][4]),
        Reasoning = NA,
        stringsAsFactors = FALSE
      )
    }
  }
  
  if (length(results) > 0) {
    return(do.call(rbind, results))
  }
  
  warning("Could not parse response: ", substr(response_text, 1, 200))
  return(data.frame())
}