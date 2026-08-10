# benchmarks/config.R

DEEPSEEK_KEY <- Sys.getenv("DEEPSEEK_API_KEY")

if (DEEPSEEK_KEY == "") warning("DEEPSEEK_API_KEY not set. DeepSeek will be skipped.")

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

MODELS <- list(
  deepseek = list(
    name = "DeepSeek",
    api_url = "https://api.deepseek.com/v1/chat/completions",
    model_id = "deepseek-chat",
    model_id_env = "DEEPSEEK_MODEL_ID",
    api_url_env = "DEEPSEEK_API_URL",
    max_tokens = 2000,
    temperature = 0,
    input_cost_per_1k = 0.00014,
    output_cost_per_1k = 0.00028,
    requires_api_key = TRUE,
    api_key_env = "DEEPSEEK_API_KEY",
    is_ollama = FALSE
  ),
  ollama = list(
    name = "Ollama",
    api_url = "http://localhost:11434/api/generate",
    api_url_env = "OLLAMA_API_URL",
    model_id = "llama3.2:latest",
    model_id_env = "OLLAMA_MODEL_ID",
    max_tokens = 2000,
    temperature = 0,
    input_cost_per_1k = 0,
    output_cost_per_1k = 0,
    requires_api_key = FALSE,
    api_key_env = character(),
    is_ollama = TRUE
  ),
  openai = list(
    name = "OpenAI",
    api_url = "https://api.openai.com/v1/responses",
    api_url_env = "OPENAI_API_URL",
    model_id = "gpt-5",
    model_id_env = "OPENAI_MODEL_ID",
    max_tokens = 6000,
    max_tokens_env = "OPENAI_MAX_OUTPUT_TOKENS",
    temperature = NULL,
    reasoning_effort = "minimal",
    reasoning_effort_env = "OPENAI_REASONING_EFFORT",
    text_verbosity = "low",
    text_verbosity_env = "OPENAI_TEXT_VERBOSITY",
    input_cost_per_1k = 0.00125,
    output_cost_per_1k = 0.01000,
    pricing_date = "2026-08-11",
    pricing_source = "https://developers.openai.com/api/docs/models/gpt-5",
    requires_api_key = TRUE,
    api_key_env = "OPENAI_API_KEY",
    is_ollama = FALSE,
    api_format = "responses"
  )
)

for (model_name in names(MODELS)) {
  api_url_override <- Sys.getenv(MODELS[[model_name]]$api_url_env %||% "", unset = "")
  model_id_override <- Sys.getenv(MODELS[[model_name]]$model_id_env %||% "", unset = "")
  max_tokens_override <- Sys.getenv(MODELS[[model_name]]$max_tokens_env %||% "", unset = "")
  reasoning_effort_override <- Sys.getenv(MODELS[[model_name]]$reasoning_effort_env %||% "", unset = "")
  text_verbosity_override <- Sys.getenv(MODELS[[model_name]]$text_verbosity_env %||% "", unset = "")

  if (nzchar(api_url_override)) {
    MODELS[[model_name]]$api_url <- api_url_override
  }
  if (nzchar(model_id_override)) {
    MODELS[[model_name]]$model_id <- model_id_override
  }
  if (nzchar(max_tokens_override)) {
    parsed_max_tokens <- suppressWarnings(as.numeric(max_tokens_override))
    if (!is.na(parsed_max_tokens) && parsed_max_tokens > 0) {
      MODELS[[model_name]]$max_tokens <- parsed_max_tokens
    }
  }
  if (nzchar(reasoning_effort_override)) {
    MODELS[[model_name]]$reasoning_effort <- reasoning_effort_override
  }
  if (nzchar(text_verbosity_override)) {
    MODELS[[model_name]]$text_verbosity <- text_verbosity_override
  }
}

DATA_DIR <- "data"
ONTOLOGY_FILE <- file.path(DATA_DIR, "cl.obo")
SCTYPE_DB <- "scType/ScTypeDB_full.xlsx"
CELLTYPIST_MODEL <- Sys.getenv("CELLTYPIST_MODEL", unset = "")

TOP_MARKERS <- 25
SEURAT_RESOLUTION <- 0.5
N_PCS <- 50
DEFAULT_BENCHMARK_REPLICATES <- 3
BENCHMARK_CACHE_VERSION <- "2026-06-12-publication-v5"
BENCHMARK_MODE <- "closed-label-marker-guided"

set.seed(42)
