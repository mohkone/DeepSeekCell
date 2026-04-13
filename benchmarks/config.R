# benchmarks/config.R
# Shared configuration for the benchmark

# API keys – set in environment before running
DEEPSEEK_KEY <- Sys.getenv("DEEPSEEK_API_KEY")
OFOX_KEY     <- Sys.getenv("OFOX_API_KEY")

if (DEEPSEEK_KEY == "" || OFOX_KEY == "") {
  warning("API keys not set. LLM methods will be skipped.")
}

MODELS <- list(
  deepseek = list(
    api_url = "https://api.deepseek.com/v1/chat/completions",
    model_id = "deepseek-chat",
    max_tokens = 2000,
    temperature = 0.1
  ),
  gpt4 = list(
    api_url = "https://api.ofox.ai/v1/chat/completions",
    model_id = "openai/gpt-4o",
    max_tokens = 2000,
    temperature = 0.1
  )
)

# Paths – relative to project root
DATA_DIR <- "data"
ONTOLOGY_FILE <- file.path(DATA_DIR, "cl.obo")
SCTYPE_DB <- "scType/ScTypeDB_full.xlsx"

set.seed(42)