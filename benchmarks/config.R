# benchmarks/config.R

DEEPSEEK_KEY <- Sys.getenv("DEEPSEEK_API_KEY")

if (DEEPSEEK_KEY == "") warning("DEEPSEEK_API_KEY not set. DeepSeek will be skipped.")

MODELS <- list(
  deepseek = list(
    name = "DeepSeek",
    api_url = "https://api.deepseek.com/v1/chat/completions",
    model_id = "deepseek-chat",
    max_tokens = 2000,
    temperature = 0,
    input_cost_per_1k = 0.00014,
    output_cost_per_1k = 0.00028
  )
)

DATA_DIR <- "data"
ONTOLOGY_FILE <- file.path(DATA_DIR, "cl.obo")
SCTYPE_DB <- "scType/ScTypeDB_full.xlsx"

TOP_MARKERS <- 25
SEURAT_RESOLUTION <- 0.5
N_PCS <- 50
DEFAULT_BENCHMARK_REPLICATES <- 3
BENCHMARK_CACHE_VERSION <- "2026-06-12-publication-v5"
BENCHMARK_MODE <- "closed-label-marker-guided"

set.seed(42)
