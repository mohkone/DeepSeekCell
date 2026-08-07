# benchmarks/external_datasets/profile_scrnaseq_candidates.R
#
# Profile candidate scRNAseq datasets before adding them to the locked external
# validation registry. This script records dimensions, metadata columns, and
# plausible truth-label columns without creating benchmark results.

script_path <- tryCatch(sys.frame(1)$ofile, error = function(e) NA_character_)
if (is.null(script_path) || is.na(script_path) || !nzchar(script_path)) {
  args0 <- commandArgs(trailingOnly = FALSE)
  file_args <- args0[grepl("^--file=", args0)]
  script_path <- if (length(file_args) > 0) sub("^--file=", "", file_args[1]) else ""
}
repo_root <- if (nzchar(script_path)) {
  normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = FALSE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}
setwd(repo_root)

candidate_calls <- list(
  KotliarovPBMC = "KotliarovPBMCData(mode = \"rna\", ensembl = FALSE, location = FALSE)",
  MairPBMC = "MairPBMCData(mode = \"rna\", ensembl = FALSE, location = FALSE)",
  BacherTCell = "BacherTCellData(filtered = TRUE, ensembl = FALSE, location = FALSE)",
  BunisHSPC = "BunisHSPCData(filtered = TRUE)",
  GiladiHSC = "GiladiHSCData(mode = \"rna\", filtered = TRUE, ensembl = FALSE, location = FALSE)",
  NestorowaHSC = "NestorowaHSCData(remove.htseq = TRUE, location = FALSE)",
  PaulHSC = "PaulHSCData(ensembl = FALSE, discard.multiple = TRUE, location = FALSE)",
  RichardTCell = "RichardTCellData(location = FALSE)",
  ReprocessedTh2 = "ReprocessedTh2Data(ensembl = FALSE, location = FALSE)",
  StoeckiusPBMC = "StoeckiusHashingData(type = \"pbmc\", ensembl = FALSE, location = FALSE)",
  SegerstolpePancreas = "SegerstolpePancreasData(ensembl = FALSE, location = FALSE)",
  LawlorPancreas = "LawlorPancreasData()",
  XinPancreas = "XinPancreasData(ensembl = FALSE, location = FALSE)",
  GrunPancreas = "GrunPancreasData(ensembl = FALSE, location = FALSE)",
  DarmanisBrain = "DarmanisBrainData(ensembl = FALSE, location = FALSE)",
  MarquesBrain = "MarquesBrainData(ensembl = FALSE, location = FALSE)",
  LaMannoBrainAdult = "LaMannoBrainData(which = \"mouse-adult\", ensembl = FALSE, location = FALSE)",
  LaMannoBrainEmbryo = "LaMannoBrainData(which = \"mouse-embryo\", ensembl = FALSE, location = FALSE)",
  RomanovBrain = "RomanovBrainData(ensembl = FALSE, location = FALSE)",
  ZhongPrefrontal = "ZhongPrefrontalData(ensembl = FALSE, location = FALSE)",
  ChenBrain = "ChenBrainData(ensembl = FALSE, location = FALSE)",
  HuCortex = "HuCortexData(mode = \"ctx\", ensembl = FALSE, location = FALSE)",
  CampbellBrain = "CampbellBrainData(ensembl = FALSE, location = FALSE)",
  NowakowskiCortex = "NowakowskiCortexData(ensembl = FALSE, location = FALSE)",
  PollenGlia = "PollenGliaData(ensembl = FALSE, location = FALSE)",
  MacoskoRetina = "MacoskoRetinaData(ensembl = FALSE, location = FALSE)",
  WuKidneyHealthy = "WuKidneyData(mode = \"healthy\", ensembl = FALSE, location = FALSE)",
  ZhaoImmuneLiver = "ZhaoImmuneLiverData(location = FALSE, filter = TRUE)",
  BachMammary = "BachMammaryData(samples = \"NP_1\", location = FALSE)",
  LedergorMyeloma = "LedergorMyelomaData(ensembl = FALSE, location = FALSE)",
  FletcherOlfactory = "FletcherOlfactoryData(filtered = TRUE, ensembl = FALSE, location = FALSE)"
)

load_candidate <- function(call_text) {
  if (!requireNamespace("scRNAseq", quietly = TRUE)) {
    stop("The scRNAseq package is required to profile candidate datasets.", call. = FALSE)
  }
  eval(parse(text = paste0("scRNAseq::", call_text)))
}

profile_candidate <- function(dataset, call_text) {
  started <- Sys.time()
  object <- tryCatch(
    suppressPackageStartupMessages(load_candidate(call_text)),
    error = function(e) e
  )
  runtime_sec <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  if (inherits(object, "error")) {
    return(data.frame(
      Dataset = dataset,
      SourceCall = call_text,
      Status = "error",
      Error = conditionMessage(object),
      NGenes = NA_integer_,
      NCells = NA_integer_,
      MetadataColumns = NA_character_,
      CandidateLabelColumns = NA_character_,
      CandidateLabelSummary = NA_character_,
      RuntimeSec = runtime_sec,
      stringsAsFactors = FALSE
    ))
  }
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("SummarizedExperiment is required to inspect scRNAseq datasets.", call. = FALSE)
  }
  metadata <- as.data.frame(SummarizedExperiment::colData(object))
  unique_counts <- vapply(
    metadata,
    function(column) length(unique(as.character(column[!is.na(column)]))),
    numeric(1)
  )
  plausible <- names(unique_counts)[
    unique_counts >= 2 &
      unique_counts <= 100 &
      !grepl(
        "id|sample|donor|individual|batch|well|lane|age|sex|barcode|percent|n_|total|sum|size|run",
        names(unique_counts),
        ignore.case = TRUE
      )
  ]
  if (length(plausible) == 0) {
    plausible <- names(unique_counts)[unique_counts >= 2 & unique_counts <= 100]
  }
  plausible <- head(plausible, 10)
  label_summary <- vapply(plausible, function(column) {
    tab <- sort(table(as.character(metadata[[column]])), decreasing = TRUE)
    paste0(
      column, " (n=", unique_counts[[column]], "): ",
      paste(names(head(tab, 8)), collapse = " | ")
    )
  }, character(1))
  data.frame(
    Dataset = dataset,
    SourceCall = call_text,
    Status = "ok",
    Error = "",
    NGenes = nrow(object),
    NCells = ncol(object),
    MetadataColumns = paste(names(metadata), collapse = ";"),
    CandidateLabelColumns = paste(plausible, collapse = ";"),
    CandidateLabelSummary = paste(label_summary, collapse = " || "),
    RuntimeSec = runtime_sec,
    stringsAsFactors = FALSE
  )
}

profile_scrnaseq_candidates <- function(output_path = file.path(
                                          "results",
                                          "external_dataset_audits",
                                          "scrnaseq_candidate_profile.csv"
                                        )) {
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  selected <- trimws(strsplit(Sys.getenv("DEEPSEEKCELL_EXTERNAL_DATASETS", unset = ""), ",", fixed = TRUE)[[1]])
  selected <- selected[nzchar(selected)]
  calls <- candidate_calls
  if (length(selected) > 0) {
    missing <- setdiff(selected, names(candidate_calls))
    if (length(missing) > 0) {
      stop(
        "DEEPSEEKCELL_EXTERNAL_DATASETS contains unknown scRNAseq candidate(s): ",
        paste(missing, collapse = ", "),
        call. = FALSE
      )
    }
    calls <- candidate_calls[selected]
  }
  rows <- lapply(names(calls), function(dataset) {
    message("Profiling ", dataset)
    row <- profile_candidate(dataset, calls[[dataset]])
    utils::write.csv(
      if (file.exists(output_path)) rbind(utils::read.csv(output_path, stringsAsFactors = FALSE), row) else row,
      output_path,
      row.names = FALSE
    )
    row
  })
  profiles <- do.call(rbind, rows)
  message("Wrote scRNAseq candidate profile to ", output_path)
  invisible(profiles)
}

args <- commandArgs(trailingOnly = TRUE)
output_path <- if (length(args) >= 1) {
  args[[1]]
} else {
  file.path("results", "external_dataset_audits", "scrnaseq_candidate_profile.csv")
}

profile_scrnaseq_candidates(output_path)
