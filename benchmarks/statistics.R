# Statistical summaries for benchmark outputs.

benchmark_pairwise_wilcoxon <- function(full,
                                        reference_method = "DeepSeek",
                                        metrics = c("ARI", "MacroF1", "Accuracy", "CladeAcc", "RuntimeSec")) {
  stopifnot(all(c("Dataset", "Replicate", "Method") %in% names(full)))

  methods <- setdiff(sort(unique(full$Method)), reference_method)
  rows <- list()

  for (metric in metrics) {
    if (!metric %in% names(full)) next

    for (method in methods) {
      pair_df <- full %>%
        dplyr::filter(Method %in% c(reference_method, method)) %>%
        dplyr::select(Dataset, Replicate, Method, Value = dplyr::all_of(metric)) %>%
        tidyr::pivot_wider(names_from = Method, values_from = Value)

      if (!all(c(reference_method, method) %in% names(pair_df))) next

      pair_df <- pair_df %>%
        dplyr::filter(!is.na(.data[[reference_method]]), !is.na(.data[[method]]))

      n_pairs <- nrow(pair_df)
      p_value <- NA_real_
      statistic <- NA_real_

      if (n_pairs >= 2 && length(unique(pair_df[[reference_method]] - pair_df[[method]])) > 1) {
        test <- tryCatch(
          stats::wilcox.test(
            pair_df[[reference_method]],
            pair_df[[method]],
            paired = TRUE,
            exact = FALSE
          ),
          error = function(e) NULL
        )

        if (!is.null(test)) {
          p_value <- unname(test$p.value)
          statistic <- unname(test$statistic)
        }
      }

      rows[[length(rows) + 1]] <- data.frame(
        Metric = metric,
        ReferenceMethod = reference_method,
        ComparatorMethod = method,
        Blocks = n_pairs,
        ReferenceMean = mean(pair_df[[reference_method]], na.rm = TRUE),
        ComparatorMean = mean(pair_df[[method]], na.rm = TRUE),
        MeanDifference = mean(pair_df[[reference_method]] - pair_df[[method]], na.rm = TRUE),
        MedianDifference = stats::median(pair_df[[reference_method]] - pair_df[[method]], na.rm = TRUE),
        WilcoxonStatistic = statistic,
        PValue = p_value,
        stringsAsFactors = FALSE
      )
    }
  }

  out <- dplyr::bind_rows(rows)
  if (nrow(out) == 0) return(out)

  out %>%
    dplyr::group_by(Metric) %>%
    dplyr::mutate(PValueBH = stats::p.adjust(PValue, method = "BH")) %>%
    dplyr::ungroup()
}

benchmark_friedman_tests <- function(full,
                                     metric_methods = list(
                                       all_datasets = c("DeepSeek", "scType", "SingleR", "scmap"),
                                       human_datasets_with_celltypist = c("DeepSeek", "scType", "SingleR", "scmap", "CellTypist")
                                     ),
                                     metrics = c("ARI", "MacroF1", "Accuracy", "CladeAcc", "RuntimeSec")) {
  rows <- list()

  for (set_name in names(metric_methods)) {
    methods <- metric_methods[[set_name]]

    for (metric in metrics) {
      if (!metric %in% names(full)) next

      wide <- full %>%
        dplyr::filter(Method %in% methods) %>%
        dplyr::select(Dataset, Replicate, Method, Value = dplyr::all_of(metric)) %>%
        dplyr::mutate(Block = paste(Dataset, Replicate, sep = "_rep")) %>%
        dplyr::select(Block, Method, Value) %>%
        tidyr::pivot_wider(names_from = Method, values_from = Value)

      if (!all(methods %in% names(wide))) next

      wide <- wide[stats::complete.cases(wide[, methods, drop = FALSE]), , drop = FALSE]
      n_blocks <- nrow(wide)
      p_value <- NA_real_
      statistic <- NA_real_

      if (n_blocks >= 2 && length(methods) >= 3) {
        mat <- as.matrix(wide[, methods, drop = FALSE])
        test <- tryCatch(stats::friedman.test(mat), error = function(e) NULL)

        if (!is.null(test)) {
          p_value <- unname(test$p.value)
          statistic <- unname(test$statistic)
        }
      }

      rows[[length(rows) + 1]] <- data.frame(
        MethodSet = set_name,
        Metric = metric,
        Methods = paste(methods, collapse = ";"),
        Blocks = n_blocks,
        FriedmanChiSquared = statistic,
        PValue = p_value,
        stringsAsFactors = FALSE
      )
    }
  }

  dplyr::bind_rows(rows)
}

benchmark_llm_stability <- function(debug_dir = file.path("results", "benchmark_debug"),
                                    method = "DeepSeek",
                                    fixed_datasets = c(
                                      "PBMC", "BaronPancreas", "MuraroPancreas",
                                      "TasicBrain", "ZeiselBrain"
                                    )) {
  files <- list.files(
    debug_dir,
    pattern = paste0("_", method, "_debug[.]csv$"),
    full.names = TRUE
  )

  if (length(files) == 0) {
    return(data.frame())
  }

  debug <- dplyr::bind_rows(lapply(files, utils::read.csv, stringsAsFactors = FALSE))
  debug <- debug %>% dplyr::filter(Dataset %in% fixed_datasets)

  if (nrow(debug) == 0) {
    return(data.frame())
  }

  cluster_stability <- debug %>%
    dplyr::group_by(Dataset, Cluster) %>%
    dplyr::summarise(
      Replicates = dplyr::n_distinct(Replicate),
      UniqueRawPredictions = dplyr::n_distinct(RawPrediction),
      UniqueHarmonisedPredictions = dplyr::n_distinct(HarmonisedPrediction),
      StableRaw = UniqueRawPredictions == 1,
      StableHarmonised = UniqueHarmonisedPredictions == 1,
      .groups = "drop"
    )

  cluster_stability %>%
    dplyr::group_by(Dataset) %>%
    dplyr::summarise(
      ComparableClusters = dplyr::n(),
      MeanReplicates = mean(Replicates),
      RawStabilityRate = mean(StableRaw),
      HarmonisedStabilityRate = mean(StableHarmonised),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      Method = method,
      StabilityScope = "fixed cluster structures only",
      .before = 1
    )
}

run_benchmark_statistical_tests <- function(full,
                                            debug_dir = file.path("results", "benchmark_debug")) {
  list(
    pairwise_wilcoxon = benchmark_pairwise_wilcoxon(full),
    friedman = benchmark_friedman_tests(full),
    llm_stability = benchmark_llm_stability(debug_dir = debug_dir)
  )
}
