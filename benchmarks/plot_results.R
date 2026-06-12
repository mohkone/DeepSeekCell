# benchmarks/plot_results.R
# Generate figures from benchmark results

suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
})

plot_benchmark <- function(results_csv = "benchmark_results.csv",
                           output_pdf = "benchmark_plot.pdf") {
  results <- read.csv(results_csv)
  df_plot <- reshape2::melt(results, id.vars = c("Dataset", "Method"), variable.name = "Metric")
  p <- ggplot(df_plot, aes(x = Method, y = value, fill = Method)) +
    geom_col(position = "dodge") +
    facet_grid(Dataset ~ Metric, scales = "free") +
    theme_minimal() +
    labs(title = "DeepSeekCell Benchmark Results (with scType + Ontology)")
  ggsave(output_pdf, p, width = 14, height = 8)
  message("Plot saved to ", output_pdf)
  invisible(p)
}


