suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

summary_path <- file.path("results", "benchmark_results_summary.csv")
if (!file.exists(summary_path)) {
  summary_path <- "benchmark_results_summary.csv"
}
summary_df <- read.csv(summary_path)

plot_df <- summary_df %>%
  select(
    Dataset,
    Method,
    ARI_mean,
    MacroF1_mean,
    Accuracy_mean,
    RuntimeSec_mean
  ) %>%
  pivot_longer(
    cols = c(ARI_mean, MacroF1_mean, Accuracy_mean, RuntimeSec_mean),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = recode(
      Metric,
      ARI_mean = "ARI",
      MacroF1_mean = "Macro-F1",
      Accuracy_mean = "Accuracy",
      RuntimeSec_mean = "Runtime"
    )
  )

p <- ggplot(plot_df, aes(x = Dataset, y = Value, fill = Method)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 2) +
  theme_minimal(base_size = 13) +
  labs(
    title = "Global Benchmark Comparison: DeepSeekCell vs SingleR vs scType",
    x = "Dataset",
    y = "Mean value"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

ggsave(file.path("results", "global_benchmark_comparison.pdf"), p, width = 14, height = 8)
ggsave(file.path("results", "global_benchmark_comparison.png"), p, width = 14, height = 8, dpi = 300)

print(p)
