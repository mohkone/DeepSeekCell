#' Visualization Functions
#' 
#' Publication-ready plotting functions for annotation results.

#' Create confidence bar plot
#' 
#' @param result Annotation result object
#' @param interactive Whether to return interactive plotly object
#' @return ggplot2 or plotly object
#' @export
plot_confidence <- function(result, interactive = FALSE) {
  
  if (!result$success || nrow(result$annotations) == 0) {
    return(.empty_plot())
  }
  
  df <- result$annotations
  df$Cluster <- factor(df$Cluster, levels = rev(unique(df$Cluster)))
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = Cluster, y = Confidence, fill = Confidence)) +
    ggplot2::geom_col(width = 0.7, alpha = 0.9) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", Confidence)), 
                       hjust = -0.1, size = 3.5) +
    ggplot2::scale_fill_gradient2(low = "#FF6B6B", mid = "#FFE66D", high = "#4ECDC4",
                                  midpoint = 0.5, limits = c(0, 1)) +
    ggplot2::coord_flip(ylim = c(0, 1)) +
    ggplot2::labs(
      title = "Cell Type Annotation Confidence",
      subtitle = sprintf("Model: %s | Tissue: %s | Species: %s",
                         result$metadata$model, 
                         result$metadata$tissue,
                         result$metadata$species),
      x = "Cluster",
      y = "Confidence Score",
      fill = "Confidence"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(color = "gray50", size = 10),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = 10),
      legend.position = "bottom"
    )
  
  if (interactive) {
    return(plotly::ggplotly(p, tooltip = c("x", "y", "fill")))
  }
  
  return(p)
}

#' Create ontology coverage plot
#' 
#' @param result Annotation result object
#' @return ggplot2 object
#' @export
plot_ontology_coverage <- function(result) {
  
  if (!result$success || !"MatchMethod" %in% colnames(result$annotations)) {
    return(.empty_plot())
  }
  
  df <- result$annotations
  match_methods <- table(df$MatchMethod)
  match_df <- data.frame(
    Method = names(match_methods),
    Count = as.numeric(match_methods)
  )
  
  ggplot2::ggplot(match_df, ggplot2::aes(x = stats::reorder(Method, -Count), 
                                         y = Count, fill = Method)) +
    ggplot2::geom_col() +
    ggplot2::geom_text(ggplot2::aes(label = Count), vjust = -0.5) +
    ggplot2::labs(
      title = "Ontology Mapping Quality",
      subtitle = "Distribution of matching methods",
      x = "Matching Method",
      y = "Number of Clusters"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")
}


#' Create confusion matrix heatmap
#' 
#' @param confusion_matrix Table from calculate_metrics()
#' @return ggplot2 object
#' @export
plot_confusion_matrix <- function(confusion_matrix) {
  
  cm_df <- as.data.frame(confusion_matrix)
  colnames(cm_df) <- c("Predicted", "True", "Count")
  
  ggplot2::ggplot(cm_df, ggplot2::aes(x = True, y = Predicted, fill = Count)) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(ggplot2::aes(label = Count), color = "white", size = 3) +
    ggplot2::scale_fill_gradient(low = "white", high = "#2c3e50") +
    ggplot2::labs(
      title = "Confusion Matrix",
      x = "True Label",
      y = "Predicted Label"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

#' Create summary dashboard
#' 
#' @param result Annotation result object
#' @return Patchwork object with multiple plots
#' @export
create_dashboard <- function(result) {
  
  p1 <- plot_confidence(result)
  p2 <- plot_ontology_coverage(result)
  
  # Combine using patchwork if available
  if (requireNamespace("patchwork", quietly = TRUE)) {
    return(p1 + p2 + patchwork::plot_layout(ncol = 1))
  }
  
  return(list(confidence_plot = p1, ontology_plot = p2))
}

.empty_plot <- function() {
  ggplot2::ggplot() + 
    ggplot2::annotate("text", x = 0.5, y = 0.5, 
                      label = "No data to display") +
    ggplot2::theme_void()
}