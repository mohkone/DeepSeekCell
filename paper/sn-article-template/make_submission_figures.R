workflow_boxes <- data.frame(
  label = c(
    "Cluster markers\n+tissue/species",
    "Marker processing\n+ prompt construction",
    "DeepSeek API\nor local Ollama",
    "Structured JSON\nannotation parsing",
    "Cell Ontology\nmapping",
    "Validation,\nvisualization, export"
  ),
  x = c(0.10, 0.26, 0.42, 0.58, 0.74, 0.90),
  y = c(0.60, 0.60, 0.60, 0.60, 0.60, 0.60),
  stringsAsFactors = FALSE
)

draw_workflow <- function() {
  par(mar = c(0, 0, 0, 0), bg = "white")
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))

  title_col <- "#2c3e50"
  accent <- "#18bc9c"
  border <- "#2c3e50"
  fill <- "#f4f8f9"

  text(
    0.5, 0.92,
    "DeepSeekCell ontology-guided marker annotation workflow",
    cex = 1.35,
    font = 2,
    col = title_col
  )

  box_w <- 0.115
  box_h <- 0.22

  for (i in seq_len(nrow(workflow_boxes))) {
    x <- workflow_boxes$x[i]
    y <- workflow_boxes$y[i]
    rect(
      x - box_w / 2,
      y - box_h / 2,
      x + box_w / 2,
      y + box_h / 2,
      col = fill,
      border = border,
      lwd = 1.6
    )
    text(x, y, workflow_boxes$label[i], cex = 0.72, col = "#1f2d3a")

    if (i < nrow(workflow_boxes)) {
      arrows(
        x + box_w / 2 + 0.006,
        y,
        workflow_boxes$x[i + 1] - box_w / 2 - 0.006,
        workflow_boxes$y[i + 1],
        length = 0.08,
        lwd = 1.6,
        col = accent
      )
    }
  }

  rect(0.08, 0.16, 0.92, 0.32, col = "#eef6f2", border = "#27ae60", lwd = 1.4)
  text(
    0.50, 0.24,
    "Outputs: cell type label, confidence, tissue consistency, mixed-cluster flag,\nmarker-based reasoning, Cell Ontology ID, match method, downloadable CSV/Excel/HTML reports",
    cex = 0.76,
    col = "#1f2d3a"
  )
}

pdf("deepseekcell_workflow.pdf", width = 11, height = 5.2, useDingbats = FALSE)
draw_workflow()
dev.off()

png("deepseekcell_workflow.png", width = 2200, height = 1040, res = 200)
draw_workflow()
dev.off()
