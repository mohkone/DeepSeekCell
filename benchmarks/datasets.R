# benchmarks/datasets.R
# Dataset loaders for PBMC, Pancreas, Brain

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratData)
  library(scRNAseq)
  library(SingleCellExperiment)
})

get_pbmc_data <- function() {
  InstallData("pbmc3k")
  pbmc <- LoadData("pbmc3k")
  pbmc <- NormalizeData(pbmc) |>
    FindVariableFeatures() |>
    ScaleData() |>
    RunPCA() |>
    FindNeighbors() |>
    FindClusters(resolution = 0.5)
  
  markers <- FindAllMarkers(pbmc, only.pos = TRUE, logfc.threshold = 0.75)
  markers <- markers[!grepl("^MT-|^RP[SL]|^MALAT", markers$gene), ]
  markers_list <- split(markers$gene, markers$cluster)
  markers_list <- lapply(markers_list, head, 25)
  
  true_labels <- pbmc$seurat_annotations
  names(true_labels) <- rownames(pbmc@meta.data)
  cluster_truth <- sapply(levels(Idents(pbmc)), function(cl) {
    cells <- WhichCells(pbmc, idents = cl)
    if (length(cells) == 0) return("Unknown")
    raw <- names(sort(table(true_labels[cells]), decreasing = TRUE))[1]
    if (grepl("CD4", raw)) return("Naive T cell")
    if (grepl("CD8", raw)) return("Cytotoxic T cell")
    if (grepl("CD14", raw)) return("Classical monocyte")
    if (grepl("FCGR3A", raw)) return("Non-classical monocyte")
    if (grepl("B", raw)) return("B cell")
    if (grepl("NK", raw)) return("Natural killer cell")
    if (grepl("DC", raw)) return("Dendritic cell")
    if (grepl("Platelet", raw)) return("Platelet")
    return("Unknown")
  })
  names(cluster_truth) <- levels(Idents(pbmc))
  cluster_truth <- cluster_truth[names(markers_list)]
  
  list(markers = markers_list, truth = cluster_truth,
       tissue = "PBMC", species = "Human", seurat_obj = pbmc)
}

get_pancreas_data <- function() {
  message("Loading pancreas data (BaronPancreasData)...")
  sce <- scRNAseq::BaronPancreasData()
  seurat_obj <- CreateSeuratObject(counts = assay(sce),
                                   meta.data = as.data.frame(colData(sce)))
  seurat_obj <- NormalizeData(seurat_obj) |>
    FindVariableFeatures() |>
    ScaleData() |>
    RunPCA() |>
    FindNeighbors() |>
    FindClusters(resolution = 0.5)
  
  markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, logfc.threshold = 0.75)
  markers <- markers[!grepl("^MT-|^RP[SL]|^MALAT", markers$gene), ]
  markers_list <- split(markers$gene, markers$cluster)
  markers_list <- lapply(markers_list, head, 25)
  
  true_labels <- seurat_obj$label
  names(true_labels) <- colnames(seurat_obj)
  cluster_truth <- sapply(levels(seurat_obj), function(cl) {
    cells <- WhichCells(seurat_obj, idents = cl)
    if (length(cells) == 0) return("Unknown")
    raw <- names(sort(table(true_labels[cells]), decreasing = TRUE))[1]
    switch(raw,
           "beta" = "Beta cell",
           "alpha" = "Alpha cell",
           "delta" = "Delta cell",
           "acinar" = "Acinar cell",
           "ductal" = "Ductal cell",
           "endothelial" = "Endothelial cell",
           "activated_stellate" = "Stellate cell",
           "quiescent_stellate" = "Stellate cell",
           "macrophage" = "Macrophage",
           "immune" = "Immune cell",
           "enteroendocrine" = "Enteroendocrine cell",
           "pericyte" = "Pericyte",
           "Unknown")
  })
  names(cluster_truth) <- levels(seurat_obj)
  cluster_truth <- cluster_truth[names(markers_list)]
  
  list(markers = markers_list, truth = cluster_truth,
       tissue = "Pancreas", species = "Human", seurat_obj = seurat_obj)
}

get_brain_data <- function() {
  message("Loading brain data (ZeiselBrainData)...")
  sce <- scRNAseq::ZeiselBrainData()
  seurat_obj <- CreateSeuratObject(counts = assay(sce),
                                   meta.data = as.data.frame(colData(sce)))
  seurat_obj <- NormalizeData(seurat_obj) |>
    FindVariableFeatures() |>
    ScaleData() |>
    RunPCA() |>
    FindNeighbors() |>
    FindClusters(resolution = 0.5)
  
  markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, logfc.threshold = 0.75)
  markers <- markers[!grepl("^MT-|^RP[SL]|^MALAT", markers$gene), ]
  markers_list <- split(markers$gene, markers$cluster)
  markers_list <- lapply(markers_list, head, 25)
  
  true_labels <- seurat_obj$level2class
  names(true_labels) <- colnames(seurat_obj)
  cluster_truth <- sapply(levels(seurat_obj), function(cl) {
    cells <- WhichCells(seurat_obj, idents = cl)
    if (length(cells) == 0) return("Unknown")
    raw <- names(sort(table(true_labels[cells]), decreasing = TRUE))[1]
    if (grepl("Pyr", raw)) return("Excitatory neuron")
    if (grepl("Int", raw)) return("Inhibitory neuron")
    if (grepl("Oligo", raw)) return("Oligodendrocyte")
    if (grepl("Astro", raw)) return("Astrocyte")
    if (grepl("Micro", raw)) return("Microglia")
    if (grepl("Endo", raw)) return("Endothelial cell")
    if (grepl("Vsmc", raw)) return("Smooth muscle cell")
    if (grepl("Peri", raw)) return("Pericyte")
    if (grepl("Epend", raw)) return("Ependymal cell")
    if (grepl("OPC", raw)) return("Oligodendrocyte precursor cell")
    return("Unknown")
  })
  names(cluster_truth) <- levels(seurat_obj)
  cluster_truth <- cluster_truth[names(markers_list)]
  
  list(markers = markers_list, truth = cluster_truth,
       tissue = "Brain", species = "Mouse", seurat_obj = seurat_obj)
}