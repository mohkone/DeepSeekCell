# benchmarks/datasets.R

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratData)
  library(scRNAseq)
  library(SingleCellExperiment)
  library(Matrix)
  library(dplyr)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

.standardize_markers <- function(markers, top_n = TOP_MARKERS) {
  
  markers <- markers[
    !grepl(
      "^MT-|^MTRNR|^RP[SL]|^MALAT|^XIST|^HB[AB]|^RP11-|^CTD-|^AC[0-9]|^AL[0-9]|^LINC",
      markers$gene,
      ignore.case = TRUE
    ),
  ]
  
  if ("p_val_adj" %in% colnames(markers)) {
    markers <- markers[markers$p_val_adj < 0.05, ]
  }
  
  markers <- markers[order(
    markers$cluster,
    -markers$pct.1,
    -markers$avg_log2FC
  ), ]
  
  markers_list <- split(markers$gene, markers$cluster)
  
  markers_list <- lapply(markers_list, function(x) {
    head(unique(x), top_n)
  })
  
  markers_list[lengths(markers_list) > 0]
}

.process_seurat <- function(seu,
                            truth_vector,
                            tissue,
                            species,
                            dataset_name,
                            resolution = SEURAT_RESOLUTION,
                            npcs = N_PCS,
                            top_n = TOP_MARKERS) {
  
  seu <- NormalizeData(seu, verbose = FALSE)
  seu <- FindVariableFeatures(seu, verbose = FALSE)
  seu <- ScaleData(seu, verbose = FALSE)
  seu <- RunPCA(seu, npcs = npcs, verbose = FALSE)
  seu <- FindNeighbors(seu, dims = seq_len(min(npcs, 30)), verbose = FALSE)
  seu <- FindClusters(seu, resolution = resolution, verbose = FALSE)
  
  markers <- FindAllMarkers(
    seu,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25,
    test.use = "wilcox",
    verbose = FALSE
  )
  
  markers_list <- .standardize_markers(markers, top_n = top_n)
  
  names(truth_vector) <- colnames(seu)
  seu$true_label <- truth_vector[colnames(seu)]
  
  cluster_truth <- sapply(names(markers_list), function(cl) {
    cells <- WhichCells(seu, idents = cl)
    labs <- truth_vector[cells]
    labs <- labs[!is.na(labs) & labs != ""]
    if (length(labs) == 0) return("Unknown")
    names(sort(table(labs), decreasing = TRUE))[1]
  })
  
  purities <- sapply(names(markers_list), function(cl) {
    cells <- WhichCells(seu, idents = cl)
    labs <- truth_vector[cells]
    labs <- labs[!is.na(labs) & labs != ""]
    if (length(labs) == 0) return(NA_real_)
    max(table(labs)) / length(labs)
  })
  
  list(
    markers = markers_list,
    truth = cluster_truth[names(markers_list)],
    tissue = tissue,
    species = species,
    dataset_name = dataset_name,
    seurat_obj = seu,
    purity = purities,
    cluster_summary = data.frame(
      Dataset = dataset_name,
      Tissue = tissue,
      Species = species,
      Cluster = names(markers_list),
      Truth = as.character(cluster_truth[names(markers_list)]),
      NMarkers = lengths(markers_list),
      ClusterPurity = as.numeric(purities[names(markers_list)]),
      stringsAsFactors = FALSE
    ),
    dataset_summary = data.frame(
      Dataset = dataset_name,
      Tissue = tissue,
      Species = species,
      NCells = ncol(seu),
      NGenes = nrow(seu),
      NClusters = length(markers_list),
      MeanClusterPurity = mean(purities, na.rm = TRUE),
      MedianClusterPurity = median(purities, na.rm = TRUE),
      MinClusterPurity = min(purities, na.rm = TRUE),
      TopMarkers = top_n,
      SeuratResolution = resolution,
      NPCs = npcs,
      stringsAsFactors = FALSE
    )
  )
}

.is_nk_label <- function(x) {
  x <- tolower(trimws(as.character(x)))
  grepl("natural killer", x) |
    grepl("(^|[^a-z0-9])[bt]?nk([^a-z0-9]|$)", x)
}

.is_unknown_label <- function(x) {
  x <- tolower(trimws(as.character(x)))
  unknown <- grepl(
    "^(unknown|unidentified|undetermined|not determined|ambiguous)(\\b|$)",
    x
  )
  unknown[is.na(unknown)] <- FALSE

  is.na(x) | x == "" | x %in% c("na", "nan") | unknown
}

.clean_label <- function(x) {
  x <- trimws(as.character(x))
  x[.is_unknown_label(x)] <- "Unknown"
  x
}

# -----------------------------------------------------------------------------
# 1. PBMC3k
# -----------------------------------------------------------------------------
get_pbmc_data <- function(seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  pbmc <- SeuratData::LoadData("pbmc3k")
  
  raw <- pbmc$seurat_annotations
  raw <- .clean_label(raw)
  
  truth <- sapply(raw, function(x) {
    x_low <- tolower(x)
    if (grepl("cd4", x_low)) return("Naive T cell")
    if (grepl("cd8", x_low)) return("Cytotoxic T cell")
    if (grepl("fcgr3a", x_low)) return("Non-classical monocyte")
    if (grepl("cd14|monocyte", x_low)) return("Classical monocyte")
    if (grepl("^b$", x_low) || grepl("b cell", x_low)) return("B cell")
    if (.is_nk_label(x_low)) return("Natural killer cell")
    if (grepl("dc|dendritic", x_low)) return("Dendritic cell")
    if (grepl("platelet", x_low)) return("Platelet")
    "Unknown"
  })
  
  .process_seurat(pbmc, truth, "PBMC", "Human", "PBMC")
}

# -----------------------------------------------------------------------------
# 2. Baron Pancreas
# -----------------------------------------------------------------------------
get_baron_pancreas_data <- function(seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  sce <- scRNAseq::BaronPancreasData()
  seu <- CreateSeuratObject(counts = assay(sce), meta.data = as.data.frame(colData(sce)))
  
  raw <- .clean_label(seu$label)
  
  truth <- dplyr::recode(
    raw,
    "beta" = "Beta cell",
    "alpha" = "Alpha cell",
    "delta" = "Delta cell",
    "gamma" = "Gamma cell",
    "epsilon" = "Epsilon cell",
    "acinar" = "Acinar cell",
    "ductal" = "Ductal cell",
    "endothelial" = "Endothelial cell",
    "activated_stellate" = "Stellate cell",
    "quiescent_stellate" = "Stellate cell",
    "macrophage" = "Macrophage",
    "immune" = "Immune cell",
    "enteroendocrine" = "Enteroendocrine cell",
    "pericyte" = "Pericyte",
    .default = "Unknown"
  )
  
  .process_seurat(seu, truth, "Pancreas", "Human", "BaronPancreas")
}

# -----------------------------------------------------------------------------
# 3. Muraro Pancreas
# -----------------------------------------------------------------------------
get_muraro_pancreas_data <- function(seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  sce <- scRNAseq::MuraroPancreasData()
  sce <- sce[, !is.na(sce$label) & sce$label != "unclear"]
  
  counts_mat <- counts(sce)
  
  gene_names <- rownames(counts_mat)
  
  gene_names <- gsub("__.*$", "", gene_names)
  gene_names <- gsub("--.*$", "", gene_names)
  
  gene_names <- make.unique(gene_names)
  
  rownames(counts_mat) <- gene_names
  
  colnames(counts_mat) <- make.unique(as.character(colnames(counts_mat)))
  
  message("Muraro cleaned genes: ", paste(head(rownames(counts_mat)), collapse = ", "))
  message("Muraro marker overlap test: ",
          sum(c("INS", "GCG", "SST", "KRT19", "PRSS1", "CPA1") %in% rownames(counts_mat)))
  
  md <- as.data.frame(colData(sce))
  rownames(md) <- colnames(counts_mat)
  
  seu <- CreateSeuratObject(
    counts = counts_mat,
    meta.data = md
  )
  
  DefaultAssay(seu) <- "RNA"
  
  raw <- .clean_label(md$label)
  
  truth <- sapply(raw, function(x) {
    x_low <- tolower(x)
    if (grepl("alpha", x_low)) return("Alpha cell")
    if (grepl("beta", x_low)) return("Beta cell")
    if (grepl("delta", x_low)) return("Delta cell")
    if (grepl("gamma|pp", x_low)) return("Gamma cell")
    if (grepl("acinar", x_low)) return("Acinar cell")
    if (grepl("duct", x_low)) return("Ductal cell")
    if (grepl("endothelial", x_low)) return("Endothelial cell")
    if (grepl("mesenchymal|stellate", x_low)) return("Stellate cell")
    "Unknown"
  })
  
  .process_seurat(seu, truth, "Pancreas", "Human", "MuraroPancreas")
}

# -----------------------------------------------------------------------------
# 4. Lawlor Pancreas
# -----------------------------------------------------------------------------
get_lawlor_pancreas_data <- function(seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  sce <- scRNAseq::LawlorPancreasData()
  sce <- sce[, !is.na(sce$celltype) & sce$celltype != ""]
  
  counts_mat <- counts(sce)
  metadata <- as.data.frame(colData(sce))
  
  if (is.null(colnames(counts_mat))) {
    colnames(counts_mat) <- paste0("cell_", seq_len(ncol(counts_mat)))
  }
  
  rownames(metadata) <- colnames(counts_mat)
  
  seu <- CreateSeuratObject(counts = counts_mat, meta.data = metadata)
  
  raw <- .clean_label(seu$celltype)
  
  truth <- sapply(raw, function(x) {
    x_low <- tolower(x)
    if (grepl("alpha", x_low)) return("Alpha cell")
    if (grepl("beta", x_low)) return("Beta cell")
    if (grepl("delta", x_low)) return("Delta cell")
    if (grepl("acinar", x_low)) return("Acinar cell")
    if (grepl("duct", x_low)) return("Ductal cell")
    if (grepl("endothelial", x_low)) return("Endothelial cell")
    if (grepl("stellate|mesenchymal", x_low)) return("Stellate cell")
    if (grepl("macrophage|immune", x_low)) return("Immune cell")
    "Unknown"
  })
  
  .process_seurat(seu, truth, "Pancreas", "Human", "LawlorPancreas")
}

# -----------------------------------------------------------------------------
# 5. Tasic Brain
# -----------------------------------------------------------------------------
get_tasic_brain_data <- function(seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  sce <- scRNAseq::TasicBrainData()
  seu <- CreateSeuratObject(counts = counts(sce), meta.data = as.data.frame(colData(sce)))
  
  raw <- .clean_label(seu$broad_type)
  
  truth <- sapply(raw, function(x) {
    x_low <- tolower(x)
    if (grepl("opc|precursor", x_low)) return("Oligodendrocyte precursor cell")
    if (grepl("glut|excit|pyr|neuron", x_low) && !grepl("gaba|inhib", x_low)) return("Excitatory neuron")
    if (grepl("gaba|inhib|interneuron", x_low)) return("Inhibitory neuron")
    if (grepl("astro", x_low)) return("Astrocyte")
    if (grepl("micro", x_low)) return("Microglia")
    if (grepl("oligo", x_low)) return("Oligodendrocyte")
    if (grepl("endo", x_low)) return("Endothelial cell")
    if (grepl("peri", x_low)) return("Pericyte")
    "Unknown"
  })
  
  .process_seurat(seu, truth, "Brain", "Mouse", "TasicBrain")
}

# -----------------------------------------------------------------------------
# 6. Zeisel Brain
# -----------------------------------------------------------------------------
get_zeisel_brain_data <- function(seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  sce <- scRNAseq::ZeiselBrainData()
  rownames(sce) <- make.names(rownames(sce))
  
  seu <- CreateSeuratObject(counts = assay(sce), meta.data = as.data.frame(colData(sce)))
  
  raw <- .clean_label(seu$level2class)
  
  truth <- sapply(raw, function(x) {
    x_low <- tolower(x)
    if (grepl("pyr|excit", x_low)) return("Excitatory neuron")
    if (grepl("int|inhib|gaba", x_low)) return("Inhibitory neuron")
    if (grepl("astro", x_low)) return("Astrocyte")
    if (grepl("^pvm|perivascular", x_low)) return("Macrophage")
    if (grepl("micro|mgl", x_low)) return("Microglia")
    if (grepl("opc", x_low)) return("Oligodendrocyte precursor cell")
    if (grepl("oligo", x_low)) return("Oligodendrocyte")
    if (grepl("endo|^vend", x_low)) return("Endothelial cell")
    if (grepl("peri", x_low)) return("Pericyte")
    if (grepl("^vsmc|smooth muscle", x_low)) return("Smooth muscle cell")
    if (grepl("epend", x_low)) return("Ependymal cell")
    "Unknown"
  })
  
  .process_seurat(seu, truth, "Brain", "Mouse", "ZeiselBrain")
}

# -----------------------------------------------------------------------------
# 7. Zilionis Lung Cancer
# -----------------------------------------------------------------------------
get_zilionis_lung_data <- function(seed = NULL, n_cells = 5000) {
  if (!is.null(seed)) set.seed(seed)
  
  sce <- scRNAseq::ZilionisLungData()
  md_all <- as.data.frame(colData(sce))
  
  label_candidates <- c(
    "Major cell type",
    "Major.cell.type",
    "major_cell_type",
    "cell_type",
    "celltype",
    "CellType",
    "label",
    "Most.likely.LM22.cell.type",
    "Most likely LM22 cell type",
    "cluster"
  )
  
  label_col <- intersect(label_candidates, colnames(md_all))[1]
  
  if (is.na(label_col)) {
    message("Available Zilionis metadata columns:")
    print(colnames(md_all))
    stop("No usable label column found for ZilionisLung.", call. = FALSE)
  }
  
  message("Zilionis label column used: ", label_col)
  annotated <- !is.na(md_all[[label_col]]) &
    nzchar(trimws(as.character(md_all[[label_col]])))

  if (sum(annotated) == 0) {
    stop("ZilionisLung label column has no annotated cells: ", label_col, call. = FALSE)
  }

  message(
    "Zilionis annotated cells retained: ",
    sum(annotated),
    " / ",
    ncol(sce)
  )

  sce <- sce[, annotated]

  if (!is.null(n_cells) && n_cells < ncol(sce)) {
    sce <- sce[, sample(colnames(sce), n_cells)]
  }

  counts_mat <- counts(sce)

  rownames(counts_mat) <- make.unique(
    gsub("_", "-", rownames(counts_mat))
  )

  colnames(counts_mat) <- make.unique(
    as.character(colnames(counts_mat))
  )

  md <- as.data.frame(colData(sce))
  rownames(md) <- colnames(counts_mat)

  stopifnot(!anyDuplicated(rownames(counts_mat)))
  stopifnot(!anyDuplicated(colnames(counts_mat)))
  stopifnot(identical(rownames(md), colnames(counts_mat)))

  seu <- CreateSeuratObject(
    counts = counts_mat,
    meta.data = md
  )
  
  raw <- .clean_label(md[[label_col]])
  
  truth <- sapply(raw, function(x) {
    x_low <- tolower(x)
    
    if (.is_unknown_label(x_low)) return("Unknown")
    if (grepl("patient[0-9]+.*specific|epithelial|alveolar|club|ciliated|tumor|malignant|cancer|type i cells|type ii cells", x_low)) return("Epithelial cell")
    if (grepl("t cell|cd4|cd8", x_low)) return("T cell")
    if (grepl("b cell", x_low)) return("B cell")
    if (grepl("plasma", x_low)) return("Plasma cell")
    if (grepl("platelet", x_low)) return("Platelet")
    if (.is_nk_label(x_low)) return("Natural killer cell")
    if (grepl("macrophage", x_low)) return("Macrophage")
    if (grepl("monocyte", x_low)) return("Monocyte")
    if (grepl("neutrophil", x_low)) return("Neutrophil")
    if (grepl("dendritic|dc", x_low)) return("Dendritic cell")
    if (grepl("mast", x_low)) return("Mast cell")
    if (grepl("fibroblast", x_low)) return("Fibroblast")
    if (grepl("endothelial", x_low)) return("Endothelial cell")
    if (grepl("smooth muscle|pericyte", x_low)) return("Smooth muscle cell")
    if (grepl("^[bt]?rbc$|erythro|red blood", x_low)) return("Erythrocyte")
    
    "Unknown"
  })
  
  message("Zilionis truth labels:")
  print(table(truth))
  
  .process_seurat(
    seu,
    truth,
    tissue = "Lung",
    species = "Human",
    dataset_name = "ZilionisLung"
  )
}

# -----------------------------------------------------------------------------
# 8. Segerstolpe Pancreas
# -----------------------------------------------------------------------------
get_segerstolpe_pancreas_data <- function(seed = NULL) {
  
  if (!is.null(seed))
    set.seed(seed)
  
  sce <- scRNAseq::SegerstolpePancreasData()
  
  md <- as.data.frame(colData(sce))
  
  label_col <- intersect(
    c("cell type", "celltype", "label"),
    colnames(md)
  )[1]
  
  if (is.na(label_col))
    stop("Could not find cell type column")
  
  sce <- sce[, !is.na(md[[label_col]]) &
               md[[label_col]] != ""]
  
  seu <- CreateSeuratObject(
    counts = counts(sce),
    meta.data = as.data.frame(colData(sce))
  )
  
  raw <- .clean_label(
    seu@meta.data[[label_col]]
  )
}
# -----------------------------------------------------------------------------
# 9. Romanov Brain
# -----------------------------------------------------------------------------
get_romanov_brain_data <- function(seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  sce <- scRNAseq::RomanovBrainData()
  seu <- CreateSeuratObject(counts = counts(sce), meta.data = as.data.frame(colData(sce)))
  
  label_col <- if ("cell.type1" %in% colnames(seu@meta.data)) "cell.type1" else colnames(seu@meta.data)[1]
  raw <- .clean_label(seu@meta.data[[label_col]])
  
  truth <- sapply(raw, function(x) {
    x_low <- tolower(x)
    if (grepl("excit|glut|pyr", x_low)) return("Excitatory neuron")
    if (grepl("inhib|gaba|interneuron", x_low)) return("Inhibitory neuron")
    if (grepl("astro", x_low)) return("Astrocyte")
    if (grepl("micro", x_low)) return("Microglia")
    if (grepl("opc|precursor", x_low)) return("Oligodendrocyte precursor cell")
    if (grepl("oligo", x_low)) return("Oligodendrocyte")
    if (grepl("endo", x_low)) return("Endothelial cell")
    if (grepl("epend", x_low)) return("Ependymal cell")
    "Unknown"
  })
  
  .process_seurat(seu, truth, "Brain", "Mouse", "RomanovBrain")
}

# -----------------------------------------------------------------------------
# 10. Tabula Muris Droplet
# -----------------------------------------------------------------------------
get_tabula_muris_droplet_data <- function(seed = NULL, tissue_filter = "Lung", n_cells = 5000) {
  if (!is.null(seed)) set.seed(seed)
  
  sce <- scRNAseq::TabulaMurisDropletData()
  cd <- as.data.frame(colData(sce))
  
  tissue_col <- intersect(c("tissue", "tissue_type", "organ"), colnames(cd))[1]
  label_col <- intersect(c("cell_ontology_class", "cell.type", "cell_type", "label"), colnames(cd))[1]
  
  if (is.na(tissue_col) || is.na(label_col)) {
    stop("Could not find tissue or label columns in TabulaMurisDropletData.", call. = FALSE)
  }
  
  keep <- cd[[tissue_col]] == tissue_filter & !is.na(cd[[label_col]])
  sce <- sce[, keep]
  
  if (!is.null(n_cells) && n_cells < ncol(sce)) {
    sce <- sce[, sample(colnames(sce), n_cells)]
  }
  
  seu <- CreateSeuratObject(counts = counts(sce), meta.data = as.data.frame(colData(sce)))
  raw <- .clean_label(seu@meta.data[[label_col]])
  
  truth <- sapply(raw, function(x) {
    x_low <- tolower(x)
    if (grepl("t cell", x_low)) return("T cell")
    if (grepl("b cell", x_low)) return("B cell")
    if (grepl("macrophage", x_low)) return("Macrophage")
    if (grepl("monocyte", x_low)) return("Monocyte")
    if (grepl("neutrophil", x_low)) return("Neutrophil")
    if (grepl("dendritic", x_low)) return("Dendritic cell")
    if (grepl("epithelial|alveolar|club|ciliated", x_low)) return("Epithelial cell")
    if (grepl("endothelial", x_low)) return("Endothelial cell")
    if (grepl("fibroblast", x_low)) return("Fibroblast")
    if (grepl("smooth muscle", x_low)) return("Smooth muscle cell")
    "Unknown"
  })
  
  .process_seurat(seu, truth, tissue_filter, "Mouse", paste0("TabulaMuris", tissue_filter))
}

# -----------------------------------------------------------------------------
# Dataset registry
# -----------------------------------------------------------------------------

cache_dataset <- function(name, seed, loader) {
  dir.create("benchmark_cache", showWarnings = FALSE)

  cache_key <- paste(
    name,
    paste0("seed", seed),
    BENCHMARK_CACHE_VERSION,
    paste0("markers", TOP_MARKERS),
    paste0("res", SEURAT_RESOLUTION),
    paste0("pcs", N_PCS),
    sep = "_"
  )
  cache_key <- gsub("[^A-Za-z0-9_.-]", "_", cache_key)
  
  cache_file <- file.path(
    "benchmark_cache",
    paste0(cache_key, ".rds")
  )
  
  if (file.exists(cache_file)) {
    message("Loading cached dataset: ", cache_file)
    return(readRDS(cache_file))
  }
  
  message("Processing dataset: ", name)
  obj <- loader()
  saveRDS(obj, cache_file)
  obj
}


load_benchmark_datasets <- function(seed = 100) {
  list(
    PBMC = cache_dataset("PBMC", seed, function() get_pbmc_data(seed)),
    
    BaronPancreas = cache_dataset("BaronPancreas", seed, function() get_baron_pancreas_data(seed)),
    
    MuraroPancreas = cache_dataset("MuraroPancreas", seed, function() get_muraro_pancreas_data(seed)),
    
    TasicBrain = cache_dataset("TasicBrain", seed, function() get_tasic_brain_data(seed)),
    
    ZeiselBrain = cache_dataset("ZeiselBrain", seed, function() get_zeisel_brain_data(seed)),
    
    ZilionisLung = cache_dataset("ZilionisLung", seed, function() get_zilionis_lung_data(seed))
  )
}
