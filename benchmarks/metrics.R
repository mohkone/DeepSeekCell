# benchmarks/metrics.R
# Evaluation metrics, ontology mapping, clade accuracy

suppressPackageStartupMessages({
  library(mclust)
  library(ontologyIndex)
})

# -----------------------------------------------------------------------------
# Cell Ontology loading (robust)
# -----------------------------------------------------------------------------

load_cell_ontology <- function(local_path = "data/cl.obo") {
  if (!file.exists(local_path)) {
    stop("Cell Ontology file not found: ", local_path)
  }
  message("Loading Cell Ontology from: ", local_path)
  ont <- try(ontologyIndex::get_ontology(local_path, extract_tags = "everything"), silent = TRUE)
  if (inherits(ont, "try-error")) stop("Failed to parse cl.obo")
  message("Ontology loaded: ", length(ont$id), " terms")
  ancestors <- lapply(ont$id, function(x) get_ancestors(ont, x))
  names(ancestors) <- ont$id
  list(ont = ont, ancestors = ancestors)
}

# -----------------------------------------------------------------------------
# Controlled vocabulary
# -----------------------------------------------------------------------------

prediction_map <- list(
  "Naive T cell" = "Naive T cell", "T cell" = "Naive T cell",
  "Monocyte" = "Classical monocyte", "Macrophage" = "Non-classical monocyte",
  "Cytotoxic T cell" = "Cytotoxic T cell", "Effector T cell" = "Cytotoxic T cell",
  "Natural Killer cell" = "Natural killer cell", "NK cell" = "Natural killer cell",
  "Dendritic cell" = "Dendritic cell", "Platelet" = "Platelet", "B cell" = "B cell",
  "Classical monocyte" = "Classical monocyte", "Non-classical monocyte" = "Non-classical monocyte",
  "Intermediate Monocyte" = "Classical monocyte",
  "Beta cell" = "Beta cell", "Alpha cell" = "Alpha cell", "Delta cell" = "Delta cell",
  "Acinar cell" = "Acinar cell", "Ductal cell" = "Ductal cell",
  "Endothelial cell" = "Endothelial cell", "Stellate cell" = "Stellate cell",
  "Enteroendocrine cell" = "Enteroendocrine cell", "Pericyte" = "Pericyte",
  "Immune cell" = "Macrophage",
  "Excitatory neuron" = "Excitatory neuron", "Neurons" = "Excitatory neuron",
  "Inhibitory neuron" = "Inhibitory neuron", "GABAergic neurons" = "Inhibitory neuron",
  "Oligodendrocyte" = "Oligodendrocyte", "Oligodendrocytes" = "Oligodendrocyte",
  "Astrocyte" = "Astrocyte", "Astrocytes" = "Astrocyte",
  "Microglia" = "Microglia", "Endothelial cells" = "Endothelial cell",
  "Oligodendrocyte precursor cells" = "Oligodendrocyte precursor cell",
  "OPC" = "Oligodendrocyte precursor cell",
  "Smooth muscle cell" = "Smooth muscle cell",
  "Vascular smooth muscle cells" = "Smooth muscle cell",
  "Ependymal cells" = "Ependymal cell", "Fibroblast" = "Fibroblast"
)

standardise_prediction <- function(cell_type) {
  if (is.na(cell_type) || cell_type == "Unknown" || cell_type == "unknown") return("Unknown")
  for (pattern in names(prediction_map)) {
    if (tolower(cell_type) == tolower(pattern)) return(prediction_map[[pattern]])
  }
  for (pattern in names(prediction_map)) {
    if (grepl(tolower(pattern), tolower(cell_type), fixed = TRUE))
      return(prediction_map[[pattern]])
  }
  "Unknown"
}

# -----------------------------------------------------------------------------
# Ontology mapping
# -----------------------------------------------------------------------------

unmatched_log <- character()

standard_to_cl <- function(cell_type, ont_data) {
  if (cell_type == "Unknown") return(NA_character_)
  norm <- tolower(cell_type)
  norm <- gsub("[^a-z]", "", norm)
  norm <- sub("s$", "", norm)
  
  cl_lookup <- c(
    "naivetcell" = "CL:0000895", "cytotoxictcell" = "CL:0000909",
    "bcell" = "CL:0000236", "classicalmonocyte" = "CL:0000870",
    "nonclassicalmonocyte" = "CL:0000871", "naturalkillercell" = "CL:0000623",
    "dendriticcell" = "CL:0000451", "platelet" = "CL:0000233",
    "macrophage" = "CL:0000235",
    "betacell" = "CL:0000169", "alphacell" = "CL:0000502",
    "deltacell" = "CL:0000503", "acinarcell" = "CL:0000168",
    "ductalcell" = "CL:0000167", "endothelialcell" = "CL:0000115",
    "stellatecell" = "CL:0000632", "enteroendocrinecell" = "CL:0000428",
    "pericyte" = "CL:0000232",
    "excitatoryneuron" = "CL:0000540", "inhibitoryneuron" = "CL:0000617",
    "oligodendrocyte" = "CL:0000128", "astrocyte" = "CL:0000127",
    "microglia" = "CL:0000129", "smoothmusclecell" = "CL:0000192",
    "ependymalcell" = "CL:0000065", "oligodendrocyteprecursorcell" = "CL:0002570",
    "fibroblast" = "CL:0000057"
  )
  
  if (norm %in% names(cl_lookup)) return(cl_lookup[[norm]])
  
  for (id in ont_data$ont$id) {
    name_norm <- tolower(ont_data$ont$name[[id]])
    name_norm <- gsub("[^a-z]", "", name_norm)
    name_norm <- sub("s$", "", name_norm)
    if (!is.na(name_norm) && name_norm == norm) return(id)
  }
  
  unmatched_log <<- c(unmatched_log, paste0(cell_type, " [", norm, "]"))
  return(NA_character_)
}

clade_accuracy <- function(pred_std, true_std, ont_data) {
  pred_cl <- sapply(pred_std, function(x) standard_to_cl(x, ont_data))
  true_cl <- sapply(true_std, function(x) standard_to_cl(x, ont_data))
  valid <- !is.na(pred_cl) & !is.na(true_cl)
  if (sum(valid) == 0) return(NA)
  pred_cl <- pred_cl[valid]; true_cl <- true_cl[valid]
  correct <- 0
  for (i in seq_along(pred_cl)) {
    if (pred_cl[i] == true_cl[i]) {
      correct <- correct + 1
    } else {
      anc_true <- ont_data$ancestors[[true_cl[i]]]
      if (pred_cl[i] %in% anc_true) correct <- correct + 1
    }
  }
  correct / length(pred_cl)
}

# -----------------------------------------------------------------------------
# Main metrics
# -----------------------------------------------------------------------------

evaluate_metrics <- function(pred, true, dataset_name, ont_data) {
  pred_std <- sapply(pred, standardise_prediction)
  true_std <- true
  valid <- pred_std != "Unknown" & true_std != "Unknown"
  pred_std <- pred_std[valid]; true_std <- true_std[valid]
  
  if (length(pred_std) < 2) {
    warning("Too few valid labels in ", dataset_name)
    return(c(ARI=NA, MacroF1=NA, Accuracy=NA, BalancedAcc=NA, CladeAcc=NA))
  }
  
  ari <- tryCatch(adjustedRandIndex(pred_std, true_std), error=function(e) NA)
  acc <- mean(pred_std == true_std)
  
  classes <- unique(true_std)
  f1s <- sapply(classes, function(cls) {
    tp <- sum(pred_std == cls & true_std == cls)
    fp <- sum(pred_std == cls & true_std != cls)
    fn <- sum(pred_std != cls & true_std == cls)
    prec <- if (tp+fp > 0) tp/(tp+fp) else 0
    rec  <- if (tp+fn > 0) tp/(tp+fn) else 0
    if (prec+rec > 0) 2*prec*rec/(prec+rec) else 0
  })
  macro_f1 <- mean(f1s)
  
  recalls <- sapply(classes, function(cls) {
    tp <- sum(pred_std == cls & true_std == cls)
    fn <- sum(pred_std != cls & true_std == cls)
    if (tp+fn > 0) tp/(tp+fn) else 0
  })
  bal_acc <- mean(recalls)
  
  clade_acc <- clade_accuracy(pred_std, true_std, ont_data)
  
  c(ARI=ari, MacroF1=macro_f1, Accuracy=acc, BalancedAcc=bal_acc, CladeAcc=clade_acc)
}