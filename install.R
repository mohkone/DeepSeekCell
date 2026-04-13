# Save the script as app.R
# Install dependencies:
install.packages(c("shiny", "ggplot2", "dplyr", "httr2", "jsonlite", 
                   "openxlsx", "ontologyIndex", "purrr", "future", 
                   "future.apply", "logger", "cachem", "DT", "shinythemes",
                   "shinycssloaders", "plotly", "stringdist"))

install.packages(c("Seurat", "SingleR", "scType", "celldex", "cluster", "mclust",
                   "ontologyIndex", "dplyr", "ggplot2", "reshape2", "httr2", "jsonlite"))

install.packages(c("SeuratData", "scRNAseq", "SingleCellExperiment"))
BiocManager::install(c("TabulaSapiens", "scMCA", "HCLR"))  # if needed