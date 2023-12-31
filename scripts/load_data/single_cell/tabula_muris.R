# Load necessary libraries
libs <- c("Seurat", "zellkonverter", "curl", "SeuratDisk",
          "TabulaMurisSenisData", "tidyverse", "Matrix", "optparse", "DESeq2")

lapply(libs, require, character.only = T)


# Tabula Muris ####
# raw fasta data available at
# https://www.ncbi.nlm.nih.gov/sra?term=SRX3607046
sn_muris <- TabulaMurisSenisDroplet(
  tissues = "Heart_and_Aorta",
  processedCounts = FALSE,
  reducedDims = TRUE,
  infoOnly = FALSE)[[1]] 

sn_muris_seurat <- sn_muris |>
  SummarizedExperiment::assay("counts")|>
  as.matrix() |>
  Seurat::CreateSeuratObject()

# Add metadata from SingleCellExperiment to Seurat
for(i in colnames(colData(sn_muris))){
  sn_muris_seurat <- Seurat::AddMetaData(object = sn_muris_seurat, col.name = i, metadata = colData(sn_muris)[[i]])
}

# Tabula Muris didn't map to mitochondrial genes
# but we need the data slot for later

sn_muris_seurat$PercentMito <- 0
sn_muris_seurat$origin <- "tabula_muris"
sn_muris_seurat$orig.ident <- paste0("tm_", sn_muris_seurat$mouse.id)
sn_muris_seurat <- subset(sn_muris_seurat, tissue_free_annotation == "Heart")

# Save the Seurat object
SaveH5Seurat(sn_muris_seurat, "data/processed/single_cell/tabula_muris")

