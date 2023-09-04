# Load necessary libraries
libs <- c("Seurat", "zellkonverter", "curl", "SeuratDisk",
          "TabulaMurisSenisData", "GEOquery", "tidyverse", "Matrix", "optparse")

lapply(libs, require, character.only = T, lib.loc = "tools/r_libs/common")


# Wu 2021 ####
# raw fasta available at:
#     https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&page_size=10&acc=SRR15248449&display=metadata

# Import data from Wu 2021 doi.org/10.1152/physiolgenomics.00016.2021
geo <- "GSM5471468"
geo_path <- "data/processed/single_cell/GSM/"
filePaths = getGEOSuppFiles(geo, baseDir = geo_path)

# Unzip files
paths <- rownames(filePaths)
for(i in paths){
  gunzip(i)
}

# Read and format as Seurat
wu.mt <- readMM(paste0(geo_path,"/GSM5471468/GSM5471468_matrix.mtx"))
cell.ids <- readLines(paste0(geo_path,"/GSM5471468/GSM5471468_barcodes.tsv"))
genes <- read.table(paste0(geo_path,"/GSM5471468/GSM5471468_genes.tsv"))
colnames(wu.mt) <- paste0(cell.ids, "-wu")
rownames(wu.mt) <- genes$V2
sn.wu <- CreateSeuratObject(wu.mt)
sn.wu$origin <- "wu"
sn.wu$orig.ident <- "wu"
sn.wu$PercentMito <- PercentageFeatureSet(sn.wu, pattern = "^mt-")

# Save the Seurat object
SaveH5Seurat(sn.wu, "data/processed/single_cell/wu")

