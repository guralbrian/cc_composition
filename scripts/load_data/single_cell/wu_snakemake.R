# Load necessary libraries
libs <- c("Seurat", "zellkonverter", "curl", "SeuratDisk",
          "TabulaMurisSenisData", "tidyverse", "Matrix", "optparse")

lapply(libs, library, character.only = T, lib.loc = "tools/r_libs/common")

# Set up command-line options
option_list <- list(
  make_option(c("-o", "--outpath"), type="character", default=NULL,
              help="Output path for the Seurat object", metavar="character")
)

# Parse command-line options
parser <- OptionParser(option_list=option_list)
opt <- parse_args(parser)

# Wu 2021 ####
# raw fasta available at:
#     https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&page_size=10&acc=SRR15248449&display=metadata

# Import data from Wu 2021 doi.org/10.1152/physiolgenomics.00016.2021
geo <- "GSM5471468"
filePaths = getGEOSuppFiles("GSM5471468")

# Unzip files
paths <- rownames(filePaths)
for(i in paths){
  gunzip(i)
}

# Read and format as Seurat
wu.mt <- readMM("GSM5471468/GSM5471468_matrix.mtx")
cell.ids <- readLines("GSM5471468/GSM5471468_barcodes.tsv")
genes <- read.table("GSM5471468/GSM5471468_genes.tsv")
colnames(wu.mt) <- paste0(cell.ids, "-wu")
rownames(wu.mt) <- genes$V2
sn.wu <- CreateSeuratObject(wu.mt)
sn.wu$origin <- "wu"
sn.wu$orig.ident <- "wu"
sn.wu$PercentMito <- PercentageFeatureSet(sn.wu, pattern = "^mt-")

# Save the Seurat object
if (!is.null(opt$outpath)) {
  SaveH5Seurat(sn.wu, opt$outpath)
} else {
  stop("No output path specified.")
}
