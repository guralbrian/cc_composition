# Load necessary libraries
libs <- c("Seurat", "zellkonverter", "curl", "SeuratDisk",
          "TabulaMurisSenisData", "tidyverse", "Matrix", "optparse")

lapply(libs, require, character.only = T)


# Martini 2019 ####

# Import data from Martini 2019 doi.org/10.1161/CIRCULATIONAHA.119.041694
# https://www.ncbi.nlm.nih.gov/sra/SRX5063185[accn]
geo <- "GSE122930"
dir <- "./data/processed/single_cell/geo"
filePaths = getGEOSuppFiles(geo, baseDir = dir)
paths <- rownames(filePaths)

# Unzip files
for(i in paths){
  if(endsWith(i, '.gz')){
    gunzip(i)
  }
}

# select files ending in .mtx or .tsv
files <- list.files(paste0(dir, "/", geo))
files.mtx <- files[endsWith(files, '.mtx')]
files.tsv <- files[endsWith(files, '.tsv')]

# group by lapply(names, strsplit( "_"), "[[", 2)

treatment <- lapply(strsplit(files.mtx, "_"), "[[", 2)
timepoint <- lapply(strsplit(files.mtx, "_"), "[[", 3)
samples.martini <- paste0(treatment,"_", timepoint)

# read in each file and merge

for(i in samples.martini){
  # List the .mtx file and read
  mtx.file <- files.mtx[grepl(i, files.mtx)]
  mtx <- readMM(paste0(dir, "/",geo, "/", mtx.file))
  
  # Get list the barcodes .tsv file
  tsv.bar <- files.tsv[grepl(i, files.tsv) & grepl("barcodes", files.tsv)]
  tsv.bar <- readLines(paste0(dir, "/",geo, "/", tsv.bar))
  
  # Get list the genes .tsv file
  tsv.genes <- files.tsv[grepl(i, files.tsv) & grepl("genes", files.tsv)]
  tsv.genes <- read.table(paste0(dir, "/",geo, "/", tsv.genes))
  
  # Format and return sparse matrix
  colnames(mtx) <- paste0(tsv.bar, "-mar")
  rownames(mtx) <- tsv.genes$V2
  assign(i, mtx)
}

# Convert to Seurats
for(i in 1:length(samples.martini)){
  seurat <- CreateSeuratObject(get(samples.martini[i]))
  seurat$martini.cond <- samples.martini[i]
  assign(samples.martini[i], seurat)
}

# Merge Seurats and save
sn.martini <- merge(get(samples.martini[1]), get(samples.martini[2])) |>
  merge(get(samples.martini[3])) |>
  merge(get(samples.martini[4]))

sn.martini$origin <- "martini"
sn.martini$orig.ident <- paste0("martini_", sn.martini$martini.cond)
sn.martini$PercentMito <- PercentageFeatureSet(sn.martini, pattern = "^mt-")

# Save the Seurat object
SaveH5Seurat(sn.martini, "data/processed/single_cell/martini")
