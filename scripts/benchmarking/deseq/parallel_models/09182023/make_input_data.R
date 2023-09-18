
# List libraries
libs <- c("Seurat", "ggplot2", "patchwork", "SeuratDisk","reshape2", "tidyverse",
          "gridExtra",  "readr", "Biobase", "ggmagnify", "stringr", 
          "ggforce", "viridis", "DESeq2", "rlist", "purrr") # list libraries here
# Require all of them
lapply(libs, require, character.only = T)

rm(libs)

sn <- LoadH5Seurat("data/processed/single_cell/merges/rau_patterson/09132023/cell_types_db_1_5.h5seurat")


# Downsample 
sn.small <- sn |>
  GetAssayData(assay = "RNA", slot = "counts") |>
  as.matrix() |>
  SampleUMI(4000, verbose = T, upsample = T)

sn.small <- CreateSeuratObject(sn.small)
meta.features <- colnames(sn@meta.data)
for(i in meta.features){
  sn.small <- AddMetaData(sn.small, sn@meta.data[[i]], col.name = i)
}

Idents(sn.small) <- Idents(sn)

library(Seurat)
library(stringi)

# Replace spaces, commas, and periods in Idents levels with underscores
new_levels <- gsub("[,\\. ]", "_", levels(Idents(sn.small)))

# Update Idents with the new levels
Idents(sn.small) <- factor(Idents(sn.small), levels = levels(Idents(sn.small)), labels = new_levels)


# Function to make a data frame of cell type ratios
simulate_ratios <- function(major_cell, major_prop, 
                            cell_types, range, step_size, replicates, noise) {
  
  # get group tiers for major cell type
  # these are the percentages of the major cell type for which groups will be centered
  major_props_groups <- seq(major_prop - range, major_prop + range, step_size)
  
  # get values for major cell type
  # each value in the major_props_groups vector is used as a mean with n replicates pulled from a normal distribution with noise
  major_props <- lapply(major_props_groups, function(x){rnorm(replicates, x, noise)}) |>
    unlist()
  
  # Set major cell type and proportion
  minor <- cell_types[!(cell_types %in% c(major_cell))]
  
  
  # Define the range for proportions
  minor_props <- rep((1-major_props)/length(minor), length(minor)) |>
    matrix(ncol = length(minor)) |>
    as.data.frame()
  colnames(minor_props) <- minor
  
  # set up props dataframe
  combos <- c(major_props) |>
    matrix(ncol = 1) |>
    as.data.frame()
  colnames(combos) <- c(major_cell)
  
  # add minor cell types to props data frame
  props <- cbind(combos, minor_props)
  props$pct.change <- as.factor(rep(major_props_groups - major_prop, each = replicates))
  props$pct.change <- relevel(props$pct.change, ref = "0")
  rownames(props) <- paste0("mix_", 1:nrow(props))
  props[,2] <- (props[,1] + rowMeans(props[,3:length(cell_types)])) / 4
  props[,3] <- (props[,1] + rowMeans(props[,4:length(cell_types)])) / 5
  props[,2:length(cell_types)] <- sapply(2:length(cell_types), function(x){props[,x] * sample(rnorm(1000, 1, noise), length(major_props))})
  
  # Normalize back to proportional sum == 1
  props[,1:length(cell_types)]  <- props[,1:length(cell_types)]/rowSums(props[,1:length(cell_types)])
  
  return(props)
}


# Function to sample cell IDs from a cell type, then aggregate their expression
AddCellProfile <- function(sn, count, cell_type){
  
  # Get the cell.ids that are the cell type of interest
  cell.ids  <- Idents(sn)[which(Idents(sn) == cell_type)] |>
    names() |>
    sample(size = count[cell_type], replace = T)
  
  # sum expression of cells with those ids
  cell.type.profile <- rowSums(sn@assays$RNA@counts[,cell.ids]) |>
    as.numeric()  
  assign(as.character(cell_type), cell.type.profile) # name the profile with the cell type
  return(get(as.character(cell_type)))
}

AggCells <- function(sn, ratios, cell_count) {
  ratio.format <- ratios[1,] |> unlist()
  names(ratio.format) <- names(ratios)
  # Ensure ratios is a named list or named vector
  if (!is.list(ratio.format) && !is.vector(ratio.format)) {
    stop("ratios must be a named list or named vector")
  }
  # Cell counts by type
  cell_counts <- round(ratio.format * cell_count)
  #names(cell_counts) <- seq(0,length(ratio.format)-1, 1)
  dat <- sapply(1:length(cell_counts), function(x)
  {
    AddCellProfile(sn, cell_counts[x], names(cell_counts)[x])
  }
  ) |>
    rowSums()
  names(dat) <- rownames(sn)
  return(dat)
}


major.cell <- "Cardiomyocytes"
major.prop <- 0.7
cell.types <- unfactor(unique(Idents(sn.small)))
range <- 0.2
step.size <- 0.01
replicates <- 5
noise <- 0.01

ratios <- simulate_ratios(major_cell = major.cell, 
                          major_prop = major.prop, 
                          cell_types = cell.types, 
                          range = range, 
                          step_size = step.size, 
                          replicates = replicates,
                          noise = noise)


start <- Sys.time()
# Function to make a list of cell type ratios
cell.target <- 5000
# take each row of cell.ratios and put it into MakePseudoBulk 
pb.add <- sapply(1:nrow(ratios), function(x){AggCells(sn.small, ratios[x,-(length(cell.types)+1)], cell.target)}) |>
  as.data.frame()
colnames(pb.add) <- rownames(ratios)

end <- Sys.time() - start


# Make clr comp
library(compositions)
clr.sample <- clr(ratios[,as.character(cell.types)])
colnames(clr.sample) <- paste0("clr.", colnames(clr.sample))

# Make ilr 
ilr.sample <- ilr(ratios[,as.character(cell.types)])
colnames(ilr.sample) <- paste0("ilr.", seq(1, ncol(ilr.sample), 1))


# Make PCA comp
pca.sample <- prcomp(t(ratios[,as.character(cell.types)]))$rotation |>
  as.data.frame()

pca.sample$pct.change <- ratios$pct.change

# add this info back to proportion dfs 
ratios <- cbind(ratios, clr.sample) |> 
  cbind(ilr.sample)  |> 
  cbind(pca.sample) 

# make df with clr, pca, and raw compositions
counts <- pb.add
total.genes <- 2000
counts.small <- counts[sample(rownames(counts), total.genes),]


# Randomly select 10% of genes to be differentially expressed
n_genes <- nrow(counts.small)
n_diff_exp_genes <- round(0.2 * n_genes)
diff_exp_genes <- sample(rownames(counts.small), n_diff_exp_genes)
stable_genes <- rownames(counts.small)[!(rownames(counts.small) %in% diff_exp_genes)]
# Define the fold-change or other modifier
fold_change <- 2  # Change this based on your specific needs
reference.group <- ratios |> 
  subset(pct.change == 0) %>%
  rownames_to_column() |> 
  pull(rowname)
# Loop over the columns (samples) in counts.small to apply the fold-change
# Skip the reference group column(s) as needed
for (col in colnames(counts.small)) {
  if (col %in% reference.group) {
    next
  }
  counts.small[diff_exp_genes, col] <- counts.small[diff_exp_genes, col] * fold_change
}

# Save necessary data
counts.dict <- data.frame(genes = row.names(counts.small),
                          diff_exp = row.names(counts.small) %in% diff_exp_genes)

write.csv(counts.dict, "results/benchmarking/deseq/comp_adj/parallel_models/09182023/input_data/counts_dict_2000g.csv", row.names = F)
write.csv(counts.small, "results/benchmarking/deseq/comp_adj/parallel_models/09182023/input_data/counts_2000g.csv")
write.csv(ratios, "results/benchmarking/deseq/comp_adj/parallel_models/09182023/input_data/ratios_2000g.csv")
