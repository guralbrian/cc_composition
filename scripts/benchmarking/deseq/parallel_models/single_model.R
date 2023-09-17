# Script to be paired with shell script to run many versions in parallel

# Load Libraries

# List libraries
libs <- c("ggplot2", "patchwork", "reshape2", "tidyverse",
          "gridExtra",  "readr", "Biobase", "ggmagnify", "stringr", 
          "ggforce", "viridis", "DESeq2", "rlist", "purrr") # list libraries here
# Require all of them
lapply(libs, require, character.only = T)

rm(libs)

# Capture command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Convert the first argument to integer (model index)
model_index <- as.integer(args[1])

# Convert the second argument to integer (chunk index)
chunk_index <- as.integer(args[2])

# Chunk size
chunk_size <- 500

# Calculate start and end row for the chunk
start_row <- (chunk_index - 1) * chunk_size + 1
end_row <- chunk_index * chunk_size

# Read and subset the dictionary and counts matrix
counts.dict <- read.csv("results/benchmarking/deseq/comp_adj/parallel_models/counts_dict_10000g.csv")[start_row:end_row,]
counts.small <- read.csv("results/benchmarking/deseq/comp_adj/parallel_models/counts_10000g.csv", row.names = 1)[start_row:end_row,]
ratios <- read.csv("results/benchmarking/deseq/comp_adj/parallel_models/ratios_10000g.csv", row.names = 1)


# Make list of models
models <- list(
  unadjusted    = ~ 0 + pct.change,
  raw_cm        = ~ 0 + pct.change + Cardiomyocytes,
  raw_cm_fb     = ~ 0 + pct.change + Cardiomyocytes + Fibroblast__Mesen_,
  raw_cm_fb_im  = ~ 0 + pct.change + Cardiomyocytes + Fibroblast__Mesen_ + Endothelial__Cor__Art_,
  raw_4         = ~ 0 + pct.change + Cardiomyocytes + Fibroblast__Mesen_ + Endothelial__Cor__Art_ + unclear,
  clr_cm        = ~ 0 + pct.change + clr.Cardiomyocytes,
  clr_cm_fb     = ~ 0 + pct.change + clr.Cardiomyocytes + clr.Fibroblast__Mesen_,
  clr_cm_fb_im  = ~ 0 + pct.change + clr.Cardiomyocytes + clr.Fibroblast__Mesen_ + clr.Endothelial__Cor__Art_,
  clr_4         = ~ 0 + pct.change + clr.Cardiomyocytes + clr.Fibroblast__Mesen_ + clr.Endothelial__Cor__Art_ + clr.unclear,
  pc1           = ~ 0 + pct.change + PC1,
  pc2           = ~ 0 + pct.change + PC1 + PC2
)

model.use <- models[model_index]
# Change formatting so that DESeq doesn't complain
sample_info <- ratios[colnames(counts.small),]
colnames(sample_info) <- colnames(sample_info) |>
  str_replace_all("`", "") |>
  str_replace_all(" ", "_")

sample_info$pct.change <- sample_info$pct.change |>
  str_replace_all("-", "_") |>
  as.factor()

# Relevel the sample info to have the no change group as the reference
sample_info$pct.change <- relevel(sample_info$pct.change, ref = "0")

# Function to run DESeq2
TestModels <- function(model){
  # Create a DESeqDataSet
  dds <- DESeqDataSetFromMatrix(
    countData = counts.small,
    colData = sample_info,
    design = model
  )
  
  # Run the DESeq 
  dds <- DESeq(dds)
  
  # Get the differential expression analysis results, contrasted against the 0 percent change group
  resultsNames(dds)
  genes <- lapply(pct.use, function(x)
  {results(dds, contrast = c("pct.change", "0", x)) |> as.data.frame()})
  return(genes)
}

# Get the list of group percents to use
pct.use <- levels(sample_info$pct.change)[-1]

# Run the analysis
deseq.results <- lapply(model.use, function(x){TestModels(x)})  

deseq.res <- deseq.results
# Iterate over the names of de.lfcse
for (sublist_name in names(deseq.res)) {
  # Get the current sublist
  sublist <- deseq.res[[sublist_name]]
  
  # Rename the elements in the sublist using pct.use
  names(sublist) <- pct.use
  
  # Update the sublist in de.lfcse
  deseq.res[[sublist_name]] <- sublist
}


# List stable and unstable genes
stable_genes <- counts.dict |> 
                subset(diff_exp == F) |> 
                pull(genes)

diff_exp_genes <- counts.dict |> 
  subset(diff_exp == T) |> 
  pull(genes)

# Function to collect values from nested list DESeq output
calculate_stats <- function(data, column_name, threshold) {
  result <- data %>%
    enframe(name = "sublist_name", value = "data") %>%
    mutate(data = map(data, ~.x %>% enframe(name = "sub_sublist_name", value = "values"))) %>%
    unnest(data) %>%
    mutate(
      mean_value = map_dbl(values, ~mean(.x[[column_name]], na.rm = TRUE)),
      std_dev = map_dbl(values, ~sd(.x[[column_name]], na.rm = TRUE)),
      false_pos = map_dbl(values, ~sum(abs(.x[[column_name]]) < threshold & rownames(.x) %in% stable_genes, na.rm = TRUE)),
      true_pos = map_dbl(values, ~sum(.x[[column_name]] < threshold & rownames(.x) %in% diff_exp_genes, na.rm = TRUE))
    ) %>%
    select(-values)
  
  return(result)
}

result <- calculate_stats(deseq.res, "padj", 0.05)

# Replace the formatting change in the percent change group labels
result$sub_sublist_name <- result$sub_sublist_name |>
  str_replace("_","-") |>
  as.numeric()

# record the actual number of DE genes
result$total_pos <- length(counts.dict$diff_exp[which(counts.dict$diff_exp == T)])

# Define the model types by adjustment used
result <- result |> 
  mutate(ModelType = factor(case_when(
    str_detect(sublist_name, "unadjusted") ~ "No adjustment",
    str_detect(sublist_name, "raw") ~ "Raw proportions",
    str_detect(sublist_name, "clr") ~ "Centered log ratios",
    str_detect(sublist_name, "ilr") ~ "Isometric log ratios",
    str_detect(sublist_name, "pc") ~ "PCA")))


# Save the output in a model-specific sub directory
mainDir <- "results/benchmarking/deseq/comp_adj/parallel_models"
subDir <- names(model.use)[1]

# Make a directory if it doesnt already exist
if (!file.exists(file.path(mainDir, subDir))){
  dir.create(file.path(mainDir, subDir))
}

# Save the results
write.csv(result, 
          file.path(mainDir, subDir, paste0(names(model.use)[1],"_",chunk_index, ".csv")), 
          row.names = F)