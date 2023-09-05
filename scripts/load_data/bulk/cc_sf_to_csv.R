# Convert CC bulk RNAseq Salmon outputs to a single counts matrix

library(tidyverse)

# List all .sf files
sf_files <- list.files("data/raw/cc_counts/06222022/quant_files/", pattern = "*.sf", full.names = TRUE)

# Read all .sf files and bind them into a single data frame
all_data <- sf_files |> 
  set_names(nm = basename(sf_files)) |> 
  map_dfr(~ read_tsv(.x) %>% select(Name, NumReads) %>% rename_with(~ .x, .cols = "NumReads"), .id = "Sample")

# Pivot to create the counts matrix
counts_matrix <- all_data |> 
  pivot_wider(names_from = Sample, values_from = NumReads) |> 
  column_to_rownames(var = "Name") |> 
  as.data.frame()

write.csv(counts_matrix, "data/raw/cc_counts/06222022/raw_matrix.csv")

