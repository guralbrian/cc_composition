# Load necessary packages
library(dplyr)
library(readr)

# Base directory where model output folders are located
base_dir <- "results/benchmarking/deseq/comp_adj/parallel_models/09182023/model_outputs/"

# List all model directories
model_dirs <- list.files(base_dir, pattern = "^[^\\.].*")

# Initialize an empty list to hold the final dataframes for each model
final_dfs <- list()

# Loop through each model directory to merge individual files
for (model_dir in model_dirs) {
  # Create the path to the model directory
  full_model_dir <- paste0(base_dir, model_dir, "/")
  
  # Get list of csv files matching the pattern
  csv_files <- list.files(full_model_dir, pattern = "_\\d+\\.csv$")
  
  # Read and merge these files
  merged_df <- bind_rows(lapply(csv_files, function(f) {
    read_csv(paste0(full_model_dir, f))
  }))
  
  # Add a new column to indicate the model
  merged_df$model_name <- model_dir
  
  # Append this dataframe to the list of final dataframes
  final_dfs[[model_dir]] <- merged_df
}

# Combine all final dataframes into a single dataframe
final_output_df <- bind_rows(final_dfs)

# Save the final merged dataframe
write.csv(final_output_df, file.path(base_dir, "merged.csv"))

