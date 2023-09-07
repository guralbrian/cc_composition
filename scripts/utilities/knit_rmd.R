library(rmarkdown)

# Get command line arguements
args <- commandArgs(trailingOnly = TRUE)

# Pull out the path
rmd_file_path <- args[1]

# Knit the file in the path
rmarkdown::render(rmd_file_path)
