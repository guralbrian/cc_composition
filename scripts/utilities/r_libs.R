# This script is meant to be run independently of the actual analysis pipeline. It just serves to load libraries.


# List libraries
libs <- c('vctrs','Seurat', 'tidyverse', 'BiocManager', 'remotes', 'curl',
          'Matrix', 'optparse')

# Install libraries through cran to the mambaforge cc_r_env R library

# Function to load libs
.installPackages <- function(x){
    utils::install.packages(x,
                 lib = "tools/r_libs/common",
                 repos = "http://cran.us.r-project.org")
}

# Apply load libs function
lapply(libs, .installPackages)

# Install packages that use BiocManager

libs <- c('zellkonverter', 'TabulaMurisSenisData', 'DESeq2')

library("BiocManager", lib.loc = "tools/r_libs/common")

# Function to load libs via Bioconductor
.installPackagesBioc <- function(x){
  BiocManager::install(x,
                       lib = "tools/r_libs/common",
                       force = T)
}


# Apply load libs function
lapply(libs, .installPackagesBioc)


# Github installs 
library("remotes", lib.loc = "tools/r_libs/common")

remotes::install_github("mojaveazure/seurat-disk",
                        lib = "tools/r_libs/common",
                        force = T)
