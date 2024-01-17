libs <- c("tidyverse",  "reshape2", "readr") # list libraries here
lapply(libs, require, character.only = T)
rm(libs)


#Load cell type fraction bulk RNAseq dataset
fractions <- read.csv("data/raw/rau_fractions/celltype_counts.csv", row.names = 1)

# Load CC data 
cc.bulk <- read.csv("data/processed/bulk/cc_counts_summed_transcripts.csv", row.names = 1)

# Pull sample IDs from colnames
colnames(cc.bulk) <- paste0("S",  readr::parse_number(colnames(cc.bulk)))


common.genes <- row.names(fractions)[row.names(fractions) %in% row.names(cc.bulk)]

bulk.all <- cbind(cc.bulk[common.genes,], fractions[common.genes,])


# Save dataset

write.csv(bulk.all, "data/processed/bulk/cc_and_fractions.csv", row.names = T)
