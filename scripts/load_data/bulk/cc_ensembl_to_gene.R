library(biomaRt)
library(dplyr)
library(tidyr)


# Load Data


counts_matrix <- read.csv("data/raw/cc_counts/06222022/raw_matrix.csv", row.names = 1)

# Initialize the biomaRt object for the organism of interest (e.g., human)
mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")

# Extract the Ensembl IDs from the counts_matrix
ensembl_ids <- rownames(counts_matrix) |> 
                as.vector()

# Fetch the corresponding gene symbols
gene_data <- getBM(attributes = c('ensembl_transcript_id_version', 'external_gene_name'), 
                   filters = 'ensembl_transcript_id_version', 
                   values = ensembl_ids, 
                   mart = mart)

# Merge the gene symbols with the counts_matrix
counts_matrix_with_gene_names <- counts_matrix |> 
  as.data.frame() |> 
  rownames_to_column("ensembl_transcript_id_version") |> 
  left_join(gene_data, by = "ensembl_transcript_id_version") |> 
  drop_na(external_gene_name) |>
  dplyr::select(-ensembl_transcript_id_version) |>
  group_by(external_gene_name) |> 
  summarise(across(colnames(counts_matrix), sum, na.rm = TRUE)) |>
  mutate(external_gene_name_compressed = make.unique(external_gene_name)) |>
  column_to_rownames(var = "external_gene_name_compressed") |>
  dplyr::select(-external_gene_name)

write.csv(counts_matrix_with_gene_names, "data/processed/bulk/cc/06222022/counts_summed_transcripts.csv")
