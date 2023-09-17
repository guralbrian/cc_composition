
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
  props[,2:length(cell.types)] <- sapply(2:length(cell.types), function(x){props[,x] * sample(rnorm(1000, 1, noise), length(major_props))})

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
major.prop <- 0.5
cell.types <- unfactor(unique(Idents(sn.small)))
range <- 0.3
step.size <- 0.01
replicates <- 5
noise <- 0.035

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
total.genes <- 10000
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



# Change formatting so that DESeq doesn't complain
sample_info <- ratios[colnames(counts),]
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
deseq.results <- lapply(models, function(x){TestModels(x)})  

save(deseq.results, file = "results/benchmarking/deseq/comp_adj/09152023_10000g.Rdata")



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
result$sub_sublist_name <- result$sub_sublist_name |>cod 
  str_replace("_","-") |>
  as.numeric()

# Define the model types by adjustment used
result <- result |> 
  mutate(ModelType = factor(case_when(
    str_detect(sublist_name, "unadjusted") ~ "No adjustment",
    str_detect(sublist_name, "raw") ~ "Raw proportions",
    str_detect(sublist_name, "clr") ~ "Centered log ratios",
    str_detect(sublist_name, "ilr") ~ "Isometric log ratios",
    str_detect(sublist_name, "pc") ~ "PCA")))


# Want to take the ratio of the false postive rate over the true positive rate

# This is an example, replace with the actual number
total_diff_exp_genes <- length(diff_exp_genes)

rates.df <- result %>%
  mutate(
    false_neg = total_diff_exp_genes - true_pos,
    true_neg = total.genes - (true_pos + false_pos + false_neg),
    Precision = true_pos / (true_pos + false_pos),
    Recall = true_pos / (true_pos + false_neg),
    "F1 score" = 2 * (Precision * Recall) / (Precision + Recall),
    Specificity = true_neg / (true_neg + false_pos)
  )
# Melt the data for plotting
rates.df <- rates.df %>%
  pivot_longer(cols = c(Precision, Recall, "F1 score", Specificity),
               names_to = "Metric",
               values_to = "Value") %>%
  group_by(Metric, sublist_name) %>%
  mutate(median_value = median(Value, na.rm = TRUE) ) %>% # 
  arrange(Metric, median_value)




pal <- c("#969696",
         "#fcc5c0", "#fa9fb5", "#f768a1", "#c51b8a",
         "#c6dbef", "#9ecae1", "#6baed6", "#3182bd",
         "#fd8d3c","#e6550d")

# Modify legend text 

models.legend <- list(
  unadjusted    = ~ "Unadjusted",
  raw_cm        = ~ "CM prop",
  raw_cm_fb     = ~ "CM prop + 1x",
  raw_cm_fb_im  = ~ "CM prop + 2x",
  raw_4         = ~ "CM prop + 3x",
  clr_cm        = ~ "CLR of CMs",
  clr_cm_fb     = ~ "CLR of CMs + 1x",
  clr_cm_fb_im  = ~ "CLR of CMs + 2x",
  clr_4         = ~ "CLR of CMs + 3x",
  pc1           = ~ "PC1",
  pc2           = ~ "PC1 + PC2"
)

p.fp <- result |> 
  subset(ModelType != "Isometric log ratios") |>
  ggplot(aes(x = sub_sublist_name, y = false_pos +0.01, color = sublist_name)) +
  geom_point(alpha = 0.4) + # smaller and more transparent points
  geom_smooth(se = F, method = "loess", size = 1.5, alpha = 0.2, span = 0.8) +
  #scale_y_continuous(trans = "log2") +
  scale_color_manual(values = pal,
                     name = "Model design",
                     labels = models.legend,
                     breaks = names(models.legend)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    text = element_text(size = 24),
    panel.grid.major = element_blank(), # remove major gridlines
    panel.grid.minor = element_blank()  # remove minor gridlines
  ) +
  #geom_magnify(from = c(-0.07, 0.07,-3,30), to = c(-0.08,0.08, 380,600), 
  #                 shape = "rect", shadow = F)  +
  labs(y= "False Positives in 10,000 genes", x = "Simulated CM proportion difference") +
  ggtitle("False Positive Rate")

p.tp <- result |> 
  subset(ModelType != "Isometric log ratios") |>
  ggplot(aes(x = sub_sublist_name, y = true_pos + 0.01, color = sublist_name, group = sublist_name)) +
  geom_point(alpha = 0.2) + # smaller and more transparent points
  geom_smooth(se = F, method = "loess", size = 1.5, alpha = 0.2, span = 0.7) +
  #scale_y_continuous(trans = "log10") +
  scale_color_manual(values = pal,
                     name = "Model design",
                     labels = models.legend,
                     breaks = names(models.legend)) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    text = element_text(size = 24),
    panel.grid.major = element_blank(), # remove major gridlines
    panel.grid.minor = element_blank(),
    legend.key.width = unit(3,"cm"),
    legend.key.height =unit(1,"cm"), # remove minor gridlines
  ) +
  #geom_magnify(from = c(-0.07, 0.07,-3,30), to = c(-0.08,0.08, 380,600), 
  #                 shape = "rect", shadow = F)  +
  labs(y= "True Positives out of 2000 DE genes", x = "Simulated CM proportion difference") +
  ggtitle("True Positive Rate") 


p.stats <- rates.df |> 
  subset(ModelType != "Isometric log ratios" &
           Metric == "F1 score") |> 
  ggplot(aes(x = factor(sublist_name, level = names(models.legend)), 
             y = Value)) +
  geom_boxplot(aes(fill = sublist_name)) +
  geom_jitter(width = 0.1, height = 0, alpha = 0.3) +
  scale_fill_manual(values = pal,
                    name = "Model design",
                    labels = models.legend,
                    breaks = names(models.legend)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    text = element_text(size = 24)
  ) +
  labs(
    title = "Model Accuracy",
    y = "F1 Score (higher is more accurate)",
    x = "Model Type"
  ) 
png(filename="results/benchmarking/deseq/comp_adj/09152023_10000g.png",
    width = 1920,
    height = 1080,
    units = "px")

design <- c("ABC")
wrap_plots(A = p.fp, B = p.tp, C = p.stats, design = design)

dev.off()

