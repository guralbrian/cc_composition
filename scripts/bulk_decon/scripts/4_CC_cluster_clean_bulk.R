# PCA on cc samples to find outliers
libs <- c("tidyverse", "ggrepel")
lapply(libs, require, character.only = T)
rm(libs)
# Load CC data 
cc.bulk <- read.csv("data/processed/bulk/cc_counts_summed_transcripts.csv", row.names = 1)
phenotypes_real <- read.csv("data/processed/bulk/cc_pheno.csv")

# Pull sample IDs from colnames
colnames(cc.bulk) <- paste0("S",  readr::parse_number(colnames(cc.bulk)))

# PCA 

cc_pca <- prcomp(cc.bulk)$rotation |> as.data.frame()
cc_pca$Adaptor <- row.names(cc_pca)
cc_pca <- cc_pca |>
  select(PC1, PC2, Adaptor) |>
  left_join(phenotypes_real)

#plot first two PCs
cc_pca |>
  as.data.frame() |>
  ggplot(aes(x = PC1, y = PC2, color = Treatment, label = Mouse)) +
  geom_point() +
  geom_text_repel(size = 4, box.padding = unit(0.35, "lines"), 
                  force = 10, segment.linetype = 2, segment.size = 0.6, color = "black", force_pull = 0.01, min.segment.length = 0)
  
# Annotate the two groups 
cc_pca <- cc_pca |> 
  mutate(cluster = factor(case_when(
    PC2 < -0.1 ~ "1",
    PC2 > -0.1 & PC1 > -0.3  ~ "2",
    PC1 < -0.3 ~ "outlier"
  )))

# Replot
my_palette <- c("#A6CEE3", "#FDBF6F", "#bd5b69")

cc_pca |>
  as.data.frame() |>
  ggplot(aes(x = PC1, y = PC2, color = cluster, label = Mouse)) +
  geom_point() +
  geom_text_repel(size = 4, box.padding = unit(0.35, "lines"), 
                  force = 10, segment.linetype = 2, segment.size = 0.6, color = "black", force_pull = 0.01, min.segment.length = 0) +
  theme(axis.text.x = element_text(color = "black", size = 10, angle =25, vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 10),
        axis.ticks = element_blank(),
        legend.position = c(0.2, 0.5),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        strip.text = element_text(size = 10),
        title = element_text(size = 10),
        legend.text = element_text(size = 10),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major = element_line(color = "grey"),
        panel.grid.minor = element_blank(),
        #legend.background = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='transparent'),
        axis.title.x = element_text(size = 10),
        plot.margin = unit(c(1,1,1,1), units = "cm"),
        axis.title.y = element_text(size = 10)) +
  labs(color = "Cluster") +
  scale_color_manual(values = my_palette)

# Save Clusters as a phenotype


write.csv(cc_pca, "data/processed/bulk/cc_pca_clusters.csv", row.names = F)


  