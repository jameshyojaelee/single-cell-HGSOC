#!/usr/bin/env Rscript

library(Seurat)
library(ggplot2)
# Read the RDS file
shah_data <- readRDS("/gpfs/commons/home/jameslee/HGSOC/data/Shah_T.rds")

# set working directory to the output folder
setwd("/gpfs/commons/home/jameslee/HGSOC/output")

# Export count matrix
count_matrix <- GetAssayData(shah_data, slot = "counts")
write.csv(as.matrix(count_matrix), "count_matrix.csv")

# Export metadata
write.csv(shah_data@meta.data, "metadata.csv")

# Export UMAP coordinates
umap_coords <- Embeddings(shah_data, reduction = "umap")
write.csv(umap_coords, "umap_coordinates.csv")

# View metadata
colnames(shah_data@meta.data)
#head(shah_data@meta.data)

# Visualize UMAP
pdf("umap_plot.pdf", width = 12, height = 8)
DimPlot(shah_data, reduction = "umap", label = TRUE, raster = FALSE)
dev.off()

pdf("umap_plot_unlabelled.pdf", width = 8, height = 7)
DimPlot(shah_data, reduction = "umap", label = FALSE, raster = FALSE) + 
  theme(legend.position = "none",
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black"))
dev.off()

# Visualize feature expression on UMAP
pdf("feature_plot.pdf")
FeaturePlot(shah_data, features = rownames(shah_data)[1:4], reduction = "umap", raster = FALSE)
dev.off()