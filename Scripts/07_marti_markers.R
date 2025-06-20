# Schwann Cell Markers UMAP Visualization - Annotated Version
# Creates dotplot and wide UMAP layout using cell type annotations

library(Seurat)
library(ggplot2)
library(patchwork)

# Load the annotated Seurat object
data.merged <- readRDS("raw_data/annotation_renamed.rds")

# Use annotated cell types
print(paste("Total cells in dataset:", ncol(data.merged)))
print(paste("Number of cell types:", length(unique(Idents(data.merged)))))
print("Cell types in dataset:")
print(table(Idents(data.merged)))

# Create the cluster UMAP plot using cell type annotations
p_cluster <- DimPlot(data.merged, reduction = "umap",
                     label = TRUE,
                     repel = TRUE) +
             ggtitle("Annotated Cell Types") +
             theme(plot.title = element_text(size = 14, face = "bold"),
                   legend.position = "right")

# Define the folder where the plots will be saved
picture_folder <- "pictures"
if (!dir.exists(picture_folder)) {
  dir.create(picture_folder, recursive = TRUE)
}

# List of Schwann cell marker genes to visualize
gene_list <- c("MS4A1", "CD8A", "CD4", "ELANE", "CD163", "CD44", "MKI67", "ACTA2")


# Create feature plots for all genes
feature_plots <- list()

for (i in 1:length(gene_list)) {
  gene <- gene_list[i]
  
  # Create feature plot for the current gene
  p_feature <- FeaturePlot(data.merged,
                           features = gene,
                           reduction = "umap",
                           min.cutoff = "q5",
                           max.cutoff = "q95",
                           order = TRUE) +
              ggtitle(paste(gene, "Expression")) +
              theme(plot.title = element_text(size = 12, face = "bold"),
                    legend.position = "right",
                    legend.key.size = unit(0.3, "cm"))
  
  # Store the plot in the list
  feature_plots[[i]] <- p_feature
}

# Create dotplot for Schwann cell markers across annotated cell types
print("Generating dotplot for Schwann cell markers across annotated cell types...")
dot_plot <- DotPlot(data.merged, 
                   features = gene_list,
                   cols = c("lightgrey", "red"),
                   dot.scale = 8,
                   scale = TRUE) +
          RotatedAxis() +
          ggtitle("Markers Across Annotated Cell Types") +
          theme(plot.title = element_text(size = 14, face = "bold"),
                axis.text.x = element_text(angle = 45, hjust = 1),
                axis.text.y = element_text(size = 10)) +
          guides(size = guide_legend(title = "% Expressing",
                                   override.aes = list(color = "black")),
                 color = guide_colorbar(title = "Avg Expression"))

# Save the dotplot
ggsave(paste0(picture_folder, "/schwann_markers/marti_markers_annotated_dotplot.pdf"), 
       dot_plot, width = 12, height = 8)

# Wide version with wider cluster map
combined_plot_wide <- p_cluster + 
                     feature_plots[[1]] + feature_plots[[2]] + feature_plots[[3]] +
                     feature_plots[[4]] + feature_plots[[5]] + feature_plots[[6]] +
                     plot_layout(ncol = 4, nrow = 2,
                               widths = c(2, 1, 1, 1),
                               design = "
                               1234
                               1567")

ggsave(paste0(picture_folder, "/schwann_markers/marti_markers_annotated_wide_UMAP.pdf"), 
       combined_plot_wide, width = 24, height = 12)

# Print summary of expression across cell types
print("Summary of Schwann marker expression across annotated cell types:")
for (gene in gene_list) {
  if (gene %in% rownames(data.merged)) {
    expr_data <- GetAssayData(data.merged, slot = "data")[gene,]
    percent_expressing <- sum(expr_data > 0) / length(expr_data) * 100
    print(paste(gene, ": ", round(percent_expressing, 2), "% of cells expressing"))
  } else {
    print(paste(gene, ": Gene not found in dataset"))
  }
}

message("Saved Schwann cell markers visualizations for annotated cell types:")
message("- marti_markers_annotated_wide_UMAP.pdf (wide UMAP layout)")
message("- marti_markers_annotated_dotplot.pdf (dotplot)")

print("Schwann cell markers annotated visualization complete!")