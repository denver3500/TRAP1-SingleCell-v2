# Schwann Cell Subclustering Analysis
# Subset Schwann cells and perform detailed clustering to identify subpopulations

library(Seurat)
library(ggplot2)
library(dplyr)

# Load the annotated Seurat object
data.merged <- readRDS("raw_data/annotation_renamed.rds")

# Check current cell type annotations
print("Current cell type counts:")
print(table(Idents(data.merged)))

# Subset for Schwann cells only
schwann_cells <- subset(data.merged, idents = "Schwann cell")

print(paste("Total Schwann cells extracted:", ncol(schwann_cells)))
print(paste("Original clusters represented:", paste(unique(schwann_cells$seurat_clusters), collapse = ", ")))

# Normalize and scale the data for subclustering
schwann_cells <- NormalizeData(schwann_cells)
schwann_cells <- FindVariableFeatures(schwann_cells, selection.method = "vst", nfeatures = 2000)
schwann_cells <- ScaleData(schwann_cells)

# Run PCA
schwann_cells <- RunPCA(schwann_cells, features = VariableFeatures(object = schwann_cells))

# Determine the number of PCs to use
ElbowPlot(schwann_cells, ndims = 50)
ggsave("pictures/subclustering/schwann_elbow_plot.pdf", width = 10, height = 6)

# Find neighbors and clusters with different resolutions
schwann_cells <- FindNeighbors(schwann_cells, dims = 1:30)

# Try multiple resolutions to find optimal subclustering
resolutions <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0)

for (res in resolutions) {
  schwann_cells <- FindClusters(schwann_cells, resolution = res)
  print(paste("Resolution", res, "produces", length(unique(Idents(schwann_cells))), "subclusters"))
}

# Set a reasonable resolution (adjust based on the output above)
schwann_cells <- FindClusters(schwann_cells, resolution = 0.3)

# Run UMAP
schwann_cells <- RunUMAP(schwann_cells, dims = 1:30)

# Create directory for subclustering results
dir.create("pictures/subclustering", recursive = TRUE, showWarnings = FALSE)
dir.create("statistics/subclustering", recursive = TRUE, showWarnings = FALSE)

# Visualize subclusters
p1 <- DimPlot(schwann_cells, reduction = "umap", label = TRUE, label.size = 6) +
      ggtitle("Schwann Cell Subclusters") +
      theme(plot.title = element_text(size = 16, face = "bold"))

p2 <- DimPlot(schwann_cells, reduction = "umap", group.by = "orig.ident") +
      ggtitle("Schwann Cell Subclusters by Sample") +
      theme(plot.title = element_text(size = 16, face = "bold"))

p3 <- DimPlot(schwann_cells, reduction = "umap", group.by = "seurat_clusters") +
      ggtitle("Original Clusters within Schwann Cells") +
      theme(plot.title = element_text(size = 16, face = "bold"))

# Save plots
ggsave("pictures/subclustering/schwann_subclusters_umap.pdf", p1, width = 10, height = 8)
ggsave("pictures/subclustering/schwann_subclusters_by_sample.pdf", p2, width = 12, height = 8)
ggsave("pictures/subclustering/schwann_original_clusters.pdf", p3, width = 10, height = 8)

# Combined plot
combined_plot <- p1 + p2 + p3
ggsave("pictures/subclustering/schwann_subclustering_combined.pdf", combined_plot, width = 18, height = 6)

# Find markers for each subcluster
schwann_markers <- FindAllMarkers(schwann_cells, 
                                  only.pos = TRUE,
                                  min.pct = 0.25,
                                  logfc.threshold = 0.25)

# Get top 10 markers per subcluster
top10_schwann_markers <- schwann_markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) %>%
    arrange(cluster, desc(avg_log2FC))

# Save markers
write.csv(schwann_markers, "statistics/subclustering/schwann_all_markers.csv", row.names = FALSE)
write.csv(top10_schwann_markers, "statistics/subclustering/schwann_top10_markers.csv", row.names = FALSE)

# Create a heatmap of top markers
top5_markers <- schwann_markers %>%
    group_by(cluster) %>%
    top_n(n = 5, wt = avg_log2FC)

heatmap_plot <- DoHeatmap(schwann_cells, features = top5_markers$gene) + 
                NoLegend() +
                theme(axis.text.y = element_text(size = 8))

ggsave("pictures/subclustering/schwann_markers_heatmap.pdf", heatmap_plot, width = 12, height = 10)

# Print summary
print("=== SCHWANN CELL SUBCLUSTERING SUMMARY ===")
print(paste("Number of Schwann cell subclusters:", length(unique(Idents(schwann_cells)))))
print("Subcluster sizes:")
print(table(Idents(schwann_cells)))

print("Distribution across original clusters:")
print(table(schwann_cells$seurat_clusters, Idents(schwann_cells)))

print("Distribution across samples:")
print(table(schwann_cells$orig.ident, Idents(schwann_cells)))

# Save the subclustered object
schwann_cells$schwann_subclusters <- Idents(schwann_cells)
saveRDS(schwann_cells, "raw_data/schwann_subclustered.rds")

print("Schwann cell subclustering analysis complete!")
print("Files saved:")
print("- raw_data/schwann_subclustered.rds (subclustered object)")
print("- pictures/subclustering/ (UMAP plots)")
print("- statistics/subclustering/ (marker genes)")


# Analyze COL6A genes in Schwann cell subclusters
col6_genes <- c("COL6A1", "COL6A2", "COL6A3")

# Check which COL6A genes are present in the dataset
available_col6_genes <- col6_genes[col6_genes %in% rownames(schwann_cells)]
print(paste("Available COL6A genes:", paste(available_col6_genes, collapse = ", ")))

if (length(available_col6_genes) > 0) {
  
  # Create feature plots for COL6A genes
  col6_feature_plots <- list()
  
  for (i in 1:length(available_col6_genes)) {
    gene <- available_col6_genes[i]
    
    p_feature <- FeaturePlot(schwann_cells,
                            features = gene,
                            reduction = "umap",
                            min.cutoff = "q5",
                            max.cutoff = "q95",
                            order = TRUE) +
                ggtitle(paste(gene, "Expression in Schwann Subclusters")) +
                theme(plot.title = element_text(size = 12, face = "bold"))
    
    col6_feature_plots[[i]] <- p_feature
  }
  
  # Create combined feature plot
  if (length(col6_feature_plots) == 3) {
    combined_col6_features <- col6_feature_plots[[1]] + col6_feature_plots[[2]] + col6_feature_plots[[3]]
  } else if (length(col6_feature_plots) == 2) {
    combined_col6_features <- col6_feature_plots[[1]] + col6_feature_plots[[2]]
  } else {
    combined_col6_features <- col6_feature_plots[[1]]
  }
  
  ggsave("pictures/subclustering/schwann_COL6A_feature_plots.pdf", 
         combined_col6_features, width = 15, height = 5)
  
  # Create dotplot for COL6A genes
  col6_dotplot <- DotPlot(schwann_cells, 
                         features = available_col6_genes,
                         cols = c("lightgrey", "red"),
                         dot.scale = 8) +
                 RotatedAxis() +
                 ggtitle("COL6A Expression Across Schwann Cell Subclusters") +
                 theme(plot.title = element_text(size = 14, face = "bold"),
                       axis.text.x = element_text(angle = 45, hjust = 1),
                       axis.text.y = element_text(size = 10))
  
  ggsave("pictures/subclustering/schwann_COL6A_dotplot.pdf", 
         col6_dotplot, width = 8, height = 6)
  
  # Create violin plots for COL6A genes
  col6_violin <- VlnPlot(schwann_cells, 
                        features = available_col6_genes,
                        ncol = length(available_col6_genes),
                        pt.size = 0.1) +
                plot_annotation(title = "COL6A Expression Distribution in Schwann Subclusters")
  
  ggsave("pictures/subclustering/schwann_COL6A_violin_plots.pdf", 
         col6_violin, width = 12, height = 6)
  
  # Calculate average expression and percentage expressing for each subcluster
  print("=== COL6A EXPRESSION ANALYSIS ===")
  
  for (gene in available_col6_genes) {
    print(paste("\n", gene, "expression by subcluster:"))
    
    # Get expression data
    expr_data <- GetAssayData(schwann_cells, slot = "data")[gene,]
    clusters <- Idents(schwann_cells)
    
    # Calculate stats per cluster
    for (cluster in sort(unique(clusters))) {
      cluster_cells <- names(clusters)[clusters == cluster]
      cluster_expr <- expr_data[cluster_cells]
      
      pct_expressing <- sum(cluster_expr > 0) / length(cluster_expr) * 100
      avg_expr <- mean(cluster_expr[cluster_expr > 0])
      
      if (is.nan(avg_expr)) avg_expr <- 0
      
      print(paste("  Cluster", cluster, ": ", 
                 round(pct_expressing, 1), "% expressing, avg expr:", 
                 round(avg_expr, 3)))
    }
  }
  
  # Create summary table
  col6_summary <- data.frame()
  
  for (gene in available_col6_genes) {
    expr_data <- GetAssayData(schwann_cells, slot = "data")[gene,]
    clusters <- Idents(schwann_cells)
    
    for (cluster in sort(unique(clusters))) {
      cluster_cells <- names(clusters)[clusters == cluster]
      cluster_expr <- expr_data[cluster_cells]
      
      pct_expressing <- sum(cluster_expr > 0) / length(cluster_expr) * 100
      avg_expr <- mean(cluster_expr[cluster_expr > 0])
      if (is.nan(avg_expr)) avg_expr <- 0
      
      col6_summary <- rbind(col6_summary, data.frame(
        Gene = gene,
        Cluster = cluster,
        Percent_Expressing = round(pct_expressing, 2),
        Average_Expression = round(avg_expr, 3),
        Cell_Count = length(cluster_cells)
      ))
    }
  }
  
  # Save summary table
  write.csv(col6_summary, "statistics/subclustering/schwann_COL6A_expression_summary.csv", row.names = FALSE)
  
  print("\nCOL6A analysis files saved:")
  print("- pictures/subclustering/schwann_COL6A_feature_plots.pdf")
  print("- pictures/subclustering/schwann_COL6A_dotplot.pdf") 
  print("- pictures/subclustering/schwann_COL6A_violin_plots.pdf")
  print("- statistics/subclustering/schwann_COL6A_expression_summary.csv")
  
} else {
  print("No COL6A genes found in the dataset")
}

# Rename Schwann cell subclusters based on functional annotations
schwann_cells <- RenameIdents(schwann_cells,
    "0" = "Repair/Regenerative Schwann Cells",
    "1" = "Invasive/Migratory Tumor Schwann Cells", 
    "2" = "Dedifferentiated Tumor Schwann Cells",
    "3" = "Myelinating Schwann Cells",
    "4" = "Immune-Activated Schwann Cells",
    "5" = "Cancer Stem-like Schwann Cells",
    "6" = "Quiescent/Resting Schwann Cells",
    "7" = "Non-Myelinating Schwann Cells",
    "8" = "Inflammatory Schwann Cells"
)

# Save renamed identities
schwann_cells$schwann_subtypes <- Idents(schwann_cells)

# Create UMAP with renamed identities
p_renamed <- DimPlot(schwann_cells, reduction = "umap", 
                     label = TRUE, label.size = 3, repel = TRUE) +
             ggtitle("Schwann Cell Functional Subtypes") +
             theme(plot.title = element_text(size = 16, face = "bold"),
                   legend.text = element_text(size = 8))

ggsave("pictures/subclustering/schwann_functional_subtypes_umap.pdf", 
       p_renamed, width = 14, height = 10)

# Save the annotated Schwann cell object
saveRDS(schwann_cells, "raw_data/schwann_functional_subtypes.rds")
print("Schwann cell functional subtypes saved to raw_data/schwann_functional_subtypes.rds")
