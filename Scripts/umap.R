# Generate annotated UMAP visualization of single-cell data

# Load required packages
library(Seurat)
library(ggplot2)

# Create output directory
dir.create("pictures/umap", recursive = TRUE, showWarnings = FALSE)

# Load the annotated Seurat object
data.merged <- readRDS("raw_data/annotation_renamed.rds")

# Check the object
print(paste("Number of cells:", ncol(data.merged)))
print(paste("Number of cell types:", length(unique(Idents(data.merged)))))
print("Cell types in dataset:")
print(table(Idents(data.merged)))

# Generate standard UMAP with annotations
p1 <- DimPlot(data.merged, 
             reduction = "umap", 
             label = TRUE, 
             label.size = 8, 
             repel = TRUE) + 
      ggtitle("Cell Type Annotations") +
      theme(plot.title = element_text(size = 16, face = "bold"))

# Save as PDF and PNG
ggsave("pictures/umap/annotated_umap.pdf", p1, width = 12, height = 10)
ggsave("pictures/umap/annotated_umap.png", p1, width = 12, height = 10, dpi = 300)

# Generate UMAP split by sample (if sample information exists)
if ("orig.ident" %in% colnames(data.merged@meta.data)) {
  p2 <- DimPlot(data.merged, 
               reduction = "umap", 
               group.by = "orig.ident", 
               label = FALSE) + 
        ggtitle("Samples") +
        theme(plot.title = element_text(size = 16, face = "bold"))
  
  ggsave("pictures/umap/samples_umap.pdf", p2, width = 12, height = 10)
  ggsave("pictures/umap/samples_umap.png", p2, width = 12, height = 10, dpi = 300)
}

print("UMAP plots generated and saved in pictures/umap/")