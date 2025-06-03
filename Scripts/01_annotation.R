library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(stringr)

# Get all folders in raw_data that end with "filtered_feature_bc_matrix"
filtered_dirs <- list.dirs(path = "raw_data", full.names = TRUE, recursive = FALSE)
filtered_dirs <- filtered_dirs[grepl("filtered_feature_bc_matrix$", filtered_dirs)]

# Create a list to store all Seurat objects
seurat_objects <- list()

# Loop through each filtered directory
for (i in seq_along(filtered_dirs)) {
    # Get the directory name
    dir_name <- basename(filtered_dirs[i])

    # Extract sample ID (assuming it's the prefix before "-filtered_feature_bc_matrix")
    sample_id <- gsub("-filtered_feature_bc_matrix$", "", dir_name)

    # Read the 10X data
    data <- Read10X(data.dir = filtered_dirs[i])

    # Create Seurat object
    seurat_obj <- CreateSeuratObject(
        counts = data,
        project = sample_id,
        min.cells = 3,
        min.features = 200
    )

    # Add sample ID as metadata
    seurat_obj$sample <- sample_id

    # Store in list
    seurat_objects[[sample_id]] <- seurat_obj
}

# Merge all Seurat objects
data.merged <- merge(seurat_objects[[1]], y = seurat_objects[2:length(seurat_objects)])


# Calculate mitochondrial percentage
data.merged[["percent.mt"]] <- PercentageFeatureSet(data.merged, pattern = "^MT")
mt_plot_by_sample <- VlnPlot(data.merged, features = "percent.mt", group.by = "sample")
ggsave("pictures/QC/qc_violin_plots.png", mt_plot_by_sample, width = 15, height = 5)

# Generate and save violin plot
vln_plot <- VlnPlot(data.merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("pictures/QC/qc_violin_plots.pdf", vln_plot, width = 15, height = 5)

# Generate scatter plots
plot1 <- FeatureScatter(data.merged, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data.merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Combine and save scatter plots
combined_plots <- CombinePlots(plots = list(plot1, plot2))
ggsave("pictures/QC/qc_scatter_plots.pdf", combined_plots, width = 12, height = 6)

data.merged <- subset(data.merged, subset = nFeature_RNA > 200 & nFeature_RNA < 9500 & percent.mt < 20)
data.merged <- NormalizeData(data.merged)
data.merged <- FindVariableFeatures(data.merged)
data.merged <- ScaleData(data.merged)
data.merged <- RunPCA(data.merged)
data.merged <- FindNeighbors(data.merged, dims = 1:30, reduction = "pca")
data.merged <- FindClusters(data.merged, resolution = 2, cluster.name = "unintegrated_clusters")
data.merged <- RunUMAP(data.merged, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

# Generate and display the first UMAP plot (by orig.ident and clusters)
p1 <- DimPlot(data.merged, reduction = "umap.unintegrated", 
              group.by = c("orig.ident", "seurat_clusters"))
ggsave("pictures/QC/umap_by_sample_and_clusters.pdf", p1, width = 12, height = 8)

# Generate and display the second UMAP plot (split by orig.ident)
p2 <- DimPlot(data.merged, reduction = "umap.unintegrated", 
              group.by = "seurat_clusters", split.by = "orig.ident")
ggsave("pictures/QC/umap_split_by_sample.pdf", p2, width = 12, height = 8)

pdf("pictures/QC/umap_plot.pdf", width = 12, height = 8)
DimPlot(data.merged, reduction = "umap.unintegrated", 
        group.by = "seurat_clusters", split.by = "orig.ident")
dev.off()

print("UMAP plots saved to pictures folder")


data.merged <- IntegrateLayers(
    object = data.merged, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
    verbose = FALSE
)
data.merged[["RNA"]] <- JoinLayers(data.merged[["RNA"]])
data.merged <- FindNeighbors(data.merged, reduction = "integrated.cca", dims = 1:30)
data.merged <- FindClusters(data.merged, resolution = 1)
data.merged <- RunUMAP(data.merged, dims = 1:30, reduction = "integrated.cca")


p3 <- DimPlot(data.merged, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"))
ggsave("pictures/QC/integrated_umap_by_sample_and_clusters.pdf", p3, width = 12, height = 8)

p4 <- DimPlot(data.merged, reduction = "umap", split.by = "orig.ident")
ggsave("pictures/QC/integrated_umap_split_by_sample.pdf", p4, width = 12, height = 8)

print("Integrated UMAP plots saved to pictures folder")
data.merged <- saveRDS(data.merged, file = "raw_data/annotation.rds")
print("Seurat object saved as annotation.rds in raw_data folder")
data.merged <- readRDS("raw_data/annotation.rds")

all.markers <- FindAllMarkers(data.merged, 
                              only.pos = TRUE,
                              min.pct = 0.25,
                              logfc.threshold = 0.25)

# Get top 30 genes for each cluster
top30_markers <- all.markers %>%
    group_by(cluster) %>%
    top_n(n = 30, wt = avg_log2FC) %>%
    arrange(cluster, desc(avg_log2FC))

# Save to CSV file
write.csv(top30_markers, "statistics/annotation/top30_markers_per_cluster.csv", row.names = FALSE)
print("Top 30 cluster markers saved to statistics folder")

data.merged <- RenameIdents(data.merged,
    "0" = "Fibroblast",
    "1" = "Fibroblast",
    "2" = "Macrophage (M2-Like)",
    "3" = "Fibroblast",
    "4" = "Osteoblast",
    "5" = "T cell",
    "6" = "Schwann cell",
    "7" = "Plasma cell",
    "8" = "Unknown",
    "9" = "Macrophage",
    "10" = "Unknown",
    "11" = "Schwann cell",
    "12" = "Schwann cell",
    "13" = "Macrophage",
    "14" = "Endothelial cell",
    "15" = "Dendritic cell",
    "16" = "Keratynocyte",
    "17" = "Pericyte",
    "18" = "NK cells",
    "19" = "Monocyte",
    "20" = "Unknown",
    "21" = "Unknown",
    "22" = "Unknown",
    "23" = "Unknown",
    "24" = "Endothelial cell",
    "25" = "Monocyte",
    "26" = "Plasmacytoid DC",
    "27" = "Endothelial cell",
    "28" = "B cells"
)
# Save the annotated Seurat object
saveRDS(data.merged, file = "raw_data/annotation_renamed.rds")