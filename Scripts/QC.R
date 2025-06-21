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


data.merged <- subset(data.merged, subset = nFeature_RNA > 200 & nFeature_RNA < 9500 & percent.mt < 20)
data.merged <- NormalizeData(data.merged)
data.merged <- FindVariableFeatures(data.merged)

plot1 <- VariableFeaturePlot(data.merged)
top10 <- head(VariableFeatures(data.merged), 10)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

ggsave("pictures/QC/variable_features.png", plot2, width = 15, height = 5)

data.merged <- ScaleData(data.merged)
data.merged <- RunPCA(data.merged)

PCA <- ElbowPlot(data.merged)
ggsave("pictures/QC/PCA_2.png", PCA, width = 15, height = 5)
