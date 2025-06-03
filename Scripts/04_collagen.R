# Collagen Gene Expression Dotplot - Annotated Version

library(Seurat)
library(ggplot2)

# Load the annotated Seurat object
data.merged <- readRDS("raw_data/annotation_renamed.rds")

# Define the folder where the plot will be saved
picture_folder <- "pictures"
if (!dir.exists(picture_folder)) {
  dir.create(picture_folder, recursive = TRUE)
}



# Collagen genes
collagen_genes <- c(
  # COL1 family
  "COL1A1", "COL1A2",
  # COL2 family
  "COL2A1",
  # COL3 family
  "COL3A1",
  # COL4 family
  "COL4A1", "COL4A2", "COL4A3", "COL4A4", "COL4A5", "COL4A6",
  # COL5 family
  "COL5A1", "COL5A2", "COL5A3",
  # COL6 family
  "COL6A1", "COL6A2", "COL6A3", "COL6A4", "COL6A5", "COL6A6",
  # COL7-28 families
  "COL7A1", "COL8A1", "COL8A2", "COL9A1", "COL9A2", "COL9A3", "COL10A1",
  "COL11A1", "COL11A2", "COL12A1", "COL13A1", "COL14A1", "COL15A1", "COL16A1", 
  "COL17A1", "COL18A1", "COL19A1", "COL20A1", "COL21A1", "COL22A1", "COL23A1",
  "COL24A1", "COL25A1", "COL26A1", "COL27A1", "COL28A1"
)

# Filter to only genes present in the dataset
available_collagen_genes <- collagen_genes[collagen_genes %in% rownames(data.merged)]
print(paste("Found", length(available_collagen_genes), "out of", length(collagen_genes), "collagen genes"))

# Create dotplot for collagen genes across annotated cell types
dot_plot <- DotPlot(data.merged, 
                   features = available_collagen_genes,
                   cols = c("lightgrey", "red"),
                   dot.scale = 6,
                   scale = FALSE) +
           RotatedAxis() +
           ggtitle("Collagen Genes Expression Across Cell Types") +
           theme(plot.title = element_text(size = 14, face = "bold"),
                 axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
                 axis.text.y = element_text(size = 10)) +
           guides(size = guide_legend(title = "% Expressing",
                                    override.aes = list(color = "black")),
                  color = guide_colorbar(title = "Avg Expression"))

# Save the dotplot
ggsave(paste0(picture_folder, "/collagens/collagen_genes_dotplot.pdf"), 
       dot_plot, width = max(16, length(available_collagen_genes) * 0.4), height = 8)

print("Collagen gene dotplot saved!")