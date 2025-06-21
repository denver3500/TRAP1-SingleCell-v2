# MPNST Tumor Marker Analysis by Cell Type and TRAP1 Expression
# Analyzes key gene sets across key cell types with TRAP1 stratification

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(tidyr)

# Load the annotated Seurat object
data.merged <- readRDS("raw_data/annotation_renamed.rds")

# Create output directories
dir.create("pictures/markers", recursive = TRUE, showWarnings = FALSE)
dir.create("statistics/markers", recursive = TRUE, showWarnings = FALSE)

# Define marker sets
# Invasion markers
invasion_markers <- c(
  # MMPs
  "MMP1", "MMP2", "MMP3", "MMP7", "MMP9", "MMP13", "MMP14",
  # TIMPs
  "TIMP1", "TIMP2", "TIMP3", "TIMP4",
  # uPA system
  "PLAU", "PLAUR",
  # Cathepsins
  "CTSB", "CTSL", "CTSD", "CTSK",
  # Other invasion markers
  "CD44", "ADAM10", "ADAM17"
)

# Metastasis markers
metastasis_markers <- c(
  # EMT markers
  "CDH1", "CDH2", "VIM", "SNAI1", "SNAI2", "ZEB1", "ZEB2", "TWIST1", "TWIST2",
  # Cell motility
  "RAC1", "RHOA", "CDC42",
  # Metastasis-specific genes
  "S100A4", "MALAT1", "MTA1",
  # Chemokine receptors
  "CXCR4", "CCR7"
)

# Angiogenesis markers
angiogenesis_markers <- c(
  # Growth factors
  "VEGFA", "VEGFB", "VEGFC", "FIGF", "PGF", "FGF2", "PDGFA", "PDGFB", "ANGPT1", "ANGPT2",
  # Receptors
  "KDR", "FLT1", "FLT4", "TEK", "PDGFRA", "PDGFRB",
  # Hypoxia-related
  "HIF1A", "EPAS1", "ARNT",
  # Endothelial markers
  "PECAM1", "CDH5", "ENG",
  # Inhibitors
  "THBS1", "THBS2", "SERPINF1"
)

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
  "COL6A1", "COL6A2", "COL6A3",
  # COL7-28 families
  "COL7A1", "COL8A1", "COL8A2", "COL9A1", "COL9A2", "COL9A3", "COL10A1",
  "COL11A1", "COL11A2", "COL12A1", "COL13A1", "COL14A1", "COL15A1", "COL16A1", 
  "COL17A1", "COL18A1", "COL19A1", "COL20A1", "COL21A1", "COL22A1", "COL23A1",
  "COL24A1", "COL25A1", "COL26A1", "COL27A1", "COL28A1"
)

# Check all marker sets for available genes
print("Checking marker availability in dataset...")
all_markers <- list(
  "Invasion" = invasion_markers,
  "Metastasis" = metastasis_markers,
  "Angiogenesis" = angiogenesis_markers,
  "Collagen" = collagen_genes
)

for (name in names(all_markers)) {
  markers <- all_markers[[name]]
  available <- markers[markers %in% rownames(data.merged)]
  print(paste(name, "markers:", length(available), "of", length(markers), "found"))
}

# Step 1: Subset for specific cell types
selected_cell_types <- c("Schwann cell", "Macrophage (M2-Like)", "Macrophage", "Fibroblast")
data.subset <- subset(data.merged, idents = selected_cell_types)

# Explicitly set identities to ensure they're maintained
Idents(data.subset) <- factor(as.character(Idents(data.subset)), levels = selected_cell_types)

print("Cell type subset summary:")
print(table(Idents(data.subset)))

# Step 2: Create a UMAP of the subset
p_subset <- DimPlot(data.subset, reduction = "umap", 
                   label = TRUE, label.size = 4, repel = TRUE) +
           ggtitle("Selected Cell Types") +
           theme(plot.title = element_text(size = 16, face = "bold"))

ggsave("pictures/markers/selected_celltypes_umap.pdf", p_subset, width = 10, height = 8)
ggsave("pictures/markers/selected_celltypes_umap.png", p_subset, width = 10, height = 8, dpi = 300)

# Step 3: Check if TRAP1 is present in the dataset
if ("TRAP1" %in% rownames(data.subset)) {
  print("TRAP1 found in dataset. Proceeding with TRAP1 high/low analysis.")
  
  # Create TRAP1 expression feature plot
  p_trap1 <- FeaturePlot(data.subset, features = "TRAP1") + 
            ggtitle("TRAP1 Expression") +
            theme(plot.title = element_text(size = 14, face = "bold"))
  
  ggsave("pictures/markers/TRAP1_expression.pdf", p_trap1, width = 8, height = 7)
  ggsave("pictures/markers/TRAP1_expression.png", p_trap1, width = 8, height = 7, dpi = 300)
  
  # Get TRAP1 expression values
  trap1_expr <- GetAssayData(data.subset, slot = "data")["TRAP1",]
  
  # Define high/low based on median expression among expressing cells
  cells_expressing <- names(trap1_expr)[trap1_expr > 0]
  if (length(cells_expressing) > 0) {
    median_expr <- median(trap1_expr[cells_expressing])
    
    # Add TRAP1 status to metadata
    data.subset$TRAP1_status <- ifelse(trap1_expr > median_expr, "TRAP1_high", "TRAP1_low")
    
    # Create violin plot of TRAP1 expression
    p_violin <- VlnPlot(data.subset, features = "TRAP1", group.by = "TRAP1_status") +
              ggtitle("TRAP1 Expression by Group") +
              theme(plot.title = element_text(size = 14, face = "bold"))
    
    ggsave("pictures/markers/TRAP1_violin_by_status.pdf", p_violin, width = 8, height = 6)
    ggsave("pictures/markers/TRAP1_violin_by_status.png", p_violin, width = 8, height = 6, dpi = 300)
    
    # Create UMAP colored by TRAP1 status
    p_trap1_status <- DimPlot(data.subset, reduction = "umap", group.by = "TRAP1_status") +
                     ggtitle("TRAP1 Status") +
                     theme(plot.title = element_text(size = 14, face = "bold"))
    
    ggsave("pictures/markers/TRAP1_status_umap.pdf", p_trap1_status, width = 8, height = 7)
    ggsave("pictures/markers/TRAP1_status_umap.png", p_trap1_status, width = 8, height = 7, dpi = 300)
    
    # UMAP split by cell type and colored by TRAP1 status
    p_trap1_by_celltype <- DimPlot(data.subset, reduction = "umap", 
                                 group.by = "TRAP1_status", 
                                 split.by = "ident") +
                         ggtitle("TRAP1 Status by Cell Type") +
                         theme(plot.title = element_text(size = 14, face = "bold"))
    
    ggsave("pictures/markers/TRAP1_status_by_celltype.pdf", p_trap1_by_celltype, width = 16, height = 6)
    ggsave("pictures/markers/TRAP1_status_by_celltype.png", p_trap1_by_celltype, width = 16, height = 6, dpi = 300)
    
    # TRAP1 high/low status available for further analysis
    trap1_analysis_possible <- TRUE
  } else {
    print("TRAP1 is not expressed in any cells in the subset. Skipping TRAP1 high/low analysis.")
    trap1_analysis_possible <- FALSE
  }
} else {
  print("TRAP1 not found in dataset. Skipping TRAP1 high/low analysis.")
  trap1_analysis_possible <- FALSE
}

# Step 4: Generate dotplots for marker sets without TRAP1 subsetting
print("Generating dotplots for all marker sets (no TRAP1 subsetting)...")

# Function to generate and save dotplots for a marker set
generate_dotplot <- function(marker_set, marker_name, data_obj, group_var = NULL) {
  # Filter for available markers
  available_markers <- marker_set[marker_set %in% rownames(data_obj)]
  
  if (length(available_markers) == 0) {
    print(paste("No", marker_name, "markers found in dataset"))
    return(NULL)
  }
  
  # Generate dotplot
  if (!is.null(group_var)) {
    p <- DotPlot(data_obj, 
                features = available_markers,
                cols = c("lightgrey", "red"),
                dot.scale = 8,
                scale = FALSE,  # Add this to prevent scaling warnings
                group.by = group_var)
  } else {
    # If no group_var provided, use current identities
    p <- DotPlot(data_obj, 
                features = available_markers,
                cols = c("lightgrey", "red"),
                dot.scale = 8,
                scale = FALSE)  # Add this to prevent scaling warnings
  }
  
  p <- p + RotatedAxis() +
          ggtitle(paste(marker_name, "Markers")) +
          theme(plot.title = element_text(size = 14, face = "bold"),
                axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
                axis.text.y = element_text(size = 10)) +
          guides(size = guide_legend(title = "% Expressing",
                                   override.aes = list(color = "black")),
                 color = guide_colorbar(title = "Avg Expression"))
  
  return(p)
}

# Generate dotplots for each marker set
for (i in seq_along(all_markers)) {
  name <- names(all_markers)[i]
  markers <- all_markers[[i]]
  
  # Generate dotplot
  p <- generate_dotplot(markers, name, data.subset)
  
  if (!is.null(p)) {
    # Save as PDF and PNG
    ggsave(paste0("pictures/markers/", tolower(name), "_markers_dotplot.pdf"), 
           p, width = max(12, length(markers) * 0.3), height = 8)
    ggsave(paste0("pictures/markers/", tolower(name), "_markers_dotplot.png"), 
           p, width = max(12, length(markers) * 0.3), height = 8, dpi = 300)
  }
}

# Step 5: If TRAP1 subsetting is possible, generate separate dotplots for each cell type
if (trap1_analysis_possible) {
  print("Generating separate dotplots for each cell type with TRAP1 high/low...")
  
  # Create a new identity column combining cell type and TRAP1 status
  data.subset$celltype_trap1 <- paste(Idents(data.subset), data.subset$TRAP1_status, sep = "_")
  
  # Save cell counts
  cell_counts <- table(data.subset$celltype_trap1)
  write.csv(as.data.frame(cell_counts), "statistics/markers/cell_counts_by_celltype_TRAP1.csv")
  
  # For each cell type, create a separate dotplot
  for (cell_type in selected_cell_types) {
    print(paste("Processing", cell_type))
    
    # Subset for the current cell type
    cells_of_type <- WhichCells(data.subset, idents = cell_type)
    
    # Skip if insufficient cells
    if (length(cells_of_type) < 20) {
      print(paste("Skipping", cell_type, "due to insufficient cells"))
      next
    }
    
    # Create a temp object with just this cell type
    temp <- subset(data.subset, cells = cells_of_type)
    
    # Set identity to TRAP1 status
    Idents(temp) <- temp$TRAP1_status
    
    # Count cells in each group
    trap1_counts <- table(Idents(temp))
    print(paste("  TRAP1_high:", trap1_counts["TRAP1_high"], "cells"))
    print(paste("  TRAP1_low:", trap1_counts["TRAP1_low"], "cells"))
    
    # Skip if either group has too few cells
    if (any(trap1_counts < 10)) {
      print(paste("  Skipping", cell_type, "due to insufficient cells in TRAP1 groups"))
      next
    }
    
    # Clean cell type name for filenames
    clean_name <- gsub(" |\\(|\\)", "_", cell_type)
    
    # Generate dotplots for each marker set for this cell type
    for (i in seq_along(all_markers)) {
      name <- names(all_markers)[i]
      markers <- all_markers[[i]]
      
      # Filter for available markers
      available_markers <- markers[markers %in% rownames(temp)]
      
      if (length(available_markers) == 0) {
        print(paste("  No", name, "markers found for", cell_type))
        next
      }
      
      # Generate dotplot for this cell type
      p <- DotPlot(temp, 
                  features = available_markers,
                  cols = c("lightgrey", "red"),
                  dot.scale = 8,
                  scale = FALSE) +
          RotatedAxis() +
          ggtitle(paste(name, "Markers in", cell_type, "by TRAP1 Status")) +
          theme(plot.title = element_text(size = 14, face = "bold"),
                axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
                axis.text.y = element_text(size = 10)) +
          guides(size = guide_legend(title = "% Expressing",
                                    override.aes = list(color = "black")),
                 color = guide_colorbar(title = "Avg Expression"))
      
      # Save as PDF and PNG
      ggsave(paste0("pictures/markers/", clean_name, "_", tolower(name), "_TRAP1_dotplot.pdf"), 
             p, width = max(10, length(available_markers) * 0.3), height = 6)
      ggsave(paste0("pictures/markers/", clean_name, "_", tolower(name), "_TRAP1_dotplot.png"), 
             p, width = max(10, length(available_markers) * 0.3), height = 6, dpi = 300)
      
      print(paste("  Generated", name, "dotplot for", cell_type))
    }
  }
  
  print("Individual cell type TRAP1 dotplots complete!")

  # Step 6: Identify statistically significant markers between TRAP1 high and low groups
  print("Identifying statistically significant markers between TRAP1 high and low groups...")
    
  # Combine all markers for statistical testing
  all_combined_markers <- unique(unlist(all_markers))
  available_combined_markers <- all_combined_markers[all_combined_markers %in% rownames(data.subset)]
    
  # For each cell type, find differentially expressed markers between TRAP1 high and low
  cell_types <- unique(Idents(data.subset))
    
  significant_markers_summary <- data.frame()
    
  for (cell_type in cell_types) {
    print(paste("Analyzing", cell_type))
    
    # Subset for the current cell type
    cells_of_type <- WhichCells(data.subset, idents = cell_type)
    
    # Get TRAP1 high and low cells of this type
    high_cells <- cells_of_type[data.subset$TRAP1_status[cells_of_type] == "TRAP1_high"]
    low_cells <- cells_of_type[data.subset$TRAP1_status[cells_of_type] == "TRAP1_low"]
    
    # Skip if insufficient cells in either group
    if (length(high_cells) < 10 || length(low_cells) < 10) {
      print(paste("Skipping", cell_type, "due to insufficient cells in TRAP1 high/low groups"))
      next
    }
    
    # Create a temporary object with TRAP1 status as identity
    temp <- subset(data.subset, cells = c(high_cells, low_cells))
    Idents(temp) <- temp$TRAP1_status
    
    # Find markers differentially expressed between TRAP1 high and low
    # Focusing only on our marker sets
    de_markers <- FindMarkers(temp, 
                             features = available_combined_markers,
                             ident.1 = "TRAP1_high", 
                             ident.2 = "TRAP1_low",
                             min.pct = 0.1,
                             test.use = "wilcox")
    
    if (nrow(de_markers) > 0) {
      # Add cell type and filter for significant markers
      de_markers$gene <- rownames(de_markers)
      de_markers$cell_type <- cell_type
      significant_de <- de_markers[de_markers$p_val_adj < 0.05,]
      
      # Only process if we have significant markers
      if (nrow(significant_de) > 0) {
        # Add marker category
        significant_de$category <- "Other"
        for (name in names(all_markers)) {
          markers <- all_markers[[name]]
          significant_de$category[significant_de$gene %in% markers] <- name
        }
        
        # Add to summary
        significant_markers_summary <- rbind(significant_markers_summary, significant_de)
      }
      
      # Save full results
      write.csv(de_markers, paste0("statistics/markers/", gsub(" |\\(|\\)", "_", cell_type), 
                                 "_TRAP1_DE_markers.csv"))
    }
  }
    
  # Save significant markers summary
  if (nrow(significant_markers_summary) > 0) {
    write.csv(significant_markers_summary, "statistics/markers/significant_markers_TRAP1_all_celltypes.csv")
    
    # Generate dotplot of significant markers per cell type
    # Group by cell type
    for (cell_type in unique(significant_markers_summary$cell_type)) {
      cell_type_markers <- significant_markers_summary[significant_markers_summary$cell_type == cell_type,]
      
      if (nrow(cell_type_markers) > 0) {
        clean_name <- gsub(" |\\(|\\)", "_", cell_type)
        sig_genes <- unique(cell_type_markers$gene)
        
        # Subset for this cell type
        cells_of_type <- WhichCells(data.subset, idents = cell_type)
        temp <- subset(data.subset, cells = cells_of_type)
        Idents(temp) <- temp$TRAP1_status
        
        # Create dotplot
        p_sig <- DotPlot(temp, 
                        features = sig_genes,
                        cols = c("lightgrey", "red"),
                        dot.scale = 8,
                        scale = FALSE) +
                RotatedAxis() +
                ggtitle(paste("Significant Markers in", cell_type, "by TRAP1 Status")) +
                theme(plot.title = element_text(size = 14, face = "bold"),
                      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
                      axis.text.y = element_text(size = 8)) +
                guides(size = guide_legend(title = "% Expressing",
                                         override.aes = list(color = "black")),
                      color = guide_colorbar(title = "Avg Expression"))
        
        # Save as PDF and PNG
        ggsave(paste0("pictures/markers/", clean_name, "_significant_markers_TRAP1_dotplot.pdf"), 
               p_sig, width = max(10, length(sig_genes) * 0.4), height = 6)
        ggsave(paste0("pictures/markers/", clean_name, "_significant_markers_TRAP1_dotplot.png"), 
               p_sig, width = max(10, length(sig_genes) * 0.4), height = 6, dpi = 300)
        
        print(paste("Generated significant marker dotplot for", cell_type))
      }
    }
  } else {
    print("No significant markers found between TRAP1 high and low groups")
  }
} else {
  print("Skipping TRAP1 high/low analysis due to insufficient data")
}

print("Analysis complete!")
print("Results saved in pictures/markers/ and statistics/markers/")