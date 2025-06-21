# Schwann Cell Subtype Gene Expression Analysis
# Analyzes collagen genes and specific gene list in functional Schwann cell subtypes

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# Load the functionally annotated Schwann cell object
schwann_cells <- readRDS("raw_data/schwann_functional_subtypes.rds")

# Create directories
dir.create("pictures/subclustering", recursive = TRUE, showWarnings = FALSE)
dir.create("statistics/subclustering", recursive = TRUE, showWarnings = FALSE)

print("Loaded Schwann cell functional subtypes:")
print(table(Idents(schwann_cells)))

# Define comprehensive collagen gene list
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

# Filter to only genes present in the Schwann cell dataset
available_collagen_genes <- collagen_genes[collagen_genes %in% rownames(schwann_cells)]
print(paste("Found", length(available_collagen_genes), "out of", length(collagen_genes), "collagen genes in Schwann cells"))

if (length(available_collagen_genes) > 0) {
  
# Create dotplot for ALL collagen genes across Schwann cell subtypes
  collagen_dotplot <- DotPlot(schwann_cells, 
                             features = available_collagen_genes,
                             cols = c("lightgrey", "red"),
                             dot.scale = 6,
                             scale = FALSE) +
                     RotatedAxis() +
                     ggtitle("Collagen Genes in Schwann Cell Subtypes") +
                     theme(plot.title = element_text(size = 14, face = "bold"),
                           axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
                           axis.text.y = element_text(size = 16)) +
                     guides(size = guide_legend(title = "% Expressing", 
                                              override.aes = list(color = "black"),
                                              title.theme = element_text(size = 16),
                                              label.theme = element_text(size = 14)),
                            color = guide_colorbar(title = "Avg Expression",
                                                 title.theme = element_text(size = 16),
                                                 label.theme = element_text(size = 14)))
  
    ggsave("pictures/subclustering/schwann_subtypes_collagen_dotplot.pdf", 
           collagen_dotplot, width = max(16, length(available_collagen_genes) * 0.4), height = 10)
  
  }
# Analyze specific gene list in Schwann cell subtypes
gene_list <- c("SORCS3", "STON2", "GPR17", "SOX2", "APOD", "CNTF")
#https://pubmed.ncbi.nlm.nih.gov/15710408/ SORCS3
#https://pubmed.ncbi.nlm.nih.gov/30518424/ STON2
#https://pmc.ncbi.nlm.nih.gov/articles/PMC9434083/ GPR17
#https://pmc.ncbi.nlm.nih.gov/articles/PMC6958085/ CXCL12
#MZB1
#SOX2 - transcription factor 
#APOD
#https://pubmed.ncbi.nlm.nih.gov/26187860/ CNTF
#CCL3L3 chemokine 


# Check which genes are present in the dataset
available_genes <- gene_list[gene_list %in% rownames(schwann_cells)]
print(paste("Available genes from list:", paste(available_genes, collapse = ", ")))
print(paste("Found", length(available_genes), "out of", length(gene_list), "genes"))

if (length(available_genes) > 0) {
  
  # Create dotplot for the gene list across Schwann cell subtypes
  gene_dotplot <- DotPlot(schwann_cells, 
                         features = available_genes,
                         cols = c("lightgrey", "red"),
                         dot.scale = 8,
                         scale = FALSE) +
                 RotatedAxis() +
                 ggtitle("Specific Gene Expression in Schwann Cell Subtypes") +
                 theme(plot.title = element_text(size = 14, face = "bold"),
                       axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
                       axis.text.y = element_text(size = 12)) +
                 guides(size = guide_legend(title = "% Expressing",
                                          override.aes = list(color = "black")),
                        color = guide_colorbar(title = "Avg Expression"))
  
  ggsave("pictures/subclustering/schwann_subtypes_specific_genes_dotplot.pdf", 
         gene_dotplot, width = 12, height = 10)
  
  # Create feature plots for the genes (limit to 6 if more available)
  feature_genes <- available_genes[1:min(6, length(available_genes))]
  
  gene_feature_plots <- list()
  for (i in 1:length(feature_genes)) {
    gene <- feature_genes[i]
    
    p_feature <- FeaturePlot(schwann_cells,
                            features = gene,
                            reduction = "umap",
                            min.cutoff = "q5",
                            max.cutoff = "q95",
                            order = TRUE) +
                ggtitle(paste(gene)) +
                theme(plot.title = element_text(size = 12, face = "bold"),
                      legend.key.size = unit(0.3, "cm"))
    
    gene_feature_plots[[i]] <- p_feature
  }
  
  # Combine feature plots
  combined_gene_features <- wrap_plots(gene_feature_plots, ncol = 3)
  
  ggsave("pictures/subclustering/schwann_subtypes_specific_genes_feature_plot.pdf", 
         combined_gene_features, width = 15, height = 10)
  
  # Create summary statistics for each gene by subtype
  print("\n=== SPECIFIC GENE EXPRESSION ANALYSIS ===")
  
  gene_summary_subtypes <- data.frame()
  
  for (gene in available_genes) {
    print(paste("\n", gene, "expression by subtype:"))
    
    expr_data <- GetAssayData(schwann_cells, slot = "data")[gene,]
    subtypes <- Idents(schwann_cells)
    
    for (subtype in levels(subtypes)) {
      subtype_cells <- names(subtypes)[subtypes == subtype]
      subtype_expr <- expr_data[subtype_cells]
      
      pct_expressing <- sum(subtype_expr > 0) / length(subtype_expr) * 100
      avg_expr <- mean(subtype_expr[subtype_expr > 0])
      if (is.nan(avg_expr)) avg_expr <- 0
      
      print(paste("  ", subtype, ": ", 
                 round(pct_expressing, 1), "% expressing, avg expr:", 
                 round(avg_expr, 3)))
      
      gene_summary_subtypes <- rbind(gene_summary_subtypes, data.frame(
        Gene = gene,
        Subtype = subtype,
        Percent_Expressing = round(pct_expressing, 2),
        Average_Expression = round(avg_expr, 3),
        Cell_Count = length(subtype_cells)
      ))
    }
  }
  
  # Save summary table
  write.csv(gene_summary_subtypes, "statistics/subclustering/schwann_subtypes_specific_genes_summary.csv", row.names = FALSE)
  
  # Find which subtypes express each gene most highly
  print("\n=== TOP EXPRESSING SUBTYPES ===")
  for (gene in available_genes) {
    gene_data <- gene_summary_subtypes[gene_summary_subtypes$Gene == gene,]
    top_subtype <- gene_data[which.max(gene_data$Average_Expression),]
    print(paste(gene, "highest in:", top_subtype$Subtype, 
               "(", top_subtype$Percent_Expressing, "% expressing, avg:", top_subtype$Average_Expression, ")"))
  }
  
  print("\nSpecific gene analysis complete!")
}

print("All gene expression analyses complete!")
print("Files saved in pictures/subclustering/ and statistics/subclustering/")