# Here we take a closer look at ILCP in the FL and compare them with our previous study
SAVE_RDS_FILE = file.path(PATH_EXPERIMENT_OUTPUT_RDSFiles, "FL.rds")
immune.combined = readRDS(SAVE_RDS_FILE)
# Loads the Seurat object saved in the RDS file.

# Comparison with our previous study
# Read gene signatures
Markers_FigS3 = read_xlsx(file.path(PATH_EXPERIMENT_REFERENCE, "Gene_Signatures/Signatures_Clusters_Chapter1.xlsx"))
# Reads an Excel file containing gene signatures from our previous study.

# Extract the best markers based on log fold-change threshold
Markers_FigS3 %>%
  group_by(cluster) %>%
  dplyr::filter(avg_logFC > LOGFC_THR) %>%
  slice_head(n = NUMBER_TOP_SCORING) %>%
  ungroup() -> top10
# Filters the gene signatures to include only those with average log fold-change above the threshold,
# and selects the top-scoring genes per cluster.

# Create a list of genes by cluster
gene_list_by_cluster <- split(top10$gene, top10$cluster)
# Splits the top genes into a named list, grouped by their clusters.

# Name the list based on the cluster numbers
names(gene_list_by_cluster) <- paste0("Cluster_", names(gene_list_by_cluster))
# Renames the clusters with a "Cluster_" prefix for clarity.

# Add module scores for each gene cluster
for (list_name in names(gene_list_by_cluster)) {
  immune.combined = AddModuleScore(
    immune.combined, 
    features = gene_list_by_cluster[list_name], 
    pool = NULL, 
    name = list_name, 
    seed = 19
  )
}
# Adds a module score to the Seurat object for each cluster of genes, indicating expression enrichment.

# Remove the trailing "1" in metadata column names
colnames(immune.combined@meta.data) <- sapply(colnames(immune.combined@meta.data), function(x) {
  if (grepl("^Cluster.*1$", x)) {
    return(gsub("1$", "", x))
  } else {
    return(x)
  }
})
# Cleans up metadata column names by removing any trailing "1" from cluster names.

# Create feature plots for selected clusters
p4 = FeaturePlot(
  immune.combined, 
  features = names(gene_list_by_cluster[1:9]), 
  pt.size = 0.3
) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) & 
  NoAxes()
# Generates a feature plot for the first 9 clusters with a reversed red-blue color gradient and no axes.

p5 = FeaturePlot(
  immune.combined, 
  features = names(gene_list_by_cluster[5:13]), 
  pt.size = 0.3
) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) & 
  NoAxes()
# Generates a feature plot for clusters 5 to 13 with the same styling.

# Save the figures
png(file = file.path(PATH_EXPERIMENT_OUTPUT_GlobalHeteroFigures, "Comparison_Chapter1", "ScoringclustersPart1.png"), 
    width = 35, height = 30, units = "cm", res = 600)
print(p4)
dev.off()
# Saves the first feature plot as a high-resolution PNG file.

pdf(file = file.path(PATH_EXPERIMENT_OUTPUT_GlobalHeteroFigures, "Comparison_Chapter1", "ScoringclustersPart1.pdf"), 
    width = 35/2.54, height = 30/2.54)
print(p4)
dev.off()
# Saves the first feature plot as a PDF file with dimensions converted to inches.

png(file = file.path(PATH_EXPERIMENT_OUTPUT_GlobalHeteroFigures, "Comparison_Chapter1", "ScoringclustersPart2.png"), 
    width = 35, height = 30, units = "cm", res = 600)
print(p5)
dev.off()
# Saves the second feature plot as a high-resolution PNG file.

pdf(file = file.path(PATH_EXPERIMENT_OUTPUT_GlobalHeteroFigures, "Comparison_Chapter1", "ScoringclustersPart2.pdf"), 
    width = 35/2.54, height = 30/2.54)
print(p5)
dev.off()
# Saves the second feature plot as a PDF file with dimensions converted to inches.





