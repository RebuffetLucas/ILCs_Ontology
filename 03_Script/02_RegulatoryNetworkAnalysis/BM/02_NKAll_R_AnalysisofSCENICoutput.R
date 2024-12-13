### Analysis of SCENIC Output for BM ###

# Load necessary packages for analysis and visualization
library(SCopeLoomR) # For loom file handling
library(AUCell) # For AUC-based analysis
library(SCENIC) # SCENIC package for regulatory network analysis
library(KernSmooth) # Kernel smoothing for visualization
library(RColorBrewer) # Color palettes for plots
library(plotly) # Interactive plots
library(BiocParallel) # Parallel processing
library(grid) # Grid-based visualization
library(ComplexHeatmap) # Heatmaps
library(data.table) # Efficient data manipulation
library(loomR) # Additional loom file handling
library(dplyr) # Data manipulation

# Load metadata for the BM dataset
Meta_info = readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/ILCs_Ontology/05_Output/02_RegulatoryNetworkAnalysis/data/BM/BM_Info.rds")
Meta_info$seurat_clusters # Inspect the Seurat clusters in the metadata

# Set the working directory and file paths
vsnDir <- "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/ILCs_Ontology/05_Output/02_RegulatoryNetworkAnalysis/BM"
scenicLoomPath <- file.path(vsnDir, "auc_mtx.loom")
motifEnrichmentFile <- file.path(vsnDir, "expr_mat.adjacencies.tsv")
file.exists(scenicLoomPath) # Check if the loom file exists
file.exists(motifEnrichmentFile) # Check if the motif enrichment file exists

# Load the initial loom file
loom <- open_loom(scenicLoomPath, mode="r+")

# Extract data from the loom file
regulonsAUC <- get_regulons_AUC(loom, column.attr.name = "RegulonsAUC")
exprMat <- get_dgem(loom) # Extract expression matrix
exprMat_log <- log2(exprMat + 1) # Log-normalization of expression matrix
regulons_incidMat <- get_regulons(loom, column.attr.name = "Regulons") # Extract regulons
regulons <- regulonsToGeneLists(regulons_incidMat) # Convert regulons to gene lists
regulonAUC <- get_regulons_AUC(loom, column.attr.name = "RegulonsAUC") # Get AUC matrix
regulonAUCTresholds <- get_regulon_thresholds(loom) # Get regulon thresholds
embeddings <- get_embeddings(loom) # Extract embeddings
close_loom(loom) # Close the loom file after extraction

# Inspect the regulons
length(regulons) # Number of regulons
head(names(regulons)) # Names of the regulons
regulonAUC # AUC matrix for the regulons

# Load motif enrichment results
motifEnrichment <- data.table::fread(motifEnrichmentFile, header = TRUE, skip = 0)
colnames(motifEnrichment)[1:2] <- c("TF", "MotifID")

# Scoring the network activity
# Match cells with clusters from metadata
cellClusterPrecise <- Meta_info[, "seurat_clusters", drop = FALSE]
Cells_to_keep <- intersect(rownames(Meta_info), colnames(regulonAUC))
cellClusterPrecise <- Meta_info[Cells_to_keep, "seurat_clusters", drop = FALSE]
regulonAUC <- regulonAUC[, Cells_to_keep]

# Calculate RSS (Regulon Specificity Score)
rss <- calcRSS(AUC = getAUC(regulonAUC), cellAnnotation = cellClusterPrecise[Cells_to_keep, "seurat_clusters"])

# Calculate average regulon activity by cell type
cellsPerCluster <- split(rownames(cellClusterPrecise), cellClusterPrecise[, "seurat_clusters"])
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(cellsPerCluster, function(cells) rowMeans(getAUC(regulonAUC)[, cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = TRUE, scale = TRUE))

# Plot the scaled regulon activity as a heatmap
options(repr.plot.width = 8, repr.plot.height = 10) # Set figure size
hm <- draw(ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name = "Regulon activity",
                                   row_names_gp = grid::gpar(fontsize = 6))) # Heatmap with custom font size
regulonOrder <- rownames(regulonActivity_byCellType_Scaled)[row_order(hm)] # Save regulon order for later

# Save the heatmap with pheatmap
pheatmap::pheatmap(regulonActivity_byCellType_Scaled,
                   color = colorRampPalette(c("blue", "white", "red"))(100),
                   breaks = seq(-3, 3, length.out = 100),
                   treeheight_row = 10, treeheight_col = 10, border_color = NA)

# Extract top regulators based on activity
topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators$CellType <- factor(as.character(topRegulators$CellType))
topRegulators <- topRegulators[which(topRegulators$RelativeActivity > 0),]
dim(topRegulators) # Dimensions of top regulators

# Visualize RSS for specific clusters
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot) # Interactive RSS plot

# Extract the top regulons for further analysis
rssPlot$df %>%
  dplyr::group_by(cellType) %>%
  top_n(n = 20, wt = RSS) -> top10Regulons

# Plot average expression of top regulons
pheatmap::pheatmap(regulonActivity_byCellType_Scaled[unique(as.character(top10Regulons$Topic)),],
                   fontsize_row = 8, fontsize_col = 20,
                   color = colorRampPalette(c("blue", "white", "red"))(100),
                   breaks = seq(-2, 2, length.out = 100),
                   treeheight_row = 10, treeheight_col = 10, border_color = "grey")

# Analyze regulators for specific clusters and visualize
RegulonsDvt <- c("Runx3(+)", "Mybl2(+)" ,"Tbx21(+)" ,"Tcf7(+)" ,"PRDM1(+)" ,"Gata3(+)","Ikzf2(+)","Bcl11b(+)","Klf6(+)","ETS1(+)","Ascl1(+)","Eomes(+)","Foxo1(+)")
RegulonsDvt <- intersect(rownames(regulonActivity_byCellType_Scaled), RegulonsDvt)
DvtRegulonsplot <- regulonActivity_byCellType_Scaled[unique(as.character(RegulonsDvt)),]

pheatmap::pheatmap(DvtRegulonsplot,
                   fontsize_row = 12, fontsize_col = 20,
                   color = colorRampPalette(c("blue", "white", "red"))(100),
                   breaks = seq(-2, 2, length.out = 100),
                   treeheight_row = 10, treeheight_col = 10, border_color = "grey")



# Plot using AUCell_plotTSNE
# Visualize the regulon activity (AUC) for Tcf7 and Eomes regulons in t-SNE space

regulonsToPlot <- "Tcf7(+)"
colorpalet = brewer.pal(n = 8, name = "RdBu")
AUCell::AUCell_plotTSNE(
  embeddings[["tsne"]], 
  exprMat_log, 
  regulonAUC[regulonsToPlot, ], 
  plots = c("AUC"), 
  cex = 0.5, 
  exprCols = colorpalet, 
  offColor = "lightgray"
)

regulonsToPlot <- "Eomes(+)"
AUCell::AUCell_plotTSNE(
  embeddings[["tsne"]], 
  exprMat_log, 
  regulonAUC[regulonsToPlot, ], 
  plots = c("AUC"), 
  cex = 0.5, 
  exprCols = c("goldenrod1", "darkorange", "brown")
)

# Plot t-SNE using ggplot for Bertrand's plot (adjust as necessary)
df <- data.frame(
  tsne_1 = embeddings[["tsne"]][, 1], 
  tsne_2 = embeddings[["tsne"]][, 2]
)
df$regulon <- regulonAUC@assays@data$AUC[regulonsToPlot, ]

ggplot(df, aes(x = tsne_1, y = tsne_2, color = regulon)) +
  geom_point() +
  scale_color_gradient2(
    low = "blue", 
    high = "darkred", 
    mid = "white", 
    midpoint = max(df$regulon) / 2
  ) +
  theme_classic()

##### The network in detail: TFs, targets, and motifs #####

# Overview of regulons
length(regulons) # Total number of regulons
sum(lengths(regulons) >= 10) # Count regulons with >= 10 genes

viewTable(cbind(nGenes = lengths(regulons)), options = list(pageLength = 10))

# Check if a specific gene is associated with a regulon
grep("Eomes", names(regulons), value = TRUE)
grep("Foxp1", names(regulons), value = TRUE)

# Inspect potential target genes for specific regulons
regulons[["Eomes(+)"]]
regulons[["Tbx21(+)"]]
regulons[["Rel(+)"]]

# Identify potential regulators for specific genes
gene <- "Cd7"
names(regulons)[which(sapply(regulons, function(x) gene %in% x))]

gene <- "Gzmb"
names(regulons)[which(sapply(regulons, function(x) gene %in% x))]

gene <- "Gzma"
names(regulons)[which(sapply(regulons, function(x) gene %in% x))]

# Analyze regulators for multiple genes
dim(regulons_incidMat) # Dimensions of the regulon incidence matrix

genes <- c("Cd7", "Gzmb", "Gzma")
incidMat_subset <- regulons_incidMat[, genes]
incidMat_subset <- incidMat_subset[rowSums(incidMat_subset) > 0, ]
incidMat_subset

# Motifs supporting the regulons
tableSubset <- motifEnrichment[TF == "Eomes"]
viewMotifs(tableSubset, colsToShow = c("logo", "NES", "TF", "Annotation"), options = list(pageLength = 5))
head(tableSubset)

# Regulon targets and motifs
regulons[c("Eomes(+)", "Tbx21(+)")]

# Average regulon activity by cluster
regulonActivity_byCellType <- sapply(
  split(rownames(cellInfo), cellInfo$CellType),
  function(cells) rowMeans(getAUC(regulonAUC)[, cells])
)
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = TRUE, scale = TRUE))

# Plot scaled regulon activity
pheatmap::pheatmap(
  regulonActivity_byCellType_Scaled,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-3, 3, length.out = 100),
  treeheight_row = 10, 
  treeheight_col = 10, 
  border_color = NA
)
