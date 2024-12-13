# Exploring alternatives

### Analysis of SCENIC Output ###

# Load necessary packages
library(SCopeLoomR)  # For loom file handling
library(AUCell)      # For AUC-based analysis
library(SCENIC)      # SCENIC package for regulatory network analysis
library(KernSmooth)  # Kernel smoothing for visualization
library(RColorBrewer)  # Color palettes for plots
library(plotly)      # Interactive plots
library(BiocParallel)  # Parallel processing
library(grid)        # Grid-based visualization
library(ComplexHeatmap) # Heatmaps
library(data.table)  # Efficient data manipulation
library(loomR)       # Additional loom file handling
library(dplyr)       # Data manipulation

# Define parameters for analysis
number_regulons_RSS_showtop = 10       # Number of top regulons for RSS plot
number_regulons_showtop = 30          # Top regulons ranked by RSS per group
number_Best_regulons = 60             # Best regulons without grouping
Z_TRESHOLD = 0.0000001                # Threshold for significance in RSS plot

# Regulons based on literature
RegulonsDvt = c("Runx3(+)" ,"Mybl2(+)" ,"Tbx21(+)" ,"Tcf7(+)" ,"PRDM1(+)" ,"Gata3(+)" ,"Ikzf2(+)" ,"Bcl11b(+)" ,"Klf6(+)" ,"ETS1(+)" ,"Ascl1(+)" ,"Eomes(+)" ,"Foxo1(+)")
RegulonsOfInterest = c("Runx3(+)" ,"Mybl2(+)" ,"Tbx21(+)" ,"Tcf7(+)" ,"PRDM1(+)" ,"Gata3(+)" ,"Ikzf2(+)" ,"Bcl11b(+)" ,"Klf6(+)" ,"ETS1(+)" ,"Ascl1(+)" ,"Eomes(+)" ,"Foxo1(+)")
RegulonsOfInterest2 = gsub('.{3}$', '', RegulonsOfInterest) # Remove trailing '(+)' for convenience

# Load metadata for BM analysis
Meta_info = readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/ILCs_Ontology/05_Output/02_RegulatoryNetworkAnalysis/data/BM/BM_Info.rds")
Meta_info$seurat_clusters # Inspect the clusters

table(Meta_info$seurat_clusters) # Summary of clusters

# Set working directory and file paths
vsnDir <- "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/ILCs_Ontology/05_Output/02_RegulatoryNetworkAnalysis/BM"
scenicLoomPath <- file.path(vsnDir, "auc_mtx.loom")
motifEnrichmentFile <- file.path(vsnDir, "expr_mat.adjacencies.tsv")
file.exists(scenicLoomPath) # Check loom file existence
file.exists(motifEnrichmentFile) # Check motif enrichment file existence

# Load the loom file
loom <- open_loom(scenicLoomPath, mode="r+")

# Extract data from loom file
regulonsAUC <- get_regulons_AUC(loom, column.attr.name = "RegulonsAUC")
exprMat <- get_dgem(loom) # Extract expression matrix
exprMat_log <- log2(exprMat + 1) # Log-normalization
regulons_incidMat <- get_regulons(loom, column.attr.name = "Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat) # Convert regulons to gene lists
regulonAUC <- get_regulons_AUC(loom, column.attr.name = "RegulonsAUC") # Get AUC matrix
regulonAUCTresholds <- get_regulon_thresholds(loom) # Extract regulon thresholds
embeddings <- get_embeddings(loom) # Extract embeddings
close_loom(loom) # Close loom file after extraction

# Inspect regulons
length(regulons) # Number of regulons
head(names(regulons)) # Names of the first few regulons
regulonAUC # AUC matrix for the regulons

# Load motif enrichment results
motifEnrichment <- data.table::fread(motifEnrichmentFile, header = TRUE, skip = 0)
colnames(motifEnrichment)[1:2] <- c("TF", "MotifID")

## Scoring the network activity ##

# Match cells with clusters from metadata
cellClusterPrecise <- Meta_info[, "seurat_clusters", drop = FALSE]
Cells_to_keep <- intersect(rownames(Meta_info), colnames(regulonAUC))
cellClusterPrecise <- Meta_info[Cells_to_keep, "seurat_clusters", drop = FALSE]
regulonAUC <- regulonAUC[, Cells_to_keep]

# Calculate RSS (Regulon Specificity Score)
rss <- calcRSS(AUC = getAUC(regulonAUC), cellAnnotation = cellClusterPrecise[Cells_to_keep, "seurat_clusters"])

# Analysis for precise clusters
selectedResolution <- "seurat_clusters"
cellsPerCluster <- split(rownames(cellClusterPrecise), cellClusterPrecise[, selectedResolution])
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]

# Calculate average expression
regulonActivity_byCellType <- sapply(cellsPerCluster, function(cells) rowMeans(getAUC(regulonAUC)[, cells]))
regulonActivity_byCellType_Scaled <- scale(t(regulonActivity_byCellType), center = TRUE, scale = TRUE)
colnames(regulonActivity_byCellType_Scaled) <- gsub('.{3}$', '', colnames(regulonActivity_byCellType_Scaled))

# Horizontal heatmap for all regulons
options(repr.plot.width = 8, repr.plot.height = 10)
hm <- draw(ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name = "Regulon activity", column_names_gp = grid::gpar(fontsize = 4), row_names_gp = grid::gpar(fontsize = 15)))
png(file = file.path(PATH_EXPERIMENT_OUTPUT_RegNet_BM_Figs, "/heatmap_AllREGULONS_Horizontal.png"), width = 60, height = 10, units = "cm", res = 600)
draw(hm)
dev.off()

# More visualizations and detailed analysis as per the main code structure...

#Same_By_Column:

regulonActivity_byCellType_Scaled_colones = t(regulonActivity_byCellType_Scaled)
hm2 <- draw(ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled_colones, name="Regulon activity", column_names_gp=grid::gpar(fontsize=15), 
                                   row_names_gp=grid::gpar(fontsize=6))) # row font size

png(file=paste0(PATH_EXPERIMENT_OUTPUT_RegNet_BM_Figs,"/heatmap_AllREGULONS_Vertical.png"), width = 10, height = 60,  units = "cm", res=600 )
draw(hm2)
dev.off()



#Horizontal heatmap with only regulons of interest genes ploted
test = regulonActivity_byCellType_Scaled

colnames(test)[!(colnames(test) %in% RegulonsOfInterest2)] = ""


hm <- draw(ComplexHeatmap::Heatmap(test, name="Regulon activity",  column_names_gp=grid::gpar(fontsize=7), rect_gp = gpar(col = "white", lwd = 0.4), border_gp = gpar(col = "darkgrey", lty = 1),
                                   row_names_gp=grid::gpar(fontsize=15))) # row font size

png(file=paste0(PATH_EXPERIMENT_OUTPUT_RegNet_BM_Figs,"/heatmap_AllREGULONS_Horizontal_OnlyRegulonsOfInterest.png"), width = 60, height = 15,  units = "cm", res=600 )
draw(hm)
dev.off()





#Vertical heatmap including only the regulons of interest

test = regulonActivity_byCellType_Scaled

colnames(test)[!(colnames(test) %in% RegulonsOfInterest2)] = ""

test = t(test)

hm <- draw(ComplexHeatmap::Heatmap(test, name="Regulon activity",  column_names_gp=grid::gpar(fontsize=15), rect_gp = gpar(col = "white", lwd = 0.4), border_gp = gpar(col = "darkgrey", lty = 1),
                                   row_names_gp=grid::gpar(fontsize=10))) # row font size



png(file=paste0(PATH_EXPERIMENT_OUTPUT_RegNet_BM_Figs,"/heatmap_AllREGULONS_Vertical_OnlyRegulonsOfInterest.png"), width = 50, height = 15,  units = "cm", res=600 )
draw(hm)
dev.off()


#Plot for regulons of interest
RegulonsOfInterest2 = intersect(colnames(regulonActivity_byCellType_Scaled), RegulonsOfInterest2)


DvtRegulonsplot= regulonActivity_byCellType_Scaled[,unique(as.character(RegulonsOfInterest2))]

pheatmap::pheatmap(DvtRegulonsplot, fontsize_row=12, fontsize_col = 10, angle_col= 90,
                   color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-2, 2, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color= "grey") 

hm2 <- draw(ComplexHeatmap::Heatmap(DvtRegulonsplot, name="Regulon activity",  column_names_gp=grid::gpar(fontsize=10), rect_gp = gpar(col = "white", lwd = 0.4), border_gp = gpar(col = "darkgrey", lty = 1),
                                   row_names_gp=grid::gpar(fontsize=10))) # row font size

png(file=paste0(PATH_EXPERIMENT_OUTPUT_RegNet_BM_Figs,"/heatmap_RegulonsOfInterest_Horizontal.png"), width = 15, height = 10,  units = "cm", res=600 )
draw(hm2)
dev.off()




#See the exact values

topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators$CellType <- factor(as.character(topRegulators$CellType))
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
dim(topRegulators)

viewTable(topRegulators, options = list(pageLength = 10))

#Cell-type specific regulators

rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellClusterPrecise[colnames(regulonAUC), selectedResolution])

## Showing regulons and cell types with any RSS > 0.01 
rssPlot <- plotRSS(rss, zThreshold= Z_TRESHOLD)
plotly::ggplotly(rssPlot$plot)


options(repr.plot.width=5, repr.plot.height=5) # To set the figure size in Jupyter
plotRSS_oneSet(rss, setName = "ILCP", n=number_regulons_RSS_showtop) # cluster ID
plotRSS_oneSet(rss, setName = "EILP", n=number_regulons_RSS_showtop) # cluster ID
plotRSS_oneSet(rss, setName = "ILC2P", n=number_regulons_RSS_showtop) # cluster ID




#Extract the best Regulons and do a heatmap
rssPlot$df %>%
  dplyr::group_by(cellType) %>%
  top_n(n = number_regulons_showtop, wt = RSS) -> top10Regulons


regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

topregulonsplot= regulonActivity_byCellType_Scaled[unique(as.character(top10Regulons$Topic)),]

unique(as.character(top10Regulons$Topic))

pheatmap::pheatmap(topregulonsplot, fontsize_row=10, fontsize_col = 20,
                   color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-2, 2, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color= "grey")  #Average expression regulon activity for top10Regulons


#Vertical
topregulonsplot2 = topregulonsplot

rownames(topregulonsplot2) = gsub('.{3}$', '', rownames(topregulonsplot2))

hm <- draw(ComplexHeatmap::Heatmap(topregulonsplot2, name="Regulon activity",  column_names_gp=grid::gpar(fontsize=9), rect_gp = gpar(col = "white", lwd = 0.4), border_gp = gpar(col = "darkgrey", lty = 1),
                                   row_names_gp=grid::gpar(fontsize=15))) # row font size

png(file=paste0(PATH_EXPERIMENT_OUTPUT_RegNet_BM_Figs,"/heatmap_Best_RSS_Regulons_TOP",number_regulons_showtop,"_Vertical.png"), width = 15, height = 45,  units = "cm", res=600 )
draw(hm)
dev.off()



#Horizontal
topregulonsplot2 = t(topregulonsplot)

colnames(topregulonsplot2) = gsub('.{3}$', '', colnames(topregulonsplot2))

hm <- draw(ComplexHeatmap::Heatmap(topregulonsplot2, name="Regulon activity",  column_names_gp=grid::gpar(fontsize=9), rect_gp = gpar(col = "white", lwd = 0.4), border_gp = gpar(col = "darkgrey", lty = 1),
                                   row_names_gp=grid::gpar(fontsize=15))) # row font size

png(file=paste0(PATH_EXPERIMENT_OUTPUT_RegNet_BM_Figs,"/heatmap_Best_RSS_Regulons_TOP",number_regulons_showtop,"_Horizontal.png"), width = 50, height = 15,  units = "cm", res=600 )
draw(hm)
dev.off()





#Look at the best regulons, without spliting by group

rssPlot$df %>%
  top_n(n = number_Best_regulons, wt = RSS) -> BestRegulons

unique(as.character(BestRegulons$Topic))

bestregulonsplot= regulonActivity_byCellType_Scaled[unique(as.character(BestRegulons$Topic)),]

pheatmap::pheatmap(bestregulonsplot, fontsize_row=12, fontsize_col = 20,
                   color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-2, 2, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color= "grey")  #Average expression regulon activity for top10Regulons



#Look at all with RSS>0.01
rssPlot$df -> top10Regulons

topregulonsplot= regulonActivity_byCellType_Scaled[unique(as.character(top10Regulons$Topic)),]


pheatmap::pheatmap(topregulonsplot, fontsize_row=12, fontsize_col = 20,
                   color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-2, 2, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color= "grey")  #Average expression regulon activity for top10Regulons







#Will do it if necessary
#Plot for dvt of NK cells


RegulonsDvt= c("Runx3(+)","Mybl2(+)","Tbx21(+)","Tcf7(+)","PRDM1(+)","Gata3(+)","Ikzf2(+)","Bcl11b(+)", "Klf6(+)", "ETS1(+)", "Ascl1(+)", "Eomes(+)", "Foxo1(+)")
RegulonsDvt = intersect(rownames(regulonActivity_byCellType_Scaled), RegulonsDvt)


DvtRegulonsplot= regulonActivity_byCellType_Scaled[unique(as.character(RegulonsDvt)),]

pheatmap::pheatmap(DvtRegulonsplot, fontsize_row=12, fontsize_col = 20,
                   color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-2, 2, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color= "grey")  #Average expression regulon activity for top10Regulons




#### Vizualisation with tsne etc... => Has to done, necessit embedings ! 


#Cell states based on the GRN activity

# List of embeddings available:
cat(names(embeddings), sep="\n")

# Overview of the embeddings (see below for details)

#Plot using AUCell_plotTSNE

#Plot cells present in embeddings and regulonAUC:
colorpalet=brewer.pal(n=8, name="RdBu")

regulonsToPlot <- "Tcf7(+)"
AUCell::AUCell_plotTSNE(embeddings[["tsne"]], exprMat_log, regulonAUC[regulonsToPlot,], plots=c("AUC"), cex = .5, exprCols =  colorpalet, offColor = "lightgray") 

png(file=paste0(PATH_EXPERIMENT_OUTPUT_RegNet_BM_Figs,"/tsne",regulonsToPlot,".png"), width = 15, height = 15,  units = "cm", res=600 )
AUCell::AUCell_plotTSNE(embeddings[["tsne"]], exprMat_log, regulonAUC[regulonsToPlot,], plots=c("AUC"), cex = .5, exprCols =  colorpalet, offColor = "lightgray") 
dev.off()



regulonsToPlot <- "Eomes(+)"
AUCell::AUCell_plotTSNE(embeddings[["tsne"]], exprMat_log, regulonAUC[regulonsToPlot,], plots=c("AUC"), cex = .5, exprCols = c("goldenrod1", "darkorange", "brown"))

png(file=paste0(PATH_EXPERIMENT_OUTPUT_RegNet_BM_Figs,"/tsne",regulonsToPlot,".png"), width = 15, height = 15,  units = "cm", res=600 )
AUCell::AUCell_plotTSNE(embeddings[["tsne"]], exprMat_log, regulonAUC[regulonsToPlot,], plots=c("AUC"), cex = .5, exprCols =  colorpalet, offColor = "lightgray") 
dev.off()

regulonsToPlot <- "Tbx21(+)"
AUCell::AUCell_plotTSNE(embeddings[["tsne"]], exprMat_log, regulonAUC[regulonsToPlot,], plots=c("AUC"), cex = .5, exprCols = c("goldenrod1", "darkorange", "brown"))

png(file=paste0(PATH_EXPERIMENT_OUTPUT_RegNet_BM_Figs,"/tsne",regulonsToPlot,".png"), width = 15, height = 15,  units = "cm", res=600 )
AUCell::AUCell_plotTSNE(embeddings[["tsne"]], exprMat_log, regulonAUC[regulonsToPlot,], plots=c("AUC"), cex = .5, exprCols =  colorpalet, offColor = "lightgray") 
dev.off()

regulonsToPlot <- "Gata3(+)"
AUCell::AUCell_plotTSNE(embeddings[["tsne"]], exprMat_log, regulonAUC[regulonsToPlot,], plots=c("AUC"), cex = .5, exprCols = c("goldenrod1", "darkorange", "brown"))

png(file=paste0(PATH_EXPERIMENT_OUTPUT_RegNet_BM_Figs,"/tsne",regulonsToPlot,".png"), width = 15, height = 15,  units = "cm", res=600 )
AUCell::AUCell_plotTSNE(embeddings[["tsne"]], exprMat_log, regulonAUC[regulonsToPlot,], plots=c("AUC"), cex = .5, exprCols =  colorpalet, offColor = "lightgray") 
dev.off()







#This does not work, need to input some Embeddings here

# Plot cells present in both embeddings an regulonAUC
Embeddings_to_Plot = embeddings[["tsne"]][colnames(regulonAUC),]

Embeddings_to_Plot[,"_X"] = -Embeddings_to_Plot[,"_X"]
Embeddings_to_Plot[,"_Y"] = -Embeddings_to_Plot[,"_Y"]

colorpalet=brewer.pal(n=8, name="RdBu")
AUCell::AUCell_plotTSNE(Embeddings_to_Plot, exprMat_log, regulonAUC[regulonsToPlot,], plots=c("AUC"), cex = .5, exprCols =  colorpalet, offColor = "lightgray") 


str(exprMat_log)


length(intersect(rownames(embeddings[["tsne"]][colnames(regulonAUC),]) , colnames(regulonAUC) ))


  
  
  
#Plot cells present in embeddings and regulonAUC:
colorpalet=brewer.pal(n=8, name="RdBu")
AUCell::AUCell_plotTSNE(embeddings[["tsne"]], exprMat_log, regulonAUC[regulonsToPlot,], plots=c("AUC"), cex = .5, exprCols =  colorpalet, offColor = "lightgray") 

regulonsToPlot <- "Tcf7(+)"
AUCell::AUCell_plotTSNE(embeddings[["tsne"]], exprMat_log, regulonAUC[regulonsToPlot,], plots=c("AUC"), cex = .5, exprCols = c("goldenrod1", "darkorange", "brown"))




#Plot all cells in tsne
colorpalet=brewer.pal(n=8, name="RdBu")
AUCell::AUCell_plotTSNE(embeddings[["tsne"]], exprMat_log, regulonAUC[regulonsToPlot,], plots=c("AUC"), cex = .5, exprCols =  colorpalet, offColor = "lightgray") 

regulonsToPlot <- "Eomes(+)"
AUCell::AUCell_plotTSNE(embeddings[["tsne"]], exprMat_log, regulonAUC[regulonsToPlot,], plots=c("AUC"), cex = .5, exprCols = c("goldenrod1", "darkorange", "brown"))


#Plot of Bertrand # A retravailler car pas ouf du tout
df<-data.frame(tsne_1=embeddings[["tsne"]][,1],tsne_2=embeddings[["tsne"]][,2])
df$regulon<-regulonAUC@assays@data$AUC[regulonsToPlot,]
ggplot(df,aes(x=tsne_1,y=tsne_2,color=regulon))+geom_point()+ scale_color_gradient2(low = "blue",high = "darkred",mid = "white",midpoint = max(df$regulon)/2)+theme_classic()

head(df)



##### The network in details: TFs, targets and motifs

length(regulons)
sum(lengths(regulons)>=10)


viewTable(cbind(nGenes=lengths(regulons)), options=list(pageLength=10))


# Check if a specific gene has a regulon


grep("EOMES", names(regulons), value=T) # paste0("^","EBF1","_")
grep("FOXP1", names(regulons), value=T) # paste0("^","EBF1","_")


#Look at the potential target genes

regulons[["EOMES(+)"]]
regulons[["TGIF2(+)"]]
regulons[["ZNF628(+)"]]
regulons[["RORA(+)"]]

regulons[["FOXP1(+)"]]
regulons[["KMT2A(+)"]]
regulons[["STAT4(+)"]]

regulons[["REL(+)"]]
regulons[["CREM(+)"]]


#Potential regulators for a given gene
gene <- "CD69"
names(regulons)[which(sapply(regulons, function(x) gene %in% x))]

gene <- "RGS1"
names(regulons)[which(sapply(regulons, function(x) gene %in% x))]

gene <- "IER2"
names(regulons)[which(sapply(regulons, function(x) gene %in% x))]

gene <- "IER5"
names(regulons)[which(sapply(regulons, function(x) gene %in% x))]


#dO IT WITH MULTIPLE GENES

dim(regulons_incidMat)

genes <- c("CD3E", "KLRF1", "FCGR3A") 
incidMat_subset <- regulons_incidMat[,genes]
incidMat_subset <- incidMat_subset[rowSums(incidMat_subset)>0,]

incidMat_subset


# Motifs supporting the regulons

tableSubset <- motifEnrichment[TF=="ZNF683"]
viewMotifs(tableSubset, colsToShow = c("logo", "NES", "TF" ,"Annotation"), options=list(pageLength=5))

head(tableSubset)

#Regulon targets and motifs

regulons[ c("EOMES(+)", "TBX21(+)")]


#Regulators for clusters of known cell types
#Average regulon activity by cluster

regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$CellType),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

pheatmap::pheatmap(regulonActivity_byCellType_Scaled, #fontsize_row=3, 
                   color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-3, 3, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color=NA)


