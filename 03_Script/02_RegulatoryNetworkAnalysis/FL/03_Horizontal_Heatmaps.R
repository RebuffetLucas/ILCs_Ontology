#Exploring alternatives

### Analysis of SCENIC Output ####

#Load packages
#For analysis
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
#For some of the plots
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(loomR)
library(dplyr)


RegulonsDvt= c("Runx3(+)","Mybl2(+)","Tbx21(+)","Tcf7(+)","PRDM1(+)","Gata3(+)","Ikzf2(+)","Bcl11b(+)", "Klf6(+)", "ETS1(+)", "Ascl1(+)", "Eomes(+)", "Foxo1(+)")
RegulonsDvt = intersect(rownames(regulonActivity_byCellType_Scaled), RegulonsDvt)

number_regulons_rss = 10
number_regulons_RSS_showtop= 10 #When you  do the RSS plot

number_regulons_showtop = 20 #When you rank them by RSS per group

number_Best_regulons = 60 #When you look at the best regulons without trying to have the same number of reg per group


Z_TRESHOLD = 0.0000001


#Lists to improve based on the litterature
RegulonsDvt= c("Runx3(+)","Mybl2(+)","Tbx21(+)","Tcf7(+)","PRDM1(+)","Gata3(+)","Ikzf2(+)","Bcl11b(+)", "Klf6(+)", "ETS1(+)", "Ascl1(+)", "Eomes(+)", "Foxo1(+)")

RegulonsOfInterest= c("Runx3(+)","Mybl2(+)","Tbx21(+)","Tcf7(+)","PRDM1(+)","Gata3(+)","Ikzf2(+)","Bcl11b(+)", "Klf6(+)", "ETS1(+)", "Ascl1(+)", "Eomes(+)", "Foxo1(+)")

RegulonsOfInterest2 = gsub('.{3}$', '', RegulonsOfInterest)

#For FL:
Meta_info=readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/ILCs_Ontology/05_Output/02_RegulatoryNetworkAnalysis/data/FL/FL_Info.rds")

Meta_info$seurat_clusters
table(Meta_info$seurat_clusters)



vsnDir <- "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/ILCs_Ontology/05_Output/02_RegulatoryNetworkAnalysis/FL"

scenicLoomPath <- paste0(vsnDir, "/auc_mtx.loom")
motifEnrichmentFile <- paste0(vsnDir, "/expr_mat.adjacencies.tsv")
file.exists(scenicLoomPath)
file.exists(motifEnrichmentFile)



#Loading the initial loom file
loom <- open_loom(scenicLoomPath, mode="r+")

regulonsAUC <- get_regulons_AUC(loom, column.attr.name = "RegulonsAUC")

#Add the results of SCENIC analysis in the loom
#read info from loom file
exprMat <- get_dgem(loom)
exprMat_log <- log2(exprMat+1) #Better with log normalization
regulons_incidMat <- get_regulons(loom, column.attr.name = "Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name = "RegulonsAUC")
regulonAUCTresholds = get_regulon_thresholds(loom)
embeddings <- get_embeddings(loom)
#cellClusters <- get_clusterings(loom)    #Indisponible sur cette analyse

close_loom(loom)

#Check before analysis

length(regulons)
head(names(regulons))
regulonAUC



#Load the motif enrichment results
motifEnrichment <- data.table::fread(motifEnrichmentFile, header=T, skip=0)
colnames(motifEnrichment)[1:2] <- c("TF", "MotifID")



##Scoring the network activity

#Regulators for known cell types or clusters
#Average Regulon Activity per cluster


cellClusterPrecise=Meta_info[,"seurat_clusters", drop=FALSE]


#For now only:
Cells_to_keep  = intersect(rownames(Meta_info),  colnames(regulonAUC) )
cellClusterPrecise=Meta_info[Cells_to_keep,"seurat_clusters", drop=FALSE]
selectedResolution <- "seurat_clusters" # select resolution

regulonAUC = regulonAUC[,Cells_to_keep] 
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellClusterPrecise[Cells_to_keep, selectedResolution])


## ANALYSIS FOR PRECISE CLUSTERS
selectedResolution <- "seurat_clusters" # select resolution


# Split the cells by cluster:
cellsPerCluster <- split(rownames(cellClusterPrecise), cellClusterPrecise[,selectedResolution]) 
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]


# Calculate average expression:
regulonActivity_byCellType <- sapply(cellsPerCluster,
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
#Horizontal heatmap with all regulons
# Scale expression:
  #regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

regulonActivity_byCellType_Scaled <- scale(t(regulonActivity_byCellType), center = T, scale=T)

#If a regulon have an std of 0 => We will have NaN on the entire row. So we identify these regulon with an std of 0:
#Identify the NaN genes
# Find rows with NaN values
col_with_nan <- apply(t(regulonActivity_byCellType_Scaled), 1, function(col) any(is.nan(col)))

# Extract row names that have NaN values
colnames_with_nan <- colnames(regulonActivity_byCellType_Scaled)[col_with_nan]

# Print the row names with NaN
print(paste0("regulons with nul std:      " , colnames_with_nan))

#And then we replace with a row of 0.
regulonActivity_byCellType_Scaled[is.nan(regulonActivity_byCellType_Scaled)] <- 0


colnames(regulonActivity_byCellType_Scaled) = gsub('.{3}$', '', colnames(regulonActivity_byCellType_Scaled))




# plot:
options(repr.plot.width=8, repr.plot.height=10) # To set the figure size



#Easy basic
hm <- draw(ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity", column_names_gp=grid::gpar(fontsize=4), 
                                   row_names_gp=grid::gpar(fontsize=15))) # row font size

regulonOrder <- rownames(regulonActivity_byCellType_Scaled)[row_order(hm)] # to save the clustered regulons for later


png(file=file.path(PATH_EXPERIMENT_OUTPUT_RegNet_FL_Figs,"heatmap_AllREGULONS_Horizontal.png"), width = 60, height = 10,  units = "cm", res=600 )
draw(hm)
dev.off()

pdf(file=file.path(PATH_EXPERIMENT_OUTPUT_RegNet_FL_Figs,"heatmap_AllREGULONS_Horizontal.pdf"), width = 60/2.54, height = 10/2.54 )
draw(hm)
dev.off()










#Same_By_Column:

regulonActivity_byCellType_Scaled_colones = t(regulonActivity_byCellType_Scaled)
hm2 <- draw(ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled_colones, name="Regulon activity", column_names_gp=grid::gpar(fontsize=15), 
                                   row_names_gp=grid::gpar(fontsize=5))) # row font size

png(file=paste0(PATH_EXPERIMENT_OUTPUT_RegNet_FL_Figs,"/heatmap_AllREGULONS_Vertical.png"), width = 10, height = 60,  units = "cm", res=600 )
draw(hm2)
dev.off()


pdf(file=file.path(PATH_EXPERIMENT_OUTPUT_RegNet_FL_Figs,"heatmap_AllREGULONS_Vertical.pdf"), width = 60/2.54, height = 10/2.54 )
draw(hm2)
dev.off()




#Horizontal heatmap with only regulons of interest genes ploted
test = regulonActivity_byCellType_Scaled

colnames(test)[!(colnames(test) %in% RegulonsOfInterest2)] = ""


hm <- draw(ComplexHeatmap::Heatmap(test, name="Regulon activity",  column_names_gp=grid::gpar(fontsize=7), rect_gp = gpar(col = "white", lwd = 0.4), border_gp = gpar(col = "darkgrey", lty = 1),
                                   row_names_gp=grid::gpar(fontsize=15))) # row font size

png(file=paste0(PATH_EXPERIMENT_OUTPUT_RegNet_FL_Figs,"/heatmap_AllREGULONS_Horizontal_OnlyRegulonsOfInterest.png"), width = 60, height = 15,  units = "cm", res=600 )
draw(hm)
dev.off()

pdf(file=file.path(PATH_EXPERIMENT_OUTPUT_RegNet_FL_Figs,"heatmap_AllREGULONS_Horizontal_OnlyRegulonsOfInterest.pdf"), width = 60/2.54, height = 10/2.54 )
draw(hm)
dev.off()




#Vertical heatmap including only the regulons of interest

test = regulonActivity_byCellType_Scaled

colnames(test)[!(colnames(test) %in% RegulonsOfInterest2)] = ""

test = t(test)

hm <- draw(ComplexHeatmap::Heatmap(test, name="Regulon activity",  column_names_gp=grid::gpar(fontsize=15), rect_gp = gpar(col = "white", lwd = 0.4), border_gp = gpar(col = "darkgrey", lty = 1),
                                   row_names_gp=grid::gpar(fontsize=7))) # row font size



png(file=paste0(PATH_EXPERIMENT_OUTPUT_RegNet_FL_Figs,"/heatmap_AllREGULONS_Vertical_OnlyRegulonsOfInterest.png"), width = 15, height = 45,  units = "cm", res=600 )
draw(hm)
dev.off()


pdf(file=file.path(PATH_EXPERIMENT_OUTPUT_RegNet_FL_Figs,"heatmap_AllREGULONS_Vertical_OnlyRegulonsOfInterest.pdf"), width = 60/2.54, height = 10/2.54 )
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

png(file=paste0(PATH_EXPERIMENT_OUTPUT_RegNet_FL_Figs,"/heatmap_RegulonsOfInterest_Horizontal.png"), width = 15, height = 10,  units = "cm", res=600 )
draw(hm2)
dev.off()

pdf(file=file.path(PATH_EXPERIMENT_OUTPUT_RegNet_FL_Figs,"heatmap_RegulonsOfInterest_Horizontal.pdf"), width = 60/2.54, height = 10/2.54 )
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
rss[is.na(rss)] <- 0 #In case there are na values
## Showing regulons and cell types with any RSS > 0.01 
rssPlot <- plotRSS(rss, zThreshold= Z_TRESHOLD)

plotly::ggplotly(rssPlot$plot)


options(repr.plot.width=5, repr.plot.height=5) # To set the figure size in Jupyter

#Print and save the RSS plots
for(  group_clust in levels(rssPlot[["df"]][["cellType"]])){
  p1 = plotRSS_oneSet(rss, setName = group_clust, n=number_regulons_rss) # cluster ID
  print(p1)
  
  png(file=paste0(PATH_EXPERIMENT_OUTPUT_RegNet_FL_Figs,"/RSS_Plots/", group_clust , ".png"), width = 40, height = 40,  units = "cm", res=600 )
  print(p1)
  dev.off()
  
  pdf(file=paste0(PATH_EXPERIMENT_OUTPUT_RegNet_FL_Figs,"/RSS_Plots/", group_clust , ".pdf"), width = 40/2.54, height = 40/2.54 )
  print(p1)
  dev.off()
  
}



#Extract the best Regulons and do a heatmap

#Do it from top 10 to top 30 5 by 5:

for (number_regulons_showtop in seq(10, 30, by = 5)){
  
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


png(file=paste0(PATH_EXPERIMENT_OUTPUT_RegNet_FL_Figs,"/heatmap_Best_RSS_Regulons_TOP",number_regulons_showtop,"_Vertical.png"), width = 15, height = 45,  units = "cm", res=600 )
draw(hm)
dev.off()


pdf(file=paste0(PATH_EXPERIMENT_OUTPUT_RegNet_FL_Figs,"/heatmap_Best_RSS_Regulons_TOP",number_regulons_showtop,"_Vertical.pdf"), width = 15/2.54, height = 45/2.54 )
draw(hm)
dev.off()



#Horizontal
topregulonsplot2 = t(topregulonsplot)

colnames(topregulonsplot2) = gsub('.{3}$', '', colnames(topregulonsplot2))

hm <- draw(ComplexHeatmap::Heatmap(topregulonsplot2, name="Regulon activity",  column_names_gp=grid::gpar(fontsize=9), rect_gp = gpar(col = "white", lwd = 0.4), border_gp = gpar(col = "darkgrey", lty = 1),
                                   row_names_gp=grid::gpar(fontsize=15))) # row font size

png(file=paste0(PATH_EXPERIMENT_OUTPUT_RegNet_FL_Figs,"/heatmap_Best_RSS_Regulons_TOP",number_regulons_showtop,"_Horizontal.png"), width = 50, height = 15,  units = "cm", res=600 )
draw(hm)
dev.off()

pdf(file=paste0(PATH_EXPERIMENT_OUTPUT_RegNet_FL_Figs,"/heatmap_Best_RSS_Regulons_TOP",number_regulons_showtop,"_Horizotal.pdf"), width = 50/2.54, height = 15/2.54 )
draw(hm)
dev.off()


}



#Look at the best regulons, without spliting by group

rssPlot$df %>%
  top_n(n = number_Best_regulons, wt = RSS) -> BestRegulons

unique(as.character(BestRegulons$Topic))

bestregulonsplot= regulonActivity_byCellType_Scaled[unique(as.character(BestRegulons$Topic)),]

pheatmap::pheatmap(bestregulonsplot, fontsize_row=12, fontsize_col = 20,
                   color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-2, 2, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color= "grey")  #Average expression regulon activity for top10Regulons

hm = pheatmap::pheatmap(bestregulonsplot, fontsize_row=12, fontsize_col = 20,
                        color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-2, 2, length.out = 100),
                        treeheight_row=10, treeheight_col=10, border_color= "grey")  #Average expression regulon activity for top10Regulons


png(file=paste0(PATH_EXPERIMENT_OUTPUT_RegNet_FL_Figs,"/heatmap_Best_RSS_Regulons_NoGrouping_Top",number_Best_regulons,"_Vertical.png"), width = 15, height = 50,  units = "cm", res=600 )
print(hm)
dev.off()




#Look at all with RSS>0.01
rssPlot$df -> top10Regulons

topregulonsplot= regulonActivity_byCellType_Scaled[unique(as.character(top10Regulons$Topic)),]


pheatmap::pheatmap(topregulonsplot, fontsize_row=12, fontsize_col = 20,
                   color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-2, 2, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color= "grey")  #Average expression regulon activity for top10Regulons







#Will do it if necessary
#Plot for dvt of NK cells




DvtRegulonsplot= regulonActivity_byCellType_Scaled[unique(as.character(RegulonsDvt)),]



hm = pheatmap::pheatmap(DvtRegulonsplot, fontsize_row=16, fontsize_col = 20,
                   color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-2, 2, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color= "grey")  #Average expression regulon activity for top10Regulons
hm


png(file=paste0(PATH_EXPERIMENT_OUTPUT_RegNet_FL_Figs,"/heatmap_Regulons_Development_Vertical.png"), width = 15, height = 30,  units = "cm", res=600 )
print(hm)
dev.off()








#### Vizualisation with umap etc... => Has to done, necessit embedings ! 


#Cell states based on the GRN activity

# List of embeddings available:
cat(names(embeddings), sep="\n")

# Overview of the embeddings (see below for details)

#Plot using AUCell_plotTSNE

#Plot cells present in embeddings and regulonAUC:
colorpalet=brewer.pal(n=8, name="RdBu")

regulonsToPlot <- "Tcf7(+)"
AUCell::AUCell_plotTSNE(embeddings[["umap"]], exprMat_log, regulonAUC[regulonsToPlot,], plots=c("AUC"), cex = .5, exprCols =  colorpalet, offColor = "lightgray") 

png(file=paste0(PATH_EXPERIMENT_OUTPUT_RegNet_FL_Figs,"/umap",regulonsToPlot,".png"), width = 15, height = 15,  units = "cm", res=600 )
AUCell::AUCell_plotTSNE(embeddings[["umap"]], exprMat_log, regulonAUC[regulonsToPlot,], plots=c("AUC"), cex = .5, exprCols =  colorpalet, offColor = "lightgray") 
dev.off()



regulonsToPlot <- "Eomes(+)"
AUCell::AUCell_plotTSNE(embeddings[["umap"]], exprMat_log, regulonAUC[regulonsToPlot,], plots=c("AUC"), cex = .5, exprCols = c("goldenrod1", "darkorange", "brown"))

png(file=paste0(PATH_EXPERIMENT_OUTPUT_RegNet_FL_Figs,"/umap",regulonsToPlot,".png"), width = 15, height = 15,  units = "cm", res=600 )
AUCell::AUCell_plotTSNE(embeddings[["umap"]], exprMat_log, regulonAUC[regulonsToPlot,], plots=c("AUC"), cex = .5, exprCols =  colorpalet, offColor = "lightgray") 
dev.off()

regulonsToPlot <- "Tbx21(+)"
AUCell::AUCell_plotTSNE(embeddings[["umap"]], exprMat_log, regulonAUC[regulonsToPlot,], plots=c("AUC"), cex = .5, exprCols = c("goldenrod1", "darkorange", "brown"))

png(file=paste0(PATH_EXPERIMENT_OUTPUT_RegNet_FL_Figs,"/umap",regulonsToPlot,".png"), width = 15, height = 15,  units = "cm", res=600 )
AUCell::AUCell_plotTSNE(embeddings[["umap"]], exprMat_log, regulonAUC[regulonsToPlot,], plots=c("AUC"), cex = .5, exprCols =  colorpalet, offColor = "lightgray") 
dev.off()

regulonsToPlot <- "Gata3(+)"
AUCell::AUCell_plotTSNE(embeddings[["umap"]], exprMat_log, regulonAUC[regulonsToPlot,], plots=c("AUC"), cex = .5, exprCols = c("goldenrod1", "darkorange", "brown"))

png(file=paste0(PATH_EXPERIMENT_OUTPUT_RegNet_FL_Figs,"/umap",regulonsToPlot,".png"), width = 15, height = 15,  units = "cm", res=600 )
AUCell::AUCell_plotTSNE(embeddings[["umap"]], exprMat_log, regulonAUC[regulonsToPlot,], plots=c("AUC"), cex = .5, exprCols =  colorpalet, offColor = "lightgray") 
dev.off()







#This does not work, need to input some Embeddings here

# Plot cells present in both embeddings an regulonAUC
Embeddings_to_Plot = embeddings[["umap"]][colnames(regulonAUC),]

Embeddings_to_Plot[,"_X"] = -Embeddings_to_Plot[,"_X"]
Embeddings_to_Plot[,"_Y"] = -Embeddings_to_Plot[,"_Y"]

colorpalet=brewer.pal(n=8, name="RdBu")
AUCell::AUCell_plotTSNE(Embeddings_to_Plot, exprMat_log, regulonAUC[regulonsToPlot,], plots=c("AUC"), cex = .5, exprCols =  colorpalet, offColor = "lightgray") 


str(exprMat_log)


length(intersect(rownames(embeddings[["umap"]][colnames(regulonAUC),]) , colnames(regulonAUC) ))


  
  
  
#Plot cells present in embeddings and regulonAUC:
colorpalet=brewer.pal(n=8, name="RdBu")
AUCell::AUCell_plotTSNE(embeddings[["umap"]], exprMat_log, regulonAUC[regulonsToPlot,], plots=c("AUC"), cex = .5, exprCols =  colorpalet, offColor = "lightgray") 

regulonsToPlot <- "Tcf7(+)"
AUCell::AUCell_plotTSNE(embeddings[["umap"]], exprMat_log, regulonAUC[regulonsToPlot,], plots=c("AUC"), cex = .5, exprCols = c("goldenrod1", "darkorange", "brown"))




#Plot all cells in tsne
colorpalet=brewer.pal(n=8, name="RdBu")
AUCell::AUCell_plotTSNE(embeddings[["umap"]], exprMat_log, regulonAUC[regulonsToPlot,], plots=c("AUC"), cex = .5, exprCols =  colorpalet, offColor = "lightgray") 

regulonsToPlot <- "Eomes(+)"
AUCell::AUCell_plotTSNE(embeddings[["umap"]], exprMat_log, regulonAUC[regulonsToPlot,], plots=c("AUC"), cex = .5, exprCols = c("goldenrod1", "darkorange", "brown"))


#Plot of Bertrand # A retravailler car pas ouf du tout
df<-data.frame(umap_1=embeddings[["umap"]][,1],umap_2=embeddings[["umap"]][,2])
df$regulon<-regulonAUC@assays@data$AUC[regulonsToPlot,]
ggplot(df,aes(x=umap_1,y=umap_2,color=regulon))+geom_point()+ scale_color_gradient2(low = "blue",high = "darkred",mid = "white",midpoint = max(df$regulon)/2)+theme_classic()

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


