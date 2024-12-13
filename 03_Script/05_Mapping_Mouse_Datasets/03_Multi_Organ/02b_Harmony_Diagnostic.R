## @knitr Harmony_Diag

## This script performs diagnostic and visualization steps for integrated datasets after batch correction using Harmony.
## It includes UMAP visualizations of Louvain clustering, organ, and organ-cluster distributions, 
## as well as analysis of cluster signatures using marker genes. The script evaluates cluster scores with imported signatures,
## visualizes score distributions using feature and violin plots, and generates bar plots showing the distribution of clusters across organs.
## Figures are saved for further analysis or presentation, providing insights into the integrated dataset's structure and biological relevance.


# Diagnostic

cat(" \n \n")
cat("## Analysis of all sample together after batch correction {.tabset .tabset-fade} \n\n")
cat(" \n \n")

#Bars for ViolinPlots
data_summary <- function(x) {
  m <- median(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}


cat(" \n \n")
cat("### Global overview ")
cat(" \n \n")

p1 = DimPlot(Merged_Seurat_Rescaled, label = TRUE) & ggtitle("Louvain clustering")
p2 = DimPlot(Merged_Seurat_Rescaled, group.by = "Organ")
p3 = DimPlot(Merged_Seurat_Rescaled, group.by = "Organ_Cluster")

cat(" \n \n")
print(p1)
cat(" \n \n")

cat(" \n \n")
print(p2)
cat(" \n \n")


cat(" \n \n")
print(p3) 
cat(" \n \n")




if (DO_SAVE_FIGURE == TRUE){
  
png(file=paste0(SAVE_FIG_PATH , "/Integrated_Analysis/UMAP_Organ_Cluster.png"), width = 20, height = 15,  units = "cm", res=600 )
print(p3)
dev.off()

pdf(file=paste0(SAVE_FIG_PATH , "/Integrated_Analysis/UMAP_Organ_Cluster.pdf"), width = 20/2.54, height = 15/2.54 )
print(p3)
dev.off()


png(file=paste0(SAVE_FIG_PATH , "/Integrated_Analysis/UMAP_Organ.png"), width = 20, height = 15,  units = "cm", res=600 )
print(p2)
dev.off()

pdf(file=paste0(SAVE_FIG_PATH , "/Integrated_Analysis/UMAP_Organ.pdf"), width = 20/2.54, height = 15/2.54 )
print(p2)
dev.off()


png(file=paste0(SAVE_FIG_PATH , "/Integrated_Analysis/UMAP_Louvain_Clust.png"), width = 20, height = 15,  units = "cm", res=600 )
print(p1)
dev.off()

pdf(file=paste0(SAVE_FIG_PATH , "/Integrated_Analysis/UMAP_Louvain_Clust.pdf"), width = 20/2.54, height = 15/2.54 )
print(p1)
dev.off()


}



cat(" \n \n")
cat("### Signatures of the clusters ")
cat(" \n \n")

All_Markers = FindAllMarkers(Merged_Seurat_Rescaled , only.pos = FINDMARKERS_ONLYPOS, method= FINDMARKERS_METHOD , min.pct =  FINDMARKERS_MINPCT, logfc.threshold = FINDMARKERS_LOGFC_THR , verbose = TRUE)

All_Markers %>%
  group_by(cluster) %>%
  slice_max(n = FINDMARKERS_SHOWTOP, order_by = avg_log2FC) -> top10

p10 = DotPlot(Merged_Seurat_Rescaled, features = unique(top10$gene) , cols = "RdBu") + NoLegend() + theme(axis.text.x = element_text(angle = 90)) &ggtitle("Spontaneous signatures")


cat(" \n \n")
print(p10)
cat(" \n \n")


#Scoring with FL signatures

#Import signatures

MONITORED_Markers = readRDS(FILE_SIGNATURES_PATH)

#Extract the best
MONITORED_Markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > LOGFC_THR) %>%
  slice_head(n = NUMBER_TOP_SCORING) %>%
  ungroup() -> top10

# Assuming your tibble is named 'top10'
gene_list_by_cluster <- split(top10$gene, top10$cluster)

# Optionally, to name the list based on the cluster numbers
names(gene_list_by_cluster) <- paste0("Cluster_", names(gene_list_by_cluster))


#Remove the + that generate problems
names(gene_list_by_cluster) <- gsub("\\+", "", names(gene_list_by_cluster))


#Add Module Score for each
for (list_name in names(gene_list_by_cluster) ){
  #print(gene_list_by_cluster[list_name])
  Merged_Seurat_Rescaled = AddModuleScore( Merged_Seurat_Rescaled , features = gene_list_by_cluster[list_name], pool= NULL ,name= list_name , seed=19) 
}

#Remove the 1 at the end of metadata:
colnames(Merged_Seurat_Rescaled@meta.data) <- sapply(colnames(Merged_Seurat_Rescaled@meta.data), function(x) {
  if (grepl("^Cluster.*1$", x)) {
    return(gsub("1$", "", x))
  } else {
    return(x)
  }
})



#Look One by One:

cat(" \n \n")
cat("### Applying the different mouse FL scores {.tabset .tabset-fade} \n\n")
cat(" \n \n")



palette_2 = hue_pal()(length(levels(Merged_Seurat_Rescaled@active.ident)))
names(palette_2) = levels(Merged_Seurat_Rescaled@active.ident)


for (names_signatures in names(gene_list_by_cluster)){
  
  cat(" \n \n")
  cat("#### ",  names_signatures, "\n")
  cat(" \n \n")
  
  p1 = FeaturePlot(Merged_Seurat_Rescaled, features = names_signatures)  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu" )) ) & NoAxes()
  
  cat(" \n \n")
  print(p1)
  cat(" \n \n")
  
  p2 = VlnPlot(Merged_Seurat_Rescaled, features = names_signatures ,sort = "decreasing", pt.size = 0, cols = palette_2) & stat_summary(fun.data=data_summary,color="black")
  
  cat(" \n \n")
  print(p2)
  cat(" \n \n")
  
  
}


#Saving Figures
if (DO_SAVE_FIGURE == TRUE){
  
  for (names_signatures in paste0(names(gene_list_by_cluster))){
    p1 = FeaturePlot(Merged_Seurat_Rescaled, features = names_signatures)  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu" )) ) & NoAxes()
    
    
    p2 = VlnPlot(Merged_Seurat_Rescaled, features = names_signatures ,sort = "decreasing", pt.size = 0, cols = palette_2) & stat_summary(fun.data=data_summary,color="black")
    
    
    p3 = VlnPlot(Merged_Seurat_Rescaled, features = names_signatures ,sort = "decreasing", pt.size = 0, group.by = "seurat_clusters" , cols = palette_2) & stat_summary(fun.data=data_summary,color="black")
    
    p4 = VlnPlot(Merged_Seurat_Rescaled, features = names_signatures ,sort = "decreasing",  group.by = "Organ"  ) & stat_summary(fun.data=data_summary,color="black")
    
    png(file=paste0( SAVE_FIG_PATH ,"/Integrated_Analysis/FeaturePlot_Score_",names_signatures,".png") , width = 30, height = 25,  units = "cm", res=600 )
    print(p1)
    dev.off()
    
    png(file=paste0( SAVE_FIG_PATH,"/Integrated_Analysis/_VlnPlot_",names_signatures,"per_clust.png") , width = 45, height = 25,  units = "cm", res=600 )
    print(p2)
    dev.off()
    
    png(file=paste0( SAVE_FIG_PATH,"/Integrated_Analysis/_VlnPlot_",names_signatures,"per_seurat_clust.png") , width = 45, height = 25,  units = "cm", res=600 )
    print(p3)
    dev.off()
    
    png(file=paste0( SAVE_FIG_PATH,"/Integrated_Analysis/_VlnPlot_",names_signatures,"per_site.png") , width = 45, height = 25,  units = "cm", res=600 )
    print(p4)
    dev.off()
    
    
    
    pdf(file=paste0( SAVE_FIG_PATH ,"/Integrated_Analysis/FeaturePlot_Score_",names_signatures,".pdf") , width = 30/2.54, height = 25/2.54  )
    print(p1)
    dev.off()
    
    pdf(file=paste0( SAVE_FIG_PATH,"/Integrated_Analysis/_VlnPlot_",names_signatures,"per_clust.pdf") , width = 45/2.54, height = 25/2.54  )
    print(p2)
    dev.off()
    
    pdf(file=paste0( SAVE_FIG_PATH,"/Integrated_Analysis/_VlnPlot_",names_signatures,"per_seurat_clust.pdf") , width = 45/2.54, height = 25/2.54  )
    print(p3)
    dev.off()
    
    pdf(file=paste0( SAVE_FIG_PATH,"/Integrated_Analysis/_VlnPlot_",names_signatures,"per_site.pdf") , width = 45/2.54, height = 25/2.54 )
    print(p4)
    dev.off()
    
  }
  
}




cat(" \n \n")
cat("### Proportions before reassignement ")
cat(" \n \n")


#Proportion old clust vs new
Clusters_Proportions = prop.table( table( Merged_Seurat_Rescaled$Organ_Cluster , Merged_Seurat_Rescaled$seurat_clusters) , margin = 1)
Clusters_Proportions = data.frame(Clusters_Proportions)

colnames(Clusters_Proportions) = c("Organ_Cluster", "seurat_clusters", "Freq")

Clusters_Proportions$seurat_clusters = factor(Clusters_Proportions$seurat_clusters)

p12 = Clusters_Proportions %>%
  ggplot(aes(x = seurat_clusters, y = Freq, group = seurat_clusters, color = Organ_Cluster)) +
  geom_boxplot(outlier.shape =  NA, color = palette_2, lwd= 0.9, alpha= 0.1) +
  geom_jitter( position = position_jitter(0.4),  aes(fill = Organ_Cluster) , colour= "black", shape = 21, inherit.aes = TRUE, size= 3) +
  theme_classic() +  theme(text = element_text(size = 15))       


#Bar Diagram of Datasets, Chemistry , Organ_Cluster:

palette_Organ = hue_pal()(length(levels(Merged_Seurat_Rescaled$Organ)))
names(palette_Organ) = levels(Merged_Seurat_Rescaled$Organ)

Merged_Seurat_Rescaled$Organ_Cluster = as.factor(Merged_Seurat_Rescaled$Organ_Cluster)
palette_Organ_Clust = hue_pal()(length(levels(Merged_Seurat_Rescaled$Organ_Cluster)))
names(palette_Organ_Clust) = levels(Merged_Seurat_Rescaled$Organ_Cluster)

Merged_Seurat_Rescaled= SetIdent(Merged_Seurat_Rescaled, value = "seurat_clusters")


#Before renaming

p5 = ggplot(Merged_Seurat_Rescaled@meta.data, aes(x=Organ_Cluster, fill= seurat_clusters)) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90))  + scale_fill_manual(values= palette_2)
p7 = ggplot(Merged_Seurat_Rescaled@meta.data, aes(x=Organ, fill=seurat_clusters)) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90))  + scale_fill_manual(values= palette_2)

cat(" \n \n")
p5 +  p7
cat(" \n \n")


if (DO_SAVE_FIGURE == TRUE){
  
  png(file=paste0(SAVE_FIG_PATH , "/Integrated_Analysis/BarPlot_Before_Renaming/Not_UseFull/BarPlot1.png"), width = 18, height = 15,  units = "cm", res=600 )
  print(p5)
  dev.off()
  
  pdf(file=paste0(SAVE_FIG_PATH , "/Integrated_Analysis/BarPlot_Before_Renaming/Not_UseFull/BarPlot1.pdf"), width = 18/2.54, height = 15/2.54 )
  print(p5)
  dev.off()
  
  
  png(file=paste0(SAVE_FIG_PATH , "/Integrated_Analysis/BarPlot_Before_Renaming/Not_UseFull/BarPlot2.png"), width = 18, height = 15,  units = "cm", res=600 )
  print(p7)
  dev.off()
  
  pdf(file=paste0(SAVE_FIG_PATH , "/Integrated_Analysis/BarPlot_Before_Renaming/Not_UseFull/BarPlot2.pdf"), width = 18/2.54, height = 15/2.54 )
  print(p7)
  dev.off()

}



p6 = ggplot(Merged_Seurat_Rescaled@meta.data, aes(x=seurat_clusters, fill= Organ_Cluster)) + geom_bar(position="fill")    + scale_fill_manual(values= palette_Organ_Clust)
p8 = ggplot(Merged_Seurat_Rescaled@meta.data, aes(x=seurat_clusters, fill=Organ)) + geom_bar(position="fill")    + scale_fill_manual(values= palette_Organ)


cat(" \n \n")
p6 + p8
cat(" \n \n")



if (DO_SAVE_FIGURE == TRUE){
  png(file=paste0(SAVE_FIG_PATH , "/Integrated_Analysis/BarPlot_Before_Renaming/UseFull/BarPlot1.png"), width = 18, height = 15,  units = "cm", res=600 )
  print(p6)
  dev.off()
  
  pdf(file=paste0(SAVE_FIG_PATH , "/Integrated_Analysis/BarPlot_Before_Renaming/UseFull/BarPlot1.pdf"), width = 18/2.54, height = 15/2.54 )
  print(p6)
  dev.off()
  
  
  png(file=paste0(SAVE_FIG_PATH , "/Integrated_Analysis/BarPlot_Before_Renaming/UseFull/BarPlot2.png"), width = 18, height = 15,  units = "cm", res=600 )
  print(p8)
  dev.off()
  
  pdf(file=paste0(SAVE_FIG_PATH , "/Integrated_Analysis/BarPlot_Before_Renaming/UseFull/BarPlot2.pdf"), width = 18/2.54, height = 15/2.54 )
  print(p8)
  dev.off()
  
  
}

