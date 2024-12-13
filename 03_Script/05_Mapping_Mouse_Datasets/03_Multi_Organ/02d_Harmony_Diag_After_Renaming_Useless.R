
#Renaming

Merged_Seurat_Rescaled$seurat_clusters2 = Merged_Seurat_Rescaled$seurat_clusters
levels(Merged_Seurat_Rescaled$seurat_clusters2) = RENAME_CLUSTERS
Merged_Seurat_Rescaled = SetIdent(Merged_Seurat_Rescaled, value= "seurat_clusters2")


#Bar Diagram of Datasets, Chemistry , Organ_Cluster:

palette_Organ = hue_pal()(length(levels(Merged_Seurat_Rescaled$Organ)))
names(palette_Organ) = levels(Merged_Seurat_Rescaled$Organ)

Merged_Seurat_Rescaled$Organ_Cluster = as.factor(Merged_Seurat_Rescaled$Organ_Cluster)
palette_Organ_Clust = hue_pal()(length(levels(Merged_Seurat_Rescaled$Organ_Cluster)))
names(palette_Organ_Clust) = levels(Merged_Seurat_Rescaled$Organ_Cluster)

table(Merged_Seurat_Rescaled$seurat_clusters2 )

#After renaming

p5 = ggplot(Merged_Seurat_Rescaled@meta.data, aes(x=Organ_Cluster, fill= seurat_clusters2)) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90))  + scale_fill_manual(values= palette_3)
p7 = ggplot(Merged_Seurat_Rescaled@meta.data, aes(x=Organ, fill=seurat_clusters2)) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90))  + scale_fill_manual(values= palette_3)

cat(" \n \n")
p5 +  p7
cat(" \n \n")


if (DO_SAVE_FIGURE == TRUE){
  
  png(file=paste0(SAVE_FIG_PATH , "/Integrated_Analysis/BarPlot_After_Renaming/UseFull/BarPlot1.png"), width = 18, height = 15,  units = "cm", res=600 )
  print(p5)
  dev.off()
  
  pdf(file=paste0(SAVE_FIG_PATH , "/Integrated_Analysis/BarPlot_After_Renaming/UseFull/BarPlot1.pdf"), width = 18/2.54, height = 15/2.54 )
  print(p5)
  dev.off()
  
  
  png(file=paste0(SAVE_FIG_PATH , "/Integrated_Analysis/BarPlot_After_Renaming/UseFull/BarPlot2.png"), width = 18, height = 15,  units = "cm", res=600 )
  print(p7)
  dev.off()
  
  pdf(file=paste0(SAVE_FIG_PATH , "/Integrated_Analysis/BarPlot_After_Renaming/UseFull/BarPlot2.pdf"), width = 18/2.54, height = 15/2.54 )
  print(p7)
  dev.off()
  
}



p6 = ggplot(Merged_Seurat_Rescaled@meta.data, aes(x=seurat_clusters2, fill= Organ_Cluster)) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90))  + scale_fill_manual(values= palette_Organ_Clust)
p8 = ggplot(Merged_Seurat_Rescaled@meta.data, aes(x=seurat_clusters2, fill=Organ)) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90))  + scale_fill_manual(values= palette_Organ)


cat(" \n \n")
p6 + p8
cat(" \n \n")



if (DO_SAVE_FIGURE == TRUE){
  png(file=paste0(SAVE_FIG_PATH , "/Integrated_Analysis/BarPlot_After_Renaming/UseFull/BarPlot1.png"), width = 18, height = 15,  units = "cm", res=600 )
  print(p6)
  dev.off()
  
  pdf(file=paste0(SAVE_FIG_PATH , "/Integrated_Analysis/BarPlot_After_Renaming/UseFull/BarPlot1.pdf"), width = 18/2.54, height = 15/2.54 )
  print(p6)
  dev.off()
  
  
  png(file=paste0(SAVE_FIG_PATH , "/Integrated_Analysis/BarPlot_After_Renaming/UseFull/BarPlot2.png"), width = 18, height = 15,  units = "cm", res=600 )
  print(p8)
  dev.off()
  
  pdf(file=paste0(SAVE_FIG_PATH , "/Integrated_Analysis/BarPlot_After_Renaming/UseFull/BarPlot2.pdf"), width = 18/2.54, height = 15/2.54 )
  print(p8)
  dev.off()
  
  
}

