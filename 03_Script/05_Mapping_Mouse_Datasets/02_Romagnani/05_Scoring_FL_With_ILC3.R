## @knitr Scoring_FL_ILC3_Roma


#Def functions
#Bars for ViolinPlots
data_summary <- function(x) {
  m <- median(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

#Score Human Lymphoid cells in different tissues



cat(" \n \n")
cat("## Score the FL cells with Chiara's ILC3 signatures")
cat(" \n \n")


#Diagnostic with their granularity

#First remove the cells not used

PBMC = readRDS(paste0(PATH_PROJECT,SAMPLE_CHAP2_FL))

MONITORED_Markers = read.csv( SIGNATURES_ROMAGNANI_ILC3, row.names = NULL)

PBMC$Cluster = as.factor(PBMC$seurat_clusters )
PBMC=SetIdent(PBMC, value= "Cluster")


palette_2 = hue_pal()(length(levels(PBMC@active.ident)))
names(palette_2) = levels(PBMC@active.ident)



#Extract the best
MONITORED_Markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_logFC > LOGFC_THR) %>%
  slice_head(n = NUMBER_TOP_SCORING_ILC3) %>%
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
  PBMC = AddModuleScore( PBMC , features = gene_list_by_cluster[list_name], pool= NULL ,name= list_name , seed=19) 
}

#Remove the 1 at the end of metadata:
colnames(PBMC@meta.data) <- sapply(colnames(PBMC@meta.data), function(x) {
  if (grepl("^Cluster.*1$", x)) {
    return(gsub("1$", "", x))
  } else {
    return(x)
  }
})



#Look One by One:

cat(" \n \n")
cat("### Applying the different mouse ILC3 scores {.tabset .tabset-fade} \n\n")
cat(" \n \n")


for (names_signatures in names(gene_list_by_cluster)){
  
  cat(" \n \n")
  cat("#### ",  names_signatures, "\n")
  cat(" \n \n")
  
  p1 = FeaturePlot(PBMC, features = names_signatures, reduction = REDUCTION_OF_INTEREST)  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu" )) ) & NoAxes()
  
  cat(" \n \n")
  print(p1)
  cat(" \n \n")
  
  p2 = VlnPlot(PBMC, features = names_signatures ,sort = "decreasing", pt.size = 0, cols = palette_2) & stat_summary(fun.data=data_summary,color="black")
  
  cat(" \n \n")
  print(p2)
  cat(" \n \n")
  
  
  
  
}




#Saving Figures

if (DO_SAVE_FIGURE == TRUE){
  
for (names_signatures in paste0(names(gene_list_by_cluster))){
  p1 = FeaturePlot(PBMC, features = names_signatures, reduction = REDUCTION_OF_INTEREST)  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu" )) ) & NoAxes()

  
  p2 = VlnPlot(PBMC, features = names_signatures ,sort = "decreasing", pt.size = 0, cols = palette_2) & stat_summary(fun.data=data_summary,color="black")

  
  p3 = VlnPlot(PBMC, features = names_signatures ,sort = "decreasing", pt.size = 0, group.by = "seurat_clusters" , cols = palette_seurat_clusters) & stat_summary(fun.data=data_summary,color="black")

  p4 = VlnPlot(PBMC, features = names_signatures ,sort = "decreasing",  group.by = "site" , cols = palette_site ) & stat_summary(fun.data=data_summary,color="black")

  png(file=paste0( SAVE_FIG_PATH ,"FeaturePlot_Score_ILC3",names_signatures,".png") , width = 30, height = 25,  units = "cm", res=600 )
  print(p1)
  dev.off()
  
  png(file=paste0( SAVE_FIG_PATH,"_VlnPlot_ILC3",names_signatures,"per_clust.png") , width = 45, height = 25,  units = "cm", res=600 )
  print(p2)
  dev.off()
  
  png(file=paste0( SAVE_FIG_PATH,"_VlnPlot_ILC3",names_signatures,"per_seurat_clust.png") , width = 45, height = 25,  units = "cm", res=600 )
  print(p3)
  dev.off()
  
  png(file=paste0( SAVE_FIG_PATH,"_VlnPlot_ILC3",names_signatures,"per_site.png") , width = 45, height = 25,  units = "cm", res=600 )
  print(p4)
  dev.off()
  
  
}
  
}


#For all cells

PBMC = readRDS(paste0(PATH_PROJECT,SAMPLE_CHAP2_FL))

MONITORED_Markers = read.csv( SIGNATURES_ROMAGNANI_All, row.names = NULL)

PBMC$Cluster = as.factor(PBMC$seurat_clusters )
PBMC=SetIdent(PBMC, value= "Cluster")


MONITORED_Markers$cluster = as.factor(MONITORED_Markers$cluster)
levels(MONITORED_Markers$cluster) = c("ILC1"   ,      "ILC2"  ,       "ILC3"    ,     "progenitors"  ,"prol.ILC1"  , "prol.ILC2"  , "prol.ILC3" ,  "transitional")
palette_2 = hue_pal()(length(levels(PBMC@active.ident)))
names(palette_2) = levels(PBMC@active.ident)



#Extract the best
MONITORED_Markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_logFC > LOGFC_THR) %>%
  slice_head(n = NUMBER_TOP_SCORING_ILC3) %>%
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
  PBMC = AddModuleScore( PBMC , features = gene_list_by_cluster[list_name], pool= NULL ,name= list_name , seed=19) 
}

#Remove the 1 at the end of metadata:
colnames(PBMC@meta.data) <- sapply(colnames(PBMC@meta.data), function(x) {
  if (grepl("^Cluster.*1$", x)) {
    return(gsub("1$", "", x))
  } else {
    return(x)
  }
})



#Look One by One:

cat(" \n \n")
cat("### Applying the different Romagnani's signatures {.tabset .tabset-fade} \n\n")
cat(" \n \n")


for (names_signatures in names(gene_list_by_cluster)){
  
  cat(" \n \n")
  cat("#### ",  names_signatures, "\n")
  cat(" \n \n")
  
  p1 = FeaturePlot(PBMC, features = names_signatures, reduction = REDUCTION_OF_INTEREST)  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu" )) ) & NoAxes()
  
  cat(" \n \n")
  print(p1)
  cat(" \n \n")
  
  p2 = VlnPlot(PBMC, features = names_signatures ,sort = "decreasing", pt.size = 0, cols = palette_2) & stat_summary(fun.data=data_summary,color="black")
  
  cat(" \n \n")
  print(p2)
  cat(" \n \n")
  
  
  
  
}




#Saving Figures

if (DO_SAVE_FIGURE == TRUE){
  
  for (names_signatures in paste0(names(gene_list_by_cluster))){
    p1 = FeaturePlot(PBMC, features = names_signatures, reduction = REDUCTION_OF_INTEREST)  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu" )) ) & NoAxes()
    
    
    p2 = VlnPlot(PBMC, features = names_signatures ,sort = "decreasing", pt.size = 0, cols = palette_2) & stat_summary(fun.data=data_summary,color="black")
    
    
    p3 = VlnPlot(PBMC, features = names_signatures ,sort = "decreasing", pt.size = 0, group.by = "seurat_clusters" , cols = palette_seurat_clusters) & stat_summary(fun.data=data_summary,color="black")
    
    p4 = VlnPlot(PBMC, features = names_signatures ,sort = "decreasing",  group.by = "site" , cols = palette_site ) & stat_summary(fun.data=data_summary,color="black")
    
    png(file=paste0( SAVE_FIG_PATH ,"FeaturePlot_Score_Allpops",names_signatures,".png") , width = 30, height = 25,  units = "cm", res=600 )
    print(p1)
    dev.off()
    
    png(file=paste0( SAVE_FIG_PATH,"_VlnPlot_Allpops",names_signatures,"per_clust.png") , width = 45, height = 25,  units = "cm", res=600 )
    print(p2)
    dev.off()
    
    png(file=paste0( SAVE_FIG_PATH,"_VlnPlot_Allpops",names_signatures,"per_seurat_clust.png") , width = 45, height = 25,  units = "cm", res=600 )
    print(p3)
    dev.off()
    
    png(file=paste0( SAVE_FIG_PATH,"_VlnPlot_Allpops",names_signatures,"per_site.png") , width = 45, height = 25,  units = "cm", res=600 )
    print(p4)
    dev.off()
    
    
  }
  
}

