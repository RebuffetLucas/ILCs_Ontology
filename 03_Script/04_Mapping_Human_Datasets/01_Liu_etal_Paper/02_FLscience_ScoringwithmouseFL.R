## @knitr Scoring_With_Mouse_FL_Signatures


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
cat("## Score Human Lymphoid cells")
cat(" \n \n")


#Diagnostic with their granularity

#First remove the cells not used

PBMC2 = readRDS(paste0(PATH_PROJECT,SAMPLE_Liu_Data_Used))

PBMC2$Cluster = as.factor(PBMC2$Cluster)
PBMC2=SetIdent(PBMC2, value= "Cluster")

PBMC2$seurat_clusters = as.factor(PBMC2$seurat_clusters)
PBMC2$site = as.factor(PBMC2$site)


palette_2 = hue_pal()(length(levels(PBMC2@active.ident)))
names(palette_2) = levels(PBMC2@active.ident)


palette_seurat_clusters = hue_pal()(length(levels(PBMC2@meta.data[["seurat_clusters"]])))
names(palette_seurat_clusters) = levels(PBMC2@meta.data[["seurat_clusters"]])


palette_site = hue_pal()(length(levels(PBMC2@meta.data[["site"]])))
names(palette_site) = levels(PBMC2@meta.data[["site"]])




p4  = DimPlot(PBMC2, reduction = "tsne", cols = palette_2 ) & ggtitle("  Liu cells Lymphoid cells only")




p5  = DimPlot(PBMC2, reduction = "tsne", group.by = "stage"  ) & ggtitle("  Liu cells Lymphoid cells stage")
p6  = DimPlot(PBMC2, reduction = "tsne", group.by = "batch"  ) & ggtitle("  Liu cells Lymphoid cells batch")
p7  = DimPlot(PBMC2, reduction = "tsne", group.by = "site"  ) & ggtitle("  Liu cells Lymphoid cells site")
p8  = DimPlot(PBMC2, reduction = "tsne", group.by = "seurat_clusters"  ) & ggtitle("  Liu cells Lymphoid cells site")


cat(" \n \n")
print(p4)
cat(" \n \n")



cat(" \n \n")
print( p5 + p6 + p7 + p8)
cat(" \n \n")


All_Markers = FindAllMarkers(PBMC2 , only.pos = FINDMARKERS_ONLYPOS, method= FINDMARKERS_METHOD , min.pct =  FINDMARKERS_MINPCT, logfc.threshold = FINDMARKERS_LOGFC_THR , verbose = TRUE)

All_Markers %>%
  group_by(cluster) %>%
  slice_max(n = FINDMARKERS_SHOWTOP, order_by = avg_log2FC) -> top10

p10 = DotPlot(PBMC2, features = unique(top10$gene) , cols = "RdBu") + NoLegend() + theme(axis.text.x = element_text(angle = 90)) &ggtitle("Spontaneous signatures")

cat(" \n \n")
print(p10)
cat(" \n \n")



cat(" \n \n")
print(p4)
cat(" \n \n")




if (DO_SAVE_FIGURE == TRUE){
  png(file=paste0( SAVE_FIG_PATH_Liu ,"DimPlot",".png") , width = 30, height = 25,  units = "cm", res=600 )
  print(p4)
  dev.off()
  
  png(file=paste0( SAVE_FIG_PATH_Liu ,"UMAP_stage_batch_site_seurat_clusters",".png") , width = 30, height = 25,  units = "cm", res=600 )
  print( p5 + p6 + p7 + p8)
  dev.off()
  
  pdf(file=paste0( SAVE_FIG_PATH_Liu ,"DimPlot",".pdf") , width = 30/2.54, height = 25/2.54 )
  print(p4)
  dev.off()
  
  pdf(file=paste0( SAVE_FIG_PATH_Liu ,"UMAP_stage_batch_site_seurat_clusters",".pdf") , width = 30/2.54, height = 25/2.54  )
  print( p5 + p6 + p7 + p8)
  dev.off()
  
  
  
  
}






#Import signatures
FILE_SIGNATURES_PATH = MOUSE_TO_HUM_FL_PATH
  
###Markers for scoring of subsets
# getting data from sheets
  sheets_names <- openxlsx::getSheetNames(FILE_SIGNATURES_PATH)


#def function

read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}




#Monitored markers
MONITORED_Markers = read_excel_allsheets(FILE_SIGNATURES_PATH)


human_gene_names_best <- lapply(MONITORED_Markers, function(df) {
  if(nrow(df) >= NUMBER_TOP_SCORING) {
    return(df$human.gene.name[1:NUMBER_TOP_SCORING])
  } else {
    return(df$human.gene.name[1:nrow(df)])
  }
})




#Remove the + that generate problems
names(human_gene_names_best) <- gsub("\\+", "", names(human_gene_names_best))



#Score the pops
for (names_Sign in names(human_gene_names_best)){
  PBMC2 = AddModuleScore( PBMC2 , features = human_gene_names_best[names_Sign] , pool= NULL ,name= names_Sign, seed=19) 
}


#Look One by One:

cat(" \n \n")
cat("### Applying the different mouse FL scores {.tabset .tabset-fade} \n\n")
cat(" \n \n")


for (names_signatures in paste0(names(human_gene_names_best)  , "1")){
  
  cat(" \n \n")
  cat("#### ",  names_signatures, "\n")
  cat(" \n \n")
  
  p1 = FeaturePlot(PBMC2, features = names_signatures, reduction = "tsne")  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu" )) ) & NoAxes()
  
  cat(" \n \n")
  print(p1)
  cat(" \n \n")
  
  p2 = VlnPlot(PBMC2, features = names_signatures ,sort = "decreasing", pt.size = 0, cols = palette_2) & stat_summary(fun.data=data_summary,color="black")
  
  cat(" \n \n")
  print(p2)
  cat(" \n \n")
  
  
  p3 = VlnPlot(PBMC2, features = names_signatures ,sort = "decreasing", pt.size = 0, group.by = "seurat_clusters" , cols = palette_seurat_clusters) & stat_summary(fun.data=data_summary,color="black")
  
  cat(" \n \n")
  print(p3)
  cat(" \n \n")
  
  p4 = VlnPlot(PBMC2, features = names_signatures ,sort = "decreasing",  group.by = "site" , cols = palette_site ) & stat_summary(fun.data=data_summary,color="black")
  
  cat(" \n \n")
  print(p4)
  cat(" \n \n")
  
  
  
}




#Saving Figures

if (DO_SAVE_FIGURE == TRUE){
  
for (names_signatures in paste0(names(human_gene_names_best)  , "1")){
  p1 = FeaturePlot(PBMC2, features = names_signatures, reduction = "tsne")  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu" )) ) & NoAxes()

  
  p2 = VlnPlot(PBMC2, features = names_signatures ,sort = "decreasing", pt.size = 0, cols = palette_2) & stat_summary(fun.data=data_summary,color="black")

  
  p3 = VlnPlot(PBMC2, features = names_signatures ,sort = "decreasing", pt.size = 0, group.by = "seurat_clusters" , cols = palette_seurat_clusters) & stat_summary(fun.data=data_summary,color="black")

  p4 = VlnPlot(PBMC2, features = names_signatures ,sort = "decreasing",  group.by = "site" , cols = palette_site ) & stat_summary(fun.data=data_summary,color="black")

  png(file=paste0( SAVE_FIG_PATH ,"/FeaturePlot_Score_",names_signatures,".png") , width = 30, height = 25,  units = "cm", res=600 )
  print(p1)
  dev.off()
  
  png(file=paste0( SAVE_FIG_PATH,"/_VlnPlot_",names_signatures,"per_clust.png") , width = 45, height = 25,  units = "cm", res=600 )
  print(p2)
  dev.off()
  
  png(file=paste0( SAVE_FIG_PATH,"/_VlnPlot_",names_signatures,"per_seurat_clust.png") , width = 45, height = 25,  units = "cm", res=600 )
  print(p3)
  dev.off()
  
  png(file=paste0( SAVE_FIG_PATH,"/_VlnPlot_",names_signatures,"per_site.png") , width = 45, height = 25,  units = "cm", res=600 )
  print(p4)
  dev.off()
  
  
  
  pdf(file=paste0( SAVE_FIG_PATH ,"/FeaturePlot_Score_",names_signatures,".pdf") , width = 30/2.54, height = 25/2.54 )
  print(p1)
  dev.off()
  
  pdf(file=paste0( SAVE_FIG_PATH,"/_VlnPlot_",names_signatures,"per_clust.pdf") , width = 45/2.54, height = 25/2.54 )
  print(p2)
  dev.off()
  
  pdf(file=paste0( SAVE_FIG_PATH,"/_VlnPlot_",names_signatures,"per_seurat_clust.pdf") , width = 45/2.54, height = 25/2.54 )
  print(p3)
  dev.off()
  
  pdf(file=paste0( SAVE_FIG_PATH,"/_VlnPlot_",names_signatures,"per_site.pdf") , width = 45/2.54, height = 25/2.54 )
  print(p4)
  dev.off()
  
  
}
  
}



