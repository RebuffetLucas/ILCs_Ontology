## @knitr SubClustand_Scoring_With_Mouse_FL_Signatures

cat(" \n \n")
cat("## Score only ILCP and Pre_ILC1 Pre_ILC2 and Pre_ILC3 ")
cat(" \n \n")


PBMC3 = readRDS(paste0(PATH_PROJECT,SAMPLE_Liu_Data_Used_ILC_Proge))

PBMC3$cluster = as.factor(PBMC3$cluster)
PBMC3=SetIdent(PBMC3, value= "cluster")

PBMC3$seurat_clusters = as.factor(PBMC3$seurat_clusters)
PBMC3$site = as.factor(PBMC3$site)


palette_2 = hue_pal()(length(levels(PBMC3@active.ident)))
names(palette_2) = levels(PBMC3@active.ident)


palette_seurat_clusters = hue_pal()(length(levels(PBMC3@meta.data[["seurat_clusters"]])))
names(palette_seurat_clusters) = levels(PBMC3@meta.data[["seurat_clusters"]])


palette_site = hue_pal()(length(levels(PBMC3@meta.data[["site"]])))
names(palette_site) = levels(PBMC3@meta.data[["site"]])


palette_cluster = hue_pal()(length(levels(PBMC3@meta.data[["cluster"]])))
names(palette_site) = levels(PBMC3@meta.data[["cluster"]])




p4  = DimPlot(PBMC3, reduction = "tsne", cols = palette_2 ) & ggtitle("  Liu cells Lymphoid cells only")




p5  = DimPlot(PBMC3, reduction = "tsne", group.by = "stage"  ) & ggtitle("  Liu cells Lymphoid cells stage")
p6  = DimPlot(PBMC3, reduction = "tsne", group.by = "batch"  ) & ggtitle("  Liu cells Lymphoid cells batch")
p7  = DimPlot(PBMC3, reduction = "tsne", group.by = "site"  ) & ggtitle("  Liu cells Lymphoid cells site")
p8  = DimPlot(PBMC3, reduction = "tsne", group.by = "seurat_clusters"  ) & ggtitle("  Liu cells Lymphoid cells site")


cat(" \n \n")
print(p4)
cat(" \n \n")



cat(" \n \n")
print( p5 + p6 + p7 + p8)
cat(" \n \n")


All_Markers = FindAllMarkers(PBMC3 , only.pos = FINDMARKERS_ONLYPOS, method= FINDMARKERS_METHOD , min.pct =  FINDMARKERS_MINPCT, logfc.threshold = FINDMARKERS_LOGFC_THR , verbose = TRUE)

All_Markers %>%
  group_by(cluster) %>%
  slice_max(n = FINDMARKERS_SHOWTOP, order_by = avg_log2FC) -> top10

p10 = DotPlot(PBMC3, features = unique(top10$gene) , cols = "RdBu") + NoLegend() + theme(axis.text.x = element_text(angle = 90)) &ggtitle("Spontaneous signatures")

cat(" \n \n")
print(p10)
cat(" \n \n")



#Barplots for proportions across tissues and stage


p5 = ggplot(PBMC3@meta.data, aes(x=site, fill= cluster)) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90))  + scale_fill_manual(values= palette_cluster)
p7 = ggplot(PBMC3@meta.data, aes(x=stage, fill= cluster)) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90))  + scale_fill_manual(values= palette_cluster)




cat(" \n \n")
print(p5 +  p7)
cat(" \n \n")


FeaturePlot(PBMC3, features = "KLRD1", reduction = "tsne")





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
  
  pdf(file=paste0( SAVE_FIG_PATH_Liu ,"UMAP_stage_batch_site_seurat_clusters",".pdf") , width = 30/2.54, height = 25/2.54 )
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
  PBMC3 = AddModuleScore( PBMC3 , features = human_gene_names_best[names_Sign] , pool= NULL ,name= names_Sign, seed=19) 
}


#Look One by One:

cat(" \n \n")
cat("### Applying the different mouse FL scores {.tabset .tabset-fade} \n\n")
cat(" \n \n")


for (names_signatures in paste0(names(human_gene_names_best)  , "1")){
  
  cat(" \n \n")
  cat("#### ",  names_signatures, "\n")
  cat(" \n \n")
  
  p1 = FeaturePlot(PBMC3, features = names_signatures, reduction = "tsne")  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu" )) ) & NoAxes()
  
  cat(" \n \n")
  print(p1)
  cat(" \n \n")
  

  p2 = VlnPlot(PBMC3, features = names_signatures ,sort = "decreasing", pt.size = 0, cols = palette_2) & stat_summary(fun.data=data_summary,color="black")
  
  cat(" \n \n")
  print(p2)
  cat(" \n \n")
  
  
  p3 = VlnPlot(PBMC3, features = names_signatures ,sort = "decreasing", pt.size = 0, group.by = "seurat_clusters" , cols = palette_seurat_clusters) & stat_summary(fun.data=data_summary,color="black")
  
  cat(" \n \n")
  print(p3)
  cat(" \n \n")
  
  p4 = VlnPlot(PBMC3, features = names_signatures ,sort = "decreasing",  group.by = "site" , cols = palette_site ) & stat_summary(fun.data=data_summary,color="black")
  
  cat(" \n \n")
  print(p4)
  cat(" \n \n")
  
  
  
  p3 = RidgePlot(PBMC3, features = names_signatures , sort= TRUE ,group.by = "seurat_clusters" , cols = palette_seurat_clusters) 
  
  cat(" \n \n")
  print(p3)
  cat(" \n \n")
  
  p4 = RidgePlot(PBMC3, features = names_signatures ,  group.by = "site" , cols = palette_site ) 
  
  cat(" \n \n")
  print(p4)
  cat(" \n \n")
  
  
  
}




#Saving Figures

if (DO_SAVE_FIGURE == TRUE){
  
  for (names_signatures in paste0(names(human_gene_names_best)  , "1")){
    p1 = FeaturePlot(PBMC3, features = names_signatures, reduction = "tsne")  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu" )) ) & NoAxes()
    
    
    p2 = VlnPlot(PBMC3, features = names_signatures ,sort = "decreasing", pt.size = 0, cols = palette_2) & stat_summary(fun.data=data_summary,color="black")
    
    
    p3 = VlnPlot(PBMC3, features = names_signatures ,sort = "decreasing", pt.size = 0, group.by = "seurat_clusters" , cols = palette_seurat_clusters) & stat_summary(fun.data=data_summary,color="black")
    
    p4 = VlnPlot(PBMC3, features = names_signatures ,sort = "decreasing",  group.by = "site" , cols = palette_site ) & stat_summary(fun.data=data_summary,color="black")
    
    png(file=paste0( SAVE_FIG_PATH ,"/Proge_Only_FeaturePlot_Score_",names_signatures,".png") , width = 30, height = 25,  units = "cm", res=600 )
    print(p1)
    dev.off()
    
    png(file=paste0( SAVE_FIG_PATH,"/Proge_Only_VlnPlot_",names_signatures,"per_clust.png") , width = 45, height = 25,  units = "cm", res=600 )
    print(p2)
    dev.off()
    
    png(file=paste0( SAVE_FIG_PATH,"/Proge_Only_VlnPlot_",names_signatures,"per_seurat_clust.png") , width = 45, height = 25,  units = "cm", res=600 )
    print(p3)
    dev.off()
    
    png(file=paste0( SAVE_FIG_PATH,"/Proge_Only_VlnPlot_",names_signatures,"per_site.png") , width = 45, height = 25,  units = "cm", res=600 )
    print(p4)
    dev.off()
    
    
    
    
    pdf(file=paste0( SAVE_FIG_PATH ,"/Proge_Only_FeaturePlot_Score_",names_signatures,".pdf") , width = 30/2.54, height = 25/2.54 )
    print(p1)
    dev.off()
    
    pdf(file=paste0( SAVE_FIG_PATH,"/Proge_Only_VlnPlot_",names_signatures,"per_clust.pdf") , width = 45/2.54, height = 25/2.54 )
    print(p2)
    dev.off()
    
    pdf(file=paste0( SAVE_FIG_PATH,"/Proge_Only_VlnPlot_",names_signatures,"per_seurat_clust.pdf") , width = 45/2.54, height = 25/2.54 )
    print(p3)
    dev.off()
    
    pdf(file=paste0( SAVE_FIG_PATH,"/Proge_Only_VlnPlot_",names_signatures,"per_site.pdf") , width = 45/2.54, height = 25/2.54 )
    print(p4)
    dev.off()
  }
  
}



