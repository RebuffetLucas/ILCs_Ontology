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
cat("## Score the cells with our signatures")
cat(" \n \n")


#Diagnostic with their granularity

#First remove the cells not used

PBMC2 = PBMC

PBMC2$Cluster = as.factor(PBMC2$Cluster)
PBMC2=SetIdent(PBMC2, value= "Cluster")


palette_2 = hue_pal()(length(levels(PBMC2@active.ident)))
names(palette_2) = levels(PBMC2@active.ident)





All_Markers = FindAllMarkers(PBMC2 , only.pos = FINDMARKERS_ONLYPOS, method= FINDMARKERS_METHOD , min.pct =  FINDMARKERS_MINPCT, logfc.threshold = FINDMARKERS_LOGFC_THR , verbose = TRUE)

All_Markers %>%
  group_by(cluster) %>%
  slice_max(n = FINDMARKERS_SHOWTOP, order_by = avg_log2FC) -> top10

p10 = DotPlot(PBMC2, features = unique(top10$gene) , cols = "RdBu") + NoLegend() + theme(axis.text.x = element_text(angle = 90)) &ggtitle("Spontaneous signatures")



p1 = DimPlot(PBMC2, group.by = "Cluster", reduction = REDUCTION_OF_INTEREST, cols = palette_2) & ggtitle(" Global Visualization") #First Vizualization




cat(" \n \n")
print(p1)
cat(" \n \n")


cat(" \n \n")
print(p10)
cat(" \n \n")





if (DO_SAVE_FIGURE == TRUE){
  png(file=paste0( SAVE_FIG_PATH ,"DotPlot_TOP_",FINDMARKERS_SHOWTOP ,".png") , width = 38, height = 12,  units = "cm", res=600 )
  print(p10)
  dev.off()

  pdf(file=paste0( SAVE_FIG_PATH ,"DotPlot_TOP_",FINDMARKERS_SHOWTOP ,".pdf") , width = 38/2.54, height = 12/2.54 )
  print(p10)
  dev.off()
}






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
  PBMC2 = AddModuleScore( PBMC2 , features = gene_list_by_cluster[list_name], pool= NULL ,name= list_name , seed=19) 
}

#Remove the 1 at the end of metadata:
colnames(PBMC2@meta.data) <- sapply(colnames(PBMC2@meta.data), function(x) {
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


for (names_signatures in names(gene_list_by_cluster)){
  
  cat(" \n \n")
  cat("#### ",  names_signatures, "\n")
  cat(" \n \n")
  
  p1 = FeaturePlot(PBMC2, features = names_signatures, reduction = REDUCTION_OF_INTEREST)  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu" )) ) & NoAxes()
  
  cat(" \n \n")
  print(p1)
  cat(" \n \n")
  
  p2 = VlnPlot(PBMC2, features = names_signatures ,sort = "decreasing", pt.size = 0, cols = palette_2) & stat_summary(fun.data=data_summary,color="black")
  
  cat(" \n \n")
  print(p2)
  cat(" \n \n")
  
  
  
  
}




#Saving Figures

if (DO_SAVE_FIGURE == TRUE){
  
for (names_signatures in paste0(names(gene_list_by_cluster))){
  p1 = FeaturePlot(PBMC2, features = names_signatures, reduction = REDUCTION_OF_INTEREST)  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu" )) ) & NoAxes()

  
  p2 = VlnPlot(PBMC2, features = names_signatures ,sort = "decreasing", pt.size = 0, cols = palette_2) & stat_summary(fun.data=data_summary,color="black")

  
  
  png(file=paste0( SAVE_FIG_PATH ,"FeaturePlot_Score_",names_signatures,".png") , width = 30, height = 25,  units = "cm", res=600 )
  print(p1)
  dev.off()
  
  png(file=paste0( SAVE_FIG_PATH,"_VlnPlot_",names_signatures,"per_clust.png") , width = 45, height = 25,  units = "cm", res=600 )
  print(p2)
  dev.off()
  
 
  

  
  pdf(file=paste0( SAVE_FIG_PATH ,"FeaturePlot_Score_",names_signatures,".pdf") , width = 30/2.54, height = 25/2.54  )
  print(p1)
  dev.off()
  
  pdf(file=paste0( SAVE_FIG_PATH,"_VlnPlot_",names_signatures,"per_clust.pdf") , width = 45/2.54, height = 25/2.54 )
  print(p2)
  dev.off()
  
  
}
  
}



