## @knitr First_Look_Human_FL_Proge
#def function

# +- std
data_summary <- function(x) {
  m <- median(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

cat(" \n \n")
cat("# Human FL Science paper {.tabset .tabset-fade} \n\n")
cat(" \n \n")


cat(" \n \n")
cat("## First Look at All Human progenitors")
cat(" \n \n")



#Load the data

FL_YSpaper = readRDS(paste0(PATH_PROJECT,SAMPLE_FL_Science_YS_Atlas, ".rds"))


p1 = DimPlot(FL_YSpaper) & ggtitle(" All Human FL cells") #First Vizualization

cat(" \n \n")
print(p1)
cat(" \n \n")
#Diagnostic with their granularity

#Quick look at the proge / progenitors like signatures

FL_YSpaper2 = subset(FL_YSpaper ,  subset= broad_cell.labels %in% NOT_PROGENITORS, invert= TRUE     )

if (DO_REMOVE_NK == TRUE){
  FL_YSpaper2 = subset(FL_YSpaper2 ,  subset= cell.labels == "NK", invert= TRUE     )
}



p2 = DimPlot(FL_YSpaper2) & ggtitle(" Subset FL cells")



cat(" \n \n")
print(p2)
cat(" \n \n")



All_Markers = FindAllMarkers(FL_YSpaper2 , only.pos = FINDMARKERS_ONLYPOS, method= FINDMARKERS_METHOD , min.pct =  FINDMARKERS_MINPCT, logfc.threshold = FINDMARKERS_LOGFC_THR , verbose = TRUE)

All_Markers %>%
  group_by(cluster) %>%
  slice_max(n = FINDMARKERS_SHOWTOP, order_by = avg_log2FC) -> top10

p3 = DotPlot(FL_YSpaper2, features = unique(top10$gene) , cols = "RdBu") + NoLegend() + theme(axis.text.x = element_text(angle = 90)) &ggtitle("Spontaneous signatures")

cat(" \n \n")
print(p3)
cat(" \n \n")


if (DO_SAVE_FIGURE == TRUE){
  
png(file=paste0( SAVE_FIG_PATH ,"/DotPlots/","TOP",FINDMARKERS_SHOWTOP,"FL_YSpaper","genesHorizontal.png") , width = 60, height = 15,  units = "cm", res=600 )
print(p3)
dev.off()

pdf(file=paste0( SAVE_FIG_PATH ,"/DotPlots/","TOP",FINDMARKERS_SHOWTOP,"FL_YSpaper","genesHorizontal.pdf") , width = 60/2.54, height = 15/2.54)
print(p3)
dev.off()


png(file=paste0( SAVE_FIG_PATH ,"/UMAP_All_Cells.png") , width = 25, height = 15,  units = "cm", res=600 )
print(p1)
dev.off()

pdf(file=paste0( SAVE_FIG_PATH ,"/UMAP_All_Cells.pdf") , width = 25/2.54, height = 15/2.54)
print(p1)
dev.off()


png(file=paste0( SAVE_FIG_PATH ,"/UMAP_after_subseting.png") , width = 20, height = 15,  units = "cm", res=600 )
print(p2)
dev.off()

pdf(file=paste0( SAVE_FIG_PATH ,"/UMAP_after_subseting.pdf") , width = 20/2.54, height = 15/2.54)
print(p2)
dev.off()


}


