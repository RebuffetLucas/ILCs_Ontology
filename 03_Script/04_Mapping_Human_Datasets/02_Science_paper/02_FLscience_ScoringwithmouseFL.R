## @knitr Scoring_With_Mouse_FL_Signatures

#Score Human FL with mouse FL signatures



cat(" \n \n")
cat("## Score Human FL Progenitors with mouse FL signatures")
cat(" \n \n")



#Bars for ViolinPlots
data_summary <- function(x) {
  m <- median(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
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




PBMC2 =  FL_YSpaper2

#Remove the + that generate problems
names(human_gene_names_best) <- gsub("\\+", "", names(human_gene_names_best))


#Global UMAP
p_all = DimPlot(PBMC2)

cat(" \n \n")
print(p_all)
cat(" \n \n")




#Score the pops
for (names_Sign in names(human_gene_names_best)){
  PBMC2 = AddModuleScore( PBMC2 , features = human_gene_names_best[names_Sign] , pool= NULL ,name= names_Sign, seed=19) 
}


#Look One by One:

cat(" \n \n")
cat("### Applying the different mouse FL scores {.tabset .tabset-fade} \n\n")
cat(" \n \n")

palette_2 = hue_pal()(length(levels(PBMC2@active.ident)))
names(palette_2) = levels(PBMC2@active.ident)



for (names_signatures in paste0(names(human_gene_names_best)  , "1")){
  
  cat(" \n \n")
  cat("#### ",  names_signatures, "\n")
  cat(" \n \n")
  
  p1 = FeaturePlot(PBMC2, features = names_signatures)  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu" )) ) & NoAxes()
  
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
  
for (names_signatures in paste0(names(human_gene_names_best)  , "1")){
  p1 = FeaturePlot(PBMC2, features = names_signatures)  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu" )) ) & NoAxes()
  print(p1)
  p2 = VlnPlot(PBMC2, features = names_signatures ,sort = "decreasing", pt.size = 0, cols = palette_2) & stat_summary(fun.data=data_summary,color="black")
  print(p2)
  
  png(file=paste0( SAVE_FIG_PATH ,"/FL_YSpaper_FeaturePlot",names_signatures,".png") , width = 30, height = 25,  units = "cm", res=600 )
  print(p1)
  dev.off()
  
  png(file=paste0( SAVE_FIG_PATH ,"/FL_YSpaper_VlnPlot",names_signatures,".png") , width = 45, height = 20,  units = "cm", res=600 )
  print(p2)
  dev.off()
  
  pdf(file=paste0( SAVE_FIG_PATH ,"/FL_YSpaper_FeaturePlot",names_signatures,".pdf") , width = 30/2.54, height = 25/2.54  )
  print(p1)
  dev.off()
  
  pdf(file=paste0( SAVE_FIG_PATH ,"/FL_YSpaper_VlnPlot",names_signatures,".pdf") , width = 45/2.54, height = 20/2.54 )
  print(p2)
  dev.off()
  
  
  
}

}


if (DO_SAVE_FIGURE == TRUE){
  
p1 = DimPlot(FL_YSpaper2) & ggtitle(" Human FL cells progenitors") & NoAxes()

png(file=paste0( SAVE_FIG_PATH ,"/UMAP_No_Axis.png") , width = 30, height = 25,  units = "cm", res=600 )
print(p1)
dev.off()

pdf(file=paste0( SAVE_FIG_PATH ,"/UMAP_No_Axis.pdf") , width = 30/2.54, height = 25/2.54  )
print(p1)
dev.off()


p2 = DimPlot(FL_YSpaper2) & ggtitle(" Human FL cells progenitors") 

png(file=paste0( SAVE_FIG_PATH ,"/UMAP_with_Axis.png") , width = 30, height = 25,  units = "cm", res=600 )
print(p2)
dev.off()

pdf(file=paste0( SAVE_FIG_PATH ,"/UMAP_with_Axis.pdf") , width = 30/2.54, height = 25/2.54  )
print(p2)
dev.off()


}



