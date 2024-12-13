## @knitr First_Look_All_Data
#def function

# +- std
data_summary <- function(x) {
  m <- median(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

cat(" \n \n")
cat("# Romagnani Data Analysis {.tabset .tabset-fade} \n\n")
cat(" \n \n")


cat(" \n \n")
cat("## First Look at All the data")
cat(" \n \n")



#Load the data
PBMC = readRDS(paste0(PATH_PROJECT,SAMPLE_ANALYSIS))
PBMC= UpdateSeuratObject(PBMC)
PBMC@meta.data[["Cluster"]] = PBMC@meta.data[[ANNOTATION_OF_INTEREST]]

PBMC = SetIdent(PBMC, value= "Cluster")
p1 = DimPlot(PBMC,  reduction = REDUCTION_OF_INTEREST) & ggtitle(" Bai et al Global Visualization") #First Vizualization




cat(" \n \n")
print(p1)
cat(" \n \n")





if (DO_SAVE_FIGURE == TRUE){
  
png(file=paste0( SAVE_FIG_PATH ,"/Global_Viz_AllData.png") , width = 30, height = 25,  units = "cm", res=600 )
print(p1)
dev.off()

pdf(file=paste0( SAVE_FIG_PATH ,"/Global_Viz_AllData.pdf") , width = 30/2.54, height = 25/2.54 )
print(p1)
dev.off()

  
}


