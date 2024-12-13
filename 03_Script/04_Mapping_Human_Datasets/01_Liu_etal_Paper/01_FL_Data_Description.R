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
cat("# Liu et al Data Analysis {.tabset .tabset-fade} \n\n")
cat(" \n \n")


cat(" \n \n")
cat("## First Look at All the data")
cat(" \n \n")



#Load the data
PBMC = readRDS(paste0(PATH_PROJECT,SAMPLE_Liu_Data_All))

PBMC = SetIdent(PBMC, value= "Cluster")
p1 = DimPlot(PBMC, group.by = "Cluster", reduction = "tsne") & ggtitle(" All Liu cells") #First Vizualization


p3 = DimPlot(PBMC, group.by = "cells_used", reduction = "tsne") & ggtitle(" All Liu cells used vs not used")  #First Vizualization



cat(" \n \n")
print(p1)
cat(" \n \n")


cat(" \n \n")
print(p3)
cat(" \n \n")



if (DO_SAVE_FIGURE == TRUE){
  
png(file=paste0( SAVE_FIG_PATH ,"/DotPlots_","_AllData_","TOP",FINDMARKERS_SHOWTOP,"_genesHorizontal.png") , width = 60, height = 15,  units = "cm", res=600 )
print(p1)
dev.off()


pdf(file=paste0( SAVE_FIG_PATH ,"/DotPlots_","_AllData_","TOP",FINDMARKERS_SHOWTOP,"_genesHorizontal.pdf") , width = 60/2.54, height = 15/2.54)
print(p1)
dev.off()

  
}


