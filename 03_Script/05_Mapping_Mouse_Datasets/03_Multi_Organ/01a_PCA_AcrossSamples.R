## @knitr PCA_Across_Tissues

#Load data
data_mouse_FL = readRDS(paste0(PATH_PROJECT,SAMPLE_CHAP2_FL))
data_mouse_BM  = readRDS(paste0(PATH_PROJECT,SAMPLE_CHAP2_BM))
data_mouse_SI = readRDS(paste0(PATH_PROJECT,SAMPLE_ROMAGNANI))

cat(" \n \n")
cat("# Multi tissue Data Analysis {.tabset .tabset-fade} \n\n")
cat(" \n \n")

cat(" \n \n")
cat("## Global overview of the three datasets")
cat(" \n \n")

       
# Have a quick look
cat(" \n \n")
print(DimPlot(data_mouse_BM))
cat(" \n \n")

cat(" \n \n")
print(DimPlot(data_mouse_FL))
cat(" \n \n")

cat(" \n \n")
print(DimPlot(data_mouse_SI, group.by = "annotation"))
cat(" \n \n")


#Data pre-treatment
data_mouse_BM$Cluster = data_mouse_BM$seurat_clusters
data_mouse_BM$Organ = "BM"

levels(data_mouse_FL$seurat_clusters) = c( "C1_ILCP" ,  "C2_Klrd1", "C3_ILC3P",  "C4_EILP",  "C5_Apoe")
data_mouse_FL$Cluster = data_mouse_FL$seurat_clusters
data_mouse_FL$Organ = "FL"



levels(data_mouse_SI$annotation) = c( "ILC2" ,  "ILC2" , "ILC1", "ILC1",  "ILC3",   "ILC3", "ILCTransitional", "ILCP")
data_mouse_SI$Cluster = data_mouse_SI$annotation 
data_mouse_SI$Organ = "SI"

data_mouse_SI@assays[["ADT"]] = NULL
data_mouse_SI@assays[["HTO"]] = NULL




############## By organ  #####################
#Merge them together
data_merged = merge(x=data_mouse_FL , y =c(data_mouse_BM , data_mouse_SI) , add.cell.ids = c("FL", "BM", "SI") )
DefaultAssay(data_merged) = "RNA"

#Abreviation
data_merged$Organ  = as.factor(data_merged$Organ)


#Triple ID
data_merged$Organ_Cluster = paste0(data_merged$Organ,"_", data_merged$Cluster)


#Find Variable features

List_Organ = SplitObject(data_merged, split.by = "Organ")

if (USE_MULTIBATCHNORM == TRUE){
  sce_Merged= as.SingleCellExperiment(data_merged)
  sce_Merged_Rescaled = multiBatchNorm(sce_Merged,  batch = sce_Merged$Organ)
  data_merged= as.Seurat(sce_Merged_Rescaled)
  List_Organ=SplitObject(data_merged, split.by = "Organ")
  
} else {for (i in 1:length(List_Organ)) {
  List_Organ[[i]] <- NormalizeData(List_Organ[[i]], verbose = FALSE) #Maybe we need to change to MultiBatchNorm Here
}
}


for (i in 1:length(List_Organ)) {
  List_Organ[[i]] <- FindVariableFeatures(List_Organ[[i]], selection.method = "vst", 
                                             nfeatures = 2000, verbose = FALSE)
}


#Variable Features
FEATURES_RETAINED = SelectIntegrationFeatures(
  List_Organ,
  nfeatures = N_FEATURES_PCA,
  assay = rep("RNA", length(List_Organ)) ,
  verbose = TRUE
)


#rm(data_mouse_BM, data_mouse_FL, data_mouse_SI)
#gc()
Merged_Seurat_Rescaled = merge(List_Organ[[1]] , y = List_Organ[-1])
#Merged_Seurat_Rescaled = merge(List_Organ[[2]] , y = c( List_Organ[[3]], List_Organ[[4]], List_Organ[[5]])) #Just a test here
Merged_Seurat_Rescaled = SetIdent(Merged_Seurat_Rescaled, value= "Organ")
Merged_Seurat_Rescaled=ScaleData(Merged_Seurat_Rescaled, features =rownames(Merged_Seurat_Rescaled) , do.scale = DATA_SCALE , do.center = DATA_CENTER )


#Look at the mean
daframe_Average = AverageExpression(
  Merged_Seurat_Rescaled,
  assays = "RNA",
  features = FEATURES_RETAINED,
  return.seurat = FALSE,
  group.by = "Organ_Cluster",
  add.ident = NULL,
  layer = "data",
  slot= "scale.data",
  verbose = TRUE
)


mat = daframe_Average$RNA
mat= t(mat)


#Save the data
if (DO_SAVE_FIGURE == TRUE){
  saveRDS(mat, paste0(PATH_ANALYSIS_OUTPUT, "/mat_all_tum_blood.rds"))
}




