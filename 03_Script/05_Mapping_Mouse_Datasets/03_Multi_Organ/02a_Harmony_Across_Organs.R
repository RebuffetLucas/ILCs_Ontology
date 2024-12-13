## @knitr Harmony_Inte

## This script performs batch correction and integration of datasets from different organs using the Harmony algorithm.
## It starts by loading Seurat objects for fetal liver (FL), bone marrow (BM), and small intestine (SI), preparing them by normalizing data,
## finding variable features, and scaling. The script merges the datasets, corrects for batch effects across organs,
## and runs dimensionality reduction (PCA and UMAP) for visualization. Clustering is performed on the integrated data
## to identify cell populations across tissues in a harmonized space.

#Load data
data_mouse_FL = readRDS(paste0(PATH_PROJECT,SAMPLE_CHAP2_FL))
data_mouse_BM  = readRDS(paste0(PATH_PROJECT,SAMPLE_CHAP2_BM))
data_mouse_SI = readRDS(paste0(PATH_PROJECT,SAMPLE_ROMAGNANI))

genes_to_keep = Reduce(intersect , list(v1 = rownames(data_mouse_FL) , v2= rownames(data_mouse_BM) , v3 = rownames(data_mouse_SI)))


cat(" \n \n")
cat("## Datasets before batch correction \n\n")
cat(" \n \n")


# Have a quick look
cat(" \n \n")
print(DimPlot(data_mouse_BM, label = TRUE))
cat(" \n \n")

cat(" \n \n")
print(DimPlot(data_mouse_FL, label = TRUE))
cat(" \n \n")

cat(" \n \n")
print(DimPlot(data_mouse_SI, group.by = "annotation", label = TRUE))
cat(" \n \n")


#Data pre-treatment
data_mouse_BM$Cluster = data_mouse_BM$seurat_clusters
data_mouse_BM$Organ = "BM"

levels(data_mouse_FL$seurat_clusters) = c( "C1_ILCP" ,  "C2_Klrd1", "C3_ILC3P",  "C4_EILP",  "C5_Apoe")
data_mouse_FL$Cluster =  data_mouse_FL$seurat_clusters
data_mouse_FL$Organ = "FL"



levels(data_mouse_SI$annotation) = c( "ILC2" ,  "ILC2" , "ILC1", "ILC1",  "ILC3",   "ILC3", "ILCTransitional", "ILCP")
data_mouse_SI$Cluster = data_mouse_SI$annotation
data_mouse_SI$Organ = "SI"

data_mouse_SI@assays[["ADT"]] = NULL
data_mouse_SI@assays[["HTO"]] = NULL




############## By organ  #####################
#Merge them together
data_merged = merge(x=data_mouse_FL , y =c(data_mouse_BM , data_mouse_SI) , add.cell.ids = c("FL", "BM", "SI") )

data_merged = subset(data_merged, features= genes_to_keep)

DefaultAssay(data_merged) = "RNA"

#Abreviation
data_merged$Organ  = as.factor(data_merged$Organ)

#Triple ID
data_merged$Organ_Cluster = paste0(data_merged$Organ,"_", data_merged$Cluster)


#Find Variable features

List_Organ = SplitObject(data_merged, split.by = "Organ")

if (USE_MULTIBATCHNORM_BEFORE_HARMONY == TRUE){
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
  nfeatures = N_INTEGRATION_FEATURES_HARMONY,
  assay = rep("RNA", length(List_Organ)) ,
  verbose = TRUE
)




Merged_Seurat_Rescaled = merge(List_Organ[[1]] , y = List_Organ[-1])
#Merged_Seurat_Rescaled = merge(List_Organ[[2]] , y = c( List_Organ[[3]], List_Organ[[4]], List_Organ[[5]])) #Just a test here
Merged_Seurat_Rescaled = SetIdent(Merged_Seurat_Rescaled, value= "Organ")

Merged_Seurat_Rescaled=ScaleData(Merged_Seurat_Rescaled, features =rownames(Merged_Seurat_Rescaled) , do.scale = DATA_SCALE , do.center = DATA_CENTER )

Merged_Seurat_Rescaled=RunPCA(Merged_Seurat_Rescaled, features= FEATURES_RETAINED)

#Dataviz PCA
#DimPlot(Merged_Seurat_Rescaled, reduction = "pca", group.by = "Organ")

#Integration
Merged_Seurat_Rescaled$Organ = as.factor(Merged_Seurat_Rescaled$Organ)
Merged_Seurat_Rescaled=Merged_Seurat_Rescaled %>%   RunHarmony(VARIABLE_HARMONY, plot_convergence = FALSE, reduction = "pca" ,  verbose = TRUE , dims.use= 1:30 )

#UMAP
Merged_Seurat_Rescaled <- RunUMAP(Merged_Seurat_Rescaled, reduction = "harmony", dims= 1:30, reduction.name = "UMAP", reduction.key = "UMAP")

#Clustering
Merged_Seurat_Rescaled <- FindNeighbors(Merged_Seurat_Rescaled, reduction = "harmony", dims= 1:30)
Merged_Seurat_Rescaled <- FindClusters(Merged_Seurat_Rescaled, resolution = FINDCLUSTERS_RESOLUTION_AFTER_HARMONY, random.seed = SEED )


